# Boundary collision models and shared execution

The generic and material paths differ only in collision-model operations:
event selection, emission velocity sampling, and reflected velocity sampling.
Particle filtering, impact reconstruction, allocation, creation, wall-charge
deposition, deletion and diagnostics all use the same driver implementation.

## Configuration

```ini
electrons.analytic_wall.ceramic_wall.material = ceramic
electrons.analytic_wall.ceramic_wall.secondary_electron_species = electrons
electrons.analytic_wall.ceramic_wall.secondary_electron_temperature_eV = 3.0
electrons.analytic_wall.ceramic_wall.deposit_wall_charge = 1
```

`ceramic` is currently the only material name; unknown names are input errors.
The target defaults to the incident species, temperature to 3 eV and deposition
to enabled. Material selection takes precedence over `behaviors`, `p_*`,
`wall_temperature` and conversion settings for that species/wall pair. Without
a material, the generic parser-based configuration is used. Without either
configuration the policy is inactive. The Hall slot example selects ceramic for
electron/ceramic interaction while its other pairs use generic policies.

The host retains `std::variant` for material configuration and converts it to a
trivially copyable `MaterialInteractionDevice`. Only collision-model functions
perform device-side material dispatch; no `std::visit` executes on the device.

## Shared preprocessing

Particle positions, velocities, normals and sampled emission arrays use
`ParticleVector`, whose three components are `amrex::ParticleReal`.
`AnalyticBoundaryPosition` is an alias of this type. Reflection, normalization
and basis transforms also operate in particle precision, without passing through
`amrex::XDim3`. The latter remains only for field-gather geometry such as the
inverse mesh spacing and mesh-box origin.

Thermal sampling scales dimensionless normal draws in particle precision. The
zero-drift normal flux uses its exact inverse CDF in that precision; this retains
the same Maxwellian flux spectrum but may consume a different RNG sequence.

`CollectOutsideParticleIndices` counts valid particles outside the wall, then
collects their original indices into an exactly sized array. This uses a reduction
and an atomic append with `amrex::For`, with no full-particle flags or offsets.
Empty tiles and zero-candidate tiles skip subsequent work. Append order is not
guaranteed, so random draws need not match the former full-particle traversal.

`PreprocessWallImpacts` computes a position at impact, corrected physical velocity,
normal and remaining time for each candidate. Coordinates are in metres,
velocities in physical m/s, and remaining time in seconds. It preserves the
existing nonrelativistic bisection and field-regather correction. Preprocessing
completes before either model is sampled or any source storage is resized.

Event arrays use compact index `k`; particle attributes use original index
`i = outside_indices[k]`. The lists are rebuilt per species/wall/tile in the
existing wall order. The source must not be reordered while consuming a list.

## One particle lifecycle

1. `SelectWallBehavior` samples one event per candidate, using either generic
   expressions or `CeramicInteraction::SelectEvent`. The selected `WallBehavior`
   is stored once and determines multiplicity through `WallEmissionCount`.
2. `CreateWallProducts` counts the actual products of each event type and uses a
   prefix sum to allocate exactly the required slots. Both models use the same
   code, including `SmartCopy`, weight inheritance, runtime-attribute initialization
   and new IDs. Source pointers are reacquired after same-species resizing.
3. `SampleWallEmission` invokes the selected collision model for the complete
   event's outgoing velocities. The common creation code writes each velocity and
   advances its position from impact through the stored remaining time. The
   chosen event is not sampled again.
4. After all products have been created and initialized, the common deposition
   kernel deposits net retained charge at the original impact point.
5. The common update kernel either samples a reflected velocity and updates the
   primary's position, or marks the primary invalid. No primary is deleted before
   its products are complete. Redistribution later removes invalid entries.
6. Existing wall-current and event-count diagnostics consume the same event list
   for both models.

There is no material-specific maximum-capacity reservation, emission append
counter, particle-write callback, deposition kernel or deletion implementation.
Creation counts, offsets and diagnostic loops all scale with candidate count.
Only collision-model operations branch on the material tag. The public behavior
vocabulary and its current maximum of two emissions are in `WallBehavior.h`.

## Preserving outgoing spectra

`SampleEmission` returns the outgoing velocities for one already-selected event.
The driver copies those samples directly; it never replaces them with a common
thermal draw, rescales their energy or independently resamples the event. Sampling
the entire event also leaves room for a material to impose correlations between
its emitted particles.

The ceramic preset retains its probability curves, with incident energy E in eV:

- Absorption: `0.5 exp(-(E/43.4592)^2)`.
- Specular reflection: `0.5 exp(-(E/30)^2)`.
- Two secondaries: `1 - exp(-(E/127.8958)^2)`.
- One secondary: the remaining probability.

`CeramicInteraction::SampleEmission` independently draws each secondary from the
existing Maxwellian flux operator using its own configured kT. SEE2 therefore
still has two independent draws. At the default kT = 3 eV the mean emitted kinetic
energy is 6 eV. There is no event energy truncation or incidence-angle-dependent
yield in this inherited ceramic preset. Reflection uses corrected impact velocity.

The generic model retains its configured thermal distribution, including the
legacy shared velocity draw for SEE2 and stored-velocity specular reflection.
These are differences in the existing model prescriptions, not in allocation,
particle writes or deposition. Reorganizing RNG calls need not preserve the old
bitwise random sequence even though the distributions are retained.

## Shared charge accounting

Both models use `DepositWallChargeToNodes` and `WallChargeGrid` from
`HallWallCharge.H`, depositing `(q_in - sum(q_out)) * weight` at the hit point.
The helper uses trilinear nodal weights, inverse cell volume and
`HostDevice::Atomic::Add`. For equally weighted electrons (e > 0):

| Event | Retained charge |
| --- | --- |
| Absorption | -e * weight |
| Specular reflection | 0 |
| SEE1 | 0 |
| SEE2 | +e * weight |

The existing `product_charge` override still applies to generic conversions.
`deposit_wall_charge = 0` disables surface deposition for either model without
changing its particle or current-diagnostic behavior. The driver obtains the
persistent field with `HallWallCharge::get`, and the existing material Poisson
path handles synchronization and source addition (currently at level 0).
Deposition uses `amrex::For` to avoid CPU SIMD independence assumptions for shared
nodes; independent sampling/creation uses the appropriate AMReX particle kernels.
