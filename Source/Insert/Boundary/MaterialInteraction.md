# Analytic-wall preprocessing, interaction handlers and finalization

The generic, ceramic and stainless-steel paths own their complete reflection/emission stage.
They share preprocessing and finalization, not an emission buffer interface.

```text
CollectOutsideParticleIndices -> PreprocessWallImpacts
    -> ProcessGenericWallPolicy OR material.Process
    -> common charge deposition -> primary invalidation -> diagnostics
```

## Configuration

```ini
electrons.analytic_wall.ceramic_wall.material = ceramic
electrons.analytic_wall.ceramic_wall.secondary_electron_species = electrons
electrons.analytic_wall.ceramic_wall.secondary_electron_temperature_eV = 3.0
electrons.analytic_wall.ceramic_wall.deposit_wall_charge = 1
```

Material selection takes precedence over the same species/wall's generic
`behaviors` and probability expressions. Without a material, generic policies
retain their existing behavior. Materials are `ceramic` and `stainless_steel`;
unknown names are input errors. The target defaults to the incident species,
ceramic temperature to 3 eV and deposition to enabled. Both material species must be
physical electrons. The Hall slot example uses ceramic for electrons at the
ceramic wall, with generic policies for the other species/wall pairs.

The host uses `std::visit` to select a material's `Process` method after shared
preprocessing. Each method owns allocation and launches its own device kernels;
only its sampler data is captured by those kernels. No variant or virtual call
runs on the device. The former common device material dispatcher is removed.

## Shared preprocessing and records

`CollectOutsideParticleIndices` counts valid particles outside the wall and
atomically appends their indices to an exactly sized array. It uses a reduction
and `amrex::For`, without full-particle flags or offsets. Empty lists skip later
work. The append does not guarantee ordering or the previous RNG sequence.

`PreprocessWallImpacts` stores position at impact, corrected physical velocity,
normal and remaining time in `WallImpact`. These calculations are unchanged.
Coordinates, velocities and directions use `ParticleVector`/`ParticleReal`;
only field-gather mesh geometry uses `XDim3`. Every downstream array is indexed
by compact event index k; source attributes use i = outside_indices[k].

Both handlers fill `WallInteractionRecord` for every candidate:

- `event`: diagnostic category.
- `keep_primary`: whether the handler retained and updated the incident particle.
- `emitted_count`: actual newly created macro-particle count, excluding the primary.
- `outgoing_charge`: total outgoing charge in C, including macro weights and any
  retained primary. Generic conversions use the configured `product_charge`
  override when present.

Handlers complete all products, IDs, runtime attributes and record writes before
returning. They do not invalidate primaries or deposit wall charge. The common
record imposes no emission-array size. Other materials can implement a different
multiplicity/energy sampler and allocation scheme without changing generic SEE.
Diagnostic categories for any new physical channels must be defined explicitly.

## Generic processing

`ProcessGenericWallPolicy` samples the configured parser probabilities once,
then calls its private `CreateGenericWallProducts` for conversion, SEE1 and SEE2.
The helper computes exact counts and offsets, expands the destination tile and
writes products using `SmartCopy`. It reacquires source data after a same-species
resize and initializes runtime attributes and IDs before returning.

The existing generic thermal distribution is preserved, including the shared
velocity draw for the two products of SEE2. Specular reflection continues to use
the stored velocity; diffuse reflection uses the configured wall temperature.
The handler performs reflection and position updates before common finalization.

## Ceramic processing

`CeramicInteraction::Process` independently owns its event-selection kernel,
emission counts, prefix offsets, destination allocation, reflection and emission
kernel. It allocates the exact combined SEE1/SEE2 product count. It calls its own
`SampleEmission` once per selected emission event and writes those velocities
directly into the destination. Only low-level copying, initialization and ID tools
are reused; it never invokes the generic product-creation helper.

The ceramic sampler retains its original probability curves (E in eV):

- Absorption: `0.5 exp(-(E/43.4592)^2)`.
- Specular reflection: `0.5 exp(-(E/30)^2)`.
- SEE2: `1 - exp(-(E/127.8958)^2)`.
- SEE1: the remaining probability.

Each secondary independently samples the same Maxwellian flux at the material's
configured kT. SEE2 retains two independent velocities. At kT = 3 eV the mean
emitted kinetic energy is 6 eV. There is no event-wise energy truncation in this
preset. Reflection uses the corrected impact velocity. The two-element velocity
array is now local to the ceramic model, not a requirement on other materials.

All container expansion occurs on the host between kernels. Both handlers retain
original source indices, reacquire pointers after resizing and avoid processing
new products as original candidates. Particle redistribution/compaction must not
occur until the current compact list has been consumed.

## Common finalization

After the selected handler returns, the driver reacquires source pointers and
uses only the completed records to perform finalization:

1. Deposit `q_in * original_weight - outgoing_charge` at the original hit point
   with `DepositWallChargeToNodes`, when deposition is enabled. Zero net charge
   is skipped. This scatter uses `amrex::For` and the existing host/device atomics.
2. Mark the primary invalid if `keep_primary` is false. All emission work is
   already complete; redistribution later removes invalid particles.
3. Accumulate current and event diagnostics. Emitted-secondary totals are reduced
   from actual `emitted_count`, rather than reconstructed from SEE1/SEE2 labels.

For equal-weight electrons the retained charge is -e*w for absorption, zero for
specular reflection/SEE1, and +e*w for SEE2. Neither deposition nor current
accounting needs to know how a material sampled its spectrum or allocated products.
The existing persistent field and material Poisson integration are unchanged.
`deposit_wall_charge = 0` disables deposition without suppressing interaction or
current diagnostics. RNG call ordering can differ after the reorganization;
the physical probability and emission distributions are retained.

## Stainless steel: Furman--Pivi SLAC 304 preset

`StainlessSteelInteraction` implements the unconditioned, etched/passivated
SLAC 304 parameter set in `不锈钢SEE模型_Furman-Pivi.md`, Tables I/II. It is a
specific surface preset, not a universal stainless-steel model.

```ini
electrons.analytic_wall.anode.material = stainless_steel
electrons.analytic_wall.anode.secondary_electron_species = electrons
electrons.analytic_wall.anode.binomial_trials = 10
electrons.analytic_wall.anode.deposit_wall_charge = 1
```

Only `binomial_trials` adjusts the multiplicity distribution; it defaults to 10
and must be in [1,10], matching the available spectral table. Fitted yield and
spectral coefficients are encoded in the device sampler. The ceramic
`secondary_electron_temperature_eV` does not control this material's spectrum.
The existing Hall example's generic anode policy is not changed automatically.

The sampler follows section 4.2's binomial alternative: first select elastic or
rediffused backscatter with probabilities delta_e and delta_r; otherwise sample
n ~ Binomial(M, delta_ts / ((1-delta_e-delta_r) M)). n=0 is absorption.
This is not the document's alternative Poisson multiplicity law. There is no
Poisson tail truncation or extrapolation of the n<=10 spectral table.

Elastic backscatter samples the finite-width (sigma=1.9 eV) Gaussian restricted
to [0,E0]. Rediffused energy is E0 U^(1/1.4). For true secondaries, the exact
conditional joint Gamma law is sampled using its independent total/fraction
factorization: S ~ Gamma(n*p_n, epsilon_n) conditional on S<=E0, with fractions
~ Dirichlet(p_n,...,p_n). A power-law rejection proposal samples the truncated
total efficiently at low E0; this is mathematically equivalent to rejecting
entire sets of independent Gamma energies. Rejection never changes the selected
multiplicity. A roundoff-only correction prevents rotated velocities from
exceeding the total energy budget. All three channels use the cosine angular
law about the normal pointing into the simulation domain.

The sampler assumes valid incident states from preprocessing and applies no
energy or angle fit-range cutoff. It evaluates the same yield formulas outside
the fitted range. The model requires A>=0 (A=0 only if delta_ts=0) and mu<M.
These conditions can fail even within the fitted range: E0=100 eV, theta=30 degrees gives mu~14.78,
invalid for M=10. Evaluation writes a device status; SelectEvent maps any
invalid status to complete absorption (zero emitted electrons), without drawing
random numbers or aborting the run. This fallback applies to the entire impact,
including its backscatter channels. It is an additional physical approximation,
not a fit-range cutoff or a probability clamp, and reduces the effective yield
for these impacts. Common finalization deposits the full incident charge when
wall-charge deposition is enabled, deletes the incident and counts an absorption.

Each emitted electron has the copied incident macro-weight and a new ID.
All channels replace the incident after emission, including backscatter.
Common deposition uses q_in*w_in - sum(q_out*w_out): absorption deposits the
incident charge, one returning electron gives zero net charge, and n true
secondaries give (n-1)*e*w. The diagnostic `removed` count includes replaced
backscattered primaries; net current is computed from charge balance.
