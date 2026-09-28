#pragma once

#include "Insert/Boundary/CeramicInteraction.h"

#include <type_traits>
#include <variant>

namespace Insert {

/** Host-side material configuration. */
using MaterialInteraction = std::variant<CeramicInteraction>;

/** Device-side collision-model dispatch. The host converts its variant once;
 * the common driver dispatches event and outgoing-state sampling on the device.
 * No material-specific container management or deposition is performed here.
 */
class MaterialInteractionDevice
{
public:
    MaterialInteractionDevice () = default;

    explicit MaterialInteractionDevice (CeramicInteraction const& ceramic)
        : m_kind(Kind::ceramic), m_ceramic(ceramic) {}

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool IsActive () const noexcept { return m_kind != Kind::none; }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    WallBehavior SelectEvent (
        ParticleVector const& velocity, amrex::RandomEngine const& engine) const noexcept
    {
        switch (m_kind) {
        case Kind::ceramic: return m_ceramic.SelectEvent(velocity, engine);
        case Kind::none: return WallBehavior::none;
        }
        return WallBehavior::none;
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::GpuArray<ParticleVector, max_wall_emissions> SampleEmission (
        WallBehavior const event, ParticleVector const& velocity,
        ParticleVector const& normal, amrex::RandomEngine const& engine) const noexcept
    {
        switch (m_kind) {
        case Kind::ceramic: return m_ceramic.SampleEmission(event, velocity, normal, engine);
        case Kind::none: return {};
        }
        return {};
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    ParticleVector SampleReflection (
        ParticleVector const& velocity, ParticleVector const& normal) const noexcept
    {
        switch (m_kind) {
        case Kind::ceramic: return m_ceramic.SampleReflection(velocity, normal);
        case Kind::none: return velocity;
        }
        return velocity;
    }

private:
    enum class Kind { none, ceramic };
    Kind m_kind = Kind::none;
    CeramicInteraction m_ceramic;
};

static_assert(std::is_trivially_copyable_v<MaterialInteractionDevice>);

} // namespace Insert
