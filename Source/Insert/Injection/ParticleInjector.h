#pragma once

#include "Insert/Injection/ParticleInjectionSource.h"

#include <AMReX_REAL.H>

#include <vector>

class WarpX;

namespace Insert {

class ParticleInjector
{
public:
    static ParticleInjector& GetInstance ();

    void ReadParameters ();
    void InitializePlasma (WarpX& warpx);
    void InjectParticles (WarpX& warpx, amrex::Real dt, int step);

private:
    ParticleInjector () = default;

    bool m_initialized = false;
    std::vector<ParticleInjectionSource> m_initial_sources;
    std::vector<ParticleInjectionSource> m_continuous_sources;
};

} // namespace Insert
