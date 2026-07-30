#pragma once

#include "Insert/Injection/HallInjectionSource.h"

#include <AMReX_REAL.H>

#include <vector>

class WarpX;

namespace Insert {

class HallInjector
{
public:
    static HallInjector& GetInstance ();

    void ReadParameters ();
    void InitializePlasma (WarpX& warpx);
    void InjectParticles (WarpX& warpx, amrex::Real dt, int step);

private:
    HallInjector () = default;

    bool m_initialized = false;
    std::vector<HallInjectionSource> m_initial_sources;
    std::vector<HallInjectionSource> m_continuous_sources;
};

} // namespace Insert
