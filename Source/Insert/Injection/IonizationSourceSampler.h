#pragma once

#include "Insert/Injection/HallPositionSampler.h"

#include <AMReX_REAL.H>
#include <AMReX_RandomEngine.H>

#include <string>
#include <vector>

namespace Insert {

class IonizationSourceSampler final : public HallPositionSampler
{
public:
    explicit IonizationSourceSampler (
        std::string source_directory = "ionization_source_fab");

    [[nodiscard]] EmissionSample
    samplePosition (amrex::RandomEngine const& engine) const override;

    [[nodiscard]] amrex::Real totalRate () const noexcept;
    [[nodiscard]] amrex::ParticleReal electronWeight () const noexcept;

private:
    [[nodiscard]] amrex::Real nodeRate (int iz, int ir) const;
    [[nodiscard]] int nodeIndex (int iz, int ir) const noexcept;

    std::string m_source_directory;
    int m_nr = 0;
    int m_nz = 0;
    amrex::Real m_x_center = 0.0;
    amrex::Real m_y_center = 0.0;
    amrex::Real m_r_min = 0.0;
    amrex::Real m_r_max = 0.0;
    amrex::Real m_z_min = 0.0;
    amrex::Real m_z_max = 0.0;
    amrex::Real m_dr = 0.0;
    amrex::Real m_dz = 0.0;
    amrex::Real m_total_rate = 0.0;
    amrex::ParticleReal m_elec_weight = 0.0;
    std::vector<amrex::Real> m_node_rate;
    std::vector<amrex::Real> m_cell_cdf;
};

} // namespace Insert
