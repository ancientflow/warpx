#include "BackgroundCoupledDensity.h"
#include "AtomDeposit.h"

/**
 * @brief init the vector of background density species
 */
void BackgroundCoupledDensity::backgroundDensityInit() {
    WarpX& warpx_instance = WarpX::GetInstance();
    MultiParticleContainer& mypc = warpx_instance.GetPartContainer();
    auto& background_species =
        mypc.GetParticleContainerFromName(m_ground_species);
    auto const flvl = background_species.finestLevel();
    m_background_density_fabs.resize(flvl + 1);
    m_background_bins.resize(flvl + 1);
    m_n_particle_in_each_cell.resize(flvl + 1);
    for (int lev = 0; lev <= flvl; lev++) {
        auto *rho = warpx_instance.m_fields.get(FieldType::rho_fp, lev);
        m_background_density_fabs[lev] =
            amrex::MultiFab(rho->boxArray(), rho->DistributionMap(), rho->nComp(),
                     rho->nGrow());
        auto geo = warpx_instance.Geom(lev);
        for (WarpXParIter pti(background_species, lev); pti.isValid(); ++pti) {
            auto& ptile = background_species.ParticlesAt(lev, pti);
            m_background_bins[lev].push_back(
                ParticleUtils::findParticlesInEachCell(geo, pti, ptile));
            long const numbins = (*m_background_bins[lev].rbegin()).numBins();
            m_n_particle_in_each_cell[lev].push_back(
                amrex::Gpu::DeviceVector<int>(numbins, 0));
        }
    }
}

/**
 * @brief update the background density data
 */
void
BackgroundCoupledDensity::backgroundDensityUpdate (
    MultiParticleContainer& mypc, amrex::ParticleReal elec_weight) {
    WarpX& warpx_instance = WarpX::GetInstance();
    auto& background_species =
        mypc.GetParticleContainerFromName(m_ground_species);
    auto const flvl = background_species.finestLevel();
    for (int lev = 0; lev <= flvl; lev++) {
        m_background_density_fabs[lev].setVal(0.0_prt);
        AtomDepostiAPI(background_species, m_background_density_fabs[lev], lev);
        auto geo = warpx_instance.Geom(lev);
        auto binIter = m_background_bins[lev].begin();
        auto npIter = m_n_particle_in_each_cell[lev].begin();
        for (WarpXParIter pti(background_species, lev); pti.isValid(); ++pti) {
            auto& ptile = background_species.ParticlesAt(lev, pti);
            (*binIter) =
                ParticleUtils::findParticlesInEachCell(geo, pti, ptile);

            const int* offsets = (*binIter).offsetsPtr();
            int const* indices = (*binIter).permutationPtr();
            int const np = ptile.numParticles();
            long const numbins = (*binIter).numBins();
            int* p_particle_num = (*npIter).dataPtr();
            auto& soa = ptile.GetStructOfArrays();
            auto& soa_arr = soa.GetRealData();
            amrex::Real const* pw = soa_arr[PIdx::w].dataPtr();

            amrex::ParallelFor(numbins, [=] AMREX_GPU_DEVICE(long ibin) {
                const int offset_start = offsets[ibin],
                          offset_end = offsets[ibin + 1];
                amrex::ParticleReal num_w = 0;
                for (int i = offset_start; i < offset_end; i++) {
                    num_w += pw[indices[i]];
                }
                p_particle_num[ibin] =
                    static_cast<int>((num_w / elec_weight) + 0.1);
            });

            ++binIter;
            ++npIter;
        }
    }
    amrex::Print() << "Update background species: " << m_ground_species << "\n";
}

/**
 * @brief delete the particles with zero weight
 */
void
BackgroundCoupledDensity::backgroundSpeciesClean (
    MultiParticleContainer& mypc) const {
    auto& background_species =
        mypc.GetParticleContainerFromName(m_ground_species);
    auto const flvl = background_species.finestLevel();
    for (int lev = 0; lev <= flvl; lev++) {
        for (WarpXParIter pti(background_species, lev); pti.isValid(); ++pti) {
            auto& ptile = background_species.ParticlesAt(lev, pti);
            long const np = ptile.numParticles();
            auto& soa = ptile.GetStructOfArrays();
            uint64_t* const AMREX_RESTRICT idcpu = soa.GetIdCPUData().data();
            auto& soa_arr = soa.GetRealData();
            amrex::Real const* pw = soa_arr[PIdx::w].dataPtr();
            ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) {
                if (std::abs(pw[ip]) < 10.0_prt) {
                    auto pidw = amrex::ParticleIDWrapper{idcpu[ip]};
                    if (pidw.is_valid()) {
                        pidw.make_invalid();
                    }
                }
            });
        }
    }
}