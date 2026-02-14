#include "BackgroundCoupledDensity.h"

/*
 * @brief deposit the atom density by WarpXParticleContainer
 */
void
AtomDepositAPI (WarpXParticleContainer& pc, amrex::MultiFab& rho,
                const int lev) {
    for (WarpXParIter pti(pc, lev); pti.isValid(); ++pti) {
        const amrex::Box& box = pti.validbox();

        auto np = pti.numParticles();
        amrex::Print() << np << " API deposit \n";
        // Extract particle data
        auto& attribs = pti.GetAttribs();
        auto& wp = attribs[PIdx::w];
        pc.DepositCharge(pti, wp, nullptr, &rho, 0, 0, np, 0, lev, lev);
    }
    // 处理边界密度
    // 尝试使用内置函数进行修改，便于进行高阶插值
    // 考虑到更多的层级，必须采用这种内置函数
    auto& warpx_instance = WarpX::GetInstance();
    /*    PEC::ApplyReflectiveBoundarytoRhofield(
            &rho, warpx_instance.field_boundary_lo,
            warpx_instance.field_boundary_hi,
       warpx_instance.particle_boundary_lo, warpx_instance.particle_boundary_hi,
       warpx_instance.Geom(lev), lev, PatchType::fine,
       warpx_instance.refRatio());*/

    amrex::Box domain = warpx_instance.Geom(lev).Domain();
    domain.surroundingNodes();
    for (amrex::MFIter mfi(rho, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const amrex::Box& box = mfi.validbox();
        const auto& rho_arr = rho[mfi].array();
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // x direction
            if (i == domain.smallEnd(0)) {
                rho_arr(i, j, k) *= 2.0_prt;
            }
            if (i == domain.bigEnd(0)) {
                rho_arr(i, j, k) *= 2.0_prt;
            }
            // y direction
            if (j == domain.smallEnd(1)) {
                rho_arr(i, j, k) *= 2.0_prt;
            }
            if (j == domain.bigEnd(1)) {
                rho_arr(i, j, k) *= 2.0_prt;
            }
            // z direction
            if (k == domain.smallEnd(2)) {
                rho_arr(i, j, k) *= 2.0_prt;
            }
            if (k == domain.bigEnd(2)) {
                rho_arr(i, j, k) *= 2.0_prt;
            }
        });
    }
}

/**
 * @brief init the vector of background density species
 */
void
BackgroundCoupledDensity::backgroundDensityInit () {
    WarpX& warpx_instance = WarpX::GetInstance();
    MultiParticleContainer& mypc = warpx_instance.GetPartContainer();
    auto& background_species =
        mypc.GetParticleContainerFromName(m_ground_species);

    // resize the vector
    auto const flvl = background_species.finestLevel();
    m_background_density_fabs.resize(flvl + 1);
    m_background_bins.resize(flvl + 1);
    m_n_particle_in_each_cell.resize(flvl + 1);

    for (int lev = 0; lev <= flvl; lev++) {
#ifndef MCC_DENSITY_MID
        auto* rho = warpx_instance.m_fields.get(FieldType::rho_fp, lev);
        m_background_density_fabs[lev] =
            amrex::MultiFab(rho->boxArray(), rho->DistributionMap(),
                            rho->nComp(), rho->nGrow());
#endif

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
#ifndef MCC_DENSITY_MID
        m_background_density_fabs[lev].setVal(0.0_prt);
        AtomDepositAPI(background_species, m_background_density_fabs[lev], lev);
#else
        m_background_density_fabs[lev] =
            background_species.GetNumberDensity(lev);
#endif
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
            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) {
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