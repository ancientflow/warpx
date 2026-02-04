#include "WarpX.H"

#include "BoundaryConditions/PML.H"
#include "Fields.H"
#ifdef WARPX_USE_FFT
#ifdef WARPX_DIM_RZ
#include "FieldSolver/SpectralSolver/SpectralSolverRZ.H"
#else
#include "FieldSolver/SpectralSolver/SpectralSolver.H"
#endif
#endif
#include "Parallelization/GuardCellManager.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/ParticleBoundaryBuffer.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "Utils/WarpXUtil.H"

#include <ablastr/utils/SignalHandling.H>
#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_BLassert.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Vector.H>

class BackgroundCoupledDensity
{
public:
    /*原则上，不需要每一个碰撞步都进行背景粒子的排列，这是
    因为背景粒子被考虑为了每一定步计算一次运动，在这没有运
    动的时间步中，背景粒子仅有权重的变化，因此所有的索引都
    是有效的。在背景粒子将要发生运动时，只要在合适的位置删
    除所有权重为0的粒子即可，同时在下一步重新计算，原则上
    这可以大大加速前期背景粒子极多而电离产物极少的时间段，
    背景粒子的排序最快是nlogn时间复杂度，依然是超线性的*/
    /* 在所有的lev上存储multifab
       在对每个lev(multifab)，存储排序粒子得到的bin，并作为
       参数传入碰撞函数
    */

    std::string m_ground_species;
    amrex::Vector<amrex::MultiFab> m_background_density_fabs;
    amrex::Vector<
        amrex::Vector<amrex::DenseBins<ParticleUtils::ParticleTileDataType>>>
        m_background_bins;
    amrex::Vector<amrex::Vector<amrex::Gpu::DeviceVector<int>>>
        m_n_particle_in_each_cell;

public:
    /**
     * @brief init the vector of background density species
     */
    void backgroundDensityInit ();

    /**
     * @brief update the background density data
     */
    void backgroundDensityUpdate (MultiParticleContainer& mypc,
                                  amrex::ParticleReal elec_weight);

    /**
     * @brief delete the particles with zero weight
     */
    void backgroudnSpeciesClean (MultiParticleContainer& mypc);
};