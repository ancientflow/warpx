#include "BackgroundCoupledDensity.h"

#include <BoundaryConditions/WarpX_PEC.H>

#if defined(MCC_DENSITY_AVERAGE_CALC) || defined(MCC_DENSITY_AVERAGE_USE)
#include <AMReX_ParmParse.H>
#include <AMReX_VisMF.H>
#endif

#ifdef MCC_DENSITY_AVERAGE_CALC
#include <AMReX_Utility.H>

#include <cctype>
#include <iomanip>
#include <sstream>
#endif

#ifdef MCC_DENSITY_AVERAGE_USE
#include <cctype>
#include <sstream>
#endif

#if defined(MCC_DENSITY_AVERAGE_CALC) || defined(MCC_DENSITY_AVERAGE_USE)
namespace {

std::string
SafeName (std::string name)
{
    for (char& c : name) {
        if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_' && c != '-') {
            c = '_';
        }
    }
    return name;
}

std::string
FabOutputParent (std::string const& base_dir, std::string const& species,
                 std::string const& quantity, int lev)
{
    std::ostringstream os;
    os << base_dir << "/" << SafeName(species) << "/" << quantity << "/lev" << lev;
    return os.str();
}

std::string
FabOutputPath (std::string const& base_dir, std::string const& species,
               std::string const& quantity, int lev)
{
    std::ostringstream os;
    os << FabOutputParent(base_dir, species, quantity, lev) << "/" << quantity;
    return os.str();
}

} // namespace
#endif

#ifdef MCC_DENSITY_AVERAGE_CALC
namespace {

std::string
FabStepOutputPath (std::string const& base_dir, std::string const& species,
                   std::string const& quantity, int lev, int step)
{
    std::ostringstream os;
    os << FabOutputParent(base_dir, species, quantity, lev) << "/step" << std::setw(8)
       << std::setfill('0') << step;
    return os.str();
}

void
CreateDirectoryTree (std::string const& dir)
{
    if (dir.empty() || dir == ".") {
        return;
    }

    amrex::Vector<std::string> paths;
    std::string current;
    auto pos = std::string::size_type{0};
    if (dir[0] == '/') {
        current = "/";
        pos = 1;
    }
    while (pos < dir.size()) {
        auto const next = dir.find('/', pos);
        auto const part = dir.substr(pos, next - pos);
        if (!part.empty()) {
            if (!current.empty() && current != "/") {
                current += "/";
            }
            current += part;
            paths.push_back(current);
        }
        if (next == std::string::npos) {
            break;
        }
        pos = next + 1;
    }

    if (amrex::ParallelDescriptor::IOProcessor()) {
        constexpr int permission_flag_rwxrxrx = 0755;
        for (auto const& path : paths) {
            if (!amrex::UtilCreateDirectory(path, permission_flag_rwxrxrx)) {
                amrex::CreateDirectoryFailed(path);
            }
        }
    }
    amrex::ParallelDescriptor::Barrier();
}

std::string
ParentPath (std::string const& path)
{
    auto const slash = path.find_last_of('/');
    if (slash == std::string::npos) {
        return ".";
    }
    if (slash == 0) {
        return "/";
    }
    return path.substr(0, slash);
}

void
WriteSingleFabMultiFab (amrex::MultiFab const& mf, std::string const& path,
                        std::string const& species, int lev)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        mf.boxArray().size() == 1,
        "Background density FAB output requires exactly one FAB for species " +
            species + " at level " + std::to_string(lev) + ".");
    CreateDirectoryTree(ParentPath(path));
    amrex::VisMF::Write(mf, path);
}

} // namespace
#endif

/**
 * @brief deposit the atom density by WarpXParticleContainer
 */
void
AtomDepositAPI (WarpXParticleContainer& pc, amrex::MultiFab& rho,
                const int lev) {
    for (WarpXParIter pti(pc, lev); pti.isValid(); ++pti) {
        const amrex::Box& box = pti.validbox();

        auto np = pti.numParticles();
        amrex::AllPrint() << "rank " << amrex::ParallelDescriptor::MyProc()
                          << ": deposit " << np << " particles" << "\n";
        // Extract particle data
        auto& attribs = pti.GetAttribs();
        auto& wp = attribs[PIdx::w];
        pc.DepositCharge(pti, wp, nullptr, &rho, 0, 0, np, 0, lev, lev);
    }
    // 处理边界密度
    // 尝试使用内置函数进行修改，便于进行高阶插值
    // 考虑到更多的层级，必须采用这种内置函数
    //考虑到PEC中的修改，这会导致边界密度被大大低估，但是一般来说这不会明细影响放电
    auto& warpx_instance = WarpX::GetInstance();
    PEC::ApplyReflectiveBoundarytoRhofield(
        &rho, WarpX::field_boundary_lo, WarpX::field_boundary_hi,
        WarpX::particle_boundary_lo, WarpX::particle_boundary_hi,
        warpx_instance.Geom(lev), lev, PatchType::fine,
        warpx_instance.refRatio());

    /*
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
        });*/
}

/**
 * @brief init the vector of background density species
 * @attention if number of background species particles at begin is zero,
 * this function will give wrong result, because findParticlesInEachCell
 * will return vacuum container
 */
void
BackgroundCoupledDensity::backgroundDensityInit () {
    WarpX& warpx_instance = WarpX::GetInstance();

#ifndef MCC_DENSITY_AVERAGE_USE
    MultiParticleContainer& mypc = warpx_instance.GetPartContainer();
    auto& background_species =
        mypc.GetParticleContainerFromName(m_ground_species);
    // resize the vector
    auto const flvl = background_species.finestLevel();
#else
    // In USE mode no neutral particle container is required;
    // the averaged density is loaded directly from file.
    int const flvl = 0;
#endif
    m_background_density_fabs.resize(flvl + 1);
#ifdef MCC_DENSITY_AVERAGE_CALC
    m_background_density_sum_fabs.resize(flvl + 1);
#endif
    m_background_bins.resize(flvl + 1);
    m_n_particle_in_each_cell.resize(flvl + 1);

    amrex::Print() << "Background density initialized with size: \n";

    for (int lev = 0; lev <= flvl; lev++) {
        auto* rho =
            warpx_instance.m_fields.get(warpx::fields::FieldType::rho_fp, lev);
        m_background_density_fabs[lev] =
            amrex::MultiFab(rho->boxArray(), rho->DistributionMap(),
                            rho->nComp(), rho->nGrow());
#ifdef MCC_DENSITY_AVERAGE_CALC
        m_background_density_sum_fabs[lev] =
            amrex::MultiFab(rho->boxArray(), rho->DistributionMap(),
                            rho->nComp(), rho->nGrow());
        m_background_density_sum_fabs[lev].setVal(0.0);
#endif

#ifndef MCC_DENSITY_AVERAGE_USE
        auto geo = warpx_instance.Geom(lev);
        const size_t box_num = warpx_instance.boxArray(lev).size();

        auto& background_bin = m_background_bins[lev];
        auto& background_np = m_n_particle_in_each_cell[lev];

        background_bin.resize(box_num);
        background_np.resize(box_num);

        for (WarpXParIter pti(background_species, lev); pti.isValid(); ++pti) {
            const int box_index = pti.index();
            auto& ptile = background_species.ParticlesAt(lev, pti);
            background_bin[box_index] =
                ParticleUtils::findParticlesInEachCell(geo, pti, ptile);
            long const numbins = background_bin[box_index].numBins();
            background_np[box_index] =
                amrex::Gpu::DeviceVector<int>(numbins, 0);
        }
#endif
    }
#ifdef MCC_DENSITY_AVERAGE_CALC
    amrex::ParmParse const pp_bg("background_density");
    pp_bg.query("average_steps_per_period", m_average_steps_per_period);
    pp_bg.query("average_periods", m_average_periods);
    pp_bg.query("raw_output_interval", m_raw_output_interval);
    pp_bg.query("output_dir", m_output_dir);
    pp_bg.query("output_fab", m_output_fab);

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_average_steps_per_period > 0,
        "background_density.average_steps_per_period must be positive.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_average_periods > 0,
        "background_density.average_periods must be positive.");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_raw_output_interval >= 0,
        "background_density.raw_output_interval must be non-negative.");

    m_average_window_steps = m_average_steps_per_period * m_average_periods;
#endif
#ifdef MCC_DENSITY_AVERAGE_USE
    amrex::ParmParse const pp_bg("background_density");
    pp_bg.query("input_fab", m_input_fab);
    if (m_input_fab.empty()) {
        m_input_fab = FabOutputPath("background_density_fab", m_ground_species, "average", 0);
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_background_density_fabs[0].boxArray().size() == 1,
        "MCC_DENSITY_AVERAGE_USE requires the current background density "
        "MultiFab to contain exactly one FAB.");
    amrex::Print() << "Reading averaged background density FAB from "
                   << m_input_fab << "\n";
    amrex::VisMF::Read(m_background_density_fabs[0], m_input_fab);
#endif
}

/**
 * @brief update the background density data
 */
void
BackgroundCoupledDensity::backgroundDensityUpdate (
    MultiParticleContainer& mypc, amrex::ParticleReal elec_weight, int step) {
    using namespace amrex::literals;
#ifndef MCC_DENSITY_AVERAGE_CALC
    static_cast<void>(step);
#endif

    WarpX& warpx_instance = WarpX::GetInstance();
    auto& background_species =
        mypc.GetParticleContainerFromName(m_ground_species);
    auto const flvl = background_species.finestLevel();
    amrex::AllPrint() << "rank " << amrex::ParallelDescriptor::MyProc()
                      << ":Start updating background density\n";
    for (int lev = 0; lev <= flvl; lev++) {
        m_background_density_fabs[lev].setVal(0.0_prt);
        AtomDepositAPI(background_species, m_background_density_fabs[lev], lev);
#ifdef MCC_DENSITY_AVERAGE_CALC
        if (m_raw_output_interval > 0 && step % m_raw_output_interval == 0) {
            WriteSingleFabMultiFab(
                m_background_density_fabs[lev],
                FabStepOutputPath(m_output_dir, m_ground_species, "raw", lev, step),
                m_ground_species, lev);
        }
        if (m_average_sample_count < m_average_window_steps) {
            amrex::MultiFab::Add(
                m_background_density_sum_fabs[lev], m_background_density_fabs[lev],
                0, 0, m_background_density_fabs[lev].nComp(),
                m_background_density_fabs[lev].nGrowVect());
        }
#endif
        auto geo = warpx_instance.Geom(lev);
        auto& background_bin = m_background_bins[lev];
        auto& background_np = m_n_particle_in_each_cell[lev];

        for (WarpXParIter pti(background_species, lev); pti.isValid(); ++pti) {
            const int box_index = pti.index();
            amrex::AllPrint()
                << "rank " << amrex::ParallelDescriptor::MyProc()
                << ": Updating lev " << lev << ": box " << box_index << "\n";

            auto& ptile = background_species.ParticlesAt(lev, pti);
            background_bin[box_index] =
                ParticleUtils::findParticlesInEachCell(geo, pti, ptile);
            int const np = ptile.numParticles();
            const int* offsets = background_bin[box_index].offsetsPtr();
            int const* indices = background_bin[box_index].permutationPtr();
            long const numbins = background_bin[box_index].numBins();
            int* p_particle_num = background_np[box_index].dataPtr();
            /*
            const int* offsets = (*binIter).offsetsPtr();
            int const* indices = (*binIter).permutationPtr();
            long const numbins = (*binIter).numBins();
            int* p_particle_num = (*npIter).dataPtr();*/
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
        }
    }
#ifdef MCC_DENSITY_AVERAGE_CALC
    if (m_average_sample_count < m_average_window_steps) {
        ++m_average_sample_count;
    }
#endif
    amrex::Print() << "rank " << amrex::ParallelDescriptor::MyProc()
                   << ": Updated background species: " << m_ground_species
                   << " end\n";
}

void
BackgroundCoupledDensity::backgroundDensityFinalize () {
#ifdef MCC_DENSITY_AVERAGE_CALC
    using namespace amrex::literals;

    if (m_average_sample_count <= 0) {
        amrex::Print() << "No background density samples collected for species "
                       << m_ground_species << "; skip averaged FAB output.\n";
        return;
    }

    auto const sample_count = static_cast<amrex::Real>(m_average_sample_count);
    amrex::Vector<amrex::MultiFab> average_fabs(m_background_density_sum_fabs.size());
    for (int lev = 0; lev < static_cast<int>(m_background_density_sum_fabs.size()); ++lev) {
        average_fabs[lev] = amrex::MultiFab(
            m_background_density_sum_fabs[lev].boxArray(),
            m_background_density_sum_fabs[lev].DistributionMap(),
            m_background_density_sum_fabs[lev].nComp(),
            m_background_density_sum_fabs[lev].nGrowVect());
        amrex::MultiFab::Copy(
            average_fabs[lev], m_background_density_sum_fabs[lev], 0, 0,
            m_background_density_sum_fabs[lev].nComp(),
            m_background_density_sum_fabs[lev].nGrowVect());
        average_fabs[lev].mult(1.0_rt / sample_count);

        std::string const output_path =
            (lev == 0 && !m_output_fab.empty())
                ? m_output_fab
                : FabOutputPath(m_output_dir, m_ground_species, "average", lev);
        WriteSingleFabMultiFab(
            average_fabs[lev], output_path, m_ground_species, lev);
    }

    amrex::Print() << "Wrote averaged background density FAB for species "
                   << m_ground_species << " from " << m_average_sample_count
                   << " samples.\n";
#endif
}

/**
 * @brief delete the particles with zero weight
 */
void
BackgroundCoupledDensity::backgroundSpeciesClean (
    MultiParticleContainer& mypc) const {
    using namespace amrex::literals;

    auto& background_species =
        mypc.GetParticleContainerFromName(m_ground_species);
    auto const flvl = background_species.finestLevel();
    for (int lev = 0; lev <= flvl; lev++) {
        for (WarpXParIter pti(background_species, lev); pti.isValid(); ++pti) {
            amrex::AllPrint()
                << "rank " << amrex::ParallelDescriptor::MyProc()
                << ": clean lev " << lev << ": box " << pti.index() << "\n";
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
