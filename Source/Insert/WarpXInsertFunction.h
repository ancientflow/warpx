#include "WarpX.H"

#include "BoundaryConditions/PML.H"
#include "Diagnostics/MultiDiagnostics.H"
#include "Diagnostics/ReducedDiags/MultiReducedDiags.H"
#include "EmbeddedBoundary/Enabled.H"
#include "Fields.H"
#include "FieldSolver/FiniteDifferenceSolver/HybridPICModel/HybridPICModel.H"
#ifdef WARPX_USE_FFT
#   ifdef WARPX_DIM_RZ
#       include "FieldSolver/SpectralSolver/SpectralSolverRZ.H"
#   else
#       include "FieldSolver/SpectralSolver/SpectralSolver.H"
#   endif
#endif
#include "Parallelization/GuardCellManager.H"
#include "Particles/MultiParticleContainer.H"
#include "Fluids/MultiFluidContainer.H"
#include "Fluids/WarpXFluidContainer.H"
#include "Particles/ParticleBoundaryBuffer.H"
#include "Python/callbacks.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXUtil.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXProfilerWrapper.H"

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

#include <algorithm>
#include <array>
#include <memory>
#include <ostream>
#include <vector>
#include "WarpXFunctionConfig.h"
#include "WarpXSimulationConfig.h"
#include "BackgroundCoupledDensity.h"

using namespace amrex;

#ifdef BENCHMARK_2D
// 电子数量计算
int
ElectronsCollection (WarpX& warpx_instance) {
    auto& boundary_buffer = warpx_instance.GetParticleBoundaryBuffer();
    int e_xlo =
            boundary_buffer.getNumParticlesInContainer("electrons", 0, false),
        // e_xhi =
        // boundary_buffer.getNumParticlesInContainer("electrons",1,false),
        xe_xlo =
            boundary_buffer.getNumParticlesInContainer("xe_ions", 0, false);
    // xe_xhi = boundary_buffer.getNumParticlesInContainer("xe_ions",1,false);
    return std::max(e_xlo - xe_xlo, 0);
}

// 粒子注入-模拟电离-阴极发射
void
ParticleInjection (WarpX& warpx_instance) {
    static const int nx = 512, ny = 256, nppc = 75;
    static const double lx = 0.025, ly = 0.0128, NPlasma = 5e16, s0 = 5.23e23,
                        const_dt = 5e-12;

    static const int once_injection = 1000; // 一次注入电子-离子对数量
    static double rest_macro_particles = 0; // 剩余应当注入粒子
    static const double one_step_real_particles =
        s0 * const_dt * 0.0128 * 2 / 3.1415926 *
        0.0075; // 每时步应当注入真实粒子
    static const Real global_weight =
        NPlasma / nx / ny * lx * ly / nppc; // 统一权重
    static const double one_step_macro_particles =
        one_step_real_particles / global_weight; // 每时步应当注入宏粒子

    int this_step_pair = 0,
        this_step_cathode_electron = ElectronsCollection(warpx_instance);
    rest_macro_particles += one_step_macro_particles;
    if (rest_macro_particles > once_injection) {
        rest_macro_particles -= once_injection;
        this_step_pair = once_injection;
    }

    auto& mypc = warpx_instance.GetPartContainer();
    auto& electons_pc =
        mypc.GetParticleContainer(mypc.getSpeciesID("electrons"));
    auto& xe_pc = mypc.GetParticleContainer(mypc.getSpeciesID("xe_ions"));

    const amrex::Vector<amrex::Vector<int>> nattr;
    amrex::RandomEngine uniform_engine(getRandState()),
        normal_engine(getRandState());

    int electron_size = this_step_pair + this_step_cathode_electron;
    std::cout << electron_size << " " << this_step_pair << " " << std::endl;
    amrex::Vector<amrex::Real> pxe(electron_size), pye(electron_size, 0),
        pze(electron_size), vxe(electron_size), vye(electron_size),
        vze(electron_size), we(electron_size, global_weight),
        pxxe(this_step_pair), pyxe(this_step_pair, 0), pzxe(this_step_pair),
        vxxe(this_step_pair), vyxe(this_step_pair), vzxe(this_step_pair),
        wxe(this_step_pair, global_weight);

    // 电子-离子对注入
    if (this_step_pair > 0) {
        // std::cout<<istart<<" "<<iend<<std::endl;
        for (size_t i = 0; i < this_step_pair; i++) {
            double r1 = amrex::Random(uniform_engine),
                   r2 = amrex::Random(uniform_engine);
            pxe[i] = 0.00625 + asin(2 * r1 - 1) / 3.14159 * 0.0075;
            pze[i] = 0.0128 * r2;
            pxxe[i] = pxe[i];
            pzxe[i] = pze[i];

            vxe[i] = amrex::RandomNormal(0, 1326232, normal_engine);
            vye[i] = amrex::RandomNormal(0, 1326232, normal_engine);
            vze[i] = amrex::RandomNormal(0, 1326232, normal_engine);
            vxxe[i] = amrex::RandomNormal(0, 606.34, normal_engine);
            vyxe[i] = amrex::RandomNormal(0, 606.34, normal_engine);
            vzxe[i] = amrex::RandomNormal(0, 606.34, normal_engine);
        }
    }

    // 阴极发射
    if (this_step_cathode_electron > 0) {
        for (size_t i = this_step_pair; i < electron_size; i++) {
            pxe[i] = 0.024;
            pze[i] = 0.0128 * amrex::Random(uniform_engine);
            vxe[i] = amrex::RandomNormal(0, 1326232, normal_engine);
            vye[i] = amrex::RandomNormal(0, 1326232, normal_engine);
            vze[i] = amrex::RandomNormal(0, 1326232, normal_engine);
        }
    }
    if (this_step_pair > 0) {
        std::cout << "add pairs number: " << this_step_pair << std::endl;
        xe_pc.AddNParticles(0, this_step_pair, pxxe, pyxe, pzxe, vxxe, vyxe,
                            vzxe, 1, {wxe}, 0, nattr, false);
    }

    // 阴极发射部分
    // 清除缓存
    if (electron_size > 0) {
        electons_pc.AddNParticles(0, electron_size, pxe, pye, pze, vxe, vye,
                                  vze, 1, {we}, 0, nattr, false);
        if (this_step_cathode_electron > 0) {
            auto& boundary_buffer = warpx_instance.GetParticleBoundaryBuffer();
            boundary_buffer.clearParticles();
        }
    }
}

// 电压修正
void
VoltageAdjustment (WarpX& warpx_instance) {
    auto phi_field = warpx_instance.m_fields.get(FieldType::phi_fp, 0);
    // fields.get(FieldType::phi_fp,0);
    Real phisum = 0;

    static MultiFab mf(phi_field->boxArray(), phi_field->DistributionMap(),
                       phi_field->nComp(), phi_field->nGrow());
    static bool DotInit = false;
    if (!DotInit) {
        mf.setVal(0);
        for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& box = mfi.tilebox();
            Array4<Real> const& arr = mf.array(mfi);
            ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j) {
                if (i == 491)
                    arr(i, j, 0) = 0.48;
                else if (i == 492)
                    arr(i, j, 0) = 0.52;
            });
        }
        DotInit = true;
    }
    phisum = MultiFab::Dot(mf, 0, *phi_field, 0, 1, 0);

    phisum /= 256;

    std::cout << "voltage adjustment: " << phisum << std::endl;

    for (MFIter mfi(*phi_field, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& box = mfi.tilebox();

        Array4<Real> const& phi = phi_field->array(mfi);
        ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j) {
            phi(i, j, 0) -= phisum * i / 491.52;
        });
    }
}
#endif


#ifdef HALL3D
// 3D阴极发射
void CathodeInjection3D ()
{
    static const double L = 0.05; // 未缩比边长
    static double Ic, dt, l_factor, elec_weight;

    static bool ifinit = false;

    if(!ifinit)
    {
        amrex::ParmParse pp_mc("my_constants");
        pp_mc.get("dt", dt);
        pp_mc.get("Ic", Ic);
        pp_mc.get("l_factor", l_factor);
        pp_mc.getWithParser("elec_weight", elec_weight);
    }

    const double num_of_real_electrons =
                     Ic * dt / l_factor / l_factor / 1.60217e-19,
                 num_of_marco_particle = num_of_real_electrons / elec_weight;

    WarpX& warpx_instance = WarpX::GetInstance();

    static double electron_rest = 0;
    const int one_step_injection = int(num_of_marco_particle) + 1;

    electron_rest += num_of_marco_particle;

    static const double r1 = 0.016 / l_factor, r2 = 0.02 / l_factor, dr = r2 - r1,
                  sumr = r1 + r2, multr = r1 * r2, z1 = 0.042 / l_factor,
                  z2 = 0.046 / l_factor, dz = z2 - z1, Pi2 = 3.1415926 * 2,
                  sigma = 592982, half_l = L / 2 / l_factor;

    double r, theta;
    static Vector<Vector<int>> nattr;
    if (electron_rest > one_step_injection) {

        electron_rest -= one_step_injection;
        //std::cout << "electron inject: " << one_step_injection << std::endl;

        Vector<ParticleReal> pz(one_step_injection), px(one_step_injection),
            py(one_step_injection), vx(one_step_injection),
            vy(one_step_injection), vz(one_step_injection),
            pw(one_step_injection, elec_weight);

        RandomEngine uniform_engine(getRandState()),
            normal_engine(getRandState());

        for (int i = 0; i < one_step_injection; i++) {
            // 柱坐标空间均匀分布
            r = Random(uniform_engine) * dr + r1;
            theta = Random(uniform_engine) * Pi2;
            pz[i] = Random(uniform_engine) * dz + z1;

            // 转化到直角坐标空间均匀分布
            r = sqrt(sumr * r - multr);
            px[i] = r * cos(theta) + half_l;
            py[i] = r * sin(theta) + half_l;

            // 速度
            vx[i] = RandomNormal(0, sigma, normal_engine);
            vy[i] = RandomNormal(0, sigma, normal_engine);
            vz[i] = RandomNormal(0, sigma, normal_engine);
        }

        auto& mypc = warpx_instance.GetPartContainer();
        auto& e_pc = mypc.GetParticleContainer(mypc.getSpeciesID("electrons"));
        e_pc.AddNParticles(0, one_step_injection, px, py, pz, vx, vy, vz, 1,
                           {pw}, 0, nattr, 0);
    }
}

// 氙气注入
void XeInjection () 
{
    static const double L = 0.05, // 未缩比边长
        NA = 6.03e23,              // 阿伏伽德罗常数
        vz0 = 250;                 // 气流整体速度
    static bool ifinit = false, ifhole = false;
    static int hole_num = 48;
    static double dt, l_factor, atom_weight, m_dot;
    static amrex::Vector<double> hole_x, hole_y;

    if (!ifinit) {
        amrex::ParmParse pp_mc("my_constants");
        pp_mc.get("dt", dt);
        pp_mc.get("l_factor", l_factor);
        pp_mc.getWithParser("xe_weight", atom_weight);
        pp_mc.get("m_dot", m_dot);
        pp_mc.query("ifhole", ifhole);
        ifinit = true;
    }
    static const double V_per_sec = 0.06 / 60. / 1e6,
                 n_per_sec = V_per_sec * 101325 / 8.314 / 273.15 * NA,
                 n_per_step = n_per_sec * dt, mxe = 2.179e-25,
                 n_marco_per_step = n_per_step / atom_weight, // 每时步应当注入宏粒子
        Pi2 = 3.1415926 * 2;

    const double n_marco_per_step_m =
        m_dot * dt / l_factor / l_factor / mxe / atom_weight;

    const double one_times_inject_particle =
        static_cast<int>(n_marco_per_step_m * 100 + 1); // 一次注入粒子数 100时步一次
    static double xe_rest = 0;

    const double rmid = (0.021 + 0.031) / 4 / l_factor,
                 rwidth = 0.001 / l_factor, r1 = rmid - rwidth,
                 r2 = rmid + rwidth, sqdr = r2 * r2 - r1 * r1, sqr1 = r1 * r1,
                 kb = 1.38064852e-23, 
                 sigmax = std::sqrt(kb * 150 / mxe),
                 sigmay = std::sqrt(kb * 150 / mxe),
                 sigmaz = std::sqrt(kb * 150 / mxe), half_L = L / 2 / l_factor,
                 sqdr_hole = rwidth * rwidth;
    double r, theta;
    static const Vector<Vector<int>> nattr;

    static bool if_hole_init = false;
    if(!if_hole_init)
    {
        double per_angle = Pi2 / hole_num, angle = 0;
        for (int i = 0; i < hole_num; i++) {
            hole_x.push_back(rmid * cos(angle));
            hole_y.push_back(rmid * sin(angle));
            angle += per_angle;
        }
        if_hole_init = true;
    }

    xe_rest += n_marco_per_step_m;

    if (xe_rest > one_times_inject_particle) {
        WarpX& warpx_instance = WarpX::GetInstance();

        xe_rest -= one_times_inject_particle;

        Vector<ParticleReal> pz(one_times_inject_particle, 0),
            px(one_times_inject_particle), py(one_times_inject_particle),
            vx(one_times_inject_particle), vy(one_times_inject_particle),
            vz(one_times_inject_particle),
            pw(one_times_inject_particle, atom_weight);

        RandomEngine uniform_engine(getRandState()),
            normal_engine(getRandState());
        int hole_start = rand() % hole_num;
        for (int i = 0; i < one_times_inject_particle; i++) {

            // 柱坐标空间均匀分布
            // 转化到直角坐标空间均匀分布
            if (!ifhole) {
                //非小孔进气
                theta = Random(uniform_engine) * Pi2;
                r = sqrt(sqdr * Random(uniform_engine) + sqr1);
                px[i] = r * cos(theta) + half_L;
                py[i] = r * sin(theta) + half_L;
            }
            else
            {
                //小孔进气
                theta = Random(uniform_engine) * Pi2;
                r = rwidth * sqrt(Random(uniform_engine)); // 小孔进气是r1 = 0, r2 = rwidth, 化简后得到该式
                px[i] = r * cos(theta) + half_L +
                        hole_x[(i + hole_start) % hole_num];
                py[i] = r * sin(theta) + half_L +
                        hole_y[(i + hole_start) % hole_num];
            }

            // 速度
            vx[i] = RandomNormal(0, sigmax, normal_engine);
            vy[i] = RandomNormal(0, sigmay, normal_engine);
            do
            {
                vz[i] = RandomNormal(0, sigmaz, normal_engine) + vz0;
            } while (vz[i] < 0);
        }

        auto& mypc = warpx_instance.GetPartContainer();
        auto& xe_pc = mypc.GetParticleContainer(mypc.getSpeciesID("xe_netural"));
        xe_pc.AddNParticles(0, one_times_inject_particle, px, py, pz, vx, vy, vz,
                           1, {pw}, 0, nattr, 0);
        amrex::Print() << "Injection Xe Atom\n";
    }
}

// 阳极电压修正
void AnodeVoltage () 
{
    WarpX& warpx_instance = WarpX::GetInstance();
    amrex::ParmParse pp_mc("my_constants");
    double ncell, voltage, L, l_factor;
    pp_mc.getWithParser("n_cell", ncell);
    pp_mc.get("L", L);
    pp_mc.get("l_factor", l_factor);
    pp_mc.query("voltage", voltage);
    double halfncell = ncell / 2;
    double sim_L = L / l_factor,gap = sim_L / ncell,
           index_sq_min =
               (0.021 / 2 / l_factor / gap) * (0.021 / 2 / l_factor / gap),
           index_sq_max =
               (0.031 / 2 / l_factor / gap) * (0.031 / 2 / l_factor / gap);

    auto phi_field = warpx_instance.m_fields.get(FieldType::phi_fp, 0);
    amrex::Box domain = warpx_instance.Geom(0).Domain();
    domain.surroundingNodes();
    for (MFIter mfi(*phi_field, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& box = mfi.tilebox();

        Array4<Real> const& phi = phi_field->array(mfi);
        if (!domain.strictly_contains(box))
        {
            ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (k == domain.smallEnd(2)) {
                    double i_ex = i - halfncell, j_ex = j - halfncell;
                    if (i_ex * i_ex + j_ex * j_ex >= index_sq_min &&
                        i_ex * i_ex + j_ex * j_ex <= index_sq_max)
                        phi(i, j, k) = voltage;
                    else
                        phi(i, j, k) = 0;
                }
            });
        }
    }
    phi_field->FillBoundary(warpx_instance.Geom(0).periodicity());
}

//叠加阳极电势
void GetPhiFromFile()
{
    static bool ifinit = false;
    WarpX& warpx_instance = WarpX::GetInstance();
    auto phi = warpx_instance.m_fields.get(FieldType::phi_fp, 0);
    static MultiFab phi_ext(phi->boxArray(), phi->DistributionMap(),
                            phi->nComp(), phi->nGrow());
    if(!ifinit)
    {
        ifinit = true;
        std::fstream filein("phiField3D.dat", std::ios::in);
        static const int ncell = 128, ndata = ncell + 1, ngap = 2;

        std::vector<Real> data(ndata * ndata * ndata);
        int ix, iy, iz;
        Real iphi;
        for (int i = 0; i < ndata; i++) {
            for (int j = 0; j < ndata; j++) {
                for (int k = 0; k < ndata; k++) {
                    filein >> ix >> iy >> iz >> iphi;
                    data[ix * ndata * ndata + iy * ndata + iz] = iphi;
                }
            }
        }

        amrex::Gpu::DeviceVector<Real> device_data(data.size());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice, data.begin(), data.end(),
                         device_data.begin());

        Real* pdevice = device_data.dataPtr();
        for (MFIter mfi(phi_ext, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            auto const& phi_arr = phi_ext.array(mfi);
            const Box& box = mfi.validbox();

            //auto lb = lbound(box);
            //auto hi = ubound(box);

            //amrex::Print() << lb.x << " " << hi.x << " " << phi_ext.nGrow() << "\n";

            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                phi_arr(i, j, k) = pdevice[i * ndata * ndata * ngap +
                                           j * ndata * ngap + k * ngap];
            });
        }
        phi_ext.FillBoundary(warpx_instance.Geom(0).periodicity());
    }
    MultiFab::Add(*phi, phi_ext, 0, 0, 1, phi->nGrow());
}

#endif

#ifdef NUMP
//粒子数监测
void ParticleNumber()
{
    WarpX& warpx_instance = WarpX::GetInstance();
    auto& mypc = warpx_instance.GetPartContainer();
    auto& xen_pc = mypc.GetParticleContainer(mypc.getSpeciesID("xe_netural"));
#if defined(HALL3D)
    auto& e_pc = mypc.GetParticleContainer(mypc.getSpeciesID("electrons"));
    auto& xe_ions_pc = mypc.GetParticleContainer(mypc.getSpeciesID("xe_ions"));
    //auto& xens_pc = mypc.GetParticleContainer(mypc.getSpeciesID("xe_netural_s"));
#endif

    amrex::Print() << "xe atom number: " << xen_pc.TotalNumberOfParticles()
                   << "\n";
#if defined(HALL3D)
    amrex::Print() << "electron number: " << e_pc.TotalNumberOfParticles()
                   << "\n";
    amrex::Print() << "xe ion number: " << xe_ions_pc.TotalNumberOfParticles()
                   << "\n";
    //amrex::Print() << "xes number: " << xens_pc.TotalNumberOfParticles()
    //               << "\n";
#endif
}

//阳极电流及全场电子计算
void AnodeCurrentCalc () 
{
    static int times = 0;
    static const int gap = 10;
    times++;

    static bool ifinit = false;
    if (times == gap) {
        // 打开输出文件
        std::fstream fileout;
        if (!ifinit) {
            fileout.open("anode_current.dat", std::ios::out);
            fileout
                << "time\txmin_electron\txmax_electron\t"
                   "ymin_electron\tymax_electron\t"
                   "zmin_electron\tzmax_electron\tzmin_ion\tanode_electron\t"
                   "anode_electron_cut\n";
            ifinit = true;
        } else {
            fileout.open("anode_current.dat", std::ios::app);
        }
        times = 0;

        // 获取粒子容器信息
        WarpX& warpx_instance = WarpX::GetInstance();
        auto& mybpc = warpx_instance.GetParticleBoundaryBuffer();

        auto pepc_xmin = mybpc.getParticleBufferPointer("electrons", 0);
        auto pepc_xmax = mybpc.getParticleBufferPointer("electrons", 1);
        auto pepc_ymin = mybpc.getParticleBufferPointer("electrons", 2);
        auto pepc_ymax = mybpc.getParticleBufferPointer("electrons", 3);
        auto pepc_zmin = mybpc.getParticleBufferPointer("electrons", 4);
        auto pepc_zmax = mybpc.getParticleBufferPointer("electrons", 5);
        auto pxepec_zmin = mybpc.getParticleBufferPointer("xe_ions", 4);

        Vector<PinnedMemoryParticleContainer*> buffers = {
            pepc_xmin, pepc_xmax, pepc_ymin,  pepc_ymax,
            pepc_zmin, pepc_zmax, pxepec_zmin};

        int data_size = buffers.size() + 2;
        amrex::Gpu::DeviceVector<ParticleReal> device_charge(data_size, 0);
        amrex::Vector<ParticleReal> host_charge(data_size, 0);

        ParticleReal* device_ptr = device_charge.dataPtr();

        ParticleReal l_factor, half_L, L, rn1, rn2;
        amrex::ParmParse pp_mc("my_constants");
        pp_mc.get("l_factor", l_factor);
        pp_mc.get("L", L);

        half_L = L / l_factor / 2;
        rn1 = 0.021 / 2 / l_factor;
        rn2 = 0.031 / 2 / l_factor;

        // 统计zmin电子权重
        for (auto pti = PinnedMemoryParticleContainer::ParIterType(*buffers[4], 0);
             pti.isValid(); ++pti) {
            auto& arr = pti.GetStructOfArrays().GetRealData();
            auto px = arr[PIdx::x].dataPtr();
            auto py = arr[PIdx::y].dataPtr();
            auto pw = arr[PIdx::w].dataPtr();
            auto pz = arr[PIdx::z].dataPtr();
            auto pvx = arr[PIdx::ux].dataPtr();
            auto pvy = arr[PIdx::uy].dataPtr();
            auto pvz = arr[PIdx::uz].dataPtr();
            int np = pti.numParticles();

            ParallelFor(np, [=] AMREX_GPU_DEVICE(long ip) {
                ParticleReal x = px[ip], y = py[ip], w = pw[ip], z = pz[ip],
                             vx = pvx[ip], vy = pvy[ip], vz = pvz[ip];
                ParticleReal dt = z / vz;
                x -= half_L;
                y -= half_L;
                ParticleReal r = sqrt(x * x + y * y);
                if (r >= rn1 && r <= rn2) {
                    amrex::Gpu::Atomic::Add(&device_ptr[data_size - 2], w);
                }
                x -= dt * vx;
                y -= dt * vy;
                r = sqrt(x * x + y * y);
                if (r >= rn1 && r <= rn2) {
                    amrex::Gpu::Atomic::Add(&device_ptr[data_size - 1], w);
                }
                amrex::Gpu::Atomic::Add(&device_ptr[4], w);
            });
        }

        //统计zmin 氙离子权重
        for (auto pti =
                 PinnedMemoryParticleContainer::ParIterType(*buffers[6], 0);
             pti.isValid(); ++pti) {
            auto& arr = pti.GetStructOfArrays().GetRealData();
            auto px = arr[PIdx::x].dataPtr();
            auto py = arr[PIdx::y].dataPtr();
            auto pw = arr[PIdx::w].dataPtr();
            auto pz = arr[PIdx::z].dataPtr();
            auto pvx = arr[PIdx::ux].dataPtr();
            auto pvy = arr[PIdx::uy].dataPtr();
            auto pvz = arr[PIdx::uz].dataPtr();
            int np = pti.numParticles();

            ParallelFor(np, [=] AMREX_GPU_DEVICE(long ip) {
                ParticleReal x = px[ip], y = py[ip], w = pw[ip], z = pz[ip],
                             vx = pvx[ip], vy = pvy[ip], vz = pvz[ip];
                ParticleReal dt = z / vz;
                x -= half_L;
                y -= half_L;
                x -= dt * vx;
                y -= dt * vy;
                ParticleReal r = sqrt(x * x + y * y);
                if (r >= rn1 && r <= rn2) {
                    amrex::Gpu::Atomic::Add(&device_ptr[6], w);
                }
            });
        }

        // 统计其余权重
        Vector<int> buffer_index = {0, 1, 2, 3, 5};
        for (int index : buffer_index) {
            for (auto pti =
                     PinnedMemoryParticleContainer::ParIterType(*buffers[index], 0);
                 pti.isValid(); ++pti) {
                auto& arr = pti.GetStructOfArrays().GetRealData();
                auto pw = arr[PIdx::w].dataPtr();
                int np = pti.numParticles();
                ParallelFor(np, [=] AMREX_GPU_DEVICE(long ip) {
                    ParticleReal w = pw[ip];
                    amrex::Gpu::Atomic::Add(&device_ptr[index], w);
                });
            }
        }

        amrex::Gpu::copy(amrex::Gpu::deviceToHost, device_charge.begin(),
                         device_charge.end(), host_charge.begin());
        // 写入文件 清理缓存
        fileout << warpx_instance.gett_new(0);
        for (int i = 0; i < host_charge.size(); i++) {
            fileout << "\t" << host_charge[i];
        }
        fileout << "\n";
        fileout.close();

        mybpc.clearParticles();
    }
}
#endif

#ifdef EXAMINE
void PhiExamine()
{
    WarpX& warpx_instance = WarpX::GetInstance();
    auto phi = warpx_instance.m_fields.get(FieldType::phi_fp, 0);
    auto step = warpx_instance.getistep(0);
    if(step!=1)
        return;
    amrex::Box domain = warpx_instance.Geom(0).Domain();
    domain.surroundingNodes();

    auto lo = domain.smallEnd();
    auto hi = domain.bigEnd();
    const int sizex = hi[0] - lo[0] + 1, sizey = hi[1] - lo[1] + 1,
              sizez = hi[2] - lo[2];
    const int datasize = sizex * sizey * sizez;

    Vector<Real> phi3d(datasize);
    amrex::Gpu::DeviceVector<Real> phi3d_device(datasize);
    Real* pphi = phi3d_device.dataPtr();

    for (MFIter mfi(*phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& box = mfi.tilebox();

        Array4<Real> const& phiarr = phi->array(mfi);
        if (!domain.strictly_contains(box)) {
            ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (k == domain.smallEnd(2)) {
                    pphi[i * sizey * sizez + j * sizez + k] = phiarr(i, j, k);
                }
            });
        }
    }

    amrex::Gpu::copy(amrex::Gpu::deviceToHost, phi3d_device.begin(),
                     phi3d_device.end(), phi3d.begin());
    std::fstream fileout("phiexam.dat", std::ios::out);

    for (int i = 0; i < sizex; i++)
    {
        for(int j = 0; j < sizey; j++)
        {
            for (int k = 0; k < sizez; k++)
            {
                fileout << i << "\t" << j << "\t" << k << "\t"
                        << phi3d[i * sizey * sizez + j * sizez + k] << "\n";
            }
        }
    }
    fileout.close();
}

void XeRhoExamine()
{
    WarpX& warpx_instance = WarpX::GetInstance();
    auto rho = warpx_instance.m_fields.get(FieldType::rho_fp, 0);
    MultiFab xe_rho(rho->boxArray(), rho->DistributionMap(), rho->nComp(),
                    rho->nGrow());
    xe_rho.setVal(0.0);
    auto& xe_pc = mypc->GetParticleContainer(mypc->getSpeciesID("xe_netural"));
    AtomDepostiAPI(xe_pc, xe_rho);
    
    auto step = warpx_instance.getistep(0);
    if(step!=1)
        return;
    amrex::Box domain = warpx_instance.Geom(0).Domain();
    domain.surroundingNodes();

    auto lo = domain.smallEnd();
    auto hi = domain.bigEnd();
    const int sizex = hi[0] - lo[0] + 1, sizey = hi[1] - lo[1] + 1,
              sizez = hi[2] - lo[2];
    const int datasize = sizex * sizey * sizez;

    Vector<Real> phi3d(datasize);
    amrex::Gpu::DeviceVector<Real> phi3d_device(datasize);
    Real* pphi = phi3d_device.dataPtr();

    for (MFIter mfi(*phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& box = mfi.tilebox();

        Array4<Real> const& phiarr = phi->array(mfi);
        if (!domain.strictly_contains(box)) {
            ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (k == domain.smallEnd(2)) {
                    pphi[i * sizey * sizez + j * sizez + k] = phiarr(i, j, k);
                }
            });
        }
    }

    amrex::Gpu::copy(amrex::Gpu::deviceToHost, phi3d_device.begin(),
                     phi3d_device.end(), phi3d.begin());
    std::fstream fileout("phiexam.dat", std::ios::out);

    for (int i = 0; i < sizex; i++)
    {
        for(int j = 0; j < sizey; j++)
        {
            for (int k = 0; k < sizez; k++)
            {
                fileout << i << "\t" << j << "\t" << k << "\t"
                        << phi3d[i * sizey * sizez + j * sizez + k] << "\n";
            }
        }
    }
    fileout.close();
}
#endif

#ifdef HALL3D_INIT
// 初始化氙气注入
void XeFastInjection () {
    double atom_weight, l_factor, m_dot;
    bool ifhole = false;
    static amrex::Vector<double> hole_x, hole_y;
    static int hole_num = 48;

    amrex::ParmParse pp_mc("my_constants");
    pp_mc.get("l_factor", l_factor);
    pp_mc.getWithParser("xe_weight", atom_weight);
    pp_mc.get("m_dot", m_dot);
    pp_mc.query("ifhole", ifhole);

    static const double L = 0.05, // 未缩比边长
        NA = 6.03e23,             // 阿伏伽德罗常数
        vz0 = 250;                  // 气流整体速度
    const double V_per_sec = 0.06 / 60. / 1e6, dt = 5.6e-10,
                 n_per_sec = V_per_sec * 101325 / 8.314 / 273.15 * NA,
                 n_per_step = n_per_sec * dt,
                 n_marco_per_step =
                     n_per_step / atom_weight, // 每时步应当注入宏粒子;
        mxe = 2.179e-25, Pi2 = 3.1415926 * 2;
    // 质量流量
    const double n_marco_per_step_m =
        m_dot * dt / mxe / l_factor / l_factor / atom_weight;

    const double one_times_inject_particle = static_cast<int>(n_marco_per_step_m + 1); // 一次注入粒子数
    static double xe_rest = 0;
    static const double rmid = (0.021 + 0.031) / 4 / l_factor,
                        rwidth = 0.001 / l_factor, r1 = rmid - rwidth,
                        r2 = rmid + rwidth, sqdr = r2 * r2 - r1 * r1,
                        sqr1 = r1 * r1, kb = 1.38064852e-23,
                        sigmax = std::sqrt(kb * 150 / mxe),
                        sigmay = std::sqrt(kb * 150 / mxe),
                        sigmaz = std::sqrt(kb * 150 / mxe),
                        half_L = L / 2 / l_factor;
    double r, theta;
    static const Vector<Vector<int>> nattr;

    static bool if_hole_init = false;
    if (!if_hole_init) {
        double per_angle = Pi2 / hole_num, angle = 0;
        for (int i = 0; i < hole_num; i++) {
            hole_x.push_back(rmid * cos(angle));
            hole_y.push_back(rmid * sin(angle));
            angle += per_angle;
        }
        if_hole_init = true;
    }

    xe_rest += n_marco_per_step_m;
    if (xe_rest > one_times_inject_particle) {
        WarpX& warpx_instance = WarpX::GetInstance();

        xe_rest -= one_times_inject_particle;

        Vector<ParticleReal> pz(one_times_inject_particle, 0),
            px(one_times_inject_particle), py(one_times_inject_particle),
            vx(one_times_inject_particle), vy(one_times_inject_particle),
            vz(one_times_inject_particle,vz0),
            pw(one_times_inject_particle, atom_weight);

        RandomEngine uniform_engine(getRandState()),
            normal_engine(getRandState());
            
        int hole_start = rand() % hole_num;
        for (int i = 0; i < one_times_inject_particle; i++) {
            // 柱坐标空间均匀分布
            // 转化到直角坐标空间均匀分布
            
            theta = Random(uniform_engine) * Pi2;
            r = sqrt(sqdr * Random(uniform_engine) + sqr1);

            if (!ifhole) {
                // 非小孔进气
                theta = Random(uniform_engine) * Pi2;
                r = sqrt(sqdr * Random(uniform_engine) + sqr1);
                px[i] = r * cos(theta) + half_L;
                py[i] = r * sin(theta) + half_L;
            } else {
                // 小孔进气
                theta = Random(uniform_engine) * Pi2;
                r = rwidth * sqrt(Random(uniform_engine)); // 小孔进气是r1 = 0, r2 = rwidth, 化简后得到该式
                px[i] = r * cos(theta) + half_L +
                        hole_x[(i + hole_start) % hole_num];
                py[i] = r * sin(theta) + half_L +
                        hole_y[(i + hole_start) % hole_num];
            }

            // 速度
            vx[i] = RandomNormal(0, sigmax, normal_engine);
            vy[i] = RandomNormal(0, sigmay, normal_engine);
            do
            {
                vz[i] = RandomNormal(0, sigmaz, normal_engine) + vz0;
            } while (vz[i] < 0);
        }

        auto& mypc = warpx_instance.GetPartContainer();
        auto& xe_pc =
            mypc.GetParticleContainer(mypc.getSpeciesID("xe_netural"));
        xe_pc.AddNParticles(0, one_times_inject_particle, px, py, pz, vx, vy,
                            vz, 1, {pw}, 0, nattr, 0);
        amrex::Print() << "Injection Xe Atom\n";
    }
}
#endif

#ifdef COLLISION_RECORD

void ShowAndWriteIonzationNum(amrex::Vector<int> num)
{
    amrex::Print() << "produce macro electron: " << num[0]
                   << "\nproduce macro xe ions: " << num[1] << "\n";
    static bool ifinit = false;
    std::fstream fileout;
    if (!ifinit) {
        fileout.open("collision_record.dat", std::ios::out);
        ifinit = true;
    } else {
        fileout.open("collision_record.dat", std::ios::app);
    }
    fileout << num[0] << "\t" << num[1] << "\n";
    fileout.close();
}
#endif

#ifndef PUSH_GAP
/**
 * 初始化各粒子推进步长
 */
void PushGapInit()
{
    const ParmParse pp_particles("particles");
    std::vector<std::string> species_names;
    pp_particles.queryarr("species_names", species_names);

    int gap;
    for (auto i : species_names) {
        ParmParse pp_species(i);
        pp_species.get("ndt", gap);
        AMREX_ASSERT_WITH_MESSAGE(gap>0,"The gap should larger than 1.");
        push_control.push_back(std::make_pair(gap, 0));
    }
}
#endif

void
DirichletPhiGuardSet () {
    WarpX& warpx_instance = WarpX::GetInstance();
    auto rho = warpx_instance.m_fields.get(FieldType::rho_fp, 0);
    auto phi = warpx_instance.m_fields.get(FieldType::phi_fp, 0);
    amrex::Box domain = warpx_instance.Geom(0).Domain();
    domain.surroundingNodes();

    const Real* dx_host = warpx_instance.Geom(0).CellSize();
    amrex::Gpu::DeviceVector<Real> dx_device(3);
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, dx_host, dx_host + 3,
                     dx_device.begin());
    Real* dx = dx_device.dataPtr();

    for (MFIter mfi(*phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        const auto& phi_arr = phi->array(mfi);
        const auto& rho_arr = rho->array(mfi);
        ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // x direction
            if (i == domain.smallEnd(0)) {
                phi_arr(i - 1, j, k) =
                    2.0_rt * phi_arr(i, j, k) - phi_arr(i + 1, j, k) -
                    rho_arr(i, j, k) * dx[0] * dx[0] / PhysConst::epsilon_0;
            }
            if (i == domain.bigEnd(0)) {
                phi_arr(i + 1, j, k) =
                    2.0_rt * phi_arr(i, j, k) - phi_arr(i - 1, j, k) -
                    rho_arr(i, j, k) * dx[0] * dx[0] / PhysConst::epsilon_0;
            }
            // y direction
            if (j == domain.smallEnd(1)) {
                phi_arr(i, j - 1, k) =
                    2.0_rt * phi_arr(i, j, k) - phi_arr(i, j + 1, k) -
                    rho_arr(i, j, k) * dx[1] * dx[1] / PhysConst::epsilon_0;
            }
            if (j == domain.bigEnd(1)) {
                phi_arr(i, j + 1, k) =
                    2.0_rt * phi_arr(i, j, k) - phi_arr(i, j - 1, k) -
                    rho_arr(i, j, k) * dx[1] * dx[1] / PhysConst::epsilon_0;
            }
            // z direction
            if (k == domain.smallEnd(2)) {
                phi_arr(i, j, k - 1) =
                    2.0_rt * phi_arr(i, j, k) - phi_arr(i, j, k + 1) -
                    rho_arr(i, j, k) * dx[2] * dx[2] / PhysConst::epsilon_0;
            }
            if (k == domain.bigEnd(2)) {
                phi_arr(i, j, k + 1) =
                    2.0_rt * phi_arr(i, j, k) - phi_arr(i, j, k - 1) -
                    rho_arr(i, j, k) * dx[2] * dx[2] / PhysConst::epsilon_0;
            }
        });
    }
}

#ifdef MCC_DENSITY
amrex::Vector<BackgroundCoupledDensity> global_background_density;
void GlobalBackgroundDensityInit()
{
    ParmParse pp_coll("collisions");
    amrex::Vector<std::string> species_names;
    pp_coll.getarr("background_species", species_names);

    global_background_density.resize(species_names.size());
    for (int i = 0; i < species_names.size(); i++) {
        global_background_density[i].m_ground_species = species_names[i];
        global_background_density[i].backgroundDensityInit();
    }
}

void
GlobalBackgroundDensityUpdate (int step) {
    amrex::ParmParse pp_mc("my_constants");
    amrex::ParticleReal elec_weight;
    pp_mc.getWithParser("elec_weight", elec_weight);

    WarpX& warpx_instance = WarpX::GetInstance();
    MultiParticleContainer& mypc = warpx_instance.GetPartContainer();

    for (int i = 0; i < global_background_density.size(); i++) {
        int ndt = mypc.GetParticleContainerFromName(
                          global_background_density[i].m_ground_species)
                      .getndt();
        if (step % ndt == 1) {
            global_background_density[i].backgroundDensityUpdate(mypc,
                                                                 elec_weight);
        }
    }
}

void
GlobalBackgroundDensityClean (int step) {
    amrex::ParmParse pp_mc("my_constants");
    amrex::ParticleReal elec_weight;
    pp_mc.getWithParser("elec_weight", elec_weight);

    WarpX& warpx_instance = WarpX::GetInstance();
    MultiParticleContainer& mypc = warpx_instance.GetPartContainer();
    /**
     * deleteInvalidParticles在处理边界条件之后被调用，push会在当前
     * 时步的碰撞之后进行，下一时步将进行背景密度的重新计算，因此这一
     * 时步的所有碰撞完成之后，进行粒子清理
     */
    for (int i = 0; i < global_background_density.size(); i++) {
        int ndt = mypc.GetParticleContainerFromName(
                          global_background_density[i].m_ground_species)
                      .getndt();
        if (step % ndt == 0) {
            global_background_density[i].backgroudnSpeciesClean(mypc);
        }
    }
}

#endif