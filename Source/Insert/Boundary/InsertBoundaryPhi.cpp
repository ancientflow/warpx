#include "InsertBoundaryPhi.h"

#include "Fields.H"
#include "Insert/Config/WarpXFunctionConfig.h"
#include "Insert/Utils/InsertUtils.h"
#include "Utils/TextMsg.H"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>
#include <AMReX_Tuple.H>

#include <cmath>
#include <iostream>

namespace Insert {

void
AnodeVoltage ()
{
#if defined(WARPX_DIM_3D)
    // Runtime switch: insert.zmin_ring_anode = 1 sets the zmin domain-face
    // nodes inside the Hall anode ring to my_constants.voltage (the rest of
    // the face keeps boundary.potential_lo_z). Skipped when the electrostatic
    // material (volume anode) path is active; see Insert::SetBoundaryPhi.
    static bool const enabled = [] () {
        amrex::ParmParse const pp_insert("insert");
        int flag = 0;
        pp_insert.query("zmin_ring_anode", flag);
        return flag != 0;
    }();
    if (!enabled) {
        return;
    }

    WarpX& warpx_instance = WarpX::GetInstance();
    auto phi_field = warpx_instance.m_fields.get(warpx::fields::FieldType::phi_fp, 0);
    auto const config = ReadHallAnodeRingConfig(warpx_instance.Geom(0));

    amrex::Box domain = warpx_instance.Geom(0).Domain();
    domain.surroundingNodes();
    int const xlo = domain.smallEnd(0);
    int const ylo = domain.smallEnd(1);
    int const zlo = domain.smallEnd(2);
    amrex::Real const problo_x = warpx_instance.Geom(0).ProbLo(0);
    amrex::Real const problo_y = warpx_instance.Geom(0).ProbLo(1);
    amrex::Real const dx = warpx_instance.Geom(0).CellSize(0);
    amrex::Real const dy = warpx_instance.Geom(0).CellSize(1);

    for (amrex::MFIter mfi(*phi_field, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const amrex::Box& box = mfi.tilebox();
        if (domain.strictly_contains(box)) {
            continue;
        }
        amrex::Array4<amrex::Real> const& phi = phi_field->array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            if (k == zlo) {
                amrex::Real const x = problo_x + amrex::Real(i - xlo) * dx;
                amrex::Real const y = problo_y + amrex::Real(j - ylo) * dy;
                if (IsHallAnodeRingHit(x, y, config)) {
                    phi(i, j, k) = config.voltage;
                }
            }
        });
    }
    phi_field->FillBoundary(warpx_instance.Geom(0).periodicity());
#endif
}

void
DirichletPhiGuardSet ()
{
#if defined(WARPX_DIM_3D)
    using namespace amrex::literals;

    // MLMG leaves Dirichlet-face ghost nodes at zero for nodal solves
    // (mlndlap_applybc only fills Neumann faces), and the collocated path
    // copies them back into phi_fp (setFinalFillBC). computeE reads those
    // ghost nodes in the centered difference at boundary nodes, so refill
    // them such that the centered difference becomes the Poisson-consistent
    // one-sided flux: phi_ghost = 2*phi_bnd - phi_int - rho_bnd*dx^2/eps0.
    if (WarpX::grid_type != ablastr::utils::enums::GridType::Collocated) {
        return;
    }

    WarpX& warpx_instance = WarpX::GetInstance();
    auto phi = warpx_instance.m_fields.get(warpx::fields::FieldType::phi_fp, 0);
    auto rho = warpx_instance.m_fields.get(warpx::fields::FieldType::rho_fp, 0);

    amrex::Box domain = warpx_instance.Geom(0).Domain();
    domain.surroundingNodes();
    int const xlo = domain.smallEnd(0);
    int const ylo = domain.smallEnd(1);
    int const zlo = domain.smallEnd(2);
    int const xhi = domain.bigEnd(0);
    int const yhi = domain.bigEnd(1);
    int const zhi = domain.bigEnd(2);

    amrex::Geometry const& geom = warpx_instance.Geom(0);
    amrex::Real const cx = geom.CellSize(0) * geom.CellSize(0) / PhysConst::epsilon_0;
    amrex::Real const cy = geom.CellSize(1) * geom.CellSize(1) / PhysConst::epsilon_0;
    amrex::Real const cz = geom.CellSize(2) * geom.CellSize(2) / PhysConst::epsilon_0;

    for (amrex::MFIter mfi(*phi, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const amrex::Box& box = mfi.validbox();
        if (domain.strictly_contains(box)) {
            continue;
        }
        amrex::Array4<amrex::Real> const& phi_arr = phi->array(mfi);
        amrex::Array4<amrex::Real const> const& rho_arr = rho->const_array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            amrex::Real const rho_bc = rho_arr(i, j, k);
            if (i == xlo) {
                phi_arr(i - 1, j, k) = 2.0_rt * phi_arr(i, j, k) -
                                       phi_arr(i + 1, j, k) - rho_bc * cx;
            }
            if (i == xhi) {
                phi_arr(i + 1, j, k) = 2.0_rt * phi_arr(i, j, k) -
                                       phi_arr(i - 1, j, k) - rho_bc * cx;
            }
            if (j == ylo) {
                phi_arr(i, j - 1, k) = 2.0_rt * phi_arr(i, j, k) -
                                       phi_arr(i, j + 1, k) - rho_bc * cy;
            }
            if (j == yhi) {
                phi_arr(i, j + 1, k) = 2.0_rt * phi_arr(i, j, k) -
                                       phi_arr(i, j - 1, k) - rho_bc * cy;
            }
            if (k == zlo) {
                    phi_arr(i, j, k - 1) = 2.0_rt * phi_arr(i, j, k) -
                                           phi_arr(i, j, k + 1) - rho_bc * cz;
                }
            if (k == zhi) {
                phi_arr(i, j, k + 1) = 2.0_rt * phi_arr(i, j, k) -
                                       phi_arr(i, j, k - 1) - rho_bc * cz;
            }
        });
    }
#endif
}

void
VoltageAdjustment ()
{
#if defined(WARPX_DIM_XZ) && defined(BENCHMARK_2D)
    WarpX& warpx_instance = WarpX::GetInstance();
    auto phi_field = warpx_instance.m_fields.get(warpx::fields::FieldType::phi_fp, 0);

    // The correction voltage is sampled at a fixed physical x position by
    // linearly interpolating between the two neighboring grid columns. This
    // avoids hard-coded column indices, which shift whenever the grid count
    // is rounded to a power of two.
    amrex::Geometry const& geom = warpx_instance.Geom(0);
    amrex::Box const& domain = geom.Domain();
    int const xlo = domain.smallEnd(0);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(!geom.isPeriodic(0) && geom.isPeriodic(1),
        "The 2D benchmark requires nonperiodic x and periodic z boundaries.");
    amrex::Real correction_x = geom.ProbHi(0) - amrex::Real(0.001);
    {
        amrex::ParmParse const pp_b2d("insert.benchmark2d");
        pp_b2d.query("correction_x", correction_x);
    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        correction_x > geom.ProbLo(0) && correction_x <= geom.ProbHi(0),
        "insert.benchmark2d.correction_x must lie above xlo and at or below xhi.");
    amrex::Real const correction_frac =
        (correction_x - geom.ProbLo(0)) / geom.CellSize(0);
    int const correction_i_lo =
        xlo + static_cast<int>(std::floor(correction_frac));
    amrex::Real const correction_w_hi =
        correction_frac - amrex::Real(correction_i_lo - xlo);
    amrex::Real const correction_w_lo = amrex::Real(1.0) - correction_w_hi;

    amrex::ReduceOps<amrex::ReduceOpSum> reduce_ops;
    amrex::ReduceData<amrex::Real> reduce_data(reduce_ops);
    for (amrex::MFIter mfi(*phi_field, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        // Use the lower node of each cell to count shared and periodic nodes
        // exactly once, independent of tiling and MPI decomposition. Include
        // the final x node as well, since only z is periodic.
        amrex::Box box = mfi.tilebox(amrex::IntVect::TheZeroVector());
        if (box.bigEnd(0) == domain.bigEnd(0)) {
            box.growHi(0, 1);
        }
        amrex::Array4<amrex::Real const> const& phi = phi_field->const_array(mfi);
        for (int correction_col = 0; correction_col < 2; ++correction_col) {
            int const correction_i = correction_i_lo + correction_col;
            amrex::Real const correction_weight =
                correction_col == 0 ? correction_w_lo : correction_w_hi;
            amrex::Box column_box = box;
            column_box.setSmall(0, correction_i);
            column_box.setBig(0, correction_i);
            if (!box.intersects(column_box)) {
                continue;
            }
            reduce_ops.eval(
                column_box, reduce_data,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) -> amrex::GpuTuple<amrex::Real> {
                    return {correction_weight * phi(i, j, k)};
                });
        }
    }

    amrex::Real phisum = amrex::get<0>(reduce_data.value());
    amrex::ParallelDescriptor::ReduceRealSum(phisum);
    phisum /= amrex::Real(domain.length(1));

    amrex::Print() << "voltage adjustment: " << phisum << std::endl;

    // Apply the same affine correction to the ghost nodes before computeE
    // takes centered differences. Iterate whole FABs without tiling so each
    // ghost is updated exactly once.
    for (amrex::MFIter mfi(*phi_field); mfi.isValid(); ++mfi) {
        amrex::Box const box = mfi.fabbox();

        amrex::Array4<amrex::Real> const& phi = phi_field->array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            phi(i, j, k) -= phisum * amrex::Real(i - xlo) / correction_frac;
        });
    }

    // MLMG leaves Dirichlet-side ghost nodes unfilled (zero) for nodal solves,
    // while computeE reads them at the boundary nodes. Refill them by quadratic
    // extrapolation so the centered difference there becomes the second-order
    // one-sided difference of the corrected potential.
    int const xhi_node = xlo + domain.length(0);
    for (amrex::MFIter mfi(*phi_field); mfi.isValid(); ++mfi) {
        amrex::Box const box = mfi.fabbox();

        amrex::Array4<amrex::Real> const& phi = phi_field->array(mfi);
        if (box.smallEnd(0) < xlo) {
            amrex::Box ghost_lo = box;
            ghost_lo.setBig(0, xlo - 1);
            amrex::ParallelFor(ghost_lo, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                amrex::Real const g = amrex::Real(xlo - i);
                phi(i, j, k) = amrex::Real(0.5) * (g + 1) * (g + 2) * phi(xlo, j, k)
                             - g * (g + 2) * phi(xlo + 1, j, k)
                             + amrex::Real(0.5) * g * (g + 1) * phi(xlo + 2, j, k);
            });
        }
        if (box.bigEnd(0) > xhi_node) {
            amrex::Box ghost_hi = box;
            ghost_hi.setSmall(0, xhi_node + 1);
            amrex::ParallelFor(ghost_hi, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                amrex::Real const g = amrex::Real(i - xhi_node);
                phi(i, j, k) = amrex::Real(0.5) * (g + 1) * (g + 2) * phi(xhi_node, j, k)
                             - g * (g + 2) * phi(xhi_node - 1, j, k)
                             + amrex::Real(0.5) * g * (g + 1) * phi(xhi_node - 2, j, k);
            });
        }
    }
#endif
}

} // namespace Insert
