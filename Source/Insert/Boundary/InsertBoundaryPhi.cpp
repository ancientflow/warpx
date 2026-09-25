#include "InsertBoundaryPhi.h"

#include "Fields.H"
#include "Insert/Config/WarpXFunctionConfig.h"
#include "Utils/TextMsg.H"
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

    // MLMG filled the ghosts for the uncorrected potential. Apply the same
    // affine correction there before computeE takes centered differences.
    // Iterate whole FABs without tiling so each ghost is updated exactly once.
    for (amrex::MFIter mfi(*phi_field); mfi.isValid(); ++mfi) {
        amrex::Box const box = mfi.fabbox();

        amrex::Array4<amrex::Real> const& phi = phi_field->array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            phi(i, j, k) -= phisum * amrex::Real(i - xlo) / correction_frac;
        });
    }
#endif
}

} // namespace Insert
