#include "InsertBoundaryPhi.h"

#include "Fields.H"
#include "Insert/Config/WarpXFunctionConfig.h"
#include "WarpX.H"

#include <AMReX_Array4.H>
#include <AMReX_GpuContainers.H>
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
    amrex::Real correction_x = geom.ProbHi(0) - amrex::Real(0.001);
    {
        amrex::ParmParse const pp_b2d("insert.benchmark2d");
        pp_b2d.query("correction_x", correction_x);
    }
    amrex::Real const correction_frac =
        (correction_x - geom.ProbLo(0)) / geom.CellSize(0);
    int const correction_i_lo =
        static_cast<int>(std::floor(correction_frac));
    amrex::Real const correction_w_hi =
        correction_frac - amrex::Real(correction_i_lo);
    amrex::Real const correction_w_lo = amrex::Real(1.0) - correction_w_hi;

    amrex::ReduceOps<amrex::ReduceOpSum> reduce_ops;
    amrex::ReduceData<amrex::Real> reduce_data(reduce_ops);
    for (amrex::MFIter mfi(*phi_field, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box const& box = mfi.tilebox();
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
    phisum /= amrex::Real(geom.Domain().length(1));

    amrex::Print() << "voltage adjustment: " << phisum << std::endl;

    for (amrex::MFIter mfi(*phi_field, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        const amrex::Box& box = mfi.tilebox();

        amrex::Array4<amrex::Real> const& phi = phi_field->array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j) {
            phi(i, j, 0) -= phisum * i / correction_frac;
        });
    }
#endif
}

} // namespace Insert
