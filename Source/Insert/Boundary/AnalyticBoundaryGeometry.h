#ifndef WARPX_INSERT_ANALYTICBOUNDARYGEOMETRY_H_
#define WARPX_INSERT_ANALYTICBOUNDARYGEOMETRY_H_

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>
#include <AMReX_Dim3.H>

#include <memory>
#include <string>

namespace Insert {

/** Particle position kept in particle precision, independently of field precision. */
struct AnalyticBoundaryPosition
{
    amrex::ParticleReal x;
    amrex::ParticleReal y;
    amrex::ParticleReal z;
};

/** \brief Geometry operator for an analytic particle boundary.
 *
 *  Stores the boundary expression and the three component functions of the
 *  normal vector field pointing into the computational domain. All four are
 *  amrex::Parser expressions in the variables (x, y, z), provided
 *  independently by the caller: the normal field is NOT derived from the
 *  boundary expression.
 *
 *  Sign convention: SignedValue(x) > 0 on the computational-domain side,
 *  < 0 on the solid side of the boundary. Only the sign and the zero
 *  crossing are meaningful; the value does not need to be a true Euclidean
 *  distance.
 *
 *  The device-side interface is DeviceView, a trivially copyable struct of
 *  amrex::ParserExecutor<3> that can be captured by value in GPU lambdas.
 *  The executors reference memory owned by this object, so the
 *  AnalyticBoundaryGeometry instance must outlive any use of its DeviceView
 *  and is therefore non-copyable and non-movable.
 */
class AnalyticBoundaryGeometry
{
public:
    struct DeviceView
    {
        amrex::ParserExecutor<3> m_signed_value;
        amrex::ParserExecutor<3> m_normal_x;
        amrex::ParserExecutor<3> m_normal_y;
        amrex::ParserExecutor<3> m_normal_z;

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal
        SignedValue (AnalyticBoundaryPosition const& x) const noexcept
        {
            // Evaluate the boundary expression F(x, y, z). By convention
            // F > 0 on the computational-domain side and F < 0 on the solid
            // side; only the sign and the zero crossing are meaningful.
            return static_cast<amrex::ParticleReal>(
                m_signed_value(x.x, x.y, x.z));
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal
        NormalX (AnalyticBoundaryPosition const& x) const noexcept
        {
            // x component of the domain-pointing normal at x.
            return static_cast<amrex::ParticleReal>(m_normal_x(x.x, x.y, x.z));
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal
        NormalY (AnalyticBoundaryPosition const& x) const noexcept
        {
            // y component of the domain-pointing normal at x.
            return static_cast<amrex::ParticleReal>(m_normal_y(x.x, x.y, x.z));
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        amrex::ParticleReal
        NormalZ (AnalyticBoundaryPosition const& x) const noexcept
        {
            // z component of the domain-pointing normal at x.
            return static_cast<amrex::ParticleReal>(m_normal_z(x.x, x.y, x.z));
        }

        /** \brief Normal vector pointing into the computational domain,
         *         assembled from the three component functions. Assumes x is
         *         on the boundary. Not necessarily normalized. */
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        AnalyticBoundaryPosition
        Normal (AnalyticBoundaryPosition const& x) const noexcept
        {
            // Assemble the full normal vector from the three independently
            // defined component functions.
            return {NormalX(x), NormalY(x), NormalZ(x)};
        }
    };

    AnalyticBoundaryGeometry (
        std::string const& signed_value_expression,
        std::string const& normal_x_expression,
        std::string const& normal_y_expression,
        std::string const& normal_z_expression);

    AnalyticBoundaryGeometry (AnalyticBoundaryGeometry const&) = delete;
    AnalyticBoundaryGeometry (AnalyticBoundaryGeometry&&) = delete;
    AnalyticBoundaryGeometry& operator= (AnalyticBoundaryGeometry const&) = delete;
    AnalyticBoundaryGeometry& operator= (AnalyticBoundaryGeometry&&) = delete;

    [[nodiscard]] DeviceView
    GetDeviceView () const noexcept { return m_view; }

private:
    std::unique_ptr<amrex::Parser> m_signed_value_parser;
    std::unique_ptr<amrex::Parser> m_normal_x_parser;
    std::unique_ptr<amrex::Parser> m_normal_y_parser;
    std::unique_ptr<amrex::Parser> m_normal_z_parser;
    DeviceView m_view;
};

} // namespace Insert

#endif // WARPX_INSERT_ANALYTICBOUNDARYGEOMETRY_H_
