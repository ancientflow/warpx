#ifndef WARPX_INSERT_HALLELECTROSTATICMATERIAL_H_
#define WARPX_INSERT_HALLELECTROSTATICMATERIAL_H_

#include <ablastr/fields/MultiFabRegister.H>

#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Parser.H>
#include <AMReX_Vector.H>
#include <AMReX_iMultiFab.H>

#include <memory>

namespace Insert {

/** Cached material data for the 3D lab-frame electrostatic Hall model. */
class HallElectrostaticMaterial
{
public:
    static HallElectrostaticMaterial& GetInstance ();

    HallElectrostaticMaterial (HallElectrostaticMaterial const&) = delete;
    HallElectrostaticMaterial (HallElectrostaticMaterial&&) = delete;
    HallElectrostaticMaterial& operator= (HallElectrostaticMaterial const&) = delete;
    HallElectrostaticMaterial& operator= (HallElectrostaticMaterial&&) = delete;

    [[nodiscard]] bool enabled ();

    void prepare (
        ablastr::fields::MultiLevelScalarField const& rho,
        ablastr::fields::MultiLevelScalarField const& phi,
        int max_level);

    [[nodiscard]] ablastr::fields::ConstMultiLevelScalarField const*
    relativePermittivity () const noexcept;

    [[nodiscard]] amrex::Vector<amrex::iMultiFab const*> const*
    anodeMasks () const noexcept;

private:
    HallElectrostaticMaterial () = default;
    ~HallElectrostaticMaterial () = default;

    void readParametersOnce ();
    [[nodiscard]] bool layoutMatches (
        ablastr::fields::MultiLevelScalarField const& phi,
        int max_level) const;
    void defineOrRebuild (
        ablastr::fields::MultiLevelScalarField const& phi,
        int max_level);

    bool m_parameters_read = false;
    bool m_enabled = false;
    bool m_defined = false;

    std::unique_ptr<amrex::Parser> m_relative_permittivity_parser;
    std::unique_ptr<amrex::Parser> m_anode_implicit_parser;
    std::unique_ptr<amrex::Parser> m_anode_potential_parser;
    amrex::ParserExecutor<3> m_relative_permittivity;
    amrex::ParserExecutor<3> m_anode_implicit;
    amrex::ParserExecutor<3> m_anode_potential;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_epsilon_r;
    amrex::Vector<std::unique_ptr<amrex::iMultiFab>> m_anode_mask;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_anode_phi;

    ablastr::fields::ConstMultiLevelScalarField m_epsilon_r_ptrs;
    amrex::Vector<amrex::iMultiFab const*> m_anode_mask_ptrs;

    amrex::Vector<amrex::Box> m_domains;
    amrex::Vector<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>> m_prob_lo;
    amrex::Vector<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>> m_prob_hi;
    amrex::Vector<amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>> m_cell_size;
};

} // namespace Insert

#endif // WARPX_INSERT_HALLELECTROSTATICMATERIAL_H_
