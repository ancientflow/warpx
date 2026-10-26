#ifndef WARPX_INSERT_ELECTROSTATIC_MATERIAL_H_
#define WARPX_INSERT_ELECTROSTATIC_MATERIAL_H_

#include <ablastr/fields/MultiFabRegister.H>

#include <AMReX_MultiFab.H>
#include <AMReX_Parser.H>
#include <AMReX_Vector.H>
#include <AMReX_iMultiFab.H>

#include <memory>
#include <string>

namespace Insert {

/** Cached material data for the 3D lab-frame electrostatic material model. */
class ElectrostaticMaterial
{
public:
    static ElectrostaticMaterial& GetInstance ();

    ElectrostaticMaterial (ElectrostaticMaterial const&) = delete;
    ElectrostaticMaterial (ElectrostaticMaterial&&) = delete;
    ElectrostaticMaterial& operator= (ElectrostaticMaterial const&) = delete;
    ElectrostaticMaterial& operator= (ElectrostaticMaterial&&) = delete;

    [[nodiscard]] bool enabled ();

    /** Initialize static material data and validate layouts without modifying rho. */
    void prepare (
        ablastr::fields::MultiLevelScalarField const& rho,
        ablastr::fields::MultiLevelScalarField const& phi,
        int max_level);

    [[nodiscard]] ablastr::fields::ConstMultiLevelScalarField const*
    relativePermittivity () const noexcept;

    [[nodiscard]] ablastr::fields::ConstMultiLevelScalarField const*
    anodePotentials () const noexcept;

    /** Refresh same-level and physical ghosts before collocated differentiation.
     * Coarse/fine ghosts are supplied by the common Poisson hierarchy.
     */
    void fillPhiGhosts (ablastr::fields::MultiLevelScalarField const& phi) const;
    void fillElectricFieldGhosts (ablastr::fields::MultiLevelVectorField const& field) const;

    [[nodiscard]] amrex::Vector<amrex::iMultiFab const*> const*
    anodeMasks () const noexcept;

private:
    ElectrostaticMaterial () = default;
    ~ElectrostaticMaterial () = default;

public:
    void readParametersOnce ();
    void initialize (
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
    std::string m_relative_permittivity_expression;
    std::string m_anode_implicit_expression;
    std::string m_anode_potential_expression;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_epsilon_r;
    amrex::Vector<std::unique_ptr<amrex::iMultiFab>> m_anode_mask;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_anode_phi;

    ablastr::fields::ConstMultiLevelScalarField m_epsilon_r_ptrs;
    ablastr::fields::ConstMultiLevelScalarField m_anode_phi_ptrs;
    amrex::Vector<amrex::iMultiFab const*> m_anode_mask_ptrs;

};

} // namespace Insert

#endif // WARPX_INSERT_ELECTROSTATIC_MATERIAL_H_
