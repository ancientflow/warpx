#include "Insert/Boundary/AnalyticBoundaryGeometry.h"

#include "Utils/Parser/ParserUtils.H"

namespace Insert {

AnalyticBoundaryGeometry::AnalyticBoundaryGeometry (
    std::string const& signed_value_expression,
    std::string const& normal_x_expression,
    std::string const& normal_y_expression,
    std::string const& normal_z_expression)
{
    // Build the four parsers from the user-provided expressions. All share
    // the variable set (x, y, z).
    m_signed_value_parser = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(signed_value_expression, {"x", "y", "z"}));
    m_normal_x_parser = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(normal_x_expression, {"x", "y", "z"}));
    m_normal_y_parser = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(normal_y_expression, {"x", "y", "z"}));
    m_normal_z_parser = std::make_unique<amrex::Parser>(
        utils::parser::makeParser(normal_z_expression, {"x", "y", "z"}));

    // Compile the device executors only after all parsers are in place. The
    // executors reference parser-owned memory, which is why this object is
    // non-copyable and must outlive any use of its DeviceView.
    m_view.m_signed_value =
        utils::parser::compileParser<3>(m_signed_value_parser.get());
    m_view.m_normal_x =
        utils::parser::compileParser<3>(m_normal_x_parser.get());
    m_view.m_normal_y =
        utils::parser::compileParser<3>(m_normal_y_parser.get());
    m_view.m_normal_z =
        utils::parser::compileParser<3>(m_normal_z_parser.get());
}

} // namespace Insert
