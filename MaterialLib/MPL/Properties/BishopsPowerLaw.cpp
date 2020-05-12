/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BishopsPowerLaw.h"

namespace MaterialPropertyLib
{
BishopsPowerLaw::BishopsPowerLaw(std::string name, double const exponent)
    : m_(exponent)
{
    name_ = std::move(name);
}

void BishopsPowerLaw::setScale(
    std::variant<Medium*, Phase*, Component*> scale_pointer)
{
    if (!std::holds_alternative<Medium*>(scale_pointer))
    {
        OGS_FATAL(
            "The property 'BishopsPowerLaw' is implemented on the 'media' "
            "scale only.");
    }
    scale_ = scale_pointer;
}

PropertyDataType BishopsPowerLaw::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    auto const S_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    return std::pow(S_L, m_);
}

PropertyDataType BishopsPowerLaw::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    (void)variable;
    assert((variable == Variable::liquid_saturation) &&
           "BishopsPowerLaw::dvalue is implemented for derivatives with "
           "respect to liquid saturation only.");

    auto const S_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    return m_ * std::pow(S_L, m_ - 1.);
}
}  // namespace MaterialPropertyLib
