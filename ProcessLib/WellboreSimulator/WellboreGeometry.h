/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <boost/math/constants/constants.hpp>

namespace ProcessLib
{
namespace WellboreSimulator
{
struct WellboreGeometry
{
    /**
     * length/depth of the wellbore
     * unit is m
     */
    double const length;

    /**
     * diameter of the wellbore
     * unit is m
     */
    double const diameter;
};
}  // namespace WellboreSimulator
}  // namespace ProcessLib
