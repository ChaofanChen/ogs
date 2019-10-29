/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <utility>

#include <Eigen/Eigen>

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "ParameterLib/Parameter.h"
#include "WellboreGeometry.h"
#include "WellboreSimulatorMaterialProperties.h"

namespace ProcessLib
{
namespace WellboreSimulator
{
struct WellboreSimulatorProcessData final
{
    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map;
    Eigen::VectorXd const specific_body_force;
    std::vector<WellboreGeometry> wellbore;
    std::vector<double> mass_flow_rate;
    bool const has_gravity;
    std::unique_ptr<WellboreSimulatorMaterialProperties> material;
};

}  // namespace WellboreSimulator
}  // namespace ProcessLib
