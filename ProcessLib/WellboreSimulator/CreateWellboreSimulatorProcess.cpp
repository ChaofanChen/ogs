/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateWellboreSimulatorProcess.h"

#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "ParameterLib/ConstantParameter.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "WellboreGeometry.h"
#include "WellboreSimulatorLocalAssemblerInterface.h"
#include "WellboreSimulatorProcess.h"
#include "WellboreSimulatorProcessData.h"

namespace ProcessLib
{
namespace WellboreSimulator
{
std::unique_ptr<Process> createWellboreSimulatorProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "WELLBORE_SIMULATOR");

    DBUG("Create WellboreSimulatorProcess.");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__WELLBORE_SIMULATOR__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;

    auto per_process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__WELLBORE_SIMULATOR__process_variables__pressure}
         "pressure",
         //! \ogs_file_param_special{prj__processes__process__WELLBORE_SIMULATOR__process_variables__enthalpy}
         "enthalpy"});
    process_variables.push_back(std::move(per_process_variables));

    // Specific body force parameter.
    Eigen::VectorXd specific_body_force;
    std::vector<double> const b =
        //! \ogs_file_param{prj__processes__process__WELLBORE_SIMULATOR__specific_body_force}
        config.getConfigParameter<std::vector<double>>("specific_body_force");
    assert(!b.empty() && b.size() < 4);
    if (b.size() < mesh.getDimension())
    {
        OGS_FATAL(
            "specific body force (gravity vector) has %d components, mesh "
            "dimension is %d",
            b.size(), mesh.getDimension());
    }
    bool const has_gravity = MathLib::toVector(b).norm() > 0;
    if (has_gravity)
    {
        specific_body_force.resize(b.size());
        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    std::vector<WellboreGeometry> wellbore_para;

    //! \ogs_file_param{prj__processes__process__WELLBORE_SIMULATOR__wellbore}
    auto const& wellbore_config = config.getConfigSubtree("wellbore");
    const auto length =
        //! \ogs_file_param{prj__processes__process__WELLBORE_SIMULATOR__wellbore__length}
        wellbore_config.getConfigParameter<double>("length");
    const auto diameter =
        //! \ogs_file_param{prj__processes__process__WELLBORE_SIMULATOR__wellbore__depth}
        wellbore_config.getConfigParameter<double>("diameter");
    WellboreGeometry const wellbore{length, diameter};
    wellbore_para.push_back(wellbore);

    std::vector<double> mass_flow_rate_data;

    //! \ogs_file_param{prj__processes__process__WELLBORE_SIMULATOR__mass_flow_rate}
    double const mass_flow_rate =
        config.getConfigParameter<double>("mass_flow_rate");
    mass_flow_rate_data.push_back(mass_flow_rate);

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);

    std::unique_ptr<WellboreSimulatorMaterialProperties> material = nullptr;

    WellboreSimulatorProcessData process_data{std::move(media_map),
                                              specific_body_force,
                                              std::move(wellbore_para),
                                              std::move(mass_flow_rate_data),
                                              has_gravity,
                                              std::move(material)};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<WellboreSimulatorProcess>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables));
}

}  // namespace WellboreSimulator
}  // namespace ProcessLib
