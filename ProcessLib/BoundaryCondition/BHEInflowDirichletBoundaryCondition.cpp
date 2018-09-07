/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BHEInflowDirichletBoundaryCondition.h"

#include <algorithm>
#include <logog/include/logog.hpp>
#include <vector>
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
void BHEInflowDirichletBoundaryCondition::getEssentialBCValues(
    const double t, GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    SpatialPosition pos;

    // auto const& bulk_node_ids_map =
    //     *_bc_mesh.getProperties().getPropertyVector<std::size_t>(
    //         "bulk_node_ids");

    bc_values.ids.clear();
    bc_values.values.clear();

    bc_values.ids.resize(_bc_values.ids.size());
    bc_values.values.resize(_bc_values.values.size());

    const size_t n_nodes = _T_out_values.size();
    double tmp_T_in(320.0); 
    for (size_t i=0; i<n_nodes; i++)
    {
        bc_values.ids[i] = _bc_values.ids[i];
        // TODO, here call the BHE functions
        // tmp_T_in = _BHE_property->;
        bc_values.values[i] = tmp_T_in;

    }

}

// update new values and corresponding indices.
void BHEInflowDirichletBoundaryCondition::preTimestep(
    const double t, const GlobalVector& x)
{
    // for each BHE, the inflow temperature is dependent on
    // the ouflow temperature of the BHE. 
    // Here the task is to get the outflow temperature and
    // save it locally
    auto const mesh_id = _bc_mesh.getID();
    auto const& nodes = _bc_mesh.getNodes();
    for (auto const* n : nodes)
    {
        std::size_t node_id = n->getID();
        auto g_idx = _bc_values.ids.at(node_id); 

        // read the T_out
        // TODO: notice that node_id + 1 is not always T_out!
        _T_out_values[node_id] = x[g_idx + 1]; 

    }
}

std::unique_ptr<BHEInflowDirichletBoundaryCondition>
createBHEInflowDirichletBoundaryCondition(
    BaseLib::ConfigTree const& config,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
    MeshLib::Mesh const& bc_mesh, int const variable_id,
    unsigned const integration_order, std::size_t const bulk_mesh_id,
    int const component_id,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters)
{
    DBUG(
        "Constructing BHEInflowDirichletBoundaryCondition from config.");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter(
        "type", "BHEInflowDirichlet");

    // find the corresponding BHE id
    auto const bhe_id = config.getConfigParameter<std::size_t>("bhe_id");

    // TODO: try to locate the BHE instance here
    // if BHE instance cannot be found,
    // then complain about it.


    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__Dirichlet__parameter}
    auto const param_name = config.getConfigParameter<std::string>("parameter");
    DBUG("Using parameter %s", param_name.c_str());

    auto& param = findParameter<double>(param_name, parameters, 1);
    return std::make_unique<BHEInflowDirichletBoundaryCondition>(
        param, dof_table_bulk, bc_mesh, variable_id, integration_order,
        bulk_mesh_id, component_id);
}

}  // namespace ProcessLib
