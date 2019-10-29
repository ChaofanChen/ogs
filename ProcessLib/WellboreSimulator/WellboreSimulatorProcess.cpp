/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "WellboreSimulatorProcess.h"

#include <cassert>

#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"

#include "WellboreSimulatorFEM-impl.h"

namespace ProcessLib
{
namespace WellboreSimulator
{
WellboreSimulatorProcess::WellboreSimulatorProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    WellboreSimulatorProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables)),
      _process_data(std::move(process_data))
{
}

void WellboreSimulatorProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    // The order of shape function can be fetched from
    // any set of the sets of process variables of the coupled processes. Here,
    // we take the one from the first process by setting process_id = 0.
    const int process_id = 0;
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

        ProcessLib::createLocalAssemblers<MonolithicWellboreSimulatorFEM>(
            mesh.getDimension(), mesh.getElements(), dof_table,
            pv.getShapeFunctionOrder(), _local_assemblers,
            mesh.isAxiallySymmetric(), integration_order, _process_data);
}

void WellboreSimulatorProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;

    DBUG("Assemble WellboreSimulator Process.");
    dof_tables.emplace_back(*_local_to_global_index_map);

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_tables, t, dt, x, xdot, process_id, M, K,
        b, _coupled_solutions);
}

void WellboreSimulatorProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
    int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
    GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian WellboreSimulator Process.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;

    dof_tables.emplace_back(std::ref(*_local_to_global_index_map));
    dof_tables.emplace_back(std::ref(*_local_to_global_index_map));

    // Call global assembler for each local assembly item.
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_tables, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac, _coupled_solutions);
}
}  // namespace WellboreSimulator
}  // namespace ProcessLib
