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

#include <array>

#include "ProcessLib/Process.h"
#include "WellboreSimulatorProcessData.h"

namespace NumLib
{
class LocalToGlobalIndexMap;
}

namespace ProcessLib
{

namespace WellboreSimulator
{
class WellboreSimulatorLocalAssemblerInterface;

class WellboreSimulatorProcess final : public Process
{
public:
    WellboreSimulatorProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        WellboreSimulatorProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables);
    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return false; }
    //! @}
    //!
private:
    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, double const dt,
                                 std::vector<GlobalVector*> const& x,
                                 std::vector<GlobalVector*> const& xdot,
                                 int const process_id, GlobalMatrix& M,
                                 GlobalMatrix& K, GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
        int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        GlobalMatrix& Jac) override;

    WellboreSimulatorProcessData _process_data;

    std::vector<std::unique_ptr<WellboreSimulatorLocalAssemblerInterface>>
        _local_assemblers;
};

}  // namespace WellboreSimulator
}  // namespace ProcessLib
