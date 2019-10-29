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

#include <Eigen/Dense>
#include <vector>

#include "WellboreSimulatorProcessData.h"

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"

#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "WellboreSimulatorFEM.h"
#include "iostream"

namespace ProcessLib
{
namespace WellboreSimulator
{
const unsigned NUM_NODAL_DOF = 2;

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class MonolithicWellboreSimulatorFEM
    : public WellboreSimulatorFEM<ShapeFunction, IntegrationMethod, GlobalDim>
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalMatrixType = typename ShapeMatricesType::template MatrixType<
        NUM_NODAL_DOF * ShapeFunction::NPOINTS,
        NUM_NODAL_DOF * ShapeFunction::NPOINTS>;
    using LocalVectorType =
        typename ShapeMatricesType::template VectorType<NUM_NODAL_DOF *
                                                        ShapeFunction::NPOINTS>;

    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;

    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;

public:
    MonolithicWellboreSimulatorFEM(
        MeshLib::Element const& element,
        std::size_t const local_matrix_size,
        bool is_axially_symmetric,
        unsigned const integration_order,
        WellboreSimulatorProcessData const& process_data)
        : WellboreSimulatorFEM<ShapeFunction, IntegrationMethod, GlobalDim>(
              element, local_matrix_size, is_axially_symmetric,
              integration_order, process_data, NUM_NODAL_DOF)
    {
        const double* coords_1 = element.getNode(0)->getCoords();
        const double* coords_2 = element.getNode(1)->getCoords();

        typename ShapeMatricesType::template VectorType<3>
            _element_direction_vector;
    }

    void assemble(double const t, double const /*dt*/,
                  std::vector<double> const& local_x,
                  std::vector<double> const& /*local_xdot*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override
    {
        auto const local_matrix_size = local_x.size();
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

        auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
            local_K_data, local_matrix_size, local_matrix_size);
        auto local_b = MathLib::createZeroedVector<LocalVectorType>(
            local_b_data, local_matrix_size);

        ParameterLib::SpatialPosition pos;
        pos.setElementID(this->_element.getID());

        double p_nodal_value = Eigen::Map<const NodalVectorType>(
                                   &local_x[pressure_index], pressure_size)
                                   .data()[0];
        double h_nodal_value = Eigen::Map<const NodalVectorType>(
                                   &local_x[enthalpy_index], enthalpy_size)
                                   .data()[0];

        auto const& process_data = this->_process_data;

        double const b = process_data.specific_body_force[0];
        double const mass_flow_rate = process_data.mass_flow_rate[0];
        double const d = process_data.wellbore[0].diameter;

        // The following coefficients are dependent on the p and h,
        // need to be calculated and updated in every element.

        double const rho_g = 15, rho_l = 821.9, h_sat_g = 2.8e6, h_sat_l = 1e6,
                     dh_sat_l_dp = 0.1, dh_sat_g_dp = 0.1, drho_g_dp = 0.01,
                     drho_l_dp = 0.01, miu_g = 0.1, miu_l = 0.0004;
        double const v_sg = 0, C_0 = 1.2, v_infin = 0.3;
        double const epsilon = 0.00984;
        // If T_f > T_ei, Q < 0. (production well)
        double const q = -0.005;

        double x = (h_nodal_value - h_sat_l) / (h_sat_g - h_sat_l);
        if (x > 1.0)
            x = 1.0;
        else if (x < 0)
            x = 0;
        auto const v_m =
            mass_flow_rate * x / rho_g + mass_flow_rate * (1 - x) / rho_l;
        auto const miu_m = miu_g * x + miu_l * (1 - x);
        double f_g = v_sg / (C_0 * v_m + v_infin);
        if (f_g > 1)
            f_g = 1;
        else if (f_g < 0)
            f_g = 0;
        auto const rho_m = rho_g * f_g + rho_l * (1 - f_g);
        auto const Re_m = d * v_m * rho_m / miu_m;
        auto const Lambda = std::pow(epsilon / d, 1.1098) / 2.8257 +
                            std::pow(7.149 / Re_m, 0.8981);
        auto const f_m = 1 / (4 * std::log(epsilon / 3.7065 * d) -
                              5.0452 * std::log(Lambda) / Re_m);

        auto const A1 =
            1 + (std::pow(mass_flow_rate, 2) *
                 ((1 / rho_g - 1 / rho_l) *
                      ((h_nodal_value - h_sat_g) /
                           std::pow((h_sat_g - h_sat_l), 2) * dh_sat_l_dp +
                       (h_sat_l - h_nodal_value) /
                           std::pow((h_sat_g - h_sat_l), 2) * dh_sat_g_dp) -
                  x * drho_g_dp / (rho_g * rho_g) +
                  (x - 1) * drho_l_dp / (rho_l * rho_l))) *
                    rho_m * (x / rho_g + (1 - x) / rho_l);
        auto const B1 = std::pow(mass_flow_rate, 2) *
                        ((1 / rho_g - 1 / rho_l) / (h_sat_g - h_sat_l) -
                         x / (rho_g * rho_g) + (x - 1) / (rho_l * rho_l)) *
                        rho_m * (x / rho_g + (1 - x) / rho_l);
        auto const C = rho_m * b + (f_m * std::pow(mass_flow_rate, 2) *
                                    (x / rho_g + (1 - x) / rho_l) *
                                    (x / rho_g + (1 - x) / rho_l) * rho_m) /
                                       (2 * d);
        auto const A2 =
            std::pow(mass_flow_rate, 2) *
            ((1 / rho_g - 1 / rho_l) *
                 ((h_nodal_value - h_sat_g) / std::pow((h_sat_g - h_sat_l), 2) *
                      dh_sat_l_dp +
                  (h_sat_l - h_nodal_value) / std::pow((h_sat_g - h_sat_l), 2) *
                      dh_sat_g_dp) -
             x * drho_g_dp / (rho_g * rho_g) +
             (x - 1) * drho_l_dp / (rho_l * rho_l)) *
            (x / rho_g + (1 - x) / rho_l);
        auto const B2 =
            1 +
            std::pow(mass_flow_rate, 2) * (x / rho_g + (1 - x) / rho_l) *
                ((1 / rho_g - 1 / rho_l) / (h_sat_g - h_sat_l) -
                 x / (rho_g * rho_g) + (x - 1) * drho_l_dp / (rho_l * rho_l));

        auto const D = q / mass_flow_rate + b;

        auto const& medium =
            *process_data.media_map->getMedium(this->_element.getID());

        Eigen::Vector2d b_pressure(C, C);
        Eigen::Vector2d b_enthalpy(D, D);

        MaterialPropertyLib::VariableArray vars;

        unsigned const n_integration_points =
            this->_integration_method.getNumberOfPoints();

        for (unsigned ip(0); ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);

            auto const& ip_data = this->_ip_data[ip];
            auto const& N = ip_data.N;
            auto const& dNdx = ip_data.dNdx;
            auto const& w = ip_data.integration_weight;

            // matrix assembly
            local_K
                .template block<pressure_size, pressure_size>(pressure_index,
                                                              pressure_index)
                .noalias() += N.transpose() * dNdx * A1 * w;
            local_K
                .template block<enthalpy_size, enthalpy_size>(pressure_index,
                                                              enthalpy_index)
                .noalias() += N.transpose() * dNdx * B1 * w;
            local_K
                .template block<pressure_size, enthalpy_size>(enthalpy_index,
                                                              pressure_index)
                .noalias() += N.transpose() * dNdx * A2 * w;
            local_K
                .template block<enthalpy_size, enthalpy_size>(enthalpy_index,
                                                              enthalpy_index)
                .noalias() += N.transpose() * dNdx * B2 * w;

            local_b.template block<2, 1>(0, 0).noalias() += b_pressure * w;
            local_b.template block<2, 1>(2, 0).noalias() += b_enthalpy * w;
        }
        // debugging
        // std::string sep = "\n----------------------------------------\n";
        // Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
        // std::cout << local_K.format(CleanFmt) << sep;
        // std::cout << local_b.format(CleanFmt) << sep;
    }

private:
    using WellboreSimulatorFEM<ShapeFunction, IntegrationMethod,
                               GlobalDim>::pressure_index;
    using WellboreSimulatorFEM<ShapeFunction, IntegrationMethod,
                               GlobalDim>::pressure_size;
    using WellboreSimulatorFEM<ShapeFunction, IntegrationMethod,
                               GlobalDim>::enthalpy_index;
    using WellboreSimulatorFEM<ShapeFunction, IntegrationMethod,
                               GlobalDim>::enthalpy_size;
};
}  // namespace WellboreSimulator
}  // namespace ProcessLib
