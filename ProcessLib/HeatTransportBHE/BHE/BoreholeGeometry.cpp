/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BoreholeGeometry.h"
#include "BaseLib/ConfigTree.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
BoreholeGeometry createBoreholeGeometry(BaseLib::ConfigTree const& config)
{
    double borehole_length(0.0), borehole_diameter(0.0);
    std::vector<double> section_length, section_diameter;

    for (
        auto section_conf :
        //! \ogs_file_param{prj__chemical_system__solution__components__component}
        config.getConfigSubtreeList("borehole_section"))
    {
        const auto borehole_length_this_section =
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__borehole__length}
            section_conf.getConfigParameter<double>("length");
        const auto borehole_diameter_this_section =
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__borehole__diameter}
            section_conf.getConfigParameter<double>("diameter");

        section_length.push_back(borehole_length_this_section);
        section_diameter.push_back(borehole_diameter_this_section);
    }

    borehole_diameter = section_diameter.front();
    borehole_length = section_length.back();

    return {borehole_length, borehole_diameter, section_length,
            section_diameter};
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
