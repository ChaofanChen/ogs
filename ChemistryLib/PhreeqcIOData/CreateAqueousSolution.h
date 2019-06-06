/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

namespace BaseLib
{
class ConfigTree;
}

namespace ChemistryLib
{
struct AqueousSolution;
}

namespace ChemistryLib
{
AqueousSolution createAqueousSolution(BaseLib::ConfigTree const& config);
}  // namespace ChemistryLib
