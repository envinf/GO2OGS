
/**
 * @file GocadNode.cpp
 * @author Thomas Fischer
 * @date 2014-06-27
 * @brief
 *
 * @copyright
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "GocadNode.h"

namespace MeshLib
{

bool operator<=(MeshLib::GocadNode const& n0, MeshLib::GocadNode const& n1)
{
	for (std::size_t k(0); k<3; k++) {
		if (n0[0] > n1[0]) return false;
		else if (n0[0] < n1[0]) return true;
		// => n0[k] == n1[k]
	}

	if (n0.getLayerTransitionIndex() > n1.getLayerTransitionIndex()) return false;
	else return true;
}

} // end namespace MeshLib

