/**
 * @file GocadNode.h
 * @author Thomas Fischer
 * @date Sep 26, 2013
 * @brief 
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#ifndef GOCADNODE_H_
#define GOCADNODE_H_

#include <vector>
#include <bitset>

// ThirdParty/logog
#include "logog/include/logog.hpp"

// MeshLib
#include "Node.h"

namespace MeshLib
{
enum class FaceIndicator : char
{
	U,
	V,
	W
};

class GocadNode : public Node
{
public:
	GocadNode(double const*const coords, std::size_t id,
	std::size_t layer_transition_idx) :
		Node(coords, id), _face_set_membership(),
		_layer_transition_idx(layer_transition_idx)
	{}

	GocadNode(GocadNode const& src) :
		Node(src.getCoords(), src._id),
		_face_indicators(src._face_indicators),
		_face_set_membership(src._face_set_membership),
		_layer_transition_idx(src._layer_transition_idx)
	{}

	void setFaceSet(std::size_t face_set_number, std::size_t face_indicator)
	{
		_face_set_membership.set(face_set_number);
		switch (face_indicator) {
		case 0:
			_face_indicators.push_back(std::pair<std::size_t, FaceIndicator> (face_set_number, FaceIndicator::U));
			break;
		case 1:
			_face_indicators.push_back(std::pair<std::size_t, FaceIndicator> (face_set_number, FaceIndicator::V));
			break;
		case 2:
			_face_indicators.push_back(std::pair<std::size_t, FaceIndicator> (face_set_number, FaceIndicator::W));
			break;
		default:
			ERR("GocadNode::setFaceSet(): unknown face indicator %d.", face_indicator);
			exit(1);
		}
	}

	/**
	 * Checks if this GocadNode is in the face set with the number
	 * face_set_number.
	 * @param face_set_number the number of the face set
	 * @return true/false
	 */
	bool isMemberOfFaceSet(std::size_t face_set_number) const
	{
		return _face_set_membership[face_set_number];
	}

	bool isMemberOfAnyFaceSet() const
	{
		return _face_set_membership.any();
	}

	void resetID(std::size_t id)
	{
		this->setID(id);
	}

	std::bitset<128> const& getFaceSetMembership() const { return _face_set_membership; }

	FaceIndicator getFaceIndicator(std::size_t face_set_number) const
	{
		auto it = _face_indicators.cbegin();
		while (it != _face_indicators.cend() && it->first != face_set_number) {
			it++;
		}
		if (it == _face_indicators.end()) {
			ERR("GocadNode %d: Could not found face indicator for face set %d", _id, face_set_number);
			abort ();
		}
		return it->second;
	}

	std::size_t getLayerTransitionIndex() const { return _layer_transition_idx; }

protected:
	friend class GocadSplitNode;
	std::vector<std::pair<std::size_t, FaceIndicator> > _face_indicators;

private:
	std::bitset<128> _face_set_membership;
	std::size_t _layer_transition_idx;
};

bool operator<=(MeshLib::GocadNode const& n0, MeshLib::GocadNode const& n1);

class GocadSplitNode : public GocadNode
{
public:
	GocadSplitNode(double const*const coords, std::size_t id,
			std::array<std::size_t,3> const& grid_coords,
			std::array<bool, 8> affected_cells,
			std::size_t layer_transition_idx) :
		GocadNode(coords, id, layer_transition_idx),
		_grid_coords(grid_coords) ,
		_affected_cells(affected_cells)
	{}

	std::array<std::size_t, 3> const& getGridCoords() const { return _grid_coords; }
	std::array<bool, 8> const& getAffectedCells() const { return _affected_cells; }
	void transmitFaceIndicators(GocadNode const& gocad_node) {
		_face_indicators = gocad_node._face_indicators;
	}
private:
	std::array<std::size_t,3> _grid_coords;
	std::array<bool, 8> _affected_cells;
};

} // end namespace MeshLib

#endif /* GOCADNODE_H_ */
