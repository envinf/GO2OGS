/**
 * @file GocadSGridReader.cpp
 * @author Thomas Fischer
 * @date Mar 7, 2013
 * @brief 
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "GocadSGridReader.h"

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// boost
#include <boost/tokenizer.hpp>

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "FileTools.h"

// GeoLib
#include "AABB.h"

// MeshLib
#include "Elements/Hex.h"
#include "Elements/Quad.h"

#include "GocadNode.h"

namespace FileIO
{

typedef boost::dynamic_bitset<> Bitset;

typedef GocadSGridReader::Region Region;
typedef GocadSGridReader::Layer Layer;
typedef GocadSGridReader::GocadProperty GocadProperty;

std::ostream& operator<<(std::ostream& os, Region const& r)
{
	return os << "(" << r.name << "|" << r.bit << ")";
}

std::ostream& operator<<(std::ostream& os, Layer const& l)
{
	std::copy(l.regions.begin(), l.regions.end(),
		std::ostream_iterator<Region>(os, " "));
	return os;
}

std::ostream& operator<<(std::ostream& os, GocadProperty const& p)
{
	return os << p._property_name << " " << p._property_id << " "
			<< p._property_data_type << " " << p._property_data_fname;
}

Region parseRegion(std::string const& line)
{
	std::istringstream iss(line);
	std::istream_iterator<std::string> it(iss);
	// Check first word is REGION or MODEL_REGION.
	if (*it != std::string("REGION") && *it != std::string("MODEL_REGION"))
	{
		ERR("Expected REGION or MODEL_REGION keyword but \"%s\" found.\n", it->c_str());
		throw std::runtime_error("In parseRegion() expected REGION or MODEL_REGION keyword not found.\n");
	}
	++it;

	Region r;
	r.name = *it;
	++it;
	r.bit = atoi(it->c_str());

	return r;
}

Layer parseLayer(std::string const& line, std::vector<Region> const& regions)
{
	std::istringstream iss(line);
	std::istream_iterator<std::string> it(iss);
	// Check first word is MODEL_LAYER.
	if (*it != std::string("MODEL_LAYER"))
	{
		ERR("Expected MODEL_LAYER keyword but \"%s\" found.\n", it->c_str());
		throw std::runtime_error("In parseRegion() expected MODEL_LAYER keyword not found.\n");
	}
	++it;

	Layer l;
	while (it != std::istream_iterator<std::string>() && *it != "END")
	{
		l.regions.push_back(
			*std::find_if(regions.begin(), regions.end(),
				[&](Region const& r) { return r.name == *it; }));
		++it;
	}

	return l;
}

GocadProperty parseGocadPropertyMetaData(std::string &line, std::istream &in, std::string const& path)
{
	boost::char_separator<char> sep("\t ");
	boost::tokenizer<boost::char_separator<char> > tokens(line, sep);
	auto tok_it(tokens.begin());
	// A property section in Gocad file starts with a line
	// PROPERTY id "property_name"
	if (tok_it->compare("PROPERTY") != 0 ) {
		ERR("Expected PROPERTY keyword but \"%s\" found.", tok_it->c_str());
		throw std::runtime_error("In parseGocadPropertyMetaData() expected PROPERTY keyword not found.");
	}
	tok_it++;

	GocadProperty prop;
	prop._property_id = BaseLib::str2number<std::size_t>(*tok_it);
	tok_it++;
	prop._property_name = *tok_it;
	BaseLib::trim(prop._property_name, '\"');

	while (getline(in, line)) {
		tokens.assign(line);

		tok_it = tokens.begin();
		// this is the last entry of the property
		if (tok_it->compare("PROP_FILE") == 0) {
			tok_it++;
			if (! prop.checkID(*tok_it)) {
				throw std::runtime_error("parseGocadPropertyMetaData(): id mismatch.");
			}
			tok_it++;
			std::string tmp(*tok_it);
			BaseLib::trim(tmp, '\"');
			prop._property_data_fname = path + tmp;
			return prop;
		}

		if (tok_it->compare("PROPERTY_CLASS") == 0) {
			tok_it++;
			if (! prop.checkID(*tok_it)) {
				throw std::runtime_error("parseGocadPropertyMetaData(): id mismatch.");
			}
			tok_it++;
			prop._property_class_name = *tok_it;
		}

		if (tok_it->compare("PROPERTY_CLASS") == 0) {
			tok_it++;
			if (! prop.checkID(*tok_it)) {
				throw std::runtime_error("parseGocadPropertyMetaData(): id mismatch.");
			}
			tok_it++;
			prop._property_class_name = *tok_it;
		}

		if (tok_it->compare("PROPERTY_SUBCLASS") == 0) {
			tok_it++;
			if (! prop.checkID(*tok_it)) {
				throw std::runtime_error("parseGocadPropertyMetaData(): id mismatch.");
			}
			tok_it++;
			if (tok_it->compare("QUANTITY") != 0) {
				ERR("Expected keyword QUANTITY but found \"%s\".", tok_it->c_str());
				throw std::runtime_error("parseGocadPropertyMetaData(): id mismatch.");
			}
			tok_it++;
			prop._property_data_type = *tok_it;
		}

		if (tok_it->compare("PROP_UNIT") == 0 ||
				tok_it->compare("PROP_ORIGINAL_UNIT") == 0) {
			tok_it++;
			if (! prop.checkID(*tok_it)) {
				throw std::runtime_error("parseGocadPropertyMetaData(): id mismatch.");
			}
			tok_it++;
			prop._property_unit = *tok_it;
		}

		if (tok_it->compare("PROP_NO_DATA_VALUE") == 0) {
			tok_it++;
			if (! prop.checkID(*tok_it)) {
				throw std::runtime_error("parseGocadPropertyMetaData(): id mismatch.");
			}
			tok_it++;
			prop._property_no_data_value = BaseLib::str2number<double>(*tok_it);
		}
	}
	return prop;
}

GocadSGridReader::GocadSGridReader(std::string const& fname) :
		_fname(fname), _path(_fname.substr(0, _fname.find_last_of("/\\") + 1)),
		_n_face_sets(0), _double_precision_binary(false), _bin_pnts_in_double_precision(false),
		_center_node(0.0, 0.0, 0.0)

{
	// check if file exists
	std::ifstream in(_fname.c_str());
	if (!in) {
		ERR("Could not open \"%s\".", _fname.c_str());
		in.close();
		return;
	}

	bool pnts_read(false);

	// read information in the stratigraphic grid file
	std::string line;
	while (std::getline(in, line)) {
		if (line.compare(0, 8, "HEADER {") == 0)
		{
			parseHeader(in);
		}
		if (line.compare(0, 7, "AXIS_N ") == 0)
		{
			parseDims(line);
		}
		else if (line.compare(0, 12, "POINTS_FILE ") == 0)
		{
			parseFileName(line, _pnts_fname);
		}
		else if (line.compare(0, 9, "PROPERTY ") == 0)
		{
			_property_meta_data_vecs.push_back(parseGocadPropertyMetaData(line, in, _path));
		}
		else if (line.compare(0, 35, "BINARY_POINTS_IN_DOUBLE_PRECISION 1") == 0) {
			_bin_pnts_in_double_precision = true;
		}
		else if (line.compare(0, 11, "FLAGS_FILE ") == 0) {
			parseFileName(line, _flags_fname);
		}
		else if (line.compare(0, 18, "REGION_FLAGS_FILE ") == 0)
		{
			parseFileName(line, _region_flags_fname);
		}
		else if (line.compare(0, 7, "REGION ") == 0 || line.compare(0, 13, "MODEL_REGION ") == 0)
		{
			regions.push_back(parseRegion(line));
		}
		else if (line.compare(0, 12, "MODEL_LAYER ") == 0)
		{
			layers.push_back(parseLayer(line, regions));
		}
		else if (line.compare(0, 24, "REGION_FLAGS_BIT_LENGTH ") == 0)
		{
			std::istringstream iss(line);
			std::istream_iterator<std::string> it(iss);
			it++;
			std::size_t bit_length = atoi(it->c_str());
			if (regions.size() != bit_length)
			{
				ERR("%d regions read but %d expected.\n", regions.size(), bit_length);
				throw std::runtime_error("Number of read regions differs from expected.\n");
			}
		}
		else if (line.compare(0, 9, "FACE_SET ") == 0) {
			// first read the points
			if (! pnts_read) {
				readNodesBinary();
				pnts_read = true;
			}
			parseFaceSet(line, in);
		}
		else
		{
			//std::cout << "Skip: \"" << line << "\"\n";
		}
	}
	std::cout << regions.size() << " regions read:\n";
	std::copy(regions.begin(), regions.end(), std::ostream_iterator<Region>(std::cout, "\t"));
	std::cout << "\n";
	std::cout << layers.size() << " layers read:\n";
	std::copy(layers.begin(), layers.end(), std::ostream_iterator<Layer>(std::cout, "\n"));

	std::cout << "meta data for " << _property_meta_data_vecs.size() << " properties read:\n";
	std::copy(_property_meta_data_vecs.begin(), _property_meta_data_vecs.end(),
			std::ostream_iterator<GocadProperty>(std::cout, "\n"));

	// if not done already read the points
	if (! pnts_read) {
		readNodesBinary();
		pnts_read = true;
	}
	readElementPropertiesBinary();
	std::vector<Bitset> region_flags = readRegionFlagsBinary();
	//mapRegionFlagsToCellProperties(region_flags);	// modifies _material_ids.

	readSplitInformation();

	GeoLib::AABB<MeshLib::Node> aabb(_nodes.begin(), _nodes.end());
	_center_node[0] = (aabb.getMaxPoint()[0] + aabb.getMinPoint()[0])/2.0;
	_center_node[1] = (aabb.getMaxPoint()[1] + aabb.getMinPoint()[1])/2.0;
	_center_node[2] = 0.0;

	GocadProperty face_set_property;
	face_set_property._property_id = 0;
	face_set_property._property_name = "CellIDs";
	face_set_property._property_class_name = "CellIDsData";
	face_set_property._property_unit = "unitless";
	face_set_property._property_data_type = "double";
	face_set_property._property_data_fname = "";
	face_set_property._property_no_data_value = -1.0;
	face_set_property._property_data.resize(_index_calculator._n_cells);
	std::iota(face_set_property._property_data.begin(), face_set_property._property_data.end(),
			0.0);

	_property_meta_data_vecs.push_back(face_set_property);
	in.close();
}

GocadSGridReader::~GocadSGridReader()
{
	for (std::size_t k(0); k<_nodes.size(); k++)
		delete _nodes[k];
	for (std::size_t k(0); k<_split_nodes.size(); k++)
		delete _split_nodes[k];
}

MeshLib::Mesh* GocadSGridReader::getMesh() const
{
	std::vector<MeshLib::Node*> nodes;
	std::transform(_nodes.cbegin(), _nodes.cend(),
		            std::back_inserter(nodes),
		            [](MeshLib::Node const*const node) {
						return new MeshLib::Node(*node);
		            });

	std::vector<MeshLib::Element*> elements;
	createElements(nodes, elements);
	applySplitInformation(nodes, elements);

	MeshLib::Node const& center(_center_node);
	INFO("translated model (-%f, -%f, -%f).", center[0], center[1], center[2]);
	std::for_each(nodes.begin(), nodes.end(),
			[&center](MeshLib::Node* node)
			{
				(*node)[0] -= center[0];
				(*node)[1] -= center[1];
			}
	);

	INFO("Creating mesh from Gocad SGrid.");
	MeshLib::Mesh *mesh (new MeshLib::Mesh("GocadSGrid", nodes, elements));
	INFO("Mesh created.");
	return mesh;
}

void GocadSGridReader::parseHeader(std::istream &in)
{
	std::string line;
	while (std::getline(in, line)) {
		if (line.compare(0, 1, "}") == 0) {
			return;
		}
		if (line.compare(0, 27, "double_precision_binary: on") == 0) {
			_double_precision_binary = true;
		}
	}
}

void GocadSGridReader::parseDims(std::string const& line)
{
	std::size_t x_dim(0), y_dim(0), z_dim(0);
	boost::tokenizer<> tok(line);
	auto it(tok.begin());
	it++; // overread token "AXIS"
	it++; // overread "N"
	std::stringstream ssx(*(it), std::stringstream::in | std::stringstream::out);
	ssx >> x_dim;
	it++;
	std::stringstream ssy(*it, std::stringstream::in | std::stringstream::out);
	ssy >> y_dim;
	it++;
	std::stringstream ssz(*it, std::stringstream::in | std::stringstream::out);
	ssz >> z_dim;
	_index_calculator = IndexCalculator(x_dim, y_dim, z_dim);
}

void GocadSGridReader::parseFileName(std::string const& line, std::string &result_string) const
{
	std::size_t beg_pos(line.find_first_of(" ") + 1);
	std::string fname(line.substr(beg_pos, line.length() - beg_pos));
	BaseLib::trim(fname, '\"');
	result_string = _path + fname;
}

/**
 * @param line input/output
 * @param in input stream containing the face set
 */
void GocadSGridReader::parseFaceSet(std::string &line, std::istream &in)
{
	// create and initialize a GocadProperty object for storing face set data
	GocadProperty face_set_property;
	face_set_property._property_id = _n_face_sets;
	face_set_property._property_name = "FaceSet";
	face_set_property._property_class_name = "FaceSetData";
	face_set_property._property_unit = "unitless";
	face_set_property._property_data_type = "double";
	face_set_property._property_data_fname = "";
	face_set_property._property_no_data_value = -1.0;
	face_set_property._property_data.resize(_index_calculator._n_cells);
	std::fill(face_set_property._property_data.begin(), face_set_property._property_data.end(),
			face_set_property._property_no_data_value);

	std::istringstream iss(line);
	std::istream_iterator<std::string> it(iss);
	// Check first word is FACE_SET
	if (*it != std::string("FACE_SET")) {
		ERR("Expected FACE_SET keyword but \"%s\" found.", it->c_str());
		throw std::runtime_error("In GocadSGridReader::parseFaceSet() expected FACE_SET keyword not found.");
	}
	++it;
	face_set_property._property_name += *it;
	++it;
	std::size_t const n_of_face_set_ids(static_cast<std::size_t>(atoi(it->c_str())));
	std::size_t face_set_id_cnt(0);

	while (getline(in, line) && face_set_id_cnt < n_of_face_set_ids) {
		boost::char_separator<char> sep("\t ");
		boost::tokenizer<boost::char_separator<char> > tokens(line, sep);

		for(auto tok_it  = tokens.begin(); tok_it != tokens.end(); ) {
			std::size_t id(static_cast<std::size_t>(atoi(tok_it->c_str())));
			tok_it++;
			std::size_t face_indicator(static_cast<std::size_t>(atoi(tok_it->c_str())));
			tok_it++;

			if (id >= _index_calculator._n_nodes) {
				ERR("Face set id %d is greater than the number of nodes (%d).", id, _index_calculator._n_nodes);
			} else {
				dynamic_cast<MeshLib::GocadNode*>(_nodes[id])->setFaceSet(_n_face_sets, face_indicator);
				std::array<std::size_t,3> c(_index_calculator.getCoordsForID(id));
				if (c[0] >= _index_calculator._x_dim-1)
					ERR("****** i coord %d to big for id %d.", c[0], id);
				if (c[1] >= _index_calculator._y_dim-1)
					ERR("****** j coord %d to big for id %d.", c[1], id);
				if (c[2] >= _index_calculator._z_dim-1)
					ERR("****** k coord %d to big for id %d.", c[2], id);
				std::size_t const cell_id( _index_calculator.getCellIdx(c[0], c[1], c[2]));
				face_set_property._property_data[cell_id] = face_indicator;
			}
			face_set_id_cnt++;
		}
	}

	if (face_set_id_cnt != n_of_face_set_ids) {
		ERR("Expected %d number of face set ids, read %d.", n_of_face_set_ids, face_set_id_cnt);
		throw std::runtime_error("Expected number of face set points does not match number of read points.");
	}
	_n_face_sets++;

	// pre condition: split nodes are read already
	for (std::size_t k(0); k<_split_nodes.size(); k++) {
		std::size_t const id (_index_calculator(_split_nodes[k]->getGridCoords()));
		_split_nodes[k]->transmitFaceIndicators(*_nodes[id]);
	}

	_property_meta_data_vecs.push_back(face_set_property);
}

// Reads given number of bits (rounded up to next byte) into a bitset.
// Used for reading region information which can be represented by some
// number of bits.
Bitset readBits(std::ifstream& in, const std::size_t bits)
{
	typedef Bitset::block_type block_t;
	std::size_t const bytes = static_cast<std::size_t>(std::ceil(bits/8.));
	std::size_t const blocks = (bytes + 1)/ sizeof(block_t);

	block_t data[blocks];
	std::fill_n(data, blocks, 0);
	in.read(reinterpret_cast<char*>(data), bytes);

	return Bitset(data, data + blocks);
}

void GocadSGridReader::readNodesBinary()
{
	std::ifstream in(_pnts_fname.c_str(), std::ios::binary);
	if (!in) {
		ERR("Could not open points file \"%s\".", _pnts_fname.c_str());
		throw std::runtime_error("Could not open points file.");
	}

	std::size_t const n = _index_calculator._n_nodes;
	_nodes.resize(n);

	double coords[3];

	std::size_t k = 0;
	while (in && k < n * 3)
	{
		if (_bin_pnts_in_double_precision) {
			coords[k % 3] = BaseLib::swapEndianness(BaseLib::readBinaryValue<double>(in));
		} else {
			coords[k % 3] = BaseLib::swapEndianness(BaseLib::readBinaryValue<float>(in));
		}
		if ((k + 1) % 3 == 0) {
			const std::size_t layer_transition_idx(_index_calculator.getCoordsForID(k/3)[2]);
			_nodes[k/3] = new MeshLib::GocadNode(coords, k/3, layer_transition_idx);
		}
		k++;
	}
	if (k != n * 3 && !in.eof())
		ERR("Read different number of points. Expected %d floats, got %d.\n", n * 3, k);
}

void GocadSGridReader::mapRegionFlagsToCellProperties(std::vector<Bitset> const& rf)
{
//	std::size_t const n = _index_calculator._n_cells;
//	_material_ids.resize(n);
//	std::fill(_material_ids.begin(), _material_ids.end(), -1);
//	// region flags are stored in each node ijk and give the region index for the
//	// ijk-th cell.
//	for (std::size_t k(0); k < _index_calculator._z_dim-1; k++) {
//		for (std::size_t j(0); j < _index_calculator._y_dim-1; j++) {
//			for (std::size_t i(0); i < _index_calculator._x_dim-1; i++) {
//				// Find layers containing regions given by bits.
//				// Run over bits and push back set bits
//				std::set<std::size_t> layers_set;
//				for (auto r = regions.begin(); r != regions.end(); ++r)
//				{
//					if (rf[_index_calculator(i,j,k)].test(r->bit))
//					{
//						// Bit is set, find a layer.
//						for (std::size_t l_id = 0; l_id < layers.size(); ++l_id)
//							if (layers[l_id].hasRegion(*r))
//								layers_set.insert(l_id);
//					}
//				}
//				if (layers_set.size() != 1)
//					ERR("A cell %d %d %d belongs to multiple (%d) layers.", i, j, k, layers_set.size());
//
//				_material_ids[_index_calculator.getCellIdx(i,j,k)] = *layers_set.begin();
//			}
//		}
//	}
}

void GocadSGridReader::readElementPropertiesBinary()
{
	for (auto prop_it(_property_meta_data_vecs.begin());
			prop_it != _property_meta_data_vecs.end();
			prop_it++) {
		std::string const& fname(prop_it->_property_data_fname);
		if (prop_it->_property_data_fname.empty()) {
			WARN("Empty filename for property %s.", prop_it->_property_name.c_str());
			continue;
		}
		std::vector<float> float_properties =
				BaseLib::readBinaryArray<float>(fname, _index_calculator._n_cells);
		std::transform(float_properties.cbegin(), float_properties.cend(), float_properties.begin(),
				[](float const& val) {
					return BaseLib::swapEndianness(val);
				});

		std::vector<double> properties;
		prop_it->_property_data.resize(float_properties.size());
		std::copy(float_properties.begin(), float_properties.end(), prop_it->_property_data.begin());
		if (prop_it->_property_data.empty()) {
			ERR("Reading of element properties file \"%s\" failed.", fname.c_str());
		}
	}
}

std::vector<int> GocadSGridReader::readFlagsBinary() const
{
	std::vector<int> result;
	if (!_double_precision_binary) {
		result = BaseLib::readBinaryArray<int32_t>(_flags_fname, _index_calculator._n_nodes);
		std::for_each(result.begin(), result.end(),
						[](int32_t& val) {
							BaseLib::swapEndianness(val);
						});
	} else {
		result = BaseLib::readBinaryArray<int>(_flags_fname, _index_calculator._n_nodes);
		std::for_each(result.begin(), result.end(),
						[](int& val) {
							BaseLib::swapEndianness(val);
						});
	}

	if (result.empty())
		ERR("Reading of flags file \"%s\" failed.", _flags_fname.c_str());

	return result;
}

std::vector<Bitset> GocadSGridReader::readRegionFlagsBinary() const
{
	std::vector<Bitset> result;

	std::ifstream in(_region_flags_fname.c_str());
	if (!in) {
		ERR("readRegionFlagsBinary(): Could not open file \"%s\" for input.\n", _region_flags_fname.c_str());
		in.close();
		return result;
	}

	std::size_t const n = _index_calculator._n_nodes;
	result.resize(n);

	std::size_t k = 0;
	while (in && k < n)
	{
		result[k++] = readBits(in, regions.size());
	}
	if (k != n && !in.eof())
		ERR("Read different number of values. Expected %d, got %d.\n", n, k);

	return result;
}

void GocadSGridReader::createElements(std::vector<MeshLib::Node*> const& nodes,
		std::vector<MeshLib::Element*> & elements) const
{
	elements.resize(_index_calculator._n_cells);
	std::array<MeshLib::Node*, 8> element_nodes;
	std::size_t cnt(0);
	for (std::size_t k(0); k < _index_calculator._z_dim-1; k++) {
		for (std::size_t j(0); j < _index_calculator._y_dim-1; j++) {
			for (std::size_t i(0); i < _index_calculator._x_dim-1; i++) {
				element_nodes[0] = nodes[_index_calculator(i,j,k)];
				element_nodes[1] = nodes[_index_calculator(i+1,j,k)];
				element_nodes[2] = nodes[_index_calculator(i+1,j+1,k)];
				element_nodes[3] = nodes[_index_calculator(i,j+1,k)];
				element_nodes[4] = nodes[_index_calculator(i,j,k+1)];
				element_nodes[5] = nodes[_index_calculator(i+1,j,k+1)];
				element_nodes[6] = nodes[_index_calculator(i+1,j+1,k+1)];
				element_nodes[7] = nodes[_index_calculator(i,j+1,k+1)];
				elements[cnt] = new MeshLib::Hex(element_nodes, _index_calculator.getCellIdx(i,j,k));
				cnt++;
			}
		}
	}
}

void GocadSGridReader::readSplitInformation()
{
	std::ifstream in(_fname.c_str());
	if (!in) {
		ERR("Could not open \"%s\".", _fname.c_str());
		in.close();
		return;
	}

	// read split information from the stratigraphic grid file
	std::string line;
	std::stringstream ss;
	while (std::getline(in, line)) {
		std::size_t pos(line.find("SPLIT "));
		if (pos != std::string::npos) {
			ss << line.substr(pos+6, line.size()-(pos+6));
			// read position in grid
			std::array<std::size_t, 3> grid_coords;
			ss >> grid_coords[0];
			ss >> grid_coords[1];
			ss >> grid_coords[2];
			// read coordinates for the split node
			double coords[3];
			ss >> coords[0];
			ss >> coords[1];
			ss >> coords[2];
			// read the id
			std::size_t id;
			ss >> id;
			// read the affected cells
			std::array<bool, 8> affected_cells;
			for (std::size_t k(0); k<affected_cells.size(); k++) {
				ss >> affected_cells[k];
			}
			const std::size_t layer_transition_index(
				_nodes[id]->getLayerTransitionIndex());
			_split_nodes.push_back(
				new MeshLib::GocadSplitNode(coords, id, grid_coords,
					affected_cells, layer_transition_index));
		}
	}
}

void GocadSGridReader::applySplitInformation(std::vector<MeshLib::Node*> &nodes,
		std::vector<MeshLib::Element*> &elements) const
{
	for (std::size_t k(0); k<_split_nodes.size(); k++) {
		std::size_t const new_node_pos(nodes.size());
		nodes.push_back(new MeshLib::Node(_split_nodes[k]->getCoords(), new_node_pos));

		// get grid coordinates
		std::array<std::size_t, 3> const& gc(_split_nodes[k]->getGridCoords());
		// get affected cells
		std::array<bool, 8> const& affected_cells(_split_nodes[k]->getAffectedCells());
		// get mesh node to substitute in elements
		MeshLib::Node const*const node2sub(nodes[_index_calculator(gc[0],gc[1],gc[2])]);

		if (affected_cells[0] && gc[0] < _index_calculator._x_dim-1
		                   && gc[1] < _index_calculator._y_dim-1
		                   && gc[2] < _index_calculator._z_dim-1) {
			const std::size_t idx(_index_calculator.getCellIdx(gc[0], gc[1], gc[2]));
			modifyElement(elements[idx], node2sub, nodes[new_node_pos]);
		}
		if (affected_cells[1] && gc[0] > 0
				&& gc[1] < _index_calculator._y_dim-1
				&& gc[2] < _index_calculator._z_dim-1) {
			const std::size_t idx(_index_calculator.getCellIdx(gc[0]-1, gc[1], gc[2]));
			modifyElement(elements[idx], node2sub, nodes[new_node_pos]);
		}
		if (affected_cells[2] && gc[1] > 0
				&& gc[0] < _index_calculator._x_dim-1
				&& gc[2] < _index_calculator._z_dim-1) {
			const std::size_t idx(_index_calculator.getCellIdx(gc[0], gc[1]-1, gc[2]));
			modifyElement(elements[idx], node2sub, nodes[new_node_pos]);
		}
		if (affected_cells[3] && gc[0] > 0 && gc[1] > 0
				&& gc[2] < _index_calculator._z_dim-1) {
			const std::size_t idx(_index_calculator.getCellIdx(gc[0]-1, gc[1]-1, gc[2]));
			modifyElement(elements[idx], node2sub, nodes[new_node_pos]);
		}
		if (affected_cells[4] && gc[2] > 0
				&& gc[0] < _index_calculator._x_dim-1
                && gc[1] < _index_calculator._y_dim-1) {
			const std::size_t idx(_index_calculator.getCellIdx(gc[0], gc[1], gc[2]-1));
			modifyElement(elements[idx], node2sub, nodes[new_node_pos]);
		}
		if (affected_cells[5] && gc[0] > 0 && gc[2] > 0
				&& gc[1] < _index_calculator._y_dim-1) {
			const std::size_t idx(_index_calculator.getCellIdx(gc[0]-1, gc[1], gc[2]-1));
			modifyElement(elements[idx], node2sub, nodes[new_node_pos]);
		}
		if (affected_cells[6] && gc[1] > 0 && gc[2] > 0
				&& gc[0] < _index_calculator._x_dim-1) {
			const std::size_t idx(_index_calculator.getCellIdx(gc[0], gc[1]-1, gc[2]-1));
			modifyElement(elements[idx], node2sub, nodes[new_node_pos]);
		}
		if (affected_cells[7] && gc[0] > 0 && gc[1] > 0 && gc[2] > 0) {
			const std::size_t idx(_index_calculator.getCellIdx(gc[0]-1, gc[1]-1, gc[2]-1));
			modifyElement(elements[idx], node2sub, nodes[new_node_pos]);
		}
	}
}

void GocadSGridReader::modifyElement(MeshLib::Element* hex,
		MeshLib::Node const* node2sub, MeshLib::Node * substitute_node) const
{
	// get the node pointers of the cell
	MeshLib::Node *const* hex_nodes(hex->getNodes());
	// search for the position the split node will be set to
	MeshLib::Node *const* node_pos(std::find(hex_nodes, hex_nodes+8, node2sub));
	// set the split node instead of the node2sub
	if (node_pos != hex_nodes+8) {
		const_cast<MeshLib::Node**>(hex_nodes)[std::distance(hex_nodes, node_pos)] = substitute_node;
	}
}

boost::optional<GocadProperty const&>
GocadSGridReader::getProperty(std::string const& name) const
{
	auto const it(std::find_if(_property_meta_data_vecs.begin(), _property_meta_data_vecs.end(),
			[&name](GocadProperty const& p) { return p._property_name.compare(name) == 0; }));
	if (it != _property_meta_data_vecs.end())
		return boost::optional<GocadProperty const&>(*it);
	else
		return boost::optional<GocadProperty const&>();
}

std::vector<std::string> GocadSGridReader::getPropertyNames() const
{
	std::vector<std::string> names;
	std::transform(_property_meta_data_vecs.begin(),
			_property_meta_data_vecs.end(),
			std::back_inserter(names),
			[](GocadProperty const& p) {
				return p._property_name;
			}
		);
	return names;
}

MeshLib::Mesh* GocadSGridReader::getFaceSetMesh(std::size_t face_set_number) const
{
	std::vector<MeshLib::Node*> face_set_nodes;
	std::vector<MeshLib::Element*> face_set_elements;

	for (std::size_t k(0); k<_nodes.size(); k++) {
		if (_nodes[k]->isMemberOfFaceSet(face_set_number)) {
			addFaceSetQuad(_nodes[k], face_set_number, face_set_nodes, face_set_elements);
		}
	}

	if (face_set_nodes.empty())
		return nullptr;

	for (std::size_t k(0); k<_split_nodes.size(); k++) {
		if (_split_nodes[k]->isMemberOfFaceSet(face_set_number)) {
			if (_split_nodes[k]->getAffectedCells()[0]) {
				addFaceSetQuad(_split_nodes[k], face_set_number, face_set_nodes, face_set_elements);
			}
		}
	}

	INFO("translated model (-%f, -%f, -%f).", _center_node[0], _center_node[1], _center_node[2]);
	MeshLib::Node const& center(_center_node);
	std::for_each(face_set_nodes.begin(), face_set_nodes.end(), [&center](MeshLib::Node* node)
	{
		(*node)[0] -= center[0];
		(*node)[1] -= center[1];
	});

	std::string mesh_name("GocadFaceSetMesh-" + std::to_string(face_set_number));
	return new MeshLib::Mesh(mesh_name, face_set_nodes, face_set_elements);
}

void GocadSGridReader::addFaceSetQuad(MeshLib::GocadNode* face_set_node, std::size_t face_set_number,
		std::vector<MeshLib::Node*> &face_set_nodes,
		std::vector<MeshLib::Element*> &face_set_elements) const
{
	std::array<MeshLib::Node*, 4> quad_nodes;
	quad_nodes[0] = new MeshLib::GocadNode(*face_set_node);
	const std::size_t id(face_set_node->getID());
	std::array<std::size_t, 3> c(_index_calculator.getCoordsForID(id));

	const MeshLib::FaceIndicator dir(face_set_node->getFaceIndicator(face_set_number));
	switch (dir) {
	case MeshLib::FaceIndicator::U:
		quad_nodes[1] = new MeshLib::GocadNode(*_nodes[_index_calculator(c[0], c[1] + 1, c[2])]);
		quad_nodes[2] = new MeshLib::GocadNode(*_nodes[_index_calculator(c[0], c[1] + 1, c[2] + 1)]);
		quad_nodes[3] = new MeshLib::GocadNode(*_nodes[_index_calculator(c[0], c[1], c[2] + 1)]);
		break;
	case MeshLib::FaceIndicator::V:
		quad_nodes[1] = new MeshLib::GocadNode(*_nodes[_index_calculator(c[0] + 1, c[1], c[2])]);
		quad_nodes[2] = new MeshLib::GocadNode(*_nodes[_index_calculator(c[0] + 1, c[1], c[2] + 1)]);
		quad_nodes[3] = new MeshLib::GocadNode(*_nodes[_index_calculator(c[0], c[1], c[2] + 1)]);
		break;
	case MeshLib::FaceIndicator::W:
		quad_nodes[1] = new MeshLib::GocadNode(*_nodes[_index_calculator(c[0] + 1, c[1], c[2])]);
		quad_nodes[2] = new MeshLib::GocadNode(*_nodes[_index_calculator(c[0] + 1, c[1] + 1, c[2])]);
		quad_nodes[3] = new MeshLib::GocadNode(*_nodes[_index_calculator(c[0], c[1] + 1, c[2])]);
		break;
	default:
		ERR("Could not create face for node with id %d.", id);
	}
	for (auto it = quad_nodes.cbegin(); it != quad_nodes.cend(); it++)
		face_set_nodes.push_back(*it);
	face_set_elements.push_back(new MeshLib::Quad(quad_nodes));
}

} // end namespace FileIO
