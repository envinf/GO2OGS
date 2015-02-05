/**
 * This file is part of "GocadSGridReader". "GocadSGridReader" is free
 * software: you can redistribute it and/or modify it under the terms of
 * the GNU General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * "GocadSGridReader" is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License with
 * "GocadSGridReader". If not, see <http://www.gnu.org/licenses/>.
 *
 * @date 2013-03-01
 * @author Thomas Fischer
 * @file main.cpp
 */

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "tclap/CmdLine.h"
#include "FileTools.h"
#include "quicksort.h"
#include "LogogSimpleFormatter.h"

// FileIO
#include "FileIO/XmlIO/Boost/BoostVtuInterface.h"

// MeshLib
#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"
#include "Elements/Quad.h"
#include "MeshSurfaceExtraction.h"

// Utils/FileConverter
#include "GocadSGridReader.h"

void regenerateFaceSetMesh(MeshLib::Mesh const& mesh,
	std::size_t face_set_number,
	std::string const& path)
{
	std::vector<MeshLib::GocadNode*> gocad_nodes;
	gocad_nodes.resize(mesh.getNNodes());
	for (std::size_t k(0); k < mesh.getNNodes(); k++) {
		gocad_nodes[k] = new MeshLib::GocadNode(*
			static_cast<MeshLib::GocadNode const*>(mesh.getNode(k)));
	}

	std::vector<std::size_t> perm;
	perm.resize(mesh.getNNodes());
	std::iota(perm.begin(), perm.end(), 0);
	BaseLib::Quicksort<MeshLib::GocadNode*>(gocad_nodes, 0, gocad_nodes.size(), perm);

	auto last = std::unique(gocad_nodes.begin(), gocad_nodes.end(),
		[](MeshLib::GocadNode* n0, MeshLib::GocadNode* n1) {
			if ((*n0) <= (*n1) && (*n1) <= (*n0))
				return true;
		}
	);
	std::size_t const n(std::distance(gocad_nodes.begin(), last));
	std::size_t const n_cols(n/14 - 1);
	for (std::size_t k(0); k<n_cols; ++k) {
		// bubble sort column according to layer_transition_idx
		for (std::size_t i(k*14); i<(k+1)*14; ++i) {
			std::size_t i_idx(gocad_nodes[i]->getLayerTransitionIndex());
			for (std::size_t j(i); j<(k+1)*14; ++j) {
				if (gocad_nodes[j]->getLayerTransitionIndex() < i_idx) {
					std::swap(gocad_nodes[i], gocad_nodes[j]);
					i_idx = gocad_nodes[i]->getLayerTransitionIndex();
				}
			}
		}
	}


	std::string sorted_nodes_fname(path+"Surfaces/SortedNodes-" +
		std::to_string(face_set_number)+".gli");
	std::ofstream os(sorted_nodes_fname);
	os << "#POINTS\n";
	for (std::size_t k(0); k<n; k++) {
		os << k << " " << *(gocad_nodes[k]) << "$NAME " << gocad_nodes[k]->getLayerTransitionIndex() << "\n";
	}
	os << "#STOP";
	os.close();

	std::vector<MeshLib::Node*> nodes;
	for (std::size_t k(0); k<n; k++) {
		nodes.push_back(new MeshLib::Node(*(gocad_nodes[k])));
	}

	std::vector<MeshLib::Element*> elements;
	std::vector<unsigned> region_properties;
	// generate quad elements
	for (std::size_t r(0); r<n_cols; ++r) {
		for (std::size_t c(0); c<13; c++) {
			std::array<MeshLib::Node*,4> quad_nodes;
			quad_nodes[0] = nodes[14*r+c];
			quad_nodes[1] = nodes[14*r+c+1];
			quad_nodes[2] = nodes[14*(r+1)+c+1];
			quad_nodes[3] = nodes[14*(r+1)+c];
			elements.push_back(new MeshLib::Quad(quad_nodes, c, 14*r+c));
			region_properties.push_back(c);
		}
	}
	MeshLib::Mesh new_mesh(mesh.getName(), nodes, elements);
	new_mesh.addPropertyVec("Regions", region_properties);

	FileIO::BoostVtuInterface vtu;
	vtu.setMesh(&new_mesh);
	// output file name
	std::string mesh_out_fname(path+"Surfaces/RegeneratedFaceSetMesh-"
		+ std::to_string(face_set_number) + ".vtu");
	INFO("Writing face set mesh \"%s\" in vtu format.", mesh_out_fname.c_str());
	vtu.writeToFile(mesh_out_fname);

	std::for_each(gocad_nodes.begin(), last,
		std::default_delete<MeshLib::GocadNode>());
}

void writeFaceSetNodesAsGLI(MeshLib::Mesh const& mesh,
	std::size_t face_set_number,
	std::string const& path)
{
	std::stringstream ss;
	ss << "#POINTS\n";
	for (std::size_t k(0); k < mesh.getNNodes(); k++) {
		MeshLib::GocadNode const& node(
			*static_cast<MeshLib::GocadNode const*>(mesh.getNode(k)));
		ss << k << " " << node << "$NAME " << node.getLayerTransitionIndex() << "\n";
	}
	ss << "#STOP";

	std::string fname(path + "Surfaces/FaceSetNodes-" + std::to_string(face_set_number) + ".gli");
	INFO("Writing nodes of face set to file \"%s\".", fname.c_str());
	std::ofstream os(fname.c_str());
	os << ss.str();
	os.close();
}


void writeFaceSetNodesAsCSV(MeshLib::Mesh const& mesh, std::size_t face_set_number, std::string const& path)
{
	std::stringstream ss;
	std::size_t cnt(0); // count face set nodes
	for (std::size_t k(0); k < mesh.getNNodes(); k++) {
		MeshLib::GocadNode* gocad_node(
				dynamic_cast<MeshLib::GocadNode*>(const_cast<MeshLib::Node*>(mesh.getNode(k))));

		bool const face_set_member(gocad_node->isMemberOfFaceSet(face_set_number));
		if (face_set_member) {
			ss << (*gocad_node)[0] << "," << (*gocad_node)[1] << "," << (*gocad_node)[2] << ",";

			switch (gocad_node->getFaceIndicator(face_set_number))
			{
			case MeshLib::FaceIndicator::U:
				ss << "0";
				break;
			case MeshLib::FaceIndicator::V:
				ss << "1";
				break;
			case MeshLib::FaceIndicator::W:
				ss << "2";
				break;
			default:
				ss << "unknown";
			}
			ss << "\n";
			cnt++;
		}
	}
	if (cnt > 0) {
		std::string fname(path + "Surfaces/FaceSetNodes-" + std::to_string(face_set_number) + ".csv");
		INFO("Writing nodes of face set to file \"%s\".", fname.c_str());
		std::ofstream os(fname.c_str());
		os << ss.str();
		os.close();
	}
}

void generateFaceSetMeshes(FileIO::GocadSGridReader const& reader, std::string const& path)
{
	for (std::size_t l(0); l<128; l++) {
		MeshLib::Mesh *face_set_mesh(reader.getFaceSetMesh(l));

		if (face_set_mesh == nullptr)
			continue;
		INFO("Creating face set mesh.");
		INFO("Face set mesh created. #nodes: %d, #elements: %d", face_set_mesh->getNNodes(),
				face_set_mesh->getNElements());


		FileIO::BoostVtuInterface vtu;
		vtu.setMesh(face_set_mesh);
		// output file name
		std::string mesh_out_fname(path+"Surfaces/FaceSetMesh-" + std::to_string(l) + ".vtu");
		INFO("Writing face set mesh \"%s\" in vtu format.", mesh_out_fname.c_str());
		vtu.writeToFile(mesh_out_fname);
		writeFaceSetNodesAsGLI(*face_set_mesh, l, path);

		regenerateFaceSetMesh(*face_set_mesh, l, path);
		delete face_set_mesh;
	}
}

void cleanUpNoDataValues(MeshLib::Mesh &mesh, double no_data_value,
		std::vector<double> const& original,
		std::vector<double> & reworked)
{
	std::vector<MeshLib::Element*> const& elements(mesh.getElements());
	for (std::size_t k(0); k<original.size(); ++k) {
		if (std::abs(original[k] - no_data_value) > std::numeric_limits<double>::epsilon()) {
			reworked[k] = original[k];
		} else {
			std::size_t const n_neighbors(elements[k]->getNNeighbors());
			double prop_value(0.0);
			std::size_t cnt(0); // count neighbors with valid property values
			// simple average
			for (std::size_t j(0); j<n_neighbors; ++j) {
				MeshLib::Element const*const j_th_neighbor(elements[k]->getNeighbor(j));
				if (j_th_neighbor != nullptr) {
					double neighbor_val(original[j_th_neighbor->getValue()]);
					if (std::abs(neighbor_val - no_data_value) > std::numeric_limits<double>::epsilon()) {
						prop_value += neighbor_val;
						cnt++;
					}
				}
			}
			prop_value /= cnt;
			reworked[k] = prop_value;
		}
	}
}

void addGocadPropertiesToMesh(FileIO::GocadSGridReader const& reader, MeshLib::Mesh &mesh)
{
	std::vector<std::string> const& prop_names(reader.getPropertyNames());
	for (auto name_it(prop_names.begin()); name_it != prop_names.end(); name_it++) {
		boost::optional<FileIO::GocadSGridReader::GocadProperty const&> prop(reader.getProperty(*name_it));
		if (prop) {
			INFO("Adding property \"%s\".", name_it->c_str());
			std::vector<double> reworked_properties;
			reworked_properties.resize((*prop)._property_data.size());
//			cleanUpNoDataValues(mesh, (*prop)._property_no_data_value, (*prop)._property_data, reworked_properties);
//			mesh.addPropertyVec(*name_it, reworked_properties);
			mesh.addPropertyVec(*name_it, (*prop)._property_data);
		}
	}
}

void writeMeshPropertiesToFile(std::string const& fname, std::string const& prop_name,
		MeshLib::Mesh const& mesh)
{
	boost::optional<std::vector<double> const&> prop_vec(mesh.getDoublePropertyVec(prop_name));
	if (!prop_vec) {
		ERR("Could not read property \"%s\" from mesh.", prop_name.c_str());
		return;
	}

	std::ofstream os(fname.c_str());
	if (!os) {
		ERR("Could not open file \"%s\".", fname.c_str());
		return;
	} else {
		INFO("Open file \"%s\" for writing.", fname.c_str());
	}

	// write header
	os << "; OpenGeoSys material property file - " << prop_name << " \n";
	os << "#MEDIUM_PROPERTIES_DISTRIBUTED\n";
	os << "$MSH_TYPE\n";
	os << " NO_PCS\n";
	os << "$MMP_TYPE\n";
	if (prop_name.compare("Porosity") == 0) {
		os << " POROSITY\n";
	}
	if (prop_name.compare("Permeability") == 0) {
		os << " PERMEABILITY\n";
	}
	os << "$DIS_TYPE\n";
	os << " ELEMENT\n";
	os << "$DATA\n";

	if (prop_name.compare("Porosity") == 0) {
		for (std::size_t k(0); k < (*prop_vec).size(); k++) {
			if ((*prop_vec)[k] == 0)
				os << k << " " << 0.0001 << "\n";
			else
				os << k << " " << (*prop_vec)[k] / 100.0 << "\n";
		}
	}
	if (prop_name.compare("Permeability") == 0) {
		const double scale(9.86923 * 1e-16 * 9.81 * 1e3 * 1e3);
		for (std::size_t k(0); k < (*prop_vec).size(); k++) {
			if ((*prop_vec)[k] == 0)
				os << k << " " << 1e-3 * scale << "\n";
			else
				os << k << " " << (*prop_vec)[k] * scale << "\n";
		}
	}

	os << "#STOP\n";
	os.close();
}

int main(int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog::Cout *logog(new logog::Cout);
	logog->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Programm reads parts of a Gocad structured mesh", ' ', "0.1");

	// Define a value argument and add it to the command line.
	// A value arg defines a flag and a type of value that it expects
	TCLAP::ValueArg<std::string> sg_file_arg("s", "sg",
		"structured grid file name", true, "", "string");

	// Add the argument sg_file_arg to the CmdLine object. The CmdLine object
	// uses this Arg to parse the command line.
	cmd.add(sg_file_arg);

	TCLAP::ValueArg<bool> mesh_output_arg("", "output-mesh", "output the mesh", false, false,
			"default false, i.e. do not output the mesh");
	cmd.add(mesh_output_arg);

	TCLAP::ValueArg<bool> face_set_arg("f", "generate-face-sets", "generate face sets", false, false,
			"default false, i.e. do not generate face sets");
	cmd.add(face_set_arg);

	TCLAP::ValueArg<bool> output_properties_arg("p", "output-properties",
			"write the properties associated with mesh elements to a file", false, false,
			"default false, i.e. do not write properties");
	cmd.add(output_properties_arg);

	cmd.parse(argc, argv);

	// read the Gocad SGrid
	INFO("Start reading Gocad SGrid.");
	FileIO::GocadSGridReader reader(sg_file_arg.getValue());
	INFO("End reading Gocad SGrid.");

	if (face_set_arg.getValue()) {
		INFO("Generating a mesh for every face set.");
		generateFaceSetMeshes(reader, BaseLib::extractPath(sg_file_arg.getValue()));
	}

	if (output_properties_arg.getValue()) {
		MeshLib::Mesh *mesh(reader.getMesh());
		INFO("Add Gocad properties to mesh.");
		addGocadPropertiesToMesh(reader, *mesh);
		INFO("Done.");
		{
			std::string prop_name("Porosity");
			std::string prop_fname(BaseLib::dropFileExtension(sg_file_arg.getValue()) +
							"-" + prop_name + ".txt");
			writeMeshPropertiesToFile(prop_fname, prop_name, *mesh);
		}

		{
			std::string prop_name("Permeability");
			std::string prop_fname(BaseLib::dropFileExtension(sg_file_arg.getValue()) +
							"-" + prop_name + ".txt");
			writeMeshPropertiesToFile(prop_fname, prop_name, *mesh);
		}
		delete mesh;
	}

	if (mesh_output_arg.getValue()) {
		MeshLib::Mesh *mesh(reader.getMesh());

		INFO("Add Gocad properties to mesh.");
		addGocadPropertiesToMesh(reader, *mesh);

		INFO("Writing mesh in vtu format.");
		FileIO::BoostVtuInterface vtu;
		vtu.setMesh(mesh);
		// output file name
		std::string mesh_out_fname(BaseLib::dropFileExtension(sg_file_arg.getValue()) + ".vtu");

		vtu.writeToFile(mesh_out_fname);

		delete mesh;
	}

	delete logog;
	delete custom_format;
	LOGOG_SHUTDOWN();

	return 0;
}
