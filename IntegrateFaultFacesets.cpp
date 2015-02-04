/*
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Copyright (C) 2014-2015 Helmholtz Centre for Environmental Research -- UFZ.

Author: Dmitri Naumov, Thomas Fischer
*/

#include <iostream>
#include <iterator>
#include <limits>

#include "tclap/CmdLine.h"

#include <vtkCell.h>
#include <vtkCellLocator.h>
#include <vtkCellData.h>
#include <vtkImageData.h>
#include <vtkTable.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkPCAStatistics.h>

#include "common.h"

void computeSearchRectangle(double p0[3], double p1[3],
	double img_origin[3], double img_spacing[3], int img_extent[6],
	int min[3], int max[3])
{
	int p0_idx[3], p1_idx[3];
	for (std::size_t i(0); i<3; ++i)
		p0_idx[i] = (p0[i] - img_origin[i])/img_spacing[i];
	for (std::size_t i(0); i<3; ++i)
		p1_idx[i] = (p1[i] - img_origin[i])/img_spacing[i];
	for (std::size_t i(0); i<3; ++i)
		min[i] = std::min(p0_idx[i], p1_idx[i]);
	for(std::size_t i(0); i<3; ++i)
		max[i] = std::max(p0_idx[i], p1_idx[i]);

	// make rectangle for search a little bit bigger
	min[0] = std::max(img_extent[0], min[0]-1);
	min[1] = std::max(img_extent[2], min[1]-1);
	max[0] = std::min(img_extent[1], max[0]+1);
	max[1] = std::min(img_extent[3], max[1]+1);
}

void markVoxelColumnsAlongFacesetsLine(
	double p0[3], double p1[3],
	int faceset_id,
	vtkSmartPointer<vtkImageData> img)
{
	int *extent(img->GetExtent());
	double *img_origin(img->GetOrigin());
	double *img_spacing(img->GetSpacing());

	int array_id;
	vtkDataArray *mat_ids = img->GetCellData()->GetArray("MaterialIDs", array_id);
	vtkDataArray *valid = img->GetCellData()->GetArray("ValidCells", array_id);

	int offset(static_cast<int>(mat_ids->GetRange()[1])+1);

	int min[3], max[3];
	computeSearchRectangle(p0, p1, img_origin, img_spacing, extent, min, max);
	int dim[3] = {extent[1]-extent[0], extent[3]-extent[2], extent[5]-extent[4]};
	int ijk[3];
	for (ijk[0] = min[0]; ijk[0]<max[0]; ++(ijk[0])) {
		for (ijk[1] = min[1]; ijk[1]<max[1]; ++(ijk[1])) {
			// compute distance of voxel center to straight line
			double const c[2] = {
				img_origin[0]+(ijk[0]+0.5)*img_spacing[0],
				img_origin[1]+(ijk[1]+0.5)*img_spacing[1]
			};

//			// this would be an option, too:
//			vtkMath::ProjectVector2D (const double a[2],
//				const double b[2],
//				double projection[2])

			double const h1(
				(c[0]-p0[0])*(p1[0]-p0[0])+(c[1]-p0[1])*(p1[1]-p0[1])
			);
			double const h2(pow((p1[0]-p0[0]),2)+pow((p1[1]-p0[1]),2));

			double const lambda (h1/h2);
			double const pc[2] = {
				p0[0] + lambda * (p1[0]-p0[0]),
				p0[1] + lambda * (p1[1]-p0[1])
			};
			double const sqr_dist(
				std::pow(pc[0]-c[0], 2) + std::pow(pc[1]-c[1],2)
			);
			if (sqr_dist < (std::pow(img_spacing[0], 2) + std::pow(img_spacing[1], 2))) {
				for (ijk[2]=0; ijk[2]<dim[2]; ++ijk[2]) {
					if (valid->GetComponent(img->ComputeCellId(ijk), 0) == 0)
						continue;
					if (mat_ids->GetComponent(img->ComputeCellId(ijk), 0) >= offset)
						continue;
					int const v(mat_ids->GetComponent(img->ComputeCellId(ijk),0));
					mat_ids->SetTuple1(img->ComputeCellId(ijk),
						v+offset*(faceset_id+1));
				}
			}
		}
	}
}

void computePCAOfMeshNodes(vtkSmartPointer<vtkUnstructuredGrid> faceset_mesh,
	std::size_t faceset_id)
{
	vtkPoints* pnts = faceset_mesh->GetPoints();

	vtkSmartPointer<vtkDoubleArray> x = vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkDoubleArray> y = vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkDoubleArray> z = vtkSmartPointer<vtkDoubleArray>::New();

	x->SetNumberOfComponents(1);
	y->SetNumberOfComponents(1);
	z->SetNumberOfComponents(1);
	x->SetName("x");
	y->SetName("y");
	z->SetName("z");
	double coords[3];
	for (vtkIdType i(0); i<pnts->GetNumberOfPoints(); ++i) {
		pnts->GetPoint(i, coords);
		x->InsertNextValue(coords[0]);
		y->InsertNextValue(coords[1]);
		z->InsertNextValue(coords[2]);
	}

	// create a table
	vtkSmartPointer<vtkTable> dataset_table = vtkSmartPointer<vtkTable>::New();
	dataset_table->AddColumn(x);
	dataset_table->AddColumn(y);
	dataset_table->AddColumn(z);

	vtkSmartPointer<vtkPCAStatistics> pca =
		vtkSmartPointer<vtkPCAStatistics>::New();
	pca->SetInputData(vtkStatisticsAlgorithm::INPUT_DATA, dataset_table);
	pca->SetColumnStatus(dataset_table->GetColumnName(0), 1);
	pca->SetColumnStatus(dataset_table->GetColumnName(1), 1);
	pca->SetColumnStatus(dataset_table->GetColumnName(2), 1);
	pca->Update();

	///////// Eigenvectors ////////////
	vtkSmartPointer<vtkDoubleArray> eigenvectors =
		vtkSmartPointer<vtkDoubleArray>::New();
	pca->GetEigenvectors(eigenvectors);
	double* evec = new double[eigenvectors->GetNumberOfComponents()];
	eigenvectors->GetTuple(2, evec);
	std::cout << faceset_id << " ";
	for(vtkIdType j = 0; j < eigenvectors->GetNumberOfComponents(); j++) {
		std::cout << evec[j] << " ";
	}
	std::cout << std::endl;
	delete evec;

}

int main(int argc, char* argv[])
{
    std::cout << std::setprecision(std::numeric_limits<long double>::digits10);

	TCLAP::CmdLine cmd("Integrate faceset info into vti", ' ', "0.1");
	TCLAP::ValueArg<std::string> vti_input_arg("i", "vti-input",
		"", true, "", "filename");
	cmd.add(vti_input_arg);
	TCLAP::ValueArg<std::string> vti_output_arg("o", "vti-output",
		"", true, "", "filename");
	cmd.add(vti_output_arg);
	TCLAP::MultiArg<std::string> faceset_input_arg("f", "faceset-input",
		"", true, "filename");
	cmd.add(faceset_input_arg);
	cmd.parse( argc, argv );

    std::string const input_faceset_fname = argv[2];

    std::cout << "Reading image: " << vti_input_arg.getValue() << ".\n";
    vtkSmartPointer<vtkImageData> img = readImage(vti_input_arg.getValue());
    std::cout << "readImage: #cells: " << img->GetNumberOfCells() << "\n";

	std::vector<std::string> const& fnames(faceset_input_arg.getValue());
	for (std::size_t faceset(0); faceset<fnames.size(); ++faceset) {
//		std::cout << "Reading faceset mesh: " << fnames[faceset] << ".\n";
		vtkSmartPointer<vtkUnstructuredGrid>
			faceset_mesh = readUGrid(fnames[faceset], false);
//		std::cout << "readUGrid: #cells: " << faceset_mesh->GetNumberOfCells() << "\n";
		for (size_t quad_id(0); quad_id<faceset_mesh->GetNumberOfCells(); ++quad_id) {
			vtkPoints *quad_pnts(faceset_mesh->GetCell(quad_id)->GetPoints());
			double p0[3], p1[3];
			quad_pnts->GetPoint(0, p0);
			quad_pnts->GetPoint(3, p1);

			markVoxelColumnsAlongFacesetsLine(p0, p1, faceset, img);
		}
		computePCAOfMeshNodes(faceset_mesh, faceset);
	}
	writeImage(img, vti_output_arg.getValue());
}

