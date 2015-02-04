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

#include <array>

#include <vtkImageData.h>
#include <vtkMath.h>
#include <vtkOBBTree.h>
#include <vtkPolyData.h>
#include <vtkSTLWriter.h>
#include <vtkSmartPointer.h>
#include <vtkThreshold.h>
#include <vtkTransform.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

typedef std::array<double, 6> Bbox;
typedef std::array<int, 3> Dimensions;

int bboxFlatness(Bbox const& bb)
{
    return (static_cast<int>((bb[1]-bb[0]) < 1e-1) +
            static_cast<int>((bb[3]-bb[2]) < 1e-1) +
            static_cast<int>((bb[5]-bb[4]) < 1e-1));
}
void
orientMesh(vtkSmartPointer<vtkPointSet> mesh,
    double const angle, double const translation[3])
{
    // Create linear transformation removing xy rotation and translation.
    vtkSmartPointer<vtkTransform> tr1 = vtkSmartPointer<vtkTransform>::New();
    vtkSmartPointer<vtkTransform> tr2 = vtkSmartPointer<vtkTransform>::New();
    tr1->Translate(translation);
    tr2->RotateZ(angle);
    tr1->Concatenate(tr2);

    vtkSmartPointer<vtkPoints> new_mesh_points = vtkSmartPointer<vtkPoints>::New();
    tr1->GetInverse()->TransformPoints(mesh->GetPoints(), new_mesh_points);
    mesh->SetPoints(new_mesh_points);
}

void
writeSTLSurface(vtkSmartPointer<vtkPolyData> const mesh, std::string const& filename)
{
    vtkSmartPointer<vtkSTLWriter> writer
        = vtkSmartPointer<vtkSTLWriter>::New();

    writer->SetInputData(mesh);
    writer->SetFileName(filename.c_str());
    writer->SetFileTypeToASCII();
    writer->Write();
}

void
writePolyData(vtkSmartPointer<vtkPolyData> const mesh, std::string const& filename)
{
    vtkSmartPointer<vtkXMLPolyDataWriter> writer
        = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

    writer->SetInputData(mesh);
    writer->SetFileName(filename.c_str());
    writer->Write();
}

void
writeMesh(vtkSmartPointer<vtkUnstructuredGrid> const mesh, std::string const& filename)
{
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer
        = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    writer->SetInputData(mesh);
    writer->SetFileName(filename.c_str());
    writer->SetCompressorTypeToNone();
    writer->SetDataModeToAscii();
    writer->Write();
}

void
writeImage(vtkImageData* const img, std::string const& filename)
{
    vtkSmartPointer<vtkXMLImageDataWriter> writer
        = vtkSmartPointer<vtkXMLImageDataWriter>::New();

    writer->SetInputData(img);
    writer->SetFileName(filename.c_str());
    writer->Write();
}

vtkSmartPointer<vtkImageData>
readImage(std::string const& filename)
{
    vtkSmartPointer<vtkXMLImageDataReader> reader
        = vtkSmartPointer<vtkXMLImageDataReader>::New();

    reader->SetFileName(filename.c_str());
    reader->Update();
    return reader->GetOutput();
}

vtkSmartPointer<vtkImageData>
createImage(Bbox const& bb, Dimensions const& dims)
{
    vtkSmartPointer<vtkImageData> img = vtkSmartPointer<vtkImageData>::New();
    img->SetOrigin(bb[0], bb[2], bb[4]);
    img->SetDimensions(dims[0], dims[1], dims[2]);
    img->SetSpacing(
            (bb[1]-bb[0])/(dims[0]-1),
            (bb[3]-bb[2])/(dims[1]-1),
            (bb[5]-bb[4])/(dims[2]-1));
    return img;
}

Bbox
getBboxWithoutMaterial12(vtkSmartPointer<vtkUnstructuredGrid> mesh)
{
    // Remove cells with materialIDs == 12.
    vtkSmartPointer<vtkThreshold> mesh_threshold = 
        vtkSmartPointer<vtkThreshold>::New();

    mesh_threshold->SetInputData(mesh);
    //mesh_threshold->SetAttributeModeToUseCellData();
    mesh_threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "reservoir_units");
    mesh_threshold->ThresholdByLower(11);
    mesh_threshold->Update();

    Bbox bbox;
    mesh_threshold->GetOutput()->GetBounds(bbox.data());
    return bbox;
}

// Project given points on the x,y,0 plane.
vtkSmartPointer<vtkPoints>
flattenPoints(vtkSmartPointer<vtkPoints> input_points)
{
    vtkSmartPointer<vtkPoints> output_points = vtkSmartPointer<vtkPoints>::New();
    for(vtkIdType pi = 0; pi < input_points->GetNumberOfPoints(); ++pi)
    {
        double p[3];
        input_points->GetPoint(pi, p);
        output_points->InsertNextPoint(p[0], p[1], 0);
    }

    return output_points;
}

vtkSmartPointer<vtkUnstructuredGrid>
readUGrid(std::string const& filename, bool const transform = true)
{
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    vtkSmartPointer<vtkUnstructuredGrid> mesh = reader->GetOutput();
    // Return reader output if no transformation is requested.
    if (!transform)
        return mesh;

    //
    // Orient mesh along cartesian coordinates for better resampling.
    //

    vtkSmartPointer<vtkPoints> points2D = flattenPoints(mesh->GetPoints());

    // Compute OBB of the 2D points.
    vtkSmartPointer<vtkOBBTree> obbtree = vtkSmartPointer<vtkOBBTree>::New();
    double corner[3], max[3], mid[3], min[3], size[3];
    obbtree->ComputeOBB(points2D, corner, max, mid, min, size);

    /*
    std::cout << "OBB of input mesh:\n"
        << " corner: " << corner[0] << " " << corner[1] << " " << corner[2]
        << " max: " << max[0] << " " << max[1] << " " << max[2]
        << " mid: " << mid[0] << " " << mid[1] << " " << mid[2]
        << " min: " << min[0] << " " << min[1] << " " << min[2]
        << " size: " << size[0] << " " << size[1] << " " << size[2]
        << std::endl;
        */

    double const angle = vtkMath::DegreesFromRadians(atan2(max[1], max[0]));
    double const translation[3] = {corner[0] + mid[0], corner[1] + mid[1], 0};

	std::ofstream os("GeometricTransformationData.bin", std::ios::binary);
	os << angle << translation[0] << translation[1] << translation[2];
	os.close();
	std::ofstream os_ascii("GeometricTransformationData.txt");
	os_ascii << angle << "\n";
	os_ascii << translation[0] << " " << translation[1] << " " << translation[2];
	os_ascii.close();
	std::cout << "GeometricTransformationData\n";
	std::cout << angle << "\n";
	std::cout << translation[0] << " " << translation[1] << " " << translation[2];

    orientMesh(mesh, angle, translation);

    return mesh;
}

vtkSmartPointer<vtkPolyData> const
readPolyData(std::string const& filename)
{
    vtkSmartPointer<vtkXMLPolyDataReader> reader
        = vtkSmartPointer<vtkXMLPolyDataReader>::New();

    reader->SetFileName(filename.c_str());
    reader->Update();

    return reader->GetOutput();
}

