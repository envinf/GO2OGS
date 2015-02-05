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

#include <vtkCell.h>
#include <vtkCellLocator.h>
#include <vtkCellData.h>
#include <vtkImageData.h>
#include <vtkCharArray.h>

#include "common.h"

void
GetCellCenter(vtkImageData* img, unsigned int const cellId, double center[3])
{
    double pcoords[3] = {0,0,0};
    double *weights = new double [img->GetMaxCellSize()];
    vtkCell* cell = img->GetCell(cellId);
    int subId = cell->GetParametricCenter(pcoords);
    cell->EvaluateLocation(subId, pcoords, center, weights);
}

vtkSmartPointer<vtkCellLocator>
createCellLocator(vtkSmartPointer<vtkPointSet> mesh)
{
    vtkSmartPointer<vtkCellLocator> locator =
        vtkSmartPointer<vtkCellLocator>::New();

    locator->SetDataSet(mesh);
    locator->SetNumberOfCellsPerNode(4);
    locator->CacheCellBoundsOff();
    locator->BuildLocator();

    return locator;
}

int
main(int argc, char* argv[])
{
    std::cout << std::setprecision(std::numeric_limits<long double>::digits10);

    if ( argc != 3 )
    {
        std::cerr << "Usage: " << argv[0] << "InputFile(vtu) OutputFile(vti).\n";
        return EXIT_FAILURE;
    }

    std::string const input_filename = argv[1];
    std::string const output_filename = argv[2];

    vtkSmartPointer<vtkUnstructuredGrid> mesh = readUGrid(input_filename, false);

    Bbox bbox = getBboxWithoutMaterial12(mesh);
    //mesh->GetBounds(bbox.data());

    std::cout << "bbox:\n";
    std::copy(bbox.begin(), bbox.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;

    int const x_cells = 100;
    int const y_cells = x_cells * (bbox[3]-bbox[2])/(bbox[1] - bbox[0]);
    int const z_cells = 150;
    vtkSmartPointer<vtkImageData> img = createImage(bbox, {x_cells, y_cells, z_cells});
    std::cout << "Input cell number: " << mesh->GetNumberOfCells() << "\n";
    std::cout << "Image dimensions: " << x_cells << " " << y_cells << " " << z_cells << "\n";

    vtkSmartPointer<vtkIntArray> material_ids =
        vtkSmartPointer<vtkIntArray>::New();
    material_ids->SetNumberOfComponents(1);
    material_ids->SetNumberOfTuples(img->GetNumberOfCells());
    material_ids->SetName("MaterialIDs");

    vtkSmartPointer<vtkCharArray> valid_cells =
        vtkSmartPointer<vtkCharArray>::New();
    valid_cells->SetNumberOfComponents(1);
    valid_cells->SetNumberOfTuples(img->GetNumberOfCells());
    valid_cells->SetName("ValidCells");

    vtkSmartPointer<vtkDataArray> mesh_material_ids = mesh->GetCellData()->GetScalars("reservoir_units");


    //
    // Locate image's cell centers in the mesh.
    //
    vtkIdType cell_id_in_mesh = 0;

    vtkSmartPointer<vtkCellLocator> locator = createCellLocator(mesh);
    std::cout << *locator;

    //vtkExecutionTimer timer;
    //timer.StartTimer();
    // Run over all image's cells.
    double center[3] = {0,0,0};
    for(vtkIdType ci = 0; ci < img->GetNumberOfCells(); ++ci)
    {
        //std::cout << ci << std::endl;
        int is_interesting = 0;

        GetCellCenter(img, ci, center);

        // Find a mesh's cell containing this center point.
        cell_id_in_mesh = locator->FindCell(center);

        // Set valid flag and material id if cell found, and non-valid otherwise.
        bool const valid = (cell_id_in_mesh == -1) ? false : true;
        valid_cells->SetValue(ci, valid);
        if (valid)
        {
            material_ids->SetValue(ci,
                mesh_material_ids->GetTuple1(cell_id_in_mesh));
        }
        else
        {
            material_ids->SetValue(ci,-1);
        }
    }
    //timer.StopTimer();
    //std::cout << "time : " << timer.GetElapsedCPUTime() << "\n";

    img->GetCellData()->AddArray(material_ids);
    img->GetCellData()->AddArray(valid_cells);

    writeImage(img, output_filename + ".vti");

    writeMesh(mesh, output_filename + ".vtu");
    return EXIT_SUCCESS;
}
