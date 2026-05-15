#!/usr/bin/env bash

# Test Workflow
#
# Author: Caio Ciardelli, Northwestern University, October 2025
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#-----------------------------------------------------------------------------------------------
#
# This script generates GMT plots of mesh data from text files containing triangle, center, and edge
# information. It creates a separate plot for each mesh number found in the available files and, if
# the maximum mesh number is 320, generates an additional plot for a nodes_spiral figure.
#
# Requirements:
# - Bash
# - GCC
# - Python 3.7 or later
# - Joint Inversion Program
# - Numpy
# - Matplotlib
# - GMT (Generic Mapping Tools) installed
# - ParaView installation with Python support
# - paraview Python package available
# - OpenMPI
#
#-----------------------------------------------------------------------------------------------

help_menu ()
{
  echo " RUN_TEST_WORKFLOW

  USAGE
    ./run_test_workflow.bash

  EXAMPLE
    ./run_test_workflow.bash

  COMMAND-LINE ARGUMENTS
    None                - The script runs the test workflow.

  DESCRIPTION
    Runs a test workflow for the Joint Inversion Package. It generates the mesh, computes seismic-wave
    phases, ray paths, and travel times using the PyTracer package, and calculates numerical travel times
    and synthetic gravity disturbances using the Joint Inversion Package. Outputs include GMT plots showing
    mesh triangles, centers, and edges, highlighting one mesh triangle and its neighbours. For a maximum
    mesh number of 320, an additional plot displays edges, spiral, and node data, saved as PNG images.
    Additionally, the mesh is visualized using ParaView for interactive 3D rendering of the generated mesh
    data.

  REQUIREMENTS
    - Bash
    - GCC
    - Python 3.7 or later
    - Joint Inversion Program
    - Numpy
    - Matplotlib
    - GMT (Generic Mapping Tools) installed
    - ParaView installation with Python support
    - paraview Python package available
    - OpenMPI"
}

# Check if help is requested
if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
  help_menu
  exit 0
fi

# Build and run the Joint Inversion Program test workflow

# Announce the start of the build process
echo "Bulding Joint Inversion Program..."

# Remove previous output files and directories to ensure a clean build
rm -fv input/1d_moho_radius.dat input/input.model
rm -fv mesh/*
rm -rfv phases_EQ001_STA01/ *.txt *.dat *.png *.vtk *.pvsm *.grd gmt.history d01a_rescaled

# Clean previous build artifacts
make clean
# Compile all source files to build the program
make all

# Grant execution permissions to utility scripts
chmod +x utils/*.bash

# Announce the start of the test workflow
echo ""
echo "Running Test Workflow..."

# Resample the 1D IASP91 model for input
echo ""
echo "Resampling 1D model IASP91..."
python utils/resample_1d_model.py -i models/iasp91.model -o input/input.model

# Compute seismic ray paths using the PyTracer package
echo ""
echo "Computing ray paths using PyTracer..."
python utils/compute_rays.py -r input/receivers.dat -s input/sources.dat

# Generate the mesh using the mesher binary
echo ""
echo "Running mesher..."
./bin/mesher 7 input/input.model

# Create GMT plots to visualize mesh triangles, centers, and edges
echo ""
echo "Plotting mesh..."
./utils/plot_mesh.bash

# Generate a ParaView state file for 3D visualization of mesh and ray paths
echo ""
echo "Creating Paraview plot to visualize mesh and ray paths..."
pvpython utils/pvsm_generator.py mesh/Mesh_VpVsRho.vtk --color-by Vp --phases-dir phases_EQ001_STA01 --phases PKIIKP PKPPm+PPcS ScSScSm-ScS SPS660-S \
         --reverse-colormap --camera-view 210,10,0,0,20 --label 'Vp [km/s]' --light-intensity 0.5 --light-azimuth 140 --light-elevation 35

# Compute travel times using the ttimes binary
echo ""
echo "Computing travel times using Ttimes..."
./bin/ttimes phases_EQ001_STA01/ > ttimes.txt

# Define parameters for parallel gravity computation
NCORES=6
DELTA=5.0
HEIGHT=255

# Compute synthetic gravity disturbances using MPI parallel execution
# This step requires around 16 GB of RAM. If you have less memory,
# you can reduce the number of cores, the number of layers, or the
# mesh refinement.
echo ""
echo "Computing synthetic gravity..."
mpiexec -n $NCORES ./bin/gravity $DELTA $HEIGHT

# Plot the computed synthetic gravity disturbances
echo ""
echo "Plotting synthetic gravity..."
./utils/plot_gravity.bash G_normal
