#!/usr/bin/env bash

# GeoTools
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
# - GMT (Generic Mapping Tools) installed
#
#-----------------------------------------------------------------------------------------------

help_menu ()
{
  echo " PLOT_MESH

  USAGE
    ./mesh_plot.bash

  EXAMPLE
    ./mesh_plot.bash

  COMMAND-LINE ARGUMENTS
    None                - The script automatically detects mesh-related files in the current directory

  DESCRIPTION
    Generates GMT plots of mesh data from text files (mesh_triangles_*.txt, mesh_center_triangle_*.txt,
    mesh_edges_*.txt). Each plot displays triangles, centers, and edges (for mesh numbers <= 20480).
    If the maximum mesh number is 320, an additional nodes_spiral plot is created with specific edges,
    spiral, and node data. Outputs are saved as PNG images.

  REQUIREMENTS
    - GMT (Generic Mapping Tools) installed"
}

# Check if help is requested
if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
  help_menu
  exit 0
fi

echo 'Checking for available mesh files...'

# List available mesh-related files
available_files=$(ls mesh_triangles_*.txt mesh_center_triangle_*.txt mesh_edges_*.txt 2>/dev/null)

# Extract unique mesh numbers from available files
mesh_numbers=$(echo "$available_files" | grep -oP '\d+' | sort -nu)

echo 'Creating figures for available mesh files...'

# Get the maximum mesh number
max_mesh_number=$(echo "$mesh_numbers" | tail -n 1)

# Loop through each mesh number and plot if corresponding files exist
for num in $mesh_numbers; do
  triangles_file="mesh_triangles_${num}.txt"
  center_file="mesh_center_triangle_${num}.txt"
  edges_file="mesh_edges_${num}.txt"

  # Check if triangles file exists and plot if it does
  if [ -f "$triangles_file" ]; then
    gmt begin mesh${num} png E300
      gmt set FONT_TITLE 20p,Helvetica
      gmt basemap -Rg -JH4.5i -Ba
      gmt coast -Rg -Dc -A10000 -W0.3
      gmt plot "$triangles_file" -Gyellow
      gmt plot "$center_file" -Gred

      # Plot edges if num is less than or equal to 20480
      if [ ${num} -le 20480 ]; then
          gmt plot "$edges_file" -W0.1,gray
      fi
    gmt end
  fi
done

# Check if the maximum mesh number is 320 and create nodes_spiral plot
if [ $max_mesh_number -eq 320 ]; then
  echo 'Creating nodes_spiral figure...'
  gmt begin nodes_spiral png E300
    gmt set FONT_TITLE 20p,Helvetica

    gmt basemap -Rg -JA280/30/4.5i -BWeSn -Byg30
    gmt coast -JA280/30/4.5i -Di -A100 -W0.2

    gmt plot mesh_edges_80.txt -W0.1,gray
    gmt plot spiral.txt -W0.1,magenta
    gmt plot nodes.txt -Sc0.1 -Ggreen -W0.1
  gmt end
fi

echo 'Figures saved!'