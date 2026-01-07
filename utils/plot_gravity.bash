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
# This script generates a GMT plot of gravitational field data from a text file.
# It creates a grid using blockmean and greenspline, applies a rescaled colormap,
# and outputs a PNG image with a colorbar. The script requires a single
# command-line argument specifying the input file (without extension).
#
# Requirements:
# - GMT (Generic Mapping Tools) installed
# - bc (basic calculator) for arithmetic operations
#
#-----------------------------------------------------------------------------------------------

help_menu ()
{
  echo " PLOT_GRAVITY

  USAGE
    ./gravity_plot.bash FILENAME

  EXAMPLE
    ./gravity_plot.bash G_normal

  COMMAND-LINE ARGUMENTS
    FILENAME              - Base name of the input data file (without .dat extension)

  DESCRIPTION
    Generates a GMT plot of gravitational field data from a text file.
    The script creates a grid using blockmean and greenspline, applies a rescaled
    colormap, and outputs a PNG image with a colorbar. The input file should
    be in a format compatible with GMT (e.g., columns of latitude, longitude,
    and gravity values).

  REQUIREMENTS
    - GMT (Generic Mapping Tools) installed
    - bc (basic calculator) for arithmetic operations"
}

# Check if help is requested or no argument is provided
if [ "$#" -ne 1 ] || [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
  help_menu
  exit 0
fi

# Define label and file names
label="@[g \left[mGal\right]@["
filename="$1".dat
grdname="$1".grd
output="$1"
filename_temp="temp.dat"
grdname_temp="temp.grd"

echo 'Creating figure...'

# Convert XYZ data to grid, trim region to avoid poles (-89/89), 
# convert back to XYZ, and interpolate using GMT surface
gmt xyz2grd $filename -R-180/180/-90/90 -I1.0 -G"$grdname_temp" -:
gmt grdcut $grdname_temp -G"$grdname_temp" -R-180/180/-89/89
gmt grd2xyz $grdname_temp > $filename_temp
gmt surface $filename_temp -R-180/180/-89/89 -I0.25 -G$grdname

# Extract min/max and adjust grid
min=$(gmt grdinfo $grdname | grep 'v_min' | cut -f3 -d' ')
max=$(gmt grdinfo $grdname | grep 'v_max' | cut -f5 -d' ')
mean=$(echo "($min + $max) / 2" | bc -l)
gmt grdmath $grdname $mean SUB = $grdname
gmt grdmath $grdname 1e5 MUL = $grdname

# Update min/max and set colorbar bounds
min=$(gmt grdinfo $grdname | grep 'v_min' | cut -f3 -d' ')
max=$(gmt grdinfo $grdname | grep 'v_max' | cut -f5 -d' ')
cbmin=$(echo "scale=15; 1.001 * $min" | bc -l)
cbmax=$(echo "scale=15; 1.001 * $max" | bc -l)

# Create colormap
gmt makecpt -Cextra/CET-D01A.cpt -T$cbmin/$cbmax -I > d01a_rescaled

# Generate figure
gmt begin $output png
  gmt set FONT_TITLE 20p,Helvetica
  gmt grdimage $grdname -JR0/12c -Rg -Bxa90fg90 -Bya30fg30 -BWeSn+t"Synthetic Gravity - ETOPO" -Cd01a_rescaled
  gmt coast -W0.2,50
  gmt colorbar -Cd01a_rescaled -Baf -DJBC+e -B+l"$label"
gmt end