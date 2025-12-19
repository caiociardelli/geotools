#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
 GeoTools

 Author: Caio Ciardelli, Northwestern University, October 2025

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along
 with this program; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

-----------------------------------------------------------------------------------------------

This script computes seismic ray paths for multiple source-receiver pairs using the SeisTracer library.
It generates a directory for each pair containing VTK files for ray paths, text files with coordinates
and branch types, phase lists, and optional PNG plots. It also saves the receiver elevation in each
directory.

Requirements:
- seistracer library
- numpy library
- os and shutil for file handling
- argparse for command-line arguments

-----------------------------------------------------------------------------------------------
"""

import os
import shutil
import numpy as np
import seistracer
import argparse

EPSILON = 1e-15  # Small value to avoid division by zero

def rThetaPhi2XYZ(r, t, p):
  # Convert spherical coordinates to Cartesian coordinates
  x = r * np.sin(t) * np.cos(p)
  y = r * np.sin(t) * np.sin(p)
  z = r * np.cos(t)
  return x, y, z

def crossProduct(x1, y1, z1, x2, y2, z2):
  # Compute the cross product of two vectors
  vec1 = np.array([x1, y1, z1])
  vec2 = np.array([x2, y2, z2])
  
  axis = np.cross(vec1, vec2)
  norm = np.linalg.norm(axis)
  
  if norm < EPSILON:
    # Handle collinear vectors by choosing a perpendicular basis
    i = np.argmin(np.abs(vec1))
    basis = np.zeros(3)
    basis[i] = 1.0
    axis = np.cross(vec1, basis)
    norm = np.linalg.norm(axis)
    
  axis /= norm
  return axis[0], axis[1], axis[2]

def vincenty(t1, p1, t2, p2):
  # Compute great-circle distance using Vincenty's formula
  sin_t1 = np.sin(t1)
  cos_t1 = np.cos(t1)
  sin_t2 = np.sin(t2)
  cos_t2 = np.cos(t2)

  sin_dp = np.sin(p2 - p1)
  cos_dp = np.cos(p2 - p1)

  return np.arctan2(np.sqrt((sin_t2 * sin_dp) ** 2
                              +(sin_t1 * cos_t2
                              - cos_t1 * sin_t2 * cos_dp) ** 2),
                     cos_t1 * cos_t2 + sin_t1 * sin_t2 * cos_dp)

def greatCircle(t1, p1, t2, p2, deltas, radii):
  # Compute points along the great circle path
  x1, y1, z1 = rThetaPhi2XYZ(1.0, t1, p1)
  x2, y2, z2 = rThetaPhi2XYZ(1.0, t2, p2)

  u, v, w = crossProduct(x1, y1, z1, x2, y2, z2)

  N = len(deltas)

  x = np.empty(N)
  y = np.empty(N)
  z = np.empty(N)

  for i in range(N):
    sinc = np.sin(deltas[i])
    cosc = np.cos(deltas[i])

    # Rotation matrix for great circle interpolation
    R = np.matrix([[cosc +(u ** 2) *(1 - cosc),
                     u * v *(1 - cosc) - w * sinc,
                     u * w *(1 - cosc) + v * sinc],
                    [v * u *(1 - cosc) + w * sinc,
                     cosc +(v ** 2) *(1 - cosc),
                     v * w *(1 - cosc) - u * sinc],
                    [w * u *(1 - cosc) - v * sinc,
                     w * v *(1 - cosc) + u * sinc,
                     cosc +(w ** 2) *(1 - cosc)]])

    p = np.matrix([[x1],
                   [y1],
                   [z1]])

    x[i], y[i], z[i] = radii[i] *(R * p).A1
  return x, y, z

def read_receivers(filename):
  # Read receiver data from file, skipping header
  receivers = []
  with open(filename, 'r') as f:
    for line in f:
      if line.startswith('#'):
        continue
      parts = line.split()
      receiver_id = parts[0]
      lat = float(parts[1])
      lon = float(parts[2])
      elev = float(parts[3])
      receivers.append((receiver_id, lat, lon, elev))
  return receivers

def read_sources(filename):
  # Read source(earthquake) data from file, skipping header
  sources = []
  with open(filename, 'r') as f:
    for line in f:
      if line.startswith('#'):
        continue
      parts = line.split()
      event_id = parts[0]
      lat = float(parts[1])
      lon = float(parts[2])
      depth = float(parts[3])
      sources.append((event_id, lat, lon, depth))
  return sources

if __name__ == '__main__':

  # Parse command-line arguments
  parser = argparse.ArgumentParser(description="Process seismic ray paths for source-receiver pairs")

  parser.add_argument("-r", "--receivers", required=True, help="Path to receivers.dat file")
  parser.add_argument("-s", "--sources", required=True, help="Path to sources.dat file")
  parser.add_argument("--no-png", action="store_true", help="Do not generate PNG plots")
  parser.add_argument("--no-vtk", action="store_true", help="Do not generate VTK files")

  args = parser.parse_args()

  # Load receiver and source data
  receivers = read_receivers(args.receivers)
  sources = read_sources(args.sources)

  # Initialize tracer object
  tracer = seistracer.Tracer()

  # Mapping for branch types
  mapping = {
      'P': 'P',
      'S': 'S',
      'p': 'P',
      's': 'S',
      'K': 'P',
      'I': 'P',
      'J': 'S'
  }

  for source in sources:

    event_id, lat_source, lon_source, source_depth = source
    print(f"Processing source {event_id}...")

    # Convert source coordinates to radians
    t1, p1 = np.radians(90.0 - lat_source), np.radians(lon_source)

    for receiver in receivers:

      receiver_id, lat_receiver, lon_receiver, receiver_elevation = receiver
      print(f"Processing receiver {receiver_id}...")

      # Convert receiver coordinates to radians
      t2, p2 = np.radians(90.0 - lat_receiver), np.radians(lon_receiver)

      # Compute great-circle distance
      gcarc = np.degrees(vincenty(t1, p1, t2, p2))

      # Compute travel times and paths
      phases = tracer.ttimesAndPaths(source_depth=source_depth, delta=gcarc, arc='minor')

      if not args.no_png:
        # Plot overall ray paths
        phases.plotPath(save_figure=True)
        shutil.move(f'ray_paths_{source_depth:.0f}_{gcarc:.0f}.png',
                    f'ray_paths_{event_id}_{receiver_id}.png')

      # Create directory for each source-receiver pair
      dir_name = f"phases_{event_id}_{receiver_id}"

      if os.path.exists(dir_name):
        shutil.rmtree(dir_name)
      os.mkdir(dir_name)

      # Prepare phase list string
      string = f'# {len(phases.phases.keys())}\n'

      for phase in phases.phases.keys():

        string += f'{phase}\n'

        # Get path and info for phase
        paths = phases.getPath([phase])
        info  = phases.getDict([phase])

        if not args.no_png:
          # Plot and save ray path figure for phase
          phases.plotPath([phase], save_figure=True)
          shutil.move(f'ray_paths_{source_depth:.0f}_{gcarc:.0f}.png',
                      f'{dir_name}/{phase}.png')

        # Adjust deltas for spherical coordinates
        deltas = np.pi / 2 - paths[phase][0]['deltas']
        radii  = paths[phase][0]['radii']
        branches = paths[phase][0]['branches']

        # Compute Cartesian coordinates along great circle
        x, y, z = greatCircle(t1, p1, t2, p2, deltas, radii)

        # Map branch types
        branches = np.vectorize(mapping.get)(branches)

        # Print phase and travel time
        print(phase, info[phase][0]['ttime'])

        N = x.size

        # Write source-receiver VTK file
        with open(f'{dir_name}/sr.vtk', 'w') as FILE:
          FILE.write('# vtk DataFile Version 3.1\n')
          FILE.write('Source and Receiver VTK file\n')
          FILE.write('ASCII\n\n')
          FILE.write('DATASET POLYDATA\n')
          FILE.write(f'POINTS 2 float\n')
          FILE.write(f'{x[0]: .10E} {y[0]: .10E} {z[0]: .10E}\n')
          FILE.write(f'{x[-1]: .10E} {y[-1]: .10E} {z[-1]: .10E}\n')

        if not args.no_vtk:
          # Write ray path VTK file
          with open(f'{dir_name}/{phase}.vtk', 'w') as FILE:
            FILE.write('# vtk DataFile Version 3.1\n')
            FILE.write(f'Ray path for {phase}\n')
            FILE.write('ASCII\n\n')
            FILE.write('DATASET POLYDATA\n')
            FILE.write(f'POINTS {N} float\n')

            for ii in range(N):
              FILE.write(f'{x[ii]: .10E} {y[ii]: .10E} {z[ii]: .10E}\n')

            FILE.write(f'\nLINES 1    {N + 1}\n')
            FILE.write(f'{N}\n')

            for ii in range(N):
              FILE.write(f'{ii}\n')

        # Write ray path text file with branches
        with open(f'{dir_name}/{phase}.txt', 'w') as FILE:    
          FILE.write(f'# {N}\n')

          for ii in range(N):
            FILE.write(f'{x[ii]: .10E} {y[ii]: .10E} {z[ii]: .10E} {branches[ii]}\n')

        # Save receiver elevation in the directory
        with open(f'{dir_name}/receiver_elevation.txt', 'w') as FILE:
          FILE.write(f'#ELEVATION {receiver_elevation} m\n')

      # Write phases list to file
      with open(f'{dir_name}/phases_list.txt', 'w') as FILE:
        FILE.write(string[:-1])