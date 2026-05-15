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

This script uses the ParaView Python API to create .pvsm files and set up
visualizations for VTK meshes. It provides both programmatic and command-line
interfaces for generating ParaView state files.

IMPORTANT: Do NOT name this script 'paraview.py' as it will cause circular import issues!
Recommended filename: 'pvsm_generator.py' or 'vtk_paraview_generator.py'

Requirements:
- ParaView installation with Python support
- paraview Python package available
"""

import os
import sys
import argparse
import math

from pathlib import Path

try:
  # Import ParaView Python modules
  import paraview
  paraview.compatibility.major = 6
  paraview.compatibility.minor = 1
  # Initialize ParaView
  from paraview.simple import *
  
  # Connect to ParaView server (required for proper initialization)
  Connect()
  
except ImportError as e:
  print("Error: ParaView Python modules not found!")
  print("Please ensure ParaView is installed with Python support.")
  print("You may need to add ParaView's Python path to your environment.")
  print("Try using 'pvpython' instead of 'python' to run this script.")
  print(f"Import error: {e}")
  sys.exit(1)
except Exception as e:
  print("Error initializing ParaView:")
  print("This script should be run with ParaView's Python interpreter (pvpython)")
  print("or with proper ParaView environment setup.")
  print(f"Initialization error: {e}")
  print("\nTry running with: pvpython pvsm_generator.py Vp.vtk")
  sys.exit(1)


class ParaViewVTKVisualizer:
  """Class to handle VTK mesh visualization using ParaView Python API."""
  
  def __init__(self):
    """Initialize the visualizer."""
    self.reader = None
    self.representation = None
    self.render_view = None

  def create_slice_visualizations(self, reader, slice_positions, representation_type, 
                                  color_by_scalar=None, reverse_colormap=False, 
                                  show_edges=False, opacity=1.0, colormap="Rainbow Uniform",
                                  font_color=[0.0, 0.0, 0.0], font_size=25, label=None):
    """
    Create slice visualizations for the VTK mesh at specified positions.
    
    Args:
        reader: ParaView reader object for the VTK file
        slice_positions (dict): Dictionary with 'x', 'y', 'z' positions for slices
        representation_type (str): Representation type ('Surface', 'Wireframe', etc.)
        color_by_scalar (str): Name of scalar field to color by
        reverse_colormap (bool): Whether to reverse the colormap
        show_edges (bool): Whether to show mesh edges
        opacity (float): Opacity value (0.0 to 1.0)
        colormap (str): Colormap to use for scalar data
        font_color (list): RGB color for colorbar label [0.0-1.0]
        font_size (int): Font size for colorbar label
        label (str): Text for colorbar label (defaults to color_by_scalar or 'VTK Mesh')
        
    Returns:
        list: List of (slice_filter, slice_representation) tuples
    """
    slice_representations = []
    
    # Set default label to color_by_scalar or 'VTK Mesh'
    if label is None:
      label = color_by_scalar if color_by_scalar else "VTK Mesh"
    
    # Create slices along X, Y, Z planes
    for axis, position in slice_positions.items():
      slice_filter = Slice(Input=reader)
      slice_filter.SliceType = 'Plane'
      
      # Set slice plane
      if axis == 'x':
        slice_filter.SliceType.Normal = [1, 0, 0]
        slice_filter.SliceType.Origin = [position, 0, 0]
      elif axis == 'y':
        slice_filter.SliceType.Normal = [0, 1, 0]
        slice_filter.SliceType.Origin = [0, position, 0]
      elif axis == 'z':
        slice_filter.SliceType.Normal = [0, 0, 1]
        slice_filter.SliceType.Origin = [0, 0, position]
      
      # Show the slice
      slice_repr = Show(slice_filter, self.render_view)
      slice_repr.Representation = representation_type
      slice_repr.Opacity = opacity
      
      if show_edges:
        slice_repr.EdgeColor = [0.0, 0.0, 0.0] # Black edges
        slice_repr.ShowEdges = 1
      
      # Color by scalar if specified
      if color_by_scalar:
        try:
          ColorBy(slice_repr, ('POINTS', color_by_scalar))
          lookup_table = GetColorTransferFunction(color_by_scalar)
          lookup_table.ApplyPreset(colormap, True)
          if reverse_colormap:
            lookup_table.InvertTransferFunction()
          slice_repr.LookupTable = lookup_table
          slice_repr.SetScalarBarVisibility(self.render_view, True)
          
          # Customize scalar bar
          scalar_bar = GetScalarBar(lookup_table, self.render_view)
          scalar_bar.Title = label
          scalar_bar.ComponentTitle = ""
          scalar_bar.TitleFontSize = font_size
          scalar_bar.LabelFontSize = font_size
          scalar_bar.TitleColor = font_color
          scalar_bar.LabelColor = font_color
          scalar_bar.Position = [0.85, 0.05] # Bottom-right
          scalar_bar.ScalarBarLength = 0.35
          scalar_bar.ScalarBarThickness = 20
          scalar_bar.HorizontalTitle = 1
          scalar_bar.Visibility = 1
          scalar_bar.AddRangeLabels = 0
          
        except Exception as e:
          print(f"Warning: Could not color by scalar '{color_by_scalar}': {e}")
      
      try:
        slice_repr.Specular = 0.0
      except:
        pass
      
      # Rename the slice for clarity
      try:
        RenameSource(f"Slice_{axis.upper()}", slice_filter)
      except:
        pass
      
      slice_representations.append((slice_filter, slice_repr))
    
    return slice_representations

  def create_basic_visualization(self, vtk_file_path, mesh_name="VTK_Mesh", color_by_scalar=None,
                                 font_color=[0.0, 0.0, 0.0], font_size=25, label=None):
    """
    Create a basic visualization for a VTK file with optional scalar coloring.
    
    Args:
        vtk_file_path (str): Path to the VTK file
        mesh_name (str): Name for the mesh object
        color_by_scalar (str): Name of scalar field to color by
        font_color (list): RGB color for colorbar label [0-1]
        font_size (int): Font size for colorbar label
        label (str): Text for colorbar label (defaults to color_by_scalar or 'VTK Mesh')
        
    Returns:
        tuple: (reader, slice_representations, render_view)
    """
    # Create render view
    render_view = GetActiveViewOrCreate('RenderView')
    render_view.Background = [1.0, 1.0, 1.0]
    render_view.UseColorPaletteForBackground = 0
    render_view.ViewSize = [1200, 900]
    
    self.render_view = render_view
    
    # Create reader
    reader = LegacyVTKReader(FileNames=[os.path.abspath(vtk_file_path)])
    reader.UpdatePipeline()
    
    try:
      RenameSource(mesh_name, reader)
    except:
      pass
    
    self.reader = reader
    
    # Get bounds for slice positions
    bounds = reader.GetDataInformation().GetBounds()
    center_x = (bounds[0] + bounds[1]) / 2.0
    center_y = (bounds[2] + bounds[3]) / 2.0
    center_z = (bounds[4] + bounds[5]) / 2.0
    slice_positions = {'x': center_x, 'y': center_y, 'z': center_z}
    
    # Create slices with colorbar configuration
    slice_representations = self.create_slice_visualizations(
      reader, slice_positions, 'Surface', color_by_scalar,
      font_color=font_color, font_size=font_size, label=label
    )
    
    return reader, slice_representations, render_view

  def create_advanced_visualization(self, vtk_file_path, mesh_name="VTK_Mesh",
                                    representation_type="Surface", color_by_scalar=None,
                                    reverse_colormap=False, show_edges=False, opacity=1.0,
                                    slice_positions=None, colormap="Rainbow Uniform",
                                    font_color=[0.0, 0.0, 0.0], font_size=25, label=None,
                                    **kwargs):
    """
    Create an advanced visualization with custom settings for a single VTK file.
    
    Args:
        vtk_file_path (str): Path to the VTK file
        mesh_name (str): Name for the mesh object
        representation_type (str): Representation type
        color_by_scalar (str): Name of scalar field to color by
        reverse_colormap (bool): Whether to reverse the colormap
        show_edges (bool): Whether to show mesh edges
        opacity (float): Opacity value (0.0 to 1.0)
        slice_positions (dict): Dictionary with 'x', 'y', 'z' positions for slices
        colormap (str): Colormap to use for scalar data
        font_color (list): RGB color for colorbar label [0-1]
        font_size (int): Font size for colorbar label
        label (str): Text for colorbar label (defaults to color_by_scalar or 'VTK Mesh')
        **kwargs: Additional keyword arguments
        
    Returns:
        tuple: (reader, slice_representations, render_view)
    """
    # Create render view
    render_view = GetActiveViewOrCreate('RenderView')
    render_view.Background = kwargs.get('background_color', [1.0, 1.0, 1.0])
    render_view.UseColorPaletteForBackground = 0
    render_view.ViewSize = kwargs.get('view_size', [1200, 900])
    
    self.render_view = render_view
    
    # Create reader
    reader = LegacyVTKReader(FileNames=[os.path.abspath(vtk_file_path)])
    reader.UpdatePipeline()
    
    try:
      RenameSource(mesh_name, reader)
    except:
      pass
    
    self.reader = reader
    
    # Get bounds for slice positions
    bounds = reader.GetDataInformation().GetBounds()
    center_x = (bounds[0] + bounds[1]) / 2.0
    center_y = (bounds[2] + bounds[3]) / 2.0
    center_z = (bounds[4] + bounds[5]) / 2.0
    
    # Use provided slice positions or default to center
    if slice_positions is None:
      slice_positions = {'x': center_x, 'y': center_y, 'z': center_z}
    
    # Create slices with colorbar configuration
    slice_representations = self.create_slice_visualizations(
      reader, slice_positions, representation_type, color_by_scalar,
      reverse_colormap, show_edges, opacity, colormap,
      font_color, font_size, label
    )
    
    return reader, slice_representations, render_view
  
  def _set_camera_view(self, render_view, camera_position=None, camera_focus=None, camera_up=None):
    """Set custom camera view."""
    camera = render_view.GetActiveCamera()
    
    if camera_position:
      camera.SetPosition(camera_position)
    if camera_focus:
      camera.SetFocalPoint(camera_focus)
    if camera_up:
      camera.SetViewUp(camera_up)
    
    render_view.StillRender()

  def set_camera_view(self, azimuth, elevation, yaw=0.0, pitch=0.0, roll=0.0, render_view=None):
    """
    Set camera view using azimuth, elevation (position), and yaw, pitch, roll (orientation).
    
    Args:
        azimuth (float): Azimuth angle in degrees (0-360, 0=North, 90=East)
        elevation (float): Elevation angle in degrees (-90 to 90, 0=horizon, 90=zenith)
        yaw (float): Yaw angle in degrees (0-360, rotation around camera's local vertical axis)
        pitch (float): Pitch angle in degrees (-90 to 90, tilt up/down around camera's local horizontal axis)
        roll (float): Roll angle in degrees (-180 to 180, rotation around view direction)
        render_view: Render view to modify (uses self.render_view if None)
    """
    if render_view is None:
      render_view = self.render_view
    
    if self.reader:
      # Get mesh center for focal point
      bounds = self.reader.GetDataInformation().GetBounds()
      center_x = (bounds[0] + bounds[1]) / 2.0
      center_y = (bounds[2] + bounds[3]) / 2.0
      center_z = (bounds[4] + bounds[5]) / 2.0
      
      # Calculate distance for camera positioning
      max_extent = max(bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4])
      distance = max_extent * 2.0
      
      # Convert angles to radians for position calculation
      azimuth_rad = math.radians(azimuth)
      elevation_rad = math.radians(elevation)
      
      # Calculate camera position using spherical coordinates
      # Azimuth: 0=+Y (North), 90=-X (East), 180=-Y (South), 270=+X (West)
      # Elevation: 0=horizontal, 90=directly above, -90=directly below
      cam_x = center_x - distance * math.cos(elevation_rad) * math.sin(azimuth_rad)
      cam_y = center_y + distance * math.cos(elevation_rad) * math.cos(azimuth_rad)
      cam_z = center_z + distance * math.sin(elevation_rad)
      
      camera_position = [cam_x, cam_y, cam_z]
      camera_focus = [center_x, center_y, center_z]
      
      # Set up vector (Z up unless looking straight down/up)
      if abs(elevation) > 85: # Near vertical views
        camera_up = [0, 1, 0] if elevation > 0 else [0, -1, 0]
      else:
        camera_up = [0, 0, 1]
      
      self._set_camera_view(render_view, camera_position, camera_focus, camera_up)
      
      # Apply orientation rotations: yaw, pitch, roll
      camera = render_view.GetActiveCamera()
      camera.Yaw(yaw) # Rotate around local vertical axis
      camera.Pitch(pitch) # Rotate around local horizontal axis
      camera.Roll(roll) # Rotate around view direction
      
      print(f"Camera set to azimuth={azimuth}°, elevation={elevation}°, yaw={yaw}°, pitch={pitch}°, roll={roll}°")

  def add_light_source(self, intensity=1.0, azimuth=0.0, elevation=45.0, render_view=None):
    """
    Add an additional light source to the visualization.
    """
    if render_view is None:
      render_view = self.render_view
    
    if render_view is None:
      raise RuntimeError("No render view available. Create a visualization first.")
    
    # Create a new light source
    light = AddLight(render_view)
    if not light:
      raise RuntimeError("Failed to create light source.")
    
    # Set intensity
    light.Intensity = max(min(intensity, 2.0), 0.0)
    
    # Calculate light position using mesh bounds (if available)
    if self.reader:
      bounds = self.reader.GetDataInformation().GetBounds()
      center_x = (bounds[0] + bounds[1]) / 2.0
      center_y = (bounds[2] + bounds[3]) / 2.0
      center_z = (bounds[4] + bounds[5]) / 2.0
      max_extent = max(bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4])
      distance = max_extent * 2.0
    else:
      center_x = center_y = center_z = 0.0
      distance = 1000.0
    
    azimuth_rad = math.radians(azimuth)
    elevation_rad = math.radians(elevation)
    
    light_x = center_x + distance * math.cos(elevation_rad) * math.sin(azimuth_rad)
    light_y = center_y + distance * math.cos(elevation_rad) * math.cos(azimuth_rad)
    light_z = center_z + distance * math.sin(elevation_rad)
    
    light.Position = [light_x, light_y, light_z]
    light.FocalPoint = [center_x, center_y, center_z]
    light.Coords = 'Scene'
    
    render_view.UseLight = 1
    render_view.OrientationAxesVisibility = 0
    render_view.StillRender()
    
    print(f"Added light source: intensity={intensity}, azimuth={azimuth}°, elevation={elevation}°")

  def _apply_scalar_coloring(self, reader, representation, color_by_scalar,
                             render_view, colormap="Rainbow Uniform",
                             reverse_colormap=False):
    """Helper method to apply scalar coloring to a representation."""
    # Get data information
    data_info = reader.GetDataInformation()
    point_data = data_info.GetPointDataInformation()
    cell_data = data_info.GetCellDataInformation()
    
    # Look for the scalar in point data first, then cell data
    scalar_found = False
    for i in range(point_data.GetNumberOfArrays()):
      array_info = point_data.GetArrayInformation(i)
      if array_info.GetName() == color_by_scalar:
        ColorBy(representation, ('POINTS', color_by_scalar))
        scalar_found = True
        break
    
    if not scalar_found:
      for i in range(cell_data.GetNumberOfArrays()):
        array_info = cell_data.GetArrayInformation(i)
        if array_info.GetName() == color_by_scalar:
          ColorBy(representation, ('CELLS', color_by_scalar))
          scalar_found = True
          break
    
    if scalar_found:
      # Show color bar
      try:
        representation.SetScalarBarVisibility(render_view, True)
        
        # Set colormap
        color_map = GetColorTransferFunction(color_by_scalar)
        color_map.ApplyPreset(colormap, True)
        # Update the representation to use the new colormap FIRST
        representation.RescaleTransferFunctionToDataRange(True, False)
        # Reverse colormap if requested (AFTER rescaling)
        if reverse_colormap:
          color_map.InvertTransferFunction()
      except Exception as e:
        print(f"Warning: Could not apply colormap '{colormap}': {e}")
    else:
      print(f"Warning: Scalar field '{color_by_scalar}' not found in the data")

  def save_state_file(self, output_path):
    """
    Save the current ParaView state to a .pvsm file.
    
    Args:
        output_path (str): Path where the .pvsm file will be saved
    """
    output_path = os.path.abspath(output_path)
    
    # Ensure directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Save state
    SaveState(output_path)
    
    print(f"ParaView state file saved: {output_path}")
    return output_path

  def save_screenshot(self, output_path, width=1200, height=900, quality=100):
    """
    Save a screenshot of the current view as PNG.
    
    Args:
        output_path (str): Path where the PNG file will be saved
        width (int): Image width in pixels
        height (int): Image height in pixels
        quality (int): Image quality (0-100, only for JPEG)
    """
    output_path = os.path.abspath(output_path)
    
    # Ensure directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Set view size for screenshot
    if self.render_view:
      self.render_view.ViewSize = [width, height]
      Render()
    
    # Save screenshot
    SaveScreenshot(output_path,
                   self.render_view,
                   ImageResolution=[width, height],
                   TransparentBackground=0)
    
    print(f"Screenshot saved: {output_path}")
    return output_path

  def save_state_and_screenshot(self, pvsm_path, png_path=None,
                                width=1200, height=900, quality=100):
    """
    Save both state file and screenshot.
    
    Args:
        pvsm_path (str): Path for .pvsm file
        png_path (str): Path for .png file (auto-generated if None)
        width (int): Screenshot width
        height (int): Screenshot height
        quality (int): Image quality
    """
    # Save state file
    self.save_state_file(pvsm_path)
    
    # Generate PNG path if not provided
    if png_path is None:
      png_path = os.path.splitext(pvsm_path)[0] + '.png'
    
    # Save screenshot
    self.save_screenshot(png_path, width, height, quality)
    
    return pvsm_path, png_path

  def create_custom_slices(self, vtk_file_path, mesh_name="VTK_Mesh",
                           slice_positions={'x': 0.0, 'y': 0.0, 'z': 0.0},
                           color_by_scalar=None):
    """
    Create three orthogonal slices at custom positions.
    
    Args:
        vtk_file_path (str): Path to the VTK file
        mesh_name (str): Name for the mesh object
        slice_positions (dict): Positions for slices {'x': val, 'y': val, 'z': val}
        color_by_scalar (str): Name of scalar field to color by
        
    Returns:
        tuple: (reader, slice_representations, render_view)
    """
    return self.create_advanced_visualization(
      vtk_file_path=vtk_file_path,
      mesh_name=mesh_name,
      slice_positions=slice_positions,
      color_by_scalar=color_by_scalar
    )

  def create_multi_representation_view(self, vtk_file_path, mesh_name="VTK_Mesh"):
    """
    Create a visualization with both the full mesh and slice representations.
    
    Args:
        vtk_file_path (str): Path to the VTK file
        mesh_name (str): Name for the mesh object
        
    Returns:
        dict: Dictionary containing all created objects
    """
    # Clear existing pipeline
    try:
      sources = GetSources()
      for source in sources.values():
        Delete(source)
    except:
      try:
        ResetSession()
      except:
        pass
    
    # Create reader
    reader = LegacyVTKReader(FileNames=[os.path.abspath(vtk_file_path)])
    reader.UpdatePipeline()
    
    try:
      RenameSource(mesh_name, reader)
    except:
      pass
    
    # Create render view
    render_view = GetActiveViewOrCreate('RenderView')
    render_view.Background = [1.0, 1.0, 1.0] # White background
    render_view.UseColorPaletteForBackground = 0 # Disable color palette for background
    render_view.ViewSize = [1200, 900]
    
    # Store references for camera operations
    self.reader = reader
    self.render_view = render_view
    
    # Get mesh bounds and center
    bounds = reader.GetDataInformation().GetBounds()
    center_x = (bounds[0] + bounds[1]) / 2.0
    center_y = (bounds[2] + bounds[3]) / 2.0
    center_z = (bounds[4] + bounds[5]) / 2.0
    
    # Create wireframe representation of full mesh (semi-transparent)
    full_mesh_repr = Show(reader, render_view)
    full_mesh_repr.Representation = 'Wireframe'
    full_mesh_repr.Opacity = 0.3
    full_mesh_repr.DiffuseColor = [0.5, 0.5, 0.5]
    
    # Create three orthogonal slices
    slice_representations = []
    
    # XY slice
    xy_slice = Slice(Input=reader)
    xy_slice.SliceType = 'Plane'
    xy_slice.SliceType.Origin = [center_x, center_y, center_z]
    xy_slice.SliceType.Normal = [0.0, 0.0, 1.0]
    xy_repr = Show(xy_slice, render_view)
    xy_repr.Representation = 'Surface'
    xy_repr.DiffuseColor = [1.0, 0.0, 0.0] # Red
    slice_representations.append(('XY_slice', xy_slice, xy_repr))
    
    # XZ slice
    xz_slice = Slice(Input=reader)
    xz_slice.SliceType = 'Plane'
    xz_slice.SliceType.Origin = [center_x, center_y, center_z]
    xz_slice.SliceType.Normal = [0.0, 1.0, 0.0]
    xz_repr = Show(xz_slice, render_view)
    xz_repr.Representation = 'Surface'
    xz_repr.DiffuseColor = [0.0, 1.0, 0.0] # Green
    slice_representations.append(('XZ_slice', xz_slice, xz_repr))
    
    # YZ slice
    yz_slice = Slice(Input=reader)
    yz_slice.SliceType = 'Plane'
    yz_slice.SliceType.Origin = [center_x, center_y, center_z]
    yz_slice.SliceType.Normal = [1.0, 0.0, 0.0]
    yz_repr = Show(yz_slice, render_view)
    yz_repr.Representation = 'Surface'
    yz_repr.DiffuseColor = [0.0, 0.0, 1.0] # Blue
    slice_representations.append(('YZ_slice', yz_slice, yz_repr))
    
    # Reset camera and render
    render_view.ResetCamera()
    Render()
    
    return {
        'reader': reader,
        'full_mesh_representation': full_mesh_repr,
        'slice_representations': slice_representations,
        'render_view': render_view
    }

  def add_source_receiver(self, sr_file_path="phases/sr.vtk", sphere_radius=None,
                          sphere_color=[1.0, 1.0, 1.0], render_view=None):
    """
    Add source and receiver points as spheres to the current visualization.
    
    Args:
        sr_file_path (str): Path to the source-receiver VTK file
        sphere_radius (float): Radius of the spheres (auto-calculated if None)
        sphere_color (list): RGB color for spheres [0-1]
        render_view: Render view to add spheres to (uses self.render_view if None)
        
    Returns:
        list: List of (reader, glyph_filter, representation) tuples
    """
    if render_view is None:
      render_view = self.render_view
    
    if render_view is None:
      raise RuntimeError("No render view available. Create a visualization first.")
    
    # Check if sr.vtk file exists
    sr_absolute_path = os.path.abspath(sr_file_path)
    if not os.path.exists(sr_absolute_path):
      print(f"Warning: Source-receiver file '{sr_file_path}' not found, skipping source/receiver visualization")
      return []
    
    sr_objects = []
    
    try:
      # Create VTK reader for source-receiver points
      sr_reader = LegacyVTKReader(FileNames=[sr_absolute_path])
      sr_reader.UpdatePipeline()
      
      # Rename the source
      try:
        RenameSource("Source_Receiver", sr_reader)
      except:
        pass
      
      # Get data bounds to calculate appropriate sphere size if not provided
      if sphere_radius is None:
        if self.reader:
          # Use main mesh bounds for scaling
          bounds = self.reader.GetDataInformation().GetBounds()
          max_extent = max(bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4])
          sphere_radius = max_extent * 0.01 # 1% of mesh extent
        else:
          # Use sr data bounds for scaling
          sr_bounds = sr_reader.GetDataInformation().GetBounds()
          max_extent = max(sr_bounds[1]-sr_bounds[0], sr_bounds[3]-sr_bounds[2], sr_bounds[5]-sr_bounds[4])
          sphere_radius = max_extent * 0.01
      
      print(f"Using sphere radius: {sphere_radius:.4f}")
      
      # Create sphere glyph filter to represent points as spheres
      glyph_filter = Glyph(Input=sr_reader, GlyphType='Sphere')
      glyph_filter.GlyphType.Radius = sphere_radius
      glyph_filter.GlyphType.ThetaResolution = 16 # Reasonable sphere resolution
      glyph_filter.GlyphType.PhiResolution = 16
      glyph_filter.ScaleFactor = 1.0
      glyph_filter.GlyphMode = 'All Points'
      
      # Create representation
      sr_repr = Show(glyph_filter, render_view)
      sr_repr.Representation = 'Surface'
      sr_repr.DiffuseColor = sphere_color
      sr_repr.Opacity = 1.0 # Fully opaque
      
      # Disable specular lighting for cleaner look
      try:
        sr_repr.Specular = 0.0
        sr_repr.Ambient = 0.2 # Add some ambient lighting
      except:
        pass
      
      sr_objects.append((sr_reader, glyph_filter, sr_repr))
      print(f"Added source/receiver points from: {sr_file_path}")
      
    except Exception as e:
      print(f"Error loading source/receiver file '{sr_file_path}': {e}")
      return []
    
    # Update the view
    if sr_objects:
      Render()
    
    return sr_objects
  
  def add_ray_paths(self, ray_phases, tube_radius=0.003, ray_color=[1.0, 1.0, 1.0], render_view=None, phases_dir="phases"):
    """
    Add seismic ray paths as tubes to the current visualization.
    
    Args:
        ray_phases (list): List of seismic phase names (e.g., ['P', 'S', 'PKIKP'])
        tube_radius (float): Radius of the ray tubes (relative to model size)
        ray_color (list): RGB color for ray tubes [0-1]
        render_view: Render view to add rays to (uses self.render_view if None)
        
    Returns:
        list: List of (phase_name, reader, tube_filter, representation) tuples
    """
    if render_view is None:
      render_view = self.render_view
    
    if render_view is None:
      raise RuntimeError("No render view available. Create a visualization first.")
    
    ray_objects = []
    
    # Get mesh bounds to scale tube radius appropriately
    if self.reader:
      bounds = self.reader.GetDataInformation().GetBounds()
      max_extent = max(bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4])
      actual_tube_radius = tube_radius * max_extent
    else:
      actual_tube_radius = tube_radius * 1000.0 # Default scaling
    
    for phase in ray_phases:
      vtk_file_path = f"{phases_dir}/{phase}.vtk"
      
      # Check if file exists
      if not os.path.exists(vtk_file_path):
        print(f"Warning: Ray file '{vtk_file_path}' not found, skipping {phase}")
        continue
      
      try:
        # Create VTK reader for ray path
        ray_reader = LegacyVTKReader(FileNames=[os.path.abspath(vtk_file_path)])
        ray_reader.UpdatePipeline()
        
        # Rename the source
        try:
          RenameSource(f"Ray_{phase}", ray_reader)
        except:
          pass
        
        # Create tube filter to make rays visible as tubes
        tube_filter = Tube(Input=ray_reader)
        tube_filter.Radius = actual_tube_radius
        tube_filter.NumberofSides = 8 # Octagonal tubes for efficiency
        tube_filter.Capping = 1 # Cap the ends
        
        # Create representation
        ray_repr = Show(tube_filter, render_view)
        ray_repr.Representation = 'Surface'
        ray_repr.DiffuseColor = ray_color
        ray_repr.Opacity = 0.8 # Slightly transparent
        
        # Disable specular lighting for cleaner look
        try:
          ray_repr.Specular = 0.0
        except:
          pass
        
        ray_objects.append((phase, ray_reader, tube_filter, ray_repr))
        print(f"Added ray path: {phase}")
        
      except Exception as e:
        print(f"Error loading ray path '{phase}': {e}")
        continue
    
    # Update the view
    if ray_objects:
      Render()
    
    return ray_objects

  def create_visualization_with_rays(self, vtk_file_path, ray_phases,
                                     mesh_name="VTK_Mesh",
                                     color_by_scalar=None,
                                     reverse_colormap=False,
                                     show_edges=False,
                                     opacity=1.0,
                                     slice_positions=None,
                                     colormap="Rainbow Uniform",
                                     font_color=[0.0, 0.0, 0.0],
                                     font_size=25,
                                     label=None,
                                     tube_radius=0.003,
                                     ray_color=[1.0, 1.0, 1.0],
                                     add_source_receiver=True,
                                     sr_file_path="phases/sr.vtk",
                                     sphere_color=[1.0, 1.0, 1.0],
                                     phases_dir="phases",
                                     **kwargs):
    """
    Create a complete visualization with mesh slices, ray paths, and source/receiver points.
    
    Args:
        vtk_file_path (str): Path to the main VTK mesh file
        ray_phases (list): List of seismic phase names for ray paths
        mesh_name (str): Name for the mesh object
        color_by_scalar (str): Name of scalar field to color mesh by
        reverse_colormap (bool): Whether to reverse the colormap
        show_edges (bool): Whether to show mesh edges
        opacity (float): Opacity value (0.0 to 1.0)
        slice_positions (dict): Positions for slices {'x': val, 'y': val, 'z': val}
        colormap (str): Colormap to use for scalar data
        font_color (list): RGB color for colorbar label [0-1]
        font_size (int): Font size for colorbar label
        label (str): Text for colorbar label (defaults to color_by_scalar or 'VTK Mesh')
        tube_radius (float): Radius of ray tubes (relative to model size)
        ray_color (list): RGB color for ray tubes [0-1]
        add_source_receiver (bool): Whether to add source/receiver spheres
        sr_file_path (str): Path to source-receiver VTK file
        sphere_color (list): RGB color for source/receiver spheres [0-1]
        **kwargs: Additional arguments passed to create_advanced_visualization
        
    Returns:
        dict: Dictionary containing all created objects
    """
    # Set default label to color_by_scalar or 'VTK Mesh'
    if label is None:
      label = color_by_scalar if color_by_scalar else "VTK Mesh"
    
    # Create the basic mesh visualization first
    reader, slice_representations, render_view = self.create_advanced_visualization(
      vtk_file_path=vtk_file_path,
      mesh_name=mesh_name,
      color_by_scalar=color_by_scalar,
      reverse_colormap=reverse_colormap,
      show_edges=show_edges,
      opacity=opacity,
      slice_positions=slice_positions,
      colormap=colormap,
      font_color=font_color,
      font_size=font_size,
      label=label,
      **kwargs # Pass representation_type and other kwargs
    )
    
    # Add ray paths
    ray_objects = self.add_ray_paths(ray_phases, tube_radius, ray_color, render_view, phases_dir=phases_dir)
    
    # Add source/receiver points if requested
    sr_objects = []
    if add_source_receiver:
      # Calculate sphere radius as twice the tube radius (scaled to model size)
      if self.reader:
        bounds = self.reader.GetDataInformation().GetBounds()
        max_extent = max(bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4])
        sphere_radius = 4.0 * tube_radius * max_extent
      else:
        sphere_radius = 4.0 * tube_radius * 1000.0
      
      sr_objects = self.add_source_receiver(sr_file_path, sphere_radius, sphere_color, render_view)
    
    return {
      'reader': reader,
      'slice_representations': slice_representations,
      'ray_objects': ray_objects,
      'sr_objects': sr_objects,
      'render_view': render_view
    }

  def set_isometric_view(self, render_view=None):
    """Set an isometric view of the mesh."""
    if render_view is None:
      render_view = self.render_view
      
    # Get data bounds to position camera appropriately
    if self.reader:
      bounds = self.reader.GetDataInformation().GetBounds()
      center_x = (bounds[0] + bounds[1]) / 2.0
      center_y = (bounds[2] + bounds[3]) / 2.0
      center_z = (bounds[4] + bounds[5]) / 2.0
      
      # Calculate distance for isometric view
      max_extent = max(bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4])
      distance = max_extent * 2.0
      
      # Isometric view position (equal angles)
      iso_pos = [center_x + distance, center_y + distance, center_z + distance]
      
      self._set_camera_view(render_view, iso_pos, [center_x, center_y, center_z], [0, 0, 1])

  def set_front_view(self, render_view=None):
    """Set front view (looking down negative Y axis)."""
    if render_view is None:
      render_view = self.render_view
      
    if self.reader:
      bounds = self.reader.GetDataInformation().GetBounds()
      center_x = (bounds[0] + bounds[1]) / 2.0
      center_y = (bounds[2] + bounds[3]) / 2.0
      center_z = (bounds[4] + bounds[5]) / 2.0
      
      max_extent = max(bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4])
      distance = max_extent * 2.0
      
      front_pos = [center_x, center_y - distance, center_z]
      self._set_camera_view(render_view, front_pos, [center_x, center_y, center_z], [0, 0, 1])

  def set_top_view(self, render_view=None):
    """Set top view (looking down negative Z axis)."""
    if render_view is None:
      render_view = self.render_view
      
    if self.reader:
      bounds = self.reader.GetDataInformation().GetBounds()
      center_x = (bounds[0] + bounds[1]) / 2.0
      center_y = (bounds[2] + bounds[3]) / 2.0
      center_z = (bounds[4] + bounds[5]) / 2.0
      
      max_extent = max(bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4])
      distance = max_extent * 2.0
      
      top_pos = [center_x, center_y, center_z + distance]
      self._set_camera_view(render_view, top_pos, [center_x, center_y, center_z], [0, 1, 0])

  def get_available_colormaps(self):
    """Get list of available colormaps in ParaView."""
    try:
      # Common ParaView colormaps
      colormaps = [
          "Rainbow Uniform",
          "Cool to Warm",
          "Cool to Warm (Extended)",
          "Viridis",
          "Plasma",
          "Inferno",
          "Magma",
          "jet",
          "hsv",
          "hot",
          "cool",
          "spring",
          "summer",
          "autumn",
          "winter",
          "bone",
          "copper",
          "pink",
          "lines",
          "colorcube",
          "prism",
          "flag",
          "Blue to Red Rainbow",
          "Red to Blue Rainbow",
          "HSV",
          "Diverging (Blue-Red)",
          "Spectral",
          "RdYlBu",
          "Turbo"
      ]
      return colormaps
    except:
      return ["Rainbow Uniform", "Cool to Warm", "Viridis", "jet"]

  def get_available_scalars(self, vtk_file_path):
    """
    Get list of available scalar fields in the VTK file.
    
    Args:
        vtk_file_path (str): Path to the VTK file
        
    Returns:
        dict: Dictionary with 'point_data' and 'cell_data' scalar arrays
    """
    # Create temporary reader
    temp_reader = LegacyVTKReader(FileNames=[os.path.abspath(vtk_file_path)])
    temp_reader.UpdatePipeline()
    
    # Get data information
    data_info = temp_reader.GetDataInformation()
    point_data = data_info.GetPointDataInformation()
    cell_data = data_info.GetCellDataInformation()
    
    # Extract array names
    point_arrays = []
    for i in range(point_data.GetNumberOfArrays()):
      array_info = point_data.GetArrayInformation(i)
      point_arrays.append(array_info.GetName())
    
    cell_arrays = []
    for i in range(cell_data.GetNumberOfArrays()):
      array_info = cell_data.GetArrayInformation(i)
      cell_arrays.append(array_info.GetName())
    
    # Clean up
    try:
      Delete(temp_reader)
    except:
      # If Delete fails, the reader will be cleaned up automatically
      pass
    
    return {
      'point_data': point_arrays,
      'cell_data': cell_arrays
    }


def apply_camera_view(visualizer, camera_view,
                      camera_azimuth=None,
                      camera_elevation=None,
                      camera_yaw=None,
                      camera_pitch=None,
                      camera_roll=None,
                      light_intensity=None,
                      light_azimuth=None,
                      light_elevation=None):
  """Apply camera view and optional light source based on arguments."""
  # Handle camera view logic (unchanged)
  camera_applied = False
  if any(param is not None for param in [camera_azimuth, camera_elevation, camera_yaw, camera_pitch, camera_roll]):
    # Use defaults for missing values
    azimuth = camera_azimuth if camera_azimuth is not None else 0.0
    elevation = camera_elevation if camera_elevation is not None else 0.0
    yaw = camera_yaw if camera_yaw is not None else 0.0
    pitch = camera_pitch if camera_pitch is not None else 0.0
    roll = camera_roll if camera_roll is not None else 0.0
    visualizer.set_camera_view(azimuth, elevation, yaw, pitch, roll)
    camera_applied = True
  elif camera_view:
    if ',' in str(camera_view):
      try:
        parts = camera_view.split(',')
        azimuth = float(parts[0].strip())
        elevation = float(parts[1].strip())
        yaw = float(parts[2].strip()) if len(parts) > 2 else 0.0
        pitch = float(parts[3].strip()) if len(parts) > 3 else 0.0
        roll = float(parts[4].strip()) if len(parts) > 4 else 0.0
        visualizer.set_camera_view(azimuth, elevation, yaw, pitch, roll)
        camera_applied = True
      except (ValueError, IndexError):
        print(f"Error: Invalid camera view format '{camera_view}'. Use 'azimuth,elevation[,yaw,pitch,roll]' or preset name.")
    elif camera_view == "isometric":
      visualizer.set_isometric_view()
      camera_applied = True
    elif camera_view == "front":
      visualizer.set_front_view()
      camera_applied = True
    elif camera_view == "top":
      visualizer.set_top_view()
      camera_applied = True
    elif camera_view != "auto":
      print(f"Warning: Unknown camera view preset '{camera_view}', using default.")
  
  # Add light source if any light parameters are specified
  if any(param is not None for param in [light_intensity, light_azimuth, light_elevation]):
    intensity = max(min(light_intensity if light_intensity is not None else 1.0, 2.0), 0.0)
    azimuth = light_azimuth if light_azimuth is not None else 0.0
    elevation = max(min(light_elevation if light_elevation is not None else 45.0, 90.0), -90.0)
    visualizer.add_light_source(intensity, azimuth, elevation)
  
  return camera_applied

def create_ray_only_visualization(ray_phases, output_pvsm_path,
                                  tube_radius=0.003,
                                  ray_color=[1.0, 1.0, 1.0],
                                  add_source_receiver=True,
                                  sr_file_path="phases/sr.vtk",
                                  sphere_color=[1.0, 1.0, 1.0],
                                  save_png=True,
                                  camera_azimuth=None,
                                  camera_elevation=None,
                                  camera_yaw=None,
                                  camera_pitch=None,
                                  camera_roll=None,
                                  light_intensity=None,
                                  light_azimuth=None,
                                  light_elevation=None,
                                  phases_dir="phases",
                                  **kwargs):
  """
  Create a visualization with only ray paths (no mesh) and optionally source/receiver points.
  
  Args:
      ray_phases (list): List of seismic phase names for ray paths
      output_pvsm_path (str): Path where the .pvsm file will be saved
      tube_radius (float): Radius of ray tubes
      ray_color (list): RGB color for ray tubes [0-1]
      add_source_receiver (bool): Whether to add source/receiver spheres
      sr_file_path (str): Path to source-receiver VTK file
      sphere_color (list): RGB color for source/receiver spheres [0-1]
      save_png (bool): Whether to also save a PNG screenshot
      camera_azimuth (float): Camera azimuth angle
      camera_elevation (float): Camera elevation angle
      camera_yaw (float): Camera yaw angle
      camera_pitch (float): Camera pitch angle
      camera_roll (float): Camera roll angle
      light_intensity (float): Light intensity
      light_azimuth (float): Light azimuth angle
      light_elevation (float): Light elevation angle
      **kwargs: Additional keyword arguments for visualization settings
  """
  visualizer = ParaViewVTKVisualizer()
  # Clear any existing pipeline
  try:
    sources = GetSources()
    for source in sources.values():
      Delete(source)
  except:
    try:
      ResetSession()
    except:
      pass
  
  # Create render view
  render_view = GetActiveViewOrCreate('RenderView')
  render_view.Background = kwargs.get('background_color', [1.0, 1.0, 1.0])
  render_view.UseColorPaletteForBackground = 0
  render_view.ViewSize = kwargs.get('view_size', [1200, 900])
  visualizer.render_view = render_view
  
  # Estimate tube radius based on first available ray path
  if not ray_phases:
    raise ValueError("No ray phases specified")
  
  # Find first available ray file to get scale
  scale_factor = 1000.0 # Default
  for phase in ray_phases:
    vtk_file_path = f"{phases_dir}/{phase}.vtk"
    if os.path.exists(vtk_file_path):
      try:
        temp_reader = LegacyVTKReader(FileNames=[os.path.abspath(vtk_file_path)])
        temp_reader.UpdatePipeline()
        bounds = temp_reader.GetDataInformation().GetBounds()
        max_extent = max(bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4])
        scale_factor = max_extent
        # Set reader for camera operations
        visualizer.reader = temp_reader
        Delete(temp_reader)
        break
      except:
        continue
  
  # Add ray paths with proper scaling
  actual_tube_radius = tube_radius * scale_factor
  ray_objects = []
  
  for phase in ray_phases:
    vtk_file_path = f"{phases_dir}/{phase}.vtk"
    
    if not os.path.exists(vtk_file_path):
      print(f"Warning: Ray file '{vtk_file_path}' not found, skipping {phase}")
      continue
    
    try:
      # Create VTK reader for ray path
      ray_reader = LegacyVTKReader(FileNames=[os.path.abspath(vtk_file_path)])
      ray_reader.UpdatePipeline()
      
      # Use first ray as reference for camera if we don't have a mesh reader
      if not visualizer.reader:
        visualizer.reader = ray_reader
      
      # Rename the source
      try:
        RenameSource(f"Ray_{phase}", ray_reader)
      except:
        pass
      
      # Create tube filter
      tube_filter = Tube(Input=ray_reader)
      tube_filter.Radius = actual_tube_radius
      tube_filter.NumberofSides = 8
      tube_filter.Capping = 1
      
      # Create representation
      ray_repr = Show(tube_filter, render_view)
      ray_repr.Representation = 'Surface'
      ray_repr.DiffuseColor = ray_color
      ray_repr.Opacity = 0.8 # Slightly transparent
      
      # Disable specular lighting for cleaner look
      try:
        ray_repr.Specular = 0.0
      except:
        pass
      
      ray_objects.append((phase, ray_reader, tube_filter, ray_repr))
      print(f"Added ray path: {phase}")
      
    except Exception as e:
      print(f"Error loading ray path '{phase}': {e}")
      continue
  
  if not ray_objects:
    raise RuntimeError("No ray paths could be loaded")
  
  # Add source/receiver points if requested
  sr_objects = []
  if add_source_receiver:
    # Calculate sphere radius as twice the tube radius (scaled)
    sphere_radius = 4.0 * actual_tube_radius
    sr_objects = visualizer.add_source_receiver(sr_file_path, sphere_radius, sphere_color, render_view)
  
  # Apply camera view if specified
  camera_applied = apply_camera_view(visualizer, None, camera_azimuth, camera_elevation,
                                     camera_yaw, camera_pitch, camera_roll,
                                     light_intensity, light_azimuth, light_elevation)
  if not camera_applied:
    render_view.ResetCamera()
  
  Render()
  
  # Save state file and optionally PNG
  if save_png:
    png_path = os.path.splitext(output_pvsm_path)[0] + '.png'
    visualizer.save_state_and_screenshot(output_pvsm_path, png_path)
    print(f"Files created:")
    print(f" State file: {output_pvsm_path}")
    print(f" Screenshot: {png_path}")
  else:
    visualizer.save_state_file(output_pvsm_path)
  
  print(f"Ray visualization created with phases: {ray_phases}")
  if add_source_receiver and sr_objects:
    print(f"Source/receiver points added from: {sr_file_path}")
  
  return {'ray_objects': ray_objects, 'sr_objects': sr_objects, 'render_view': render_view}

def create_mesh_and_rays_visualization(vtk_file_path, ray_phases, output_pvsm_path,
                                       mesh_name="VTK_Mesh",
                                       color_by_scalar=None,
                                       reverse_colormap=False,
                                       tube_radius=0.003,
                                       ray_color=[1.0, 1.0, 1.0],
                                       add_source_receiver=True,
                                       sr_file_path="phases/sr.vtk",
                                       sphere_color=[1.0, 1.0, 1.0],
                                       save_png=True,
                                       camera_azimuth=None,
                                       camera_elevation=None,
                                       camera_yaw=None,
                                       camera_pitch=None,
                                       camera_roll=None,
                                       light_intensity=None,
                                       light_azimuth=None,
                                       light_elevation=None,
                                       font_color=[0.0, 0.0, 0.0],
                                       font_size=25,
                                       label=None,
                                       representation_type="Surface",
                                       phases_dir="phases",
                                       **kwargs):
  """
  Create a visualization with both mesh slices, ray paths, and source/receiver points.
  
  Args:
      vtk_file_path (str): Path to the main VTK mesh file
      ray_phases (list): List of seismic phase names for ray paths
      output_pvsm_path (str): Path where the .pvsm file will be saved
      mesh_name (str): Name for the mesh object
      color_by_scalar (str): Name of scalar field to color mesh by
      reverse_colormap (bool): Whether to reverse the colormap
      tube_radius (float): Radius of ray tubes (relative to model size)
      ray_color (list): RGB color for ray tubes [0-1]
      add_source_receiver (bool): Whether to add source/receiver spheres
      sr_file_path (str): Path to source-receiver VTK file
      sphere_color (list): RGB color for source/receiver spheres [0-1]
      save_png (bool): Whether to also save a PNG screenshot
      camera_azimuth (float): Camera azimuth angle
      camera_elevation (float): Camera elevation angle
      camera_yaw (float): Camera yaw angle
      camera_pitch (float): Camera pitch angle
      camera_roll (float): Camera roll angle
      light_intensity (float): Light intensity
      light_azimuth (float): Light azimuth angle
      light_elevation (float): Light elevation angle
      font_color (list): RGB color for colorbar label [0-1]
      font_size (int): Font size for colorbar label
      label (str): Text for colorbar label (defaults to color_by_scalar or 'VTK Mesh')
      representation_type (str): Representation type for slices (default: 'Surface')
      **kwargs: Additional arguments for mesh visualization
  """
  visualizer = ParaViewVTKVisualizer()
  # Set default view size
  width, height = kwargs.get('view_size', [1200, 900])
  # Set default label to color_by_scalar or 'VTK Mesh'
  if label is None:
    label = color_by_scalar if color_by_scalar else "VTK Mesh"
  
  # Create complete visualization
  result = visualizer.create_visualization_with_rays(
    vtk_file_path=vtk_file_path,
    ray_phases=ray_phases,
    mesh_name=mesh_name,
    color_by_scalar=color_by_scalar,
    reverse_colormap=reverse_colormap,
    tube_radius=tube_radius,
    ray_color=ray_color,
    add_source_receiver=add_source_receiver,
    sr_file_path=sr_file_path,
    sphere_color=sphere_color,
    font_color=font_color,
    font_size=font_size,
    label=label,
    representation_type=representation_type,
    phases_dir=phases_dir,
    **kwargs
  )
  # Apply camera view if specified
  camera_applied = apply_camera_view(visualizer, None, camera_azimuth, camera_elevation,
                                     camera_yaw, camera_pitch, camera_roll,
                                     light_intensity, light_azimuth, light_elevation)
  if not camera_applied:
    visualizer.render_view.ResetCamera()
  # Render after camera changes
  Render()
  # Save state file and optionally PNG
  if save_png:
    png_path = os.path.splitext(output_pvsm_path)[0] + '.png'
    visualizer.save_state_and_screenshot(output_pvsm_path, png_path, width, height)
    print(f"Files created:")
    print(f" State file: {output_pvsm_path}")
    print(f" Screenshot: {png_path}")
  else:
    visualizer.save_state_file(output_pvsm_path)
  
  print(f"Mesh and ray visualization created:")
  print(f" Mesh file: {os.path.abspath(vtk_file_path)}")
  print(f" Ray phases: {ray_phases}")
  if color_by_scalar:
    print(f" Colored by: {color_by_scalar}")
    print(f" Colorbar label: {label}")
  if add_source_receiver and result.get('sr_objects'):
    print(f" Source/receiver points added from: {sr_file_path}")
  
  # Prevent Segmentation fault when exiting pvpython 6.1
  try:
    Disconnect()
    ResetSession()
  except:
    pass

  return result

def create_pvsm_file(vtk_file_path, output_pvsm_path,
                     mesh_name="VTK_Mesh",
                     save_png=True,
                     color_by_scalar=None,
                     camera_azimuth=None,
                     camera_elevation=None,
                     camera_yaw=None,
                     camera_pitch=None,
                     camera_roll=None,
                     light_intensity=None,
                     light_azimuth=None,
                     light_elevation=None,
                     font_color=[0.0, 0.0, 0.0],
                     font_size=25,
                     label=None):
  """
  Create a basic ParaView state file for VTK mesh visualization.
  
  Args:
      vtk_file_path (str): Path to the VTK file
      output_pvsm_path (str): Path where the .pvsm file will be saved
      mesh_name (str): Name for the mesh object
      save_png (bool): Whether to also save a PNG screenshot
      color_by_scalar (str): Name of scalar field to color by
      camera_azimuth (float): Camera azimuth angle
      camera_elevation (float): Camera elevation angle
      camera_yaw (float): Camera yaw angle
      camera_pitch (float): Camera pitch angle
      camera_roll (float): Camera roll angle
      light_intensity (float): Light intensity
      light_azimuth (float): Light azimuth angle
      light_elevation (float): Light elevation angle
      font_color (list): RGB color for colorbar label [0-1]
      font_size (int): Font size for colorbar label
      label (str): Text for colorbar label (defaults to color_by_scalar or 'VTK Mesh')
  """
  visualizer = ParaViewVTKVisualizer()
  
  # Set default label to color_by_scalar or 'VTK Mesh'
  if label is None:
    label = color_by_scalar if color_by_scalar else "VTK Mesh"
  
  # Create visualization
  reader, slice_representations, render_view = visualizer.create_basic_visualization(
    vtk_file_path, mesh_name, color_by_scalar, font_color=font_color, font_size=font_size, label=label)
  
  # Apply camera view if specified
  camera_applied = apply_camera_view(visualizer, None, camera_azimuth, camera_elevation,
                                     camera_yaw, camera_pitch, camera_roll,
                                     light_intensity, light_azimuth, light_elevation)
  if not camera_applied:
    render_view.ResetCamera()
  Render()
  # Save state file and optionally PNG
  if save_png:
    png_path = os.path.splitext(output_pvsm_path)[0] + '.png'
    pvsm_path, png_path = visualizer.save_state_and_screenshot(output_pvsm_path, png_path)
    print(f"Files created:")
    print(f" State file: {pvsm_path}")
    print(f" Screenshot: {png_path}")
  else:
    visualizer.save_state_file(output_pvsm_path)
  
  print(f"Basic visualization created for: {os.path.abspath(vtk_file_path)}")
  if color_by_scalar:
    print(f"Colorbar label: {label}")
  
  return {'reader': reader, 'slice_representations': slice_representations, 'render_view': render_view}

def create_advanced_pvsm_file(vtk_file_path, output_pvsm_path,
                              mesh_name="VTK_Mesh",
                              representation="Surface",
                              color_by_scalar=None,
                              reverse_colormap=False,
                              show_edges=False,
                              opacity=1.0,
                              save_png=True,
                              camera_azimuth=None,
                              camera_elevation=None,
                              camera_yaw=None,
                              camera_pitch=None,
                              camera_roll=None,
                              light_intensity=None,
                              light_azimuth=None,
                              light_elevation=None,
                              font_color=[0.0, 0.0, 0.0],
                              font_size=25,
                              label=None,
                              **kwargs):
  """
  Create an advanced ParaView state file with custom visualization settings.
  
  Args:
      vtk_file_path (str): Path to the VTK file
      output_pvsm_path (str): Path where the .pvsm file will be saved
      mesh_name (str): Name for the mesh object
      representation (str): Representation type
      color_by_scalar (str): Name of scalar field to color by
      reverse_colormap (bool): Whether to reverse the colormap
      show_edges (bool): Whether to show mesh edges
      opacity (float): Opacity value (0.0 to 1.0)
      save_png (bool): Whether to also save a PNG screenshot
      camera_azimuth (float): Camera azimuth angle
      camera_elevation (float): Camera elevation angle
      camera_yaw (float): Camera yaw angle
      camera_pitch (float): Camera pitch angle
      camera_roll (float): Camera roll angle
      light_intensity (float): Light intensity
      light_azimuth (float): Light azimuth angle
      light_elevation (float): Light elevation angle
      font_color (list): RGB color for colorbar label [0-1]
      font_size (int): Font size for colorbar label
      label (str): Text for colorbar label (defaults to color_by_scalar or 'VTK Mesh')
      **kwargs: Additional keyword arguments passed to create_advanced_visualization
  """
  visualizer = ParaViewVTKVisualizer()
  
  # Set default label to color_by_scalar or 'VTK Mesh'
  if label is None:
    label = color_by_scalar if color_by_scalar else "VTK Mesh"
  
  # Create advanced visualization with all keyword arguments
  reader, slice_representations, render_view = visualizer.create_advanced_visualization(
    vtk_file_path=vtk_file_path,
    mesh_name=mesh_name,
    representation_type=representation,
    color_by_scalar=color_by_scalar,
    reverse_colormap=reverse_colormap,
    show_edges=show_edges,
    opacity=opacity,
    font_color=font_color,
    font_size=font_size,
    label=label,
    **kwargs
  )
  # Apply camera view if specified
  camera_applied = apply_camera_view(visualizer, None, camera_azimuth, camera_elevation,
                                     camera_yaw, camera_pitch, camera_roll,
                                     light_intensity, light_azimuth, light_elevation)
  if not camera_applied:
    render_view.ResetCamera()
  Render()
  # Save state file and optionally PNG
  if save_png:
    png_path = os.path.splitext(output_pvsm_path)[0] + '.png'
    pvsm_path, png_path = visualizer.save_state_and_screenshot(output_pvsm_path, png_path)
    print(f"Files created:")
    print(f" State file: {pvsm_path}")
    print(f" Screenshot: {png_path}")
  else:
    visualizer.save_state_file(output_pvsm_path)
  
  print(f"Advanced visualization created:")
  print(f" File: {os.path.abspath(vtk_file_path)}")
  print(f" Representation: {representation}")
  print(f" Opacity: {opacity}")
  print(f" Show edges: {show_edges}")
  if color_by_scalar:
    print(f" Colored by: {color_by_scalar}")
    print(f" Colorbar label: {label}")
  
  return {'reader': reader, 'slice_representations': slice_representations, 'render_view': render_view}

def main():
  """Command line interface for the ParaView PVSM generator."""
  parser = argparse.ArgumentParser(
      description="Generate ParaView state files for VTK meshes using ParaView Python API",
      epilog="Examples:\n"
             " # Basic mesh visualization\n"
             " python pvsm_generator.py mesh.vtk\n"
             " python pvsm_generator.py mesh.vtk -o viz.pvsm --representation Wireframe\n"
             " python pvsm_generator.py mesh.vtk --color-by pressure --edges\n"
             "\n"
             " # Ray paths only (with source/receiver)\n"
             " python pvsm_generator.py --rays-only --phases P S PKIKP\n"
             " python pvsm_generator.py --rays-only --phases P S --tube-radius 0.005 --ray-color 0.8 0.2 0.2\n"
             " python pvsm_generator.py --rays-only --phases P S --no-source-receiver\n"
             " python pvsm_generator.py --rays-only --phases-dir phases_EQ001_STA01 --phases P S PKIKP\n"
             "\n"
             " # Mesh with ray paths and source/receiver\n"
             " python pvsm_generator.py mesh.vtk --phases P S PKIKP SKS\n"
             " python pvsm_generator.py mesh.vtv --phases P S --color-by Vp --tube-radius 0.001\n"
             " python pvsm_generator.py mesh.vtv --phases P S --sphere-color 0.8 0.0 0.0\n"
             "\n"
             " # Camera control examples\n"
             " python pvsm_generator.py mesh.vtk --camera-view isometric\n"
             " python pvsm_generator.py mesh.vtk --camera-view 210,10,0,0,20 # azimuth,elevation,yaw,pitch,roll\n"
             " python pvsm_generator.py mesh.vtk --camera-azimuth 50 --camera-elevation 30 --camera-yaw 90 --camera-pitch 45 --camera-roll 15\n"
             " python pvsm_generator.py mesh.vtk --camera-view 180,-15,0,0,0 # South view, 15° below horizon\n"
             "\n"
             " # Light source examples\n"
             " python pvsm_generator.py mesh.vtk --light-intensity 0.5 --light-azimuth 140 --light-elevation 35\n"
             " python pvsm_generator.py mesh.vtv --phases P S --light-intensity 0.8 --light-azimuth 270 --light-elevation 30\n"
             "\n"
             " # Colorbar customization examples\n"
             " python pvsm_generator.py mesh.vtk --color-by Vp --font-color 1.0 0.0 0.0 --font-size 24 --label 'Seismic Velocity'\n"
             " python pvsm_generator.py mesh.vtk --color-by Vp --font-color 0.0 0.0 1.0\n",
      formatter_class=argparse.RawDescriptionHelpFormatter
  )
  
  parser.add_argument("vtk_file", nargs='?', help="Path to the VTK file (optional if using --rays-only)")
  parser.add_argument("-o", "--output", help="Output PVSM file path",
                      default="Mesh_Visualization.pvsm")
  parser.add_argument("-n", "--name", help="Name for the mesh", default="VTK_Mesh")
  parser.add_argument("-r", "--representation",
                      choices=["Surface", "Wireframe", "Points", "Surface With Edges"],
                      default="Surface", help="Representation type")
  parser.add_argument("--color-by", dest="color_by_scalar", default="Vp",
                      help="Color by scalar field (use --list-scalars to see available)")
  parser.add_argument("--edges", action="store_true", help="Show mesh edges")
  parser.add_argument("--opacity", type=float, default=1.0,
                      help="Opacity (0.0 to 1.0)")
  parser.add_argument("--list-scalars", action="store_true",
                      help="List available scalar fields and exit")
  parser.add_argument("--advanced", action="store_true",
                      help="Use advanced configuration options")
  # Ray path arguments
  parser.add_argument("--phases", "--rays", nargs='+',
                      help="List of seismic phases to plot as ray paths (e.g., P S PKIKP SKS)")
  parser.add_argument("--rays-only", action="store_true",
                      help="Plot only ray paths without mesh")
  parser.add_argument("--tube-radius", type=float, default=0.003,
                      help="Radius of ray tubes relative to model size (default: 0.003)")
  parser.add_argument("--ray-color", nargs=3, type=float,
                      default=[1.0, 1.0, 1.0], metavar=('R', 'G', 'B'),
                      help="RGB color for ray tubes (0.0-1.0, default: 1.0 1.0 1.0)")
  # Source/receiver arguments
  parser.add_argument("--sr-file", default="phases/sr.vtk",
                      help="Path to source-receiver VTK file (default: phases/sr.vtk)")
  parser.add_argument("--sphere-color", nargs=3, type=float,
                      default=[1.0, 1.0, 1.0], metavar=('R', 'G', 'B'),
                      help="RGB color for source/receiver spheres (0.0-1.0, default: 1.0 1.0 1.0)")
  parser.add_argument("--no-source-receiver", action="store_true",
                      help="Don't add source/receiver spheres when using --phases")
  parser.add_argument("--phases-dir", default="phases",
                      help="Directory containing phase VTK files (default: phases)")
  
  # Colorbar label arguments
  parser.add_argument("--font-color", nargs=3, type=float,
                      default=[0.0, 0.0, 0.0], metavar=('R', 'G', 'B'),
                      help="RGB color for colorbar label (0.0-1.0, default: 0.0 0.0 0.0)")
  parser.add_argument("--font-size", type=int, default=25,
                      help="Font size for colorbar label (default: 25)")
  parser.add_argument("--label", help="Text for colorbar label (defaults to scalar field or 'VTK Mesh')")
  
  # PNG arguments
  parser.add_argument("--no-png", action="store_true",
                      help="Don't save PNG screenshot (only save .pvsm file)")
  parser.add_argument("--png-size", nargs=2, type=int, default=[1200, 900],
                      metavar=('WIDTH', 'HEIGHT'),
                      help="PNG screenshot size in pixels (default: 1200 900)")
  parser.add_argument("--png-only", action="store_true",
                      help="Only save PNG screenshot (don't save .pvsm file)")
  
  # Slice arguments
  parser.add_argument("--slice-positions", nargs=3, type=float,
                      metavar=('X', 'Y', 'Z'),
                      help="Custom slice positions (default: mesh center)")
  parser.add_argument("--show-full-mesh", action="store_true",
                      help="Show wireframe of full mesh along with slices")
  
  # Colormap arguments
  parser.add_argument("--list-colormaps", action="store_true",
                      help="List available colormaps and exit")
  parser.add_argument("--colormap", default="Rainbow Uniform",
                      help="Colormap to use for scalar data (default: Rainbow Uniform)")
  parser.add_argument("--reverse-colormap", action="store_true",
                      help="Reverse the colormap")
  
  # Camera arguments
  parser.add_argument("--camera-view",
                      help="Camera view: preset name (isometric/front/top/auto) or 'azimuth,elevation,yaw,pitch,roll' in degrees (default: isometric)")
  parser.add_argument("--camera-azimuth", type=float,
                      help="Camera azimuth angle in degrees (0-360)")
  parser.add_argument("--camera-elevation", type=float,
                      help="Camera elevation angle in degrees (-90 to 90)")
  parser.add_argument("--camera-yaw", type=float,
                      help="Camera yaw angle in degrees (0-360)")
  parser.add_argument("--camera-pitch", type=float,
                      help="Camera pitch angle in degrees (-90 to 90)")
  parser.add_argument("--camera-roll", type=float,
                      help="Camera roll angle in degrees (-180 to 180)")
  
  # Light source arguments
  parser.add_argument("--light-intensity", type=float, default=1.0,
                      help="Additional light intensity (0.0-2.0, default: 1.0)")
  parser.add_argument("--light-azimuth", type=float,
                      help="Additional light azimuth angle in degrees (0-360)")
  parser.add_argument("--light-elevation", type=float,
                      help="Additional light elevation angle in degrees (-90 to 90)")
  
  args = parser.parse_args()
  
  # Validate input arguments
  if not args.rays_only and not args.vtk_file:
    print("Error: VTK file is required unless using --rays-only")
    return 1
    
  if args.vtk_file and not os.path.exists(args.vtk_file):
    print(f"Error: VTK file '{args.vtk_file}' does not exist")
    return 1
  
  # Validate ray color
  if any(c < 0.0 or c > 1.0 for c in args.ray_color):
    print("Error: Ray color values must be between 0.0 and 1.0")
    return 1
  
  # Validate sphere color
  if any(c < 0.0 or c > 1.0 for c in args.sphere_color):
    print("Error: Sphere color values must be between 0.0 and 1.0")
    return 1
  
  # Validate font color
  if any(c < 0.0 or c > 1.0 for c in args.font_color):
    print("Error: Font color values must be between 0.0 and 1.0")
    return 1
  
  # Validate font size
  if args.font_size < 1:
    print("Error: Font size must be positive")
    return 1
  
  # List colormaps and exit if requested
  if args.list_colormaps:
    try:
      visualizer = ParaViewVTKVisualizer()
      colormaps = visualizer.get_available_colormaps()
      print("Available colormaps:")
      for i, cmap in enumerate(colormaps, 1):
        print(f" {i:2d}. {cmap}")
      return 0
    except Exception as e:
      print(f"Error reading VTK file: {e}")
      return 1
  # List scalars and exit if requested
  if args.list_scalars:
    if not args.vtk_file:
      print("Error: VTK file required for --list-scalars")
      return 1
    try:
      list_scalars(args.vtk_file)
      return 0
    except Exception as e:
      print(f"Error reading VTK file: {e}")
      return 1
  # Validate opacity
  if not (0.0 <= args.opacity <= 1.0):
    print("Error: Opacity must be between 0.0 and 1.0")
    return 1
  
  # Validate tube radius
  if args.tube_radius <= 0:
    print("Error: Tube radius must be positive")
    return 1
  # Validate light parameters
  if args.light_intensity is not None and not (0.0 <= args.light_intensity <= 2.0):
    print("Error: Light intensity must be between 0.0 and 2.0")
    return 1
  if args.light_azimuth is not None and not (0 <= args.light_azimuth <= 360):
    print("Error: Light azimuth must be between 0 and 360 degrees")
    return 1
  if args.light_elevation is not None and not (-90 <= args.light_elevation <= 90):
    print("Error: Light elevation must be between -90 and 90 degrees")
    return 1
  
  # Parse camera view string format if needed
  parsed_azimuth = args.camera_azimuth
  parsed_elevation = args.camera_elevation
  parsed_yaw = args.camera_yaw
  parsed_pitch = args.camera_pitch
  parsed_roll = args.camera_roll
  # Parse comma-separated format if not provided
  if args.camera_view and ',' in str(args.camera_view):
    try:
      parts = args.camera_view.split(',')
      if len(parts) < 2:
        raise ValueError
      parsed_azimuth = float(parts[0].strip())
      parsed_elevation = float(parts[1].strip())
      parsed_yaw = float(parts[2].strip()) if len(parts) > 2 else None
      parsed_pitch = float(parts[3].strip()) if len(parts) > 3 else None
      parsed_roll = float(parts[4].strip()) if len(parts) > 4 else None
    except (ValueError, IndexError):
      print("Error: Camera view format should be 'azimuth,elevation[,yaw,pitch,roll]' with valid numbers")
      return 1
  elif args.camera_view and args.camera_view not in ["isometric", "front", "top", "auto"]:
    print(f"Error: Unknown camera view preset '{args.camera_view}'. Use 'isometric', 'front', 'top', 'auto', or 'azimuth,elevation,yaw,pitch,roll'")
    return 1
  # Validate parsed camera angles
  if parsed_azimuth is not None and not (0 <= parsed_azimuth <= 360):
    print("Error: Camera azimuth must be between 0 and 360 degrees")
    return 1
  if parsed_elevation is not None and not (-90 <= parsed_elevation <= 90):
    print("Error: Camera elevation must be between -90 and 90 degrees")
    return 1
  if parsed_yaw is not None and not (0 <= parsed_yaw <= 360):
    print("Error: Camera yaw must be between 0 and 360 degrees")
    return 1
  if parsed_pitch is not None and not (-90 <= parsed_pitch <= 90):
    print("Error: Camera pitch must be between -90 and 90 degrees")
    return 1
  if parsed_roll is not None and not (-180 <= parsed_roll <= 180):
    print("Error: Camera roll must be between -180 and 180 degrees")
    return 1
  
  try:
    # Determine what to save
    save_png = not args.no_png
    save_pvsm = not args.png_only
    
    # Determine whether to add source/receiver
    add_source_receiver = not args.no_source_receiver and bool(args.phases)
    
    # Handle rays-only mode
    if args.rays_only:
      if not args.phases:
        print("Error: --phases argument is required when using --rays-only")
        return 1
      
      # Generate output filename if not specified
      if args.output == "Mesh_Visualization.pvsm":
        args.output = f"rays_{'_'.join(args.phases)}.pvsm"
      
      create_ray_only_visualization(
        ray_phases=args.phases,
        output_pvsm_path=args.output,
        tube_radius=args.tube_radius,
        ray_color=args.ray_color,
        add_source_receiver=add_source_receiver,
        sr_file_path=args.sr_file,
        sphere_color=args.sphere_color,
        save_png=save_png,
        camera_azimuth=parsed_azimuth,
        camera_elevation=parsed_elevation,
        camera_yaw=parsed_yaw,
        camera_pitch=parsed_pitch,
        camera_roll=parsed_roll,
        light_intensity=args.light_intensity,
        light_azimuth=args.light_azimuth,
        light_elevation=args.light_elevation,
        phases_dir=args.phases_dir,
        view_size=args.png_size,
        background_color=[1.0, 1.0, 1.0]
      )
      return 0
    
    # Handle mesh with optional rays
    if args.phases:
      # Create visualization with both mesh and rays
      slice_positions = None
      if args.slice_positions:
        slice_positions = {
          'x': args.slice_positions[0],
          'y': args.slice_positions[1],
          'z': args.slice_positions[2]
        }
      
      sr_file = args.sr_file
      if sr_file == "phases/sr.vtk":
          sr_file = f"{args.phases_dir}/sr.vtk"
      
      create_mesh_and_rays_visualization(
        vtk_file_path=args.vtk_file,
        ray_phases=args.phases,
        output_pvsm_path=args.output,
        mesh_name=args.name,
        color_by_scalar=args.color_by_scalar,
        tube_radius=args.tube_radius,
        ray_color=args.ray_color,
        add_source_receiver=add_source_receiver,
        sr_file_path=sr_file,
        sphere_color=args.sphere_color,
        save_png=save_png,
        camera_azimuth=parsed_azimuth,
        camera_elevation=parsed_elevation,
        camera_yaw=parsed_yaw,
        camera_pitch=parsed_pitch,
        camera_roll=parsed_roll,
        light_intensity=args.light_intensity,
        light_azimuth=args.light_azimuth,
        light_elevation=args.light_elevation,
        font_color=args.font_color,
        font_size=args.font_size,
        label=args.label,
        representation_type=args.representation,
        show_edges=args.edges,
        opacity=args.opacity,
        slice_positions=slice_positions,
        colormap=args.colormap,
        reverse_colormap=args.reverse_colormap,
        phases_dir=args.phases_dir,
        view_size=args.png_size,
        background_color=[1.0, 1.0, 1.0]
      )
      return 0
    
    # Handle mesh-only visualizations
    if args.png_only:
      # Only save PNG
      visualizer = ParaViewVTKVisualizer()
      
      # Prepare slice positions if provided
      slice_positions = None
      if args.slice_positions:
        slice_positions = {
          'x': args.slice_positions[0],
          'y': args.slice_positions[1],
          'z': args.slice_positions[2]
        }
      
      # Create appropriate visualization
      if args.show_full_mesh:
        # Multi-representation view with full mesh + slices
        result = visualizer.create_multi_representation_view(args.vtk_file, args.name)
        # Set the reader for camera operations
      elif args.advanced or args.color_by_scalar or args.edges or args.opacity != 1.0 or slice_positions or args.reverse_colormap:
        visualizer.create_advanced_visualization(
          vtk_file_path=args.vtk_file,
          mesh_name=args.name,
          representation_type=args.representation,
          color_by_scalar=args.color_by_scalar,
          show_edges=args.edges,
          opacity=args.opacity,
          slice_positions=slice_positions,
          colormap=args.colormap,
          reverse_colormap=args.reverse_colormap,
          font_color=args.font_color,
          font_size=args.font_size,
          label=args.label
        )
      else:
        visualizer.create_basic_visualization(
          vtk_file_path=args.vtk_file,
          mesh_name=args.name,
          color_by_scalar=args.color_by_scalar,
          font_color=args.font_color,
          font_size=args.font_size,
          label=args.label
        )
      
      # Apply camera view if specified
      if any(param is not None for param in [parsed_azimuth, parsed_elevation, parsed_yaw, parsed_pitch, parsed_roll]):
        apply_camera_view(visualizer, None, parsed_azimuth, parsed_elevation, parsed_yaw, parsed_pitch, parsed_roll)
      elif args.camera_view:
        apply_camera_view(visualizer, args.camera_view, None, None, None, None, None)
      else:
        visualizer.render_view.ResetCamera()

      # Essential render after camera changes
      Render()
      
      # Generate PNG filename from VTK filename if output not specified
      if args.output == "Mesh_Visualization.pvsm":
        png_path = os.path.splitext(args.vtk_file)[0] + '_slices.png'
      else:
        png_path = os.path.splitext(args.output)[0] + '.png'
      
      visualizer.save_screenshot(png_path, args.png_size[0], args.png_size[1])
      print(f"PNG screenshot saved: {png_path}")
      
    else:
      # Save PVSM and optionally PNG
      slice_positions = None
      if args.slice_positions:
        slice_positions = {
          'x': args.slice_positions[0],
          'y': args.slice_positions[1],
          'z': args.slice_positions[2]
        }
      
      if args.show_full_mesh:
        # Create multi-representation view
        visualizer = ParaViewVTKVisualizer()
        visualizer.create_multi_representation_view(args.vtk_file, args.name)

        # Apply camera view if specified
        if any(param is not None for param in [parsed_azimuth, parsed_elevation, parsed_yaw, parsed_pitch, parsed_roll]):
          apply_camera_view(visualizer, None, parsed_azimuth, parsed_elevation, parsed_yaw, parsed_pitch, parsed_roll)
        elif args.camera_view:
          apply_camera_view(visualizer, args.camera_view, None, None, None, None, None)
        else:
          visualizer.render_view.ResetCamera()

        # Essential render after camera changes
        Render()
        
        if save_png:
          png_path = os.path.splitext(args.output)[0] + '.png'
          visualizer.save_state_and_screenshot(args.output, png_path, 
                                               args.png_size[0], args.png_size[1])
        else:
          visualizer.save_state_file(args.output)
          
      elif args.advanced or args.color_by_scalar or args.edges or args.opacity != 1.0 or slice_positions:
        # Use create_advanced_pvsm_file
        create_advanced_pvsm_file(
          vtk_file_path=args.vtk_file,
          output_pvsm_path=args.output,
          mesh_name=args.name,
          representation=args.representation,
          color_by_scalar=args.color_by_scalar,
          camera_azimuth=parsed_azimuth,
          camera_elevation=parsed_elevation,
          camera_yaw=parsed_yaw,
          camera_pitch=parsed_pitch,
          camera_roll=parsed_roll,
          light_intensity=args.light_intensity,
          light_azimuth=args.light_azimuth,
          light_elevation=args.light_elevation,
          font_color=args.font_color,
          font_size=args.font_size,
          label=args.label,
          show_edges=args.edges,
          opacity=args.opacity,
          save_png=save_png,
          slice_positions=slice_positions,
          colormap=args.colormap,
          reverse_colormap=args.reverse_colormap
        )
      else:
        create_pvsm_file(
          vtk_file_path=args.vtk_file,
          output_pvsm_path=args.output,
          mesh_name=args.name,
          save_png=save_png,
          color_by_scalar=args.color_by_scalar,
          camera_azimuth=parsed_azimuth,
          camera_elevation=parsed_elevation,
          camera_yaw=parsed_yaw,
          camera_pitch=parsed_pitch,
          camera_roll=parsed_roll,
          light_intensity=args.light_intensity,
          light_azimuth=args.light_azimuth,
          light_elevation=args.light_elevation,
          font_color=args.font_color,
          font_size=args.font_size,
          label=args.label
        )
      
      return 0
      
  except Exception as e:
    print(f"Error creating visualization: {e}")
    return 1

if __name__ == "__main__":
  # Check if running in demo mode (no arguments)
  if len(sys.argv) == 1:
    print("ParaView Python API PVSM Generator")
    print("===================================")
    print()
    print("This script uses the ParaView Python API to create .pvsm state files.")
    print("It supports mesh visualization, seismic ray path plotting, and source/receiver points.")
    print()
    print("Requirements:")
    print("- ParaView installation with Python support")
    print("- ParaView Python modules in your Python path")
    print()
    print("Example usage:")
    print()
    print("Mesh visualization:")
    print("  python pvsm_generator.py mesh.vtk")
    print("  python pvsm_generator.py mesh.vtk --list-scalars")
    print("  python pvsm_generator.py mesh.vtv --color-by pressure --representation Surface")
    print("  python pvsm_generator.py mesh.vtk --advanced --edges --opacity 0.8")
    print()
    print("Ray paths only (with automatic source/receiver detection):")
    print("  python pvsm_generator.py --rays-only --phases P S PKIKP")
    print("  python pvsm_generator.py --rays-only --phases P S --tube-radius 0.005")
    print("  python pvsm_generator.py --rays-only --phases SKS --ray-color 0.8 0.2 0.2")
    print("  python pvsm_generator.py --rays-only --phases P S --no-source-receiver")
    print("  python pvsm_generator.py --rays-only --phases-dir phases_EQ001_STA01 --phases P S PKIKP")
    print()
    print("Mesh with ray paths and source/receiver points:")
    print("  python pvsm_generator.py mesh.vtk --phases P S PKIKP SKS")
    print("  python pvsm_generator.py mesh.vtv --phases P S --color-by Vp --tube-radius 0.001")
    print("  python pvsm_generator.py mesh.vtv --phases P S --sphere-color 0.8 0.0 0.0")
    print("  python pvsm_generator.py mesh.vtv --phases P S --sr-file sources_receivers.vtk")
    print()
    print("Camera control examples:")
    print("  python pvsm_generator.py mesh.vtk --camera-view isometric")
    print("  python pvsm_generator.py mesh.vtk --camera-view 210,10,0,0,20  # azimuth,elevation,yaw,pitch,roll")
    print("  python pvsm_generator.py mesh.vtk --camera-azimuth 50 --camera-elevation 30 --camera-yaw 90 --camera-pitch 45 --camera-roll 15")
    print("  python pvsm_generator.py mesh.vtk --camera-view 180,-15,0,0,0  # South view, 15° below horizon")
    print()
    print("Light source examples:")
    print("  python pvsm_generator.py mesh.vtk --light-intensity 0.5 --light-azimuth 140 --light-elevation 35")
    print("  python pvsm_generator.py mesh.vtv --phases P S --light-intensity 0.8 --light-azimuth 270 --light-elevation 30")
    print()
    print("Colorbar customization examples:")
    print("  python pvsm_generator.py mesh.vtk --color-by Vp --font-color 1.0 0.0 0.0 --font-size 24 --label 'Seismic Velocity'")
    print("  python pvsm_generator.py mesh.vtk --color-by Vp --font-color 0.0 0.0 1.0")
    print()
    print("Complete command line usage example:")
    print("  python pvsm_generator.py Mesh_VpVsRho.vtk --color-by Vp --phases-dir phases_EQ001_STA01 --phases PKIIKP PKPPm+PPcS ScSScSm-ScS SPS660-S --reverse-colormap" \
             " --camera-view 210,10,0,0,20 --label 'Vp [km/s]' --light-intensity 0.5 --light-azimuth 140 --light-elevation 35")
    print()
    print("Programmatic usage:")
    print("  from pvsm_generator import (create_pvsm_file, create_advanced_pvsm_file,")
    print("                              create_ray_only_visualization,")
    print("                              create_mesh_and_rays_visualization)")
    print("  create_pvsm_file('mesh.vtv', 'output.pvsm')")
    print("  create_ray_only_visualization(['P', 'S'], 'rays.pvsm')")
    print("  create_mesh_and_rays_visualization('mesh.vtv', ['P', 'S'], 'combined.pvsm')")
    print()
    print("Notes:")
    print("- Ray path VTK files should be named after their phases (e.g., P.vtk, S.vtk)")
    print("- Source/receiver file should be named 'sr.vtk' by default")
    print("- Source/receiver spheres are automatically added when using --phases (disable with --no-source-receiver)")
    print("- Sphere radius is automatically set to twice the tube radius")
    print("Run with --help for full command line options.")
  else:
    sys.exit(main())