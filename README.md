# scale_region

A Python tool for scaling regional MPAS meshes

## Overview

The `scale_region.py` tool is used to scale the size of elements in a regional
MPAS mesh, so that the resulting Voronoi cells and Delaunay triangles are some
factor of their original size. Typically, a coarser regional mesh is scaled to
higher resolution, with a corresponding decrease in the area covered by the
mesh.

Scaling is performed on a stereographic projection plane, and because mesh cells
move between areas of different map scale factor on the plane, the scaled MPAS
mesh will not be exactly the specified factor smaller (or larger) than the
original mesh in all areas.

## Requirements

The `scale_region.py` tool requires:
- Python 3.11+
- netCDF4
- NumPy

## Usage

Five command-line arguments are required by the `scale_region.py` tool, and
these are summarized in the help message provided when `scale_region.py` is run
with the `-h` / `--help` option. Conceptually, these arguments serve three
purposes:

1. To provide the name of the input regional mesh to be scaled, as well as to
specify the name of the scaled output mesh to be created.

2. To specify the factor by which the input regional mesh is to be scaled. Real
values larger than 1.0 correspond to a reduction in mesh size (i.e., an increase
in resolution); for example, a scaling factor of 4.0 reduces the size of all
cells in the mesh by a factor of 4.0.

3. To provide the (lat, lon) coordinates of the tangent point of the
stereographic projection. Generally, the tangent point should be chosen to lie
close to the geographic center of the regional mesh.
