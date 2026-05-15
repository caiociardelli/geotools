# GeoTools Package

The _GeoTools_ Package is a geosciences software package for computing seismic ray paths, travel times, and synthetic gravity disturbances in 3D Earth models, incorporating surface topography and Moho topography. Written in C, Python, and Bash, it leverages parallel processing via OpenMPI for high-performance computations. The package is under active development, with planned features including 3D ray tracing, Earth's ellipticity, and inversion capabilities for 3D velocity and density distributions as well as Moho topography.

* **Author**: Caio Ciardelli, Northwestern University
* **Supervisor**: Prof. Suzan van der Lee
* **License**: GNU General Public License (GPL) version 3 or later

---

## Features

* **Seismic Ray Tracing**: Computes ray paths for source-receiver pairs using the `SeisTracer` library.
* **Travel Time Calculation**: Calculates seismic travel times for specified phases using a spherical mesh, accounting for Moho topography by locally deforming the ray path as an approximation.
* **Synthetic Gravity**: Computes gravitational field components in parallel, accounting for surface topography and Moho topography.
* **Mesh Generation**: Creates a spherical mesh based on icosahedron refinement, with optional ETOPO surface topography and CRUST1.0 Moho topography integration.
* **Visualization**: Generates GMT plots (2D) and ParaView state files (3D) for mesh, rays, and gravity fields.
* **Modular Workflow**: Automated test workflow for model resampling, ray tracing, meshing, and visualization.

---

## Installation

The _GeoTools_ Package requires no installation beyond setting up dependencies and compiling the C binaries. Follow these steps:
**1. Obtain the Package**
Extract the package directory to your desired location (e.g., `~/geotools`).
**2. Install Dependencies**

Ensure the following are installed:
* **Bash**: For running scripts (e.g., `run_test_workflow.bash`).
* **GCC**: For compiling C binaries (`mesher`, `gravity`, `ttimes`).
* **Python 3.7+**: For scripts (`compute_rays.py`, `resample_1d_model.py`, `pvsm_generator.py`).
* **Python Libraries**: Install via pip:
    ```bash
    pip install numpy matplotlib paraview
    ```
    *Note: `SeisTracer` may require additional setup; see its documentation.*
* **GMT (Generic Mapping Tools)**: For 2D plotting (`plot_mesh.bash`, `plot_gravity.bash`).
* **ParaView**: With Python support for 3D visualization (`pvsm_generator.py`).
* **OpenMPI**: For parallel gravity computation (`gravity.c`).
* **bc**: For arithmetic in `plot_gravity.bash`.
**3. Compile the Package**

Navigate to the package directory and compile the C binaries using the provided `Makefile`:
```bash
cd geotools
make all
```
This creates a `bin/` directory with `mesher`, `gravity`, and `ttimes`. To clean previous builds:
```bash
make clean
```
**4. Set Permissions**
Ensure utility scripts are executable:
```bash
chmod +x utils/*.bash
```

-----

## Directory Structure
```
geotools/
├── docs/
│ ├── Manual.pdf # Manual in PDF
│ └── Manual.tex # Manual in TEX
├── extra/
│ ├── CET-D01A.cpt # GMT color palette
│ ├── crust1_moho.dat # CRUST1.0 Moho data
│ └── earth_relief_15m.topo # ETOPO topography data
├── input/
│ ├── receivers.dat # Receiver locations (ID, lat, lon, elev)
│ └── sources.dat # Source locations (ID, lat, lon, depth)
├── LICENSE # GPL v3 license
├── Makefile # Compilation instructions
├── mesh/ # Output directory for mesh files
│ ├── facets.bin
│ ├── Mesh_VpVsRho.vtk
│ ├── r.bin
│ ├── rho.bin
│ ├── vertices.bin
│ ├── vp.bin
│ └── vs.bin
├── models/ # 1D Earth models
│ ├── ak135f.model
│ ├── iasp91.model
│ ├── prem.model
│ ├── rem1d.model
│ └── stw105.model
├── README.md
├── run_test_workflow.bash # Main workflow script
├── setup/
│ └── constants.h # Configuration constants for C Packages
├── src/
│ ├── gravity.c # Gravity computation
│ ├── headers/ # Header files (coordinates.h, exmath.h, etc.)
│ ├── mesher.c # Mesh generation
│ ├── shared/ # Shared utilities (coordinates.c, io.c, etc.)
│ └── ttimes.c # Travel time computation
└── utils/
    ├── compute_rays.py        # Ray path computation
    ├── plot_gravity.bash # Gravity plotting
    ├── plot_mesh.bash # Mesh plotting
    ├── pvsm_generator.py      # ParaView visualization
    └── resample_1d_model.py   # Model resampling
```

-----

## Usage

The main entry point is `run_test_workflow.bash`, which automates the full forward modeling process. Run it from the package root directory:
```bash
./run_test_workflow.bash
```

This script:
1. Cleans previous outputs (e.g., `.txt`, `.dat`, `.png`, `.vtk`, `.pvsm`, `phases_EQ001_STA01/`).
2. Resamples the IASP91 model to `input/input.model`.
3. Computes ray paths for sample source/receiver pairs, outputting to `phases_EQ001_STA01/`.
4. Generates a mesh with refinement level 7.
5. Creates GMT plots (e.g., `mesh320.png`, `nodes_spiral.png`).
6. Generates a ParaView state file for 3D visualization.
7. Computes travel times, outputting to `ttimes.txt`.
8. Computes synthetic gravity in parallel (6 cores, grid spacing 5°, height 255 km), outputting `G_x.dat`, `G_y.dat`, `G_z.dat`, `G_normal.dat`.
9. Plots gravity as `G_normal.png`.

For help:
```bash
./run_test_workflow.bash --help
```

### Individual Components
  * **`resample_1d_model.py`**: Resamples 1D Earth models (e.g., `models/iasp91.model`).
    ```bash
    python utils/resample_1d_model.py -i models/iasp91.model -o input/input.model
    ```
  * **`compute_rays.py`**: Computes seismic ray paths using `SeisTracer`.
    ```bash
    python utils/compute_rays.py -r input/receivers.dat -s input/sources.dat
    ```
  * **`mesher`**: Generates a spherical mesh with specified refinement, incorporating surface and Moho topography if enabled in `constants.h`.
    ```bash
    ./bin/mesher 7 input/input.model
    ```
  * **`ttimes`**: Computes travel times for phases in a directory, accounting for Moho topography if incorporated.
    ```bash
    ./bin/ttimes phases_EQ001_STA01/ > ttimes.txt
    ```
  * **`gravity`**: Computes gravity field in parallel, affected by Moho topography if incorporated.
    ```bash
    mpiexec -n 6 ./bin/gravity 5.0 255
    ```
  * **`plot_mesh.bash`**: Plots mesh data as PNGs.
    ```bash
    ./utils/plot_mesh.bash
    ```
  * **`plot_gravity.bash`**: Plots gravity data.
    ```bash
    ./utils/plot_gravity.bash G_normal
    ```
  * **`pvsm_generator.py`**: Generates ParaView state files.
    ```bash
    pvpython utils/pvsm_generator.py mesh/Mesh_VpVsRho.vtk --color-by Vp --phases-dir phases_EQ001_STA01 --phases PKIIKP PKPPm+PPcS ScSScSm-ScS SPS660-S
    ```

-----

## Velocity Models

The `models/` directory includes standard 1D Earth models:
  * `iasp91.model` (default, provided by Suzan van der Lee)
  * `ak135f.model`
  * `prem.model`
  * `rem1d.model`
  * `stw105.model`
Users can add custom models in the same format (depth, Vp, Vs, rho, with 4 header lines).

-----

## Development

The package is under active development with planned features:
  * 3D ray tracing.
  * Earth's ellipticity corrections.
  * Inversion for 3D velocity/density and Moho topography.
Development is hosted on GitHub in the [caio.ciardelli/geotools repository](<https://github.com/caiociardelli/geotools>).
Contributions are welcome; please contact the author for collaboration details.

-----

## Contact
For questions, suggestions, or bug reports, contact:

**Caio Ciardelli**
*Email*: `<caio.ciardelli@gmail.com>`
Northwestern University, Department of Earth and Planetary Sciences

-----

## Acknowledgments
  * Prof. Suzan van der Lee for supervision and providing the IASP91 model.
  * Northwestern University for support.
  * Open-source communities for GMT, ParaView, OpenMPI, and SeisTracer.
