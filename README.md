# Topology Optimization for Geometrical Nonlinearity Structure

The project modifies the original 230-line MATLAB code for **3D geometrically nonlinear topology optimization using the SIMP method** (from the paper: "MATLAB implementations for 3D geometrically nonlinear topology optimization: 230-line code for SIMP method and 280-line code for MMB method" by Yanfang Zhao, Guoqiang Guo, Wen, Wenjing Zuo, 2023) to also support **2D models**.

## How to Use

1. Open MATLAB and navigate to the repository folder.
2. Set the dimension mode in the main script
3. Adjust parameter:
    - Volume fraction ('volfrac')
    - Penalization ('penal')
    - Filter radius ('rmin')
    - Load and boundary conditions
4. Run the main file '2D_GNTO_SIMP.m'

Example for 2D cantilever beam:
```matlab
nelx = 120; nely = 40; nelz = 1;        % 2D mode
volfrac = 0.3; penal = 3; rmin = 1.5;