# ElectronCloud

Visualise Molecules and Molecular Orbitals (MOs) with Python.

[<img width="548" height="132" alt="image" src="https://github.com/user-attachments/assets/143345a6-e07b-45e5-836a-d8a55c1332c4" />](https://aleksandraleszczyk.github.io/ElectronCloud/examples/ethylene_mo6_rhf.html)

Architecture
------------
**BasisGTO**

Parses a shell-based GTO basis into a flat list of Cartesian AO
functions and pre-computes (norm × contraction) coefficients so that
evaluation on a grid is a single vectorised loop.

**plot_molecular_orbital(Basis, mo_coeffs, mo_index, ...)**
1. Builds a regular 3-D grid around the molecule.
2. Evaluates every AO on the grid, then contracts them with the MO coefficients to obtain ψ(r).
3. Renders:
- Blue isosurface  at  ψ = +isovalue
- Gold isosurface  at  ψ = −isovalue
- Grey cylinder    for every inferred bond
- CPK sphere       for every atom (Scatter3d)

Installation
------------

    >> git clone git@github.com:AleksandraLeszczyk/ElectronCloud.git
    >> cd ElectronCloud
    >> pip install .
