# ElectronCloud

Visualise Gaussian-Type Orbital wave-functions interactively with Plotly.

Architecture
------------
BasisGTO
    Parses a shell-based GTO basis into a flat list of Cartesian AO
    functions and pre-computes (norm × contraction) coefficients so that
    evaluation on a grid is a single vectorised loop.

plot_molecular_orbital(Basis, mo_coeffs, n, ...)
    1. Builds a regular 3-D grid around the molecule.
    2. Evaluates every AO on the grid, then contracts them with the MO
       coefficients to obtain ψ(r).
    3. Renders:
         • Blue  isosurface  at  ψ = +isovalue
         • Gold  isosurface  at  ψ = −isovalue
         • Grey  cylinder    for every inferred bond
         • CPK sphere        for every atom (Scatter3d)