This repo contains jupyter notebooks implementing solutions to selected problems from Taras Gerya's *Introduction to Numerical Geodynamic Modeling* in the Julia programming language.

These exercises were assigned to students as part of a graduate level GEL298 course at UC Davis, Winter 2022.

### Chapter-01-02-Divergence-EOS.ipynb
1. Computation of material properties (density, thermal expansivity, compressibility) using finite differences.
2. Plot a vector velocity field and its divergence in 2D.

### Chapter03-PoissonEquation.ipynb     
1. Solve the Poisson equation in 1D and compare the numerical solution with an analytic solution.
2. Solve the poisson equation in 2D cartesian coordinates using a 5-point stencil.

### Chapter05-StokesEquations.ipynb 
- Solve the Stokes equations with constant viscosity using a vorticity-stream function approach in 2D Cartesian coordinates.

### Chapter07-NumericalStokes.ipynb
- Solve the Stokes equations with variable viscosity in 2D Cartesian coordinates.

### Chapter18-PoissonMultigrid.ipynb
- Solution of a Poisson problem using multigrid. The implementation here loosely follows what is described in Gerya's book.
