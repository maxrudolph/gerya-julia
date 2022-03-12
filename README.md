This repo contains jupyter notebooks implementing solutions to selected problems from Taras Gerya's *Introduction to Numerical Geodynamic Modeling* in the Julia programming language.

These exercises were assigned to students as part of a graduate level GEL298 course at UC Davis, Winter 2022.

### Chapter-01-02-Divergence-EOS.ipynb
1. Computation of material properties (density, thermal expansivity, compressibility) using finite differences.
2. Plot a vector velocity field and its divergence in 2D.

### Chapter03-PoissonEquation.ipynb     
1. Solve the Poisson equation in 1D and compare the numerical solution with an analytic solution.
2. Solve the poisson equation in 2D cartesian coordinates using a 5-point stencil.
3. Use ML and GMRES to solve the poisson equation with a geometric multigrid preconditioner.

### Chapter05-StokesEquations.ipynb 
- Solve the Stokes equations with constant viscosity using a vorticity-stream function approach in 2D Cartesian coordinates.

### Chapter07-NumericalStokes.ipynb
- Solve the Stokes equations with variable viscosity in 2D Cartesian coordinates.

### Chapter08-Stokes-Continuity-Markers.ipynb
- Solve the Stokes equations with variable viscosity and advect density and viscosity using markers.

### Chapter10-HeatEquation.ipynb
- Solve the heat conduction equation.

### Chapter11-I2VIS.ipynb
- Solve the coupled Stokes and energy equations in 2D. 
- Advect temperature using markers.
- Use a cell-centered discretization for the temperature equation.
- This functionality is similar to the code I2VIS (Gerya and Yuen, 2003).

### Chapter13-I2ELVIS.ipynb
- Include the ability to model the deformation of a viscoelastic material.
- This functionality is similar to the code I2ELVIS (Gerya and Yuen, 2007).
- This notebook hasn't been finished/

### Chapter18-PoissonMultigrid.ipynb
- Solution of a Poisson problem using multigrid. The implementation here loosely follows what is described in Gerya's book.

## Additional Notebooks:
### M04-YinYang.ipynb
- Solution of a Poisson problem on a spherical surface using an overset grid, called the yin-yang grid, described in Kageyama and Sato (2004). The point of this problem was to think about how to set up the constraint equations and indexing for a code like StagYY. It is guaranteed to be far from optimal and the disretization doesn't appear to be well-suited to the use of iterative solvers.
