# Numerical Solution for the Nonlinear Schrodinger Equation


This is the the code for the final project for ESAM 446 with Prof. Daniel Lecoanet. This project was prepared by Divjyot Singh and Gonzalo Ferrandez Quinto.

# Running the code

To regenerate the simulations presented in our report `Numerical Solution for the Nonlinear Schrodinger
Equation.pdf`, simply run the scripts in `SchrodingerPDE.ipynb` using a Python Kernel. This Jupyter notebook includes the code for all of our simulations and visualizations

# Changing Simulation Parameters

Most simulation parameters can be conveniently changed in `SchrodingerPDE.ipynb` by creating an instance of the `equations.SchrodingerBCNonLinear` class with desired parameters and taking the desired number of time steps.

We assume a square uniform domain, but the values and resolution of this domain can be adjusted by changing the variables `grid_x`, `grid_y`, and `resolution`.

To change the ratio of the timestep size to grad spacing $\alpha = dt/dx$, the variable `alpha` can be changed.

To monitor the simulation time, the instance variable `t` can be accessed and compared against desired values. Or more simply, the variable `T` in `SchrodingerPDE.ipynb` can be used to set the final time.

To change the initial conditions, the definition of `IC` can be changed.

Changing the boundary conditions is a little less straightforward. To use different boundary conditions, go to `equations.py` and create a copy of the `SchrodingerBCNonLinear` and change the function `BC` in the `Reaction` sub-class to fit your needs. An example of this can be seen in the class `SchrodingerBCLinearSlit`. However, it must be noted that any boundary conditions will be set in addition to particle-in-box condition which is set implictly in the diffusion timestepper.
