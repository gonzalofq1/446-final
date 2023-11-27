# Numerical Solution for the Nonlinear Schrodinger Equation


This is the the code for the final project for ESAM 446 with Prof. Daniel Lecoanet. This project was prepared by Divjyot Singh and Gonzalo Ferrandez Quinto.

# Running the code

To regenerate the simulations presented in our report (Numerical Solution for the Nonlinear Schrodinger
Equation.pdf), simply run the scripts in 'SchrodingerPDE.ipynb' using a Python Kernel. This Jupyter notebook includes the code for all of our simulations and visualizations

# Changing Domain, Resolution, and Initial Conditions

Most simulation parameters can be conveniently changed in SchrodingerPDE.ipynb by creating an instance of the equations.SchrodingerBCNonLinear class with desired parameters and taking the desired number of time steps. To monitor the simulation, the instance variable $t$
To use different initial conditions, go to SchrodingerPDE.ipynb

# Changing Boundary Conditions

To use different boundary conditions, go to equations.py and create a copy of the SchrodingerBCNonLinear and change the function BC(X) in the Advection sub-class to fit your needs. An example of this can be seen in the class SchrodingerBCLinearSlit.
