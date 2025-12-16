# nomad-heatsink
This repository is used to determine the steady-state heat transfer for USC LPL's Nomad Engine. This code accounts for Nomad's chamber profile and copper heatsink parameters to determine temperature across the heatsink and chamber survivability for short-duration fires. 

## Pre-Requisites
This repository was intended with a functional C++ compiler over version 17. 

You must have the pyvista library installed in your local environment to plot the results. Ensure that you are in the correct environment, and then to install pyvista, run the following in terminal:

conda install pyvista

## Repository Structure
The main file used for this project is 'heatsink2D.cpp'. To compile this file, run the following:

g++ ./heatsink2D.cpp -o heatsink2D

To run the executable, use this command. This should save the results in an output file called 'field.vti'

./heatsink2D

To view the results, run this command:

python show_field.py

## Numerical Method & Discussion
This project uses the Gauss-Seidel solver to determine a steady-state temperature for the heatsink after a given number of iterations. The five-point stencil is used as in class, but now includes the convective heat transfer coefficient for copper. The mesh is a rectangle for ease of use, and includes a portion of hot gas in the chamber and shows the chamber and nozzle profile as designed. 

### Boundary Constraints
For the left and right boundaries, a zero neumann BC was applied as the goal of this simulation was to determine whether the selected heatsink thickness will survive. 

For the upper and lower boundaries, a Robin BC for convective heat transfer was applied. This is a mixed BC that is more applicable to this scenario compared to the Neumann or Dirchlet assumptions. The upper boundary was assumed to be ambient gas at 30K, and the lower boundary was assumed to be combustion gases of kerosene and liquid oxygen at 3670K. 

### Numerical Method
The Gauss-Seidel solver is a method used for steady-state calculations and should be used in conjunction with a convergence check to ensure the error decreases over iterations. While I took this route for this solution, one thing I realized is that with combustion and heat transfer over rockets, a steady-state solution is very difficult to reach. With the convergence check in place, the error would never decrease to below the threshold we discussed in class, even after 100,000 iterations. The new goal of this project should be to determine maximum firing time of the copper heatsink, which would have benefited from a FTCS and Crank-Nicholson solver instead that has defined time steps.