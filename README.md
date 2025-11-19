# EDO-Solver
A differential equation numerical method comparison


This is an R App that uses Shiny t implement different numerical methods such as Euler, RK4 and Dormand Prince (RK45).

The app displays two different differential equations modelling real life phenomena: 1. The absorption of a pharmaceutical, 2. The absorption of Caffeine.

The model includes a comparison to the _exact_ analytical solution, and it also compares the library use of **deSolve**. 


It includes a graph of the steps, the visible table of the steps and the different table's solution, plots of each method and plots of the absolute error (comparing with the exact analytical solution).

Finally, it also displays the implementation _in code_ of the numerical methods (not a library) and briefly explains how are they used.
