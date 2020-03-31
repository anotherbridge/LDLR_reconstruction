# Author
André Breuer

# Information
The solver uses a semi-discrete method to solve the 1D-Euler-equations. 
To approximate the flux function the following numerical fluxes are implemented
* van-Leer flux vector splitting (option: vanLeer),
* Lax-Friedrichs flux (option: LF) and
* Harten-Lax-van-Leer (option: HLL).

The reconstruction of the values on the cell interfaces is calculated
via the third order Local-Double-Logarithmic-Reconstruction (LDLR). 

Finally the ODE in time is solved by using a Strong-Stability-Preserving
(SSP) Runge-Kutta method of third order to obtain an overall third order 
numerical solution.
 
Calculations are done on an equidistant grid. 
Initial conditions need to be given in the form of densitiy, velocity and
pressure.
Two test cases are implemented. One is the Sod problem. The other one is
a smooth problem for which there exists an exact solution which could be
used for testing accuracy. Other test cases can be implemented.
Periodic boundary conditions can used.

# Literature
* **[1]** Artebrant, R. and Schroll, H.-J. "Limiter-Free Third Order Logarithmic Reconstruction"
          SIAM J. Sci. Comput. 28.2 (2006): 359-381. DOI: 10.1137/040620187  
* **[2]** Gottlieb, S. and Shu, C.-W. and Tadmor, E. "Strong Stability-Preserving High-Order Time Discretization Methods"
          SIAM Rev. 43.1 (2001). DOI: 10.1137/S003614450036757X  
* **[3]** Van Leer, B. (1997) "Flux-vector Splitting for the Euler Equation" in 
          M. Y. Hussaini, B. van Leer, J. van Rosendale (ed.) *Upwind and High-Resolution Schemes*.
          Springer, Berlin Heidelberg, pp. 80-89. DOI: 10.1007/978-3-642-60543-7_5
    
