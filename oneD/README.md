# Information
The solver uses a semi-discrete method to solve the 1D-Euler-equations. 
To approximate the flux function the following numerical fluxes are implemented
* van-Leer flux vector splitting (option: vanLeer),
* Lax-Friedrichs flux (option: LF) and
* Harten-Lax-van-Leer flux (option: HLL).

For the reconstruction of the cell interface values several options are
available:
* no reconstruction (option: none),
* min-mod-limiter (option: minMod),
* van-Leer-limiter (option: vanLeer),
* local double logarithmic reconstruction (option: LDLR).

Finally the ODE in time is solved by using a Strong-Stability-Preserving
(SSP) Runge-Kutta method of third order to obtain an overall third order 
numerical solution.
 
Calculations are done on an equidistant grid. 
Initial conditions need to be given in the form of densitiy, velocity and
pressure.

Implemented boundary conditions (BCs):
* transmissive boundaries (option: transmissive)
* periodic BCs (option: periodic)
* transmissive left boundary, reflective right boundary (option: reflectiveRight)
* both boundaries reflective (option: reflectiveFull)
If no BC is given to the solver class transmissive boundaries are set as default.

Test cases available in *main.m*:
* Sod problem (testCase = 1)
* Smooth problem (testCase = 2)
* Lax problem with reflective right boundary (testCase = 3)
* Lax problem with both boundaries reflective (testCase = 4)
* Shu-Osher shock-acoustic problem (testCase = 5)



# Literature
* **[1]** Artebrant, R. and Schroll, H.-J. "Limiter-Free Third Order Logarithmic Reconstruction"
          SIAM J. Sci. Comput. 28.2 (2006): 359-381. DOI: 10.1137/040620187  
* **[2]** Gottlieb, S. and Shu, C.-W. and Tadmor, E. "Strong Stability-Preserving High-Order Time Discretization Methods"
          SIAM Rev. 43.1 (2001). DOI: 10.1137/S003614450036757X  
* **[3]** Van Leer, B. (1997) "Flux-vector Splitting for the Euler Equation" in 
          M. Y. Hussaini, B. van Leer, J. van Rosendale (ed.) *Upwind and High-Resolution Schemes*.
          Springer, Berlin Heidelberg, pp. 80-89. DOI: 10.1007/978-3-642-60543-7_5
    
