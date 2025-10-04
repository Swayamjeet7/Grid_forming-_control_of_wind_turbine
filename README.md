This project gives an idea to connect a PMSG based wind turbine 
with a pre-existing grid through a LCL-filter.

The file `grid.slx` is a simulink model containing the grid 
connected DC source, where parameters of LCL-filter are 
determined by direct formulae and the PWM signals are generated 
using PLL and d-q control method. 
For calculating the parameters of filter, the `lcl_param.m` 
function is used.

File `L_filter.slx` is the model of system containing L-filter.
`L_filter_code.m` is the code where parameters of L-filter are 
evaluated using EHD method and then, are optimized.

File `LCL_filter_code.m` contains the code to determine the 
optimized parameters of the LCL-filter using EHD method.