function y = DetectMinima(x,thrsa,pfa)
%DetectMinima: Detect minima in the input signal x
%
%   y = DetectMinima(x,thrs,pf)
%
%   x      Input signal       
%   thrs   Threshold value (optional)
%   pf     Plot flag: 0=none (default), 1=screen
%
%   y      Location (index) of the mimima  
%
%   Finds the location (index) of all the minima in the 
%   input signal x and returns them in a sorted vector (MIN). 
%   If the thrs input is used, the funtion returns those peaks that
%   are smaller than thrs. 
%
%   Example: Detect all the minima in the ICP signal provided with the 
%   bsp toolbox data (ICP.m) and plot the results.
%
%      load ICP; 
%      DetectMinima(icp);
%
%   Version 1.0 MA
%
%   See also DetectMaxima, ManualDetector, and PressureDetector.
