function y = DetectMaxima(x,thrsa,pfa)
%DetectMaxima: Detect maxima (peaks) in the input signal x
%
%   y = DetectMaxima(x, thrs, pfa)
%
%   x      Input signal       
%   thrs   Threshold value (optional)
%   pf     Plot flag: 0=none (default), 1=screen
%  
%   y      Location (index) of all maxima in x
%
%   Finds the location (index) of all the maxima (peaks) in the 
%   input signal x and returns them in a sorted vector (y). 
%   If the thrs input is used, the funtion returns those peaks that
%   are greater than thrs. 
%
%   Example: Detect all the peaks in the ICP signal provided with the 
%   bsp toolbox data (ICP.m) and plot the results.
%
%      load ICP; 
%      DetectMaxima(icp);
%
%   Version 1.0 MA
%
%   See also DetectMinima, ManualDetector, and PressureDetect.
