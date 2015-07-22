function [y] = ComparisonSelection(x, wla, da, pfa)
%ComparisonSelection: Comparison and Selection Filter 
%
%   [y] = ComparisonSelection(x,w,d,pf)
%
%   x       Input signal
%   wl      Width of the sliding window in samples (odd integer). 
%           Default = ceiling of length(x)/100.
%   d       Distance in rank from the median that the output value is
%           to be taken from. Default = w/10.
%   pf      Plot format: 0=none (default), 1=screen.
%
%   y       Filtered signal
%
%   Filters the signal using a Comparison and Selection Filter. A 
%   window of width w is placed at the beginning of the input vector. 
%   The input values within the window are sorted in ascending order, 
%   and the mean and median of the values within the window are 
%   calculated. If the mean is greater than the median, the value 
%   which is d places before the central value (the median value)is
%   stored in the output vector y. Conversely, if the median is 
%   greater then the mean, the value which is d places after the 
%   central value is stored in the output vector y. If both mean and 
%   median are equal, the central value is stored into the output 
%   vector y. The same procedure is repeated as the window slides 
%   through the data, advancing one sample at a time.
%
%   Example: Filter the nonlinear filters test signal using a  
%   Comparison Selection Filter with window length w = 31 and distance
%   value d = 3 samples, and plot the results.
%
%      load NFSignal.mat
%      [y] = ComparisonSelection(x, 31, 3, 1);
%
%   Astola, J. and Kuosmanen, P., "Fundamentals of Nonlinear Digital 
%   Filtering," CRC Press, pp.107, 1997.
%
%   Version 1.02 CC
%
%   See also SelectiveAverage and SelectiveMedian.
