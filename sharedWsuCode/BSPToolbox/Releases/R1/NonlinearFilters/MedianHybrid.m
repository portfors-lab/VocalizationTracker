function [y] = MedianHybrid(x, wla, pfa)
%MedianHybrid:  Median Hybrid Filter 
%
%   [y] = MedianHybrid(x,w,pf)
%
%   x       Input signal
%   w       Width of the sliding window in samples (odd integer).
%           Default=11.
%   pf      Plot format: 0=none (default), 1=screen.
%
%   y       Filtered signal
%
%   Filters the signal using a Median Hybrid Filter. A window of width
%   w is placed at the beginning of the input vector. The window is
%   divided into three parts: the value in the central position, the 
%   values before the central value, and the values after the central
%   value. The mean of the values before the central value is computed
%   and stored in mean1, and the mean of the values after the central
%   value is computed and stored in mean2. Then, the median of mean1,
%   the central value, and mean2 are computed and stored into the output 
%   vector y. The same procedure is repeated as the window slides through 
%   the data, advancing one sample at a time.
%
%   Example: Filter the nonlinear filters test signal using a Median
%   Hybrid Filter with w = 31, and plot the results.
%
%      load NFSignal.mat;
%      [y] = MedianHybrid(x, 31, 1);
%
%   Astola, J. and Kuosmanen, P., "Fundamentals of Nonlinear Digital 
%   Filtering," CRC Press, pp.97-104, 1997.
%
%   Version 1.00 CC
%
%   See also MedianFilter and ComparisonSelection. 
