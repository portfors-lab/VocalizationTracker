function [y] = MedianFilter(x, wla, pfa)
%MedianFilter:  Median Filter 
%
%   [y] = MedianFilter(x,w,pf)
%
%   x       Input signal
%   wl      Width of the sliding window in samples (odd integer).
%           Default=21.
%   pf      Plot format: 0=none (default), 1=screen.
%
%   y       Filtered signal
%
%   Filters the signal using a Median Filter. A window of width w
%   is placed at the beginning of the input vector. The median of
%   the values within the window is calculated and stored in the 
%   output vector y. The same procedure is repeated as the window 
%   slides through the data, advancing one sample at a time.
%
%   Example: Filter the nonlinear filters test signal using a Median 
%   Filter with  wl = 31, and plot the results.
%
%      load NFSignal.mat;
%      [y] = MedianFilter(x, 31, 1);
%
%   Astola, J. and Kuosmanen, P., "Fundamentals of Nonlinear Digital 
%   Filtering," CRC Press, pp.52-55, 1997.
%
%   Version 1.00 CC
%
%   See also WeightedMedian, RankOrder, and WeightedOrderStatistics. 
