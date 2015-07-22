function [y] = ModifiedTrimmedMean(x, wla, qa, pfa)
%ModifiedTrimmedMean: Modified Trimmed Mean Filter
%
%   [y] = ModifiedTrimmedMean(x,wl,q,pf)
%
%   x       Input signal
%   wl      Width of the sliding window in samples (odd integer).
%           Default=31.
%   q       Maximum distance allowed between the median and each of
%           the data points within a window so that the data point
%           will not be trimmed when computing the output. Default=0.2.
%   pf      Plot format: 0=none (default), 1=screen.
%
%   y       Filtered signal
%
%   Filters the signal using a Modified Trimmed Mean Filter. A window
%   of width wl is placed at the beginning of the input vector. The
%   input values within the window are sorted in ascending order, and
%   all the values whose distance from the median is greater than q are 
%   trimmed. The remaining values are averaged, and their mean is stored
%   in the output vector y. The same procedure is repeated as the window 
%   slides through the data, advancing one sample at a time.
%
%   Example: Filter the nonlinear filters test signal using a Modified
%   Trimmed Mean Filter with q = 0.2 and wl = 21, and plot the results.
%
%      load NFSignal.mat;
%      [y] = ModifiedTrimmedMean(x, 21, 0.2, 1);
%
%   Astola, J. and Kuosmanen, P., "Fundamentals of Nonlinear Digital 
%   Filtering," CRC Press, p.56, 1997.
%
%   Version 1.00 CC
%
%   See also TrimmedMean, WinsorizedTrimmedMean, and DWModifiedTrimmedMean.
