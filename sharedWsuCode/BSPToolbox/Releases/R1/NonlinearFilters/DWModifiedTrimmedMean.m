function [y] = DWModifiedTrimmedMean(x, w1a, w2a, qa, pfa)
%DWModifiedTrimmedMean: Double-Window Modified Trimmed Mean Filter
%
%   [y] = DWModifiedTrimmedMean(x,w1,w2,q,pf)
%
%   x       Input signal
%   w1      Width of the median window in samples (odd integer).
%           Default=11.
%   w2      Width of the mean window in samples (odd integer)
%           (w2 must be >= w1). Default=31.
%   q       Maximum distance allowed (positive integer). Default=0.2
%   pf      Plot format: 0=none (default), 1=screen.
%
%   y       Filtered signal
%
%   Filters the signal using a Double-Window Modified Trimmed Mean Filter. 
%   A window of width w2 is placed at the beginning of the input vector. 
%   The input values within the window are sorted in ascending order, and
%   a window of width w1 (w1 <= w2) is placed in the center of the sorted 
%   vector. The median of these w1 values is calculated. All the values in   
%   the initial window w2 whose distance from the median is greater than q 
%   are trimmed. The remaining values are averaged, and their mean is stored
%   in the output vector y. The same procedure is repeated as the window 
%   slides through the data, advancing one sample at a time.
%
%   Example: Filter the nonlinear filters test signal using a Double-Window
%   Modified Trimmed Mean Filter with q = 0.2, w1 = 11 and w2 = 21, and plot
%   the results.
%
%      load NFSignal.mat;
%      [y] = DWModifiedTrimmedMean(x, 11, 21, 0.2, 1);
%
%   Astola, J. and Kuosmanen, P., "Fundamentals of Nonlinear Digital 
%   Filtering," CRC Press, p.56, 1997.
%
%   Version 1.02 CC
%
%   See also TrimmedMean, WinsorizedTrimmedMean, and ModifiedTrimmedMean.
