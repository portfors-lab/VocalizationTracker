function [y] = SelectiveAverage(x, wa, pfa)
%SelectiveAverage: Selective Average Filter 
%
%   [y] = SelectiveAverage(x,w,pf)
%
%   x       Input signal
%   w       Width of the sliding window in samples (odd integer).
%           Default=11.
%   pf      Plot format: 0=none (default), 1=screen.
%
%   y       Filtered signal
%
%   Filters the signal using a Selective Average Filter. A window
%   of width w is placed at the beginning of the input vector. The
%   input values within the window are divided into three groups:
%   the central value, the values before the central value (lower 
%   values), and the values after the central value (upper values). 
%   The mean of the lower values and the mean of the upper values 
%   are then computed and compared to the central value. The one whose
%   value is closest to the central value is stored in the output 
%   vector y. If both of them are at the same distance from the central 
%   value, the mean of the whole window is stored in the output vector
%   y. The same procedure is repeated as the window slides through the
%   data, advancing one sample at a time.
%
%   Example: Filter the nonlinear filters test signal using a Selective 
%   Average Filter with w = 31, and plot the results.
%
%      load NFSignal.mat;
%      [y] = SelectiveAverage(x, 31, 1);
%
%   Astola, J. and Kuosmanen, P., "Fundamentals of Nonlinear Digital 
%   Filtering," CRC Press, pp.108-109, 1997.
%
%   Version 1.02 CC
%
%   See also SelectiveMedian and ComparisonSelection.
