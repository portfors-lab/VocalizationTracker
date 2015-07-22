function [y] = RankOrder(x, wla, ra, pfa)
%RankOrder: Rank Order Filter 
%
%   [y] = RankOrder(x,wl,r,pf)
%
%   x       Input signal
%   wl      Width of the sliding window in samples (odd integer)
%   r       Rank value (number between 0 and 100)
%   pf      Plot format: 0=none (default), 1=screen.
%
%   y       Filtered signal
%
%   Filters the signal using a Rank Order Filter. A window of width
%   wl is placed at the beginning of the input vector. The value 
%   associated with the rth percentile is computed and stored in the 
%   output vector y. The same procedure is repeated as the window 
%   slides through the data, advancing one sample at a time.
%
%   Example: Filter the nonlinear filters test signal using a Rank
%   Order Filter with r = 75 and wl = 31, and plot the results.
%
%      load NFSignal.m
%      [y] = RankOrder(x, 31, 75, 1);
% 
%   Astola, J. and Kuosmanen, P., "Fundamentals of Nonlinear Digital 
%   Filtering," CRC Press, pp.78-81, 1997.
%
%   Version 1.00 CC
%
%   See also WeightedOrderStatistics and MedianFilter. 
