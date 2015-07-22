function [y] = WeightedOrderStatistics(x, wla, aa, ra, pfa)
%WeightedOrderStatistics: Weighted Order Statistics Filter 
%
%   [y] = WeightedOrderStatistics(x,wl,a,r,pf)
%
%   x       Input signal
%   wl      Width of the sliding window in samples (odd integer).
%           Default=11.
%   a       Weight vector (must have the same length as wl). Default=
%           vector of ones.
%   r       Rank value (between 0 and 100). Default=50.
%   pf      Plot format: 0=none (default), 1=screen.
%
%   y       Filtered signal
%
%   Filters the signal using a Weighted Order Statistics Filter. A 
%   window of width wl is placed at the beginning of the input vector. 
%   The values within the window are weighted by computing a(i)<>x(i)
%   for values of i from 1 to N, where a(i)<>x(i) is a vector which
%   contains the value x(i) a(i) times (for instance, 3<>2 = [2 2 2])
%   The vectors obtained from the different i values are concatenated
%   and the value associated with the rth percentile of the new vector
%   is computed and stored in the output vector y. The same procedure 
%   is repeated as the window slides through the data, advancing one 
%   sample at a time. 
%
%   Example: Filter the nonlinear filters test signal using a Weighted 
%   Order Statistics Filter with wl = 5, a = [1 1 2 3 2 1 1], and
%   r = 25, and plot the results.
%
%      load NFSignal.mat;
%      [y] = WeightedOrderStatistics(x, 7, [1 1 2 3 2 1 1], 25, 1);
%
%   Astola, J. and Kuosmanen, P., "Fundamentals of Nonlinear Digital 
%   Filtering," CRC Press, pp.81-87, 1997.
%
%   Version 1.00 CC
%
%   See also WeightedMedian and RankOrder.
