function [y] = WeightedMedian(x, wla, aa, pfa)
%WeightedMedian: Weighted Median Filter 
%
%   [y] = WeightedMedian(x,wl,a,pf)
%
%   x       Input signal
%   wl      Width of the sliding window in samples (odd integer).
%           Default=11.
%   a       Weight vector (must have the same length as wl). 
%           Default=vector of ones.
%   pf      Plot format: 0=none (default), 1=screen.
%
%   y       Filtered signal
%
%   Filters the signal using a Weighted Median Filter. A window
%   of width wl is placed at the beginning of the input vector. The
%   values within the window are weighted by computing a(i)<>x(i)
%   for values of i from 1 to N, where a(i)<>x(i) is a vector which
%   contains the value x(i) a(i) times (for instance, 3<>2 = [2 2 2])
%   The vectors obtained from the different i values are concatenated
%   and the median of the new vector is computed and stored in the 
%   output vector y. The same procedure is repeated as the window 
%   slides through the data, advancing one sample at a time.
%
%   Example: Filter the nonlinear filters test signal using a Weighted 
%   Median Filter with wl = 7 and a = [0 1 1 2 5 2 1 1 0], and plot the 
%   results.
%
%      load NFSignal.mat;
%      [y] = WeightedMedian(x, 31, [0 1 1 2 5 2 1 1 0], 1);
%
%   Astola, J. and Kuosmanen, P., "Fundamentals of Nonlinear Digital 
%   Filtering," CRC Press, pp.73-77, 1997.
%
%   Version 1.00 CC
%
%   See also MedianFilter, WeightedOrderStatistics, and RankOrder.
