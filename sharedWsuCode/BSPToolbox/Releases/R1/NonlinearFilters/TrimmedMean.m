function [y] = TrimmedMean(x, wla, alphaa, pfa)
%TrimmedMean: Trimmed Mean Filter
%
%   [y] = TrimmedMean(x,wl,alpha,pf)
%
%   x       Input signal
%   wl      Width of the sliding window in samples (odd integer).
%           Default=31.
%   alpha   Proportion of values to be trimmed at both ends of the 
%           window (0 <= alpha <=0.48). Default=0.20.
%   pf      Plot format: 0=none (default), 1=screen.
%
%   y       Filtered signal
%
%   Filters the signal using an alpha-Trimmed Mean Filter. A window
%   of width wl is placed at the beginning of the input vector. The
%   input values within the window are sorted in ascending order, and
%   a number of alpha*wl values are trimmed at each end. The remaining
%   values are averaged, and their mean is stored in the output vector 
%   y. The same procedure is repeated as the window slides through the
%   data, advancing one sample at a time.
%
%   Example: Filter the nonlinear filters test signal using a Trimmed 
%   Mean Filter with alpha = 0.15 and wl = 31, and plot the results.
%
%     load NFSignal.mat;
%     [y] = TrimmedMean(x, 31, 0.15, 1);
%
%   Astola, J. and Kuosmanen, P., "Fundamentals of Nonlinear Digital 
%   Filtering," CRC Press, pp.52-55, 1997.
%
%   Version 1.00 CC
%
%   See also WinsorizedTrimmedMean, ModifiedTrimmedMean, and 
%   DWModifiedTrimmedMean.
