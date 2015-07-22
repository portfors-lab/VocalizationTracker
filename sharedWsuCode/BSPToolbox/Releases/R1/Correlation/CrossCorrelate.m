function [cc,lg] = CrossCorrelate(x1,x2,lag,method,pfa)
%CrossCorrelate: Estimate the cross-correlation of two signals.
%
%   [cc,lg] = CrossCorrelate(x1,x2,lag,method,pf)
%
%   x1       Input Signal 1
%   x2       Input Signal 2
%   lag      Length of lag (samples).  Default=length(x).
%   method   'biased', 'unbiased' (default), 'fast'.
%   pf       Plot format:  0=none (default), 1=screen
%
%   cc       Estimated cross-correlation
%   lg      X-axis lags
%
%   CrossCorrelate is a function that gives a normalized correlation
%   of two equal length random stationary signals, x1 and x2.  The 
%   cross-correlation is defined as:
%
%      cc(m) = sum((x1(m) - mean(x1)).*(x2(m+nx) - mean(x2)))
%
%   where m is the correlation index and nx is the length of the
%   signal.  This figure is divided by the var(x1).*var(x2) so that
%   -1<=cc<=+1.  The range of the output is specified by the number
%   of lags chosen, and will be from -lag:+lag.
%
%   The method chosen determines how the output will be normalized.
%   The 'unbiased' method is the most accurate, but slowest method,
%   and uses a normalization factor of 1/(nx-abs(m)).  The 'biased'
%   method is also slow, but will give a positive definite output and
%   uses a normalization factor of 1./nx.  The 'fast' method is the
%   least accurate method.  It calculates the cross-correlation by 
%   taking the inverse fft of (X1 .* conj(X2)).
%
%   If no output argument is specified, the default will plot to the 
%   screen.  
%
%   Example:  Calculate the cross-correlation of an ABP data segment 
%   and an ICP data segment with a lag of 500, a biased output, and 
%   output plotted to the screen.
%
%      load ABPICP.mat
%      x1 = abp(1:1000);
%      x2 = icp(1:1000);
%      [cc lg] = CrossCorrelate(x1,x2,500,'Biased',1);  
%
%   Priestley, M., "Spectral Analysis and Time Series," Academic 
%   Press Limited, pp. 330-333, 1981.
%
%   Version 1.00 LJ
%
%   See Also Autocorrelate, CrossCovariance, AutoCovariance, and 
%   CrossCorrelogram.
