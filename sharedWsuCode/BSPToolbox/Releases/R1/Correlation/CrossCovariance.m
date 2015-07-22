function [cc,lg] = CrossCovariance(x1,x2,lag,method,pfa)
%CrossCovariance: Estimate the cross-covariance of two signals.
%
%   [cc,lg] = CrossCovariance(x1,x2,lag,method,pf)
%
%   x1       Input signal 1
%   x2       Input signal 2
%   lag      Length of lag (samples).  Default=length(x).
%   method   'biased', 'unbiased' (default), 'fast'.
%   pf       Plot format:  0=none (default), 1=screen.
%
%   cc       Estimated cross-covariance
%   lg       X-axis lags
%
%   CrossCovariance is a function that gives an unnormalized 
%   covariance of two equal length random stationary signals, x1 
%   and x2.  The cross-covariance is defined as
%
%      cc(m) = sum((x1(m) - mean(x1)).*(x2(m+nx) - mean(x2)))
%
%   where m is the covariance index and nx is the length of
%   signal.  The range of the output is specified by the number
%   of lags chosen, and will be from -lag:+lag.
%
%   The method chosen determines how the output will be normalized.
%   The 'unbiased' method is the most accurate, but slowest method,
%   and uses a normalization factor of 1/(nx-abs(m)).  The 'biased'
%   method is also slow, but will give a positive definite output and
%   uses a normalization factor of 1./nx.  The 'fast' method is the
%   least accurate method.  It calculates the cross-covariance by 
%   taking the inverse fft of (X1 .* conj(X2)).  If no output 
%   argument is specified, the default will plot to the screen.  
%
%   Example:  Calculate the cross-covariance of an ABP signal and an
%   ICP signal with a lag of 1000 and a biased output printed to the 
%   screen.
%
%      load ABPICP.mat
%      x1 = abp(1:2000);
%      x2 = icp(1:2000);
%      [cc,lg] = CrossCovariance(x1,x2,1000,'biased',1);
%
%   Challis, R. E., and Kitney, R. I., "Biomedical Signal Processing
%   (in four parts), Part 1, Time-domain Methods", Med. & Biol. Eng.
%   & Comput., 1990, 28, pp. 509-524.
%
%   Version 1.00 LJ
%
%   See Also AutoCov, Autocorrelate, and CrossCorrelate.
