function [ac,lg] = AutoCovariance(x,lag,method,pfa)
%AutoCovariance: Estimates the signal autocovariance.
%
%   [ac,lg] = AutoCovariance(x,lag,method,pf)
%
%   x        Input Signal
%   lag      Length of lag (samples).  Default=length(x).
%   method   'biased', 'unbiased' (default), 'fast'.
%   pf       Plot flag:  0=none (default), 1=screen.
%
%   ac       Estimated autocovariance
%   lg       X-axis lag
%
%   AutoCovariance is a function that gives an unnormalized 
%   covariance of a random stationary signal with itself.  
%   The auto-covariance is defined as
%
%        ac(m) = sum((x(m) - mean(x)).*(x(m+nx) - mean(x)))
%
%   where m is the covariance index and nx is the the length of
%   signal x.  The output is specified by the number of lags chosen,
%   so that the auto-covariance will have a maximum value at lag 0 
%   and will extend out to the number of lags chosen.
%
%   The method chosen determines how the output will be normalized.
%   The 'unbiased' method is the most accurate, but slowest method,
%   and uses a normalization factor of 1/(nx-abs(m)).  The 'biased'
%   method is also slow, but will give a positive definite output and
%   uses a normalization factor of 1./nx.  The 'fast' method is the
%   least accurate method.  It calculates the autocovariance by 
%   taking the inverse fft of the squared power spectral density.
%   If no output argument is specified, the default will plot to the 
%   screen.  
%
%   Example: Calculate the autocovariance of an ABP data segment with
%   a lag of 1000 and a biased output and plot the results to the 
%   screen.
%
%      load ABPICP.mat
%      x = abp(1:2000);
%      [ac,lg] = Autocovariance(x,1000,'biased',1);
%
%   Shumway, Robert and Stoffer, David, "Time Series Analysis and Its
%   Applications," Springer, pp.15-37, 2000.
%
%   Version 1.00 LJ
%
%   See Also CrossCovariance, Autocorrelate, and CrossCorrelate.
