function [r,lg,e] = Autocorrelate(x,lag,method,pfa)
%Autocorrelate:  Estimates the signal autocorrelation. 
%
%   [r,lg,e] = Autocorrelate(x,lag,method,pf)
%   
%   x        Input Signal
%   lag      Length of lag (samples).  Default=length(x).
%   method   'biased', 'unbiased' (default), 'fast'.
%   pf       Plot flag: 0=none (default), 1=screen.
%
%   r        Estimated autocorrelation
%   lg       X-axis lags
%   e        Variance of estimated autocorrelation
%                       
%   Autocorrelate is a function that gives a normalized correlation
%   of a random stationary signal with itself.  The autocorrelation
%   is defined as
%
%      r(m) = sum((x(m) - mean(x)).*(x(m+nx) - mean(x)))./var(x)
%
%   where m is the correlation index and nx is the the length of
%   signal x.  The output is specified by the number of lags chosen,
%   so that the autocorrelation will be a value of 1 at lag 0 and
%   will extend out to the number of lags chosen.
%
%   The method chosen determines how the output will be normalized.
%   The 'unbiased' method is the most accurate, but slowest method,
%   and uses a normalization factor of 1/(nx-abs(m)).  The 'biased'
%   method is also slow, but will give a positive definite output and
%   uses a normalization factor of 1./nx.  The 'fast' method is the
%   least accurate method.  It calculates the autocorrelation by 
%   taking the inverse fft of the squared power spectral density.
%
%   If no output argument is specified, the default will plot to the 
%   screen.  If the 'e' input argument is specified, the variance of
%   the autocorrelation will also be plotted.
%
%   Example: Calculate the autocorrelation of an ABP data segment 
%   with a lag of 500 samples and a biased output plotted to the 
%   screen.
%
%      load ABPICP.mat
%      x = abp(1:1000);
%      [r,lg] = autocorrelate(x,500,'Biased',1);  
%
%   Priestley, M., "Spectral Analysis and Time Series," Academic Press
%   Limited, pp. 330-333, 1981.
%  
%   Version 1.00 LJ
%
%   See also AutoCorrelogram, CrossCorrelate, CrossCovariance, 
%   AutoCovariance, and CrossCorrelogram.
