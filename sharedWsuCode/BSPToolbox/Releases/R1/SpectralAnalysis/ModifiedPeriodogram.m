function [p,f] = ModifiedPeriodogram(x,fsa,wna,nfa,pfa)
%ModifiedPeriodogram: Estimate PSD using a modified periodogram.
%
%   [p,f] = ModifiedPeriodogram(x,fs,wn,nf,pf);
%
%   x    Input signal.
%   fs   Sample rate (Hz). Default = 1 Hz.
%   wn   Window to use. Default = Blackman. Must be same length as x.
%   nf   Number of frequencies to evaluate.
%        Default = max(128,round(length(x)/2)).
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   p    Power spectral density. 
%   f    Frequencies at which p is estimated (Hz).
%
%   Estimates the power spectral density (PSD) of an input signal 
%   using a Modified Periodogram. This method estimates the PSD by 
%   calculating the FFT of the windowed signal. The window 
%   multiplication in the time domain is equivalent to convolution 
%   (a form of smoothing) in the frequency domain and trades variance 
%   of the PSD estimate for increased bias. Note that, like the 
%   periodogram, this is not a consistent estimator and the 
%   variability of the estimate cannot be controlled by any of the 
%   user-specified parameters.
%
%   Since the PSD is an estimate of power, multiplication by the 
%   window in the time domain is equivalent to convolution by the
%   square of the window in the frequency domain. The mean is removed 
%   from the signal prior to estimation to prevent an impulse at 0 Hz 
%   from dominating the smoothed estimate.
%
%   The estimated PSD is scaled such that Parseval's relation is 
%   approximately approximately satisfied: 
%                          +pi
%      var(x) ~= inv(2*pi) int p(w) dw ~= sum(p)/length(p).
%                          -pi
%
%   Example: Estimate the PSD of an electrocardiogram signal and 
%   plot the results. 
%
%      load NOISYECG.mat;
%      x  = decimate(noisyecg(1:50e3),5);
%      fs = fs/5;
%      ModifiedPeriodogram(x,fs,[],5000);
%
%   M. Hayes, Statistical Digital Signal Processing and Modeling. 
%   New York: John Wiley & Sons, 1996, pp. 408-412.
%
%   J. G. Proakis, C. M. Rader, F. Ling, C. L. Nikias, M. Moonen, 
%   and I. K. Proudler, Algorithms for Statistical Signal Processing.
%   Saddle River, NJ: Prentice Hall, 2002, pp. 447-448.
%
%   Version 1.00 JM
%
%   See also SPECTRUM, WINDOW, and SpectralAnalysis.
