function [p,f] = Welch(x,fsa,wla,ola,nfa,pfa)
%Welch: Estimate PSD using the Welch's method.
%
%   [p,f] = Welch(x,fs,wl,ol,nf,pf)
%
%   x    Input signal.
%   fs   Sample rate (Hz). Default = 1 Hz.
%   wl   Length of window to use (sec). Default = (max. possible).
%        If a vector, specifies actual window.
%   ol   Percent window overlap.  Default = 50%.
%   nf   Number of frequencies to evaluate.
%        Default = max(128,round(wl/2)).
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   p    Power spectral density. 
%   f    Frequencies at which p is estimated (Hz).
%
%   Estimates the power spectral density (PSD) of an input signal 
%   using Welch's method. This method estimates the PSD by 
%   calculating the FFT of overlapping windowed segments of the 
%   signal. The window multiplication in the time domain is 
%   equivalent to convolution (a form of smoothing) in the frequency 
%   domain and trades variance of the PSD estimate for increased 
%   bias. 
%
%   Since the PSD is an estimate of power, multiplication by the 
%   window in the time domain is equivalent to convolution by the
%   square of the window in the frequency domain. Thus, it may be
%   more intuitive to multiply the signal by the square root of
%   the popular windows.
%
%   The mean is removed from the signal prior to estimation to 
%   prevent an impulse at 0 Hz from dominating the smoothed estimate.
%   If only the window length is specified, the square root of the 
%   Blackman window is used. 
%
%   The estimated PSD is scaled such that the estimate is 
%   asymptotically unbiased and Parseval's relation is approximately 
%   satisfied:
%                          +pi
%      var(x) ~= inv(2*pi) int p(w) dw ~= sum(p)/length(p).
%                          -pi
%
%   Example: Estimate the PSD of an electrocardiogram signal and 
%   plot the results. Include at least 5000 points in the PSD 
%   estimate and use a window length of 10 seconds.
%
%      load NOISYECG.mat;
%      x  = decimate(noisyecg(1:50e3),5);
%      fs = fs/5;
%      Welch(x,fs,10,[],5000);
%
%   M. Hayes, Statistical Digital Signal Processing and Modeling. 
%   New York: John Wiley & Sons, 1996, pp. 415-420.
%
%   J. G. Proakis, C. M. Rader, F. Ling, C. L. Nikias, M. Moonen, 
%   and I. K. Proudler, Algorithms for Statistical Signal Processing.
%   Saddle River, NJ: Prentice Hall, 2002, pp. 447-449.
%
%   Version 1.00 JM
%
%   See also SPECTRUM, WINDOW, and SpectralAnalysis.
