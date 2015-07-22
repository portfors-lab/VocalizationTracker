function [p,f] = BlackmanTukey(x,fsa,wla,nfa,pfa)
%BlackmanTukey: Estimate PSD using the Blackman-Tukey method.
%
%   [p,f] = BlackmanTukey(x,fs,wl,nf,pf);        
%
%   x    Input signal.
%   fs   Sample rate (Hz). Default = 1 Hz.
%   wl   Length of window to use (sec). Default = (max. possible).
%        If a vector, specifies actual window.
%   nf   Number of frequencies to evaluate.
%        Default = max(128,round(wl/2)).
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   p    Power spectral density. 
%   f    Frequencies at which p is estimated (Hz).
%
%   Estimates the power spectral density (PSD) of an input signal 
%   using the Blackman-Tukey method. This method estimates the
%   autocorrelation sequence, multiplies the autocorrelation sequence
%   by a window, and then calculates the Fourier transform of the
%   windowed autocorrelation. The window multiplication in the time
%   domain is equivalent to convolution (a form of smoothing) in the
%   frequency domain and trades variance of the PSD estimate for
%   increased bias. 
%
%   This implementation uses an unbiased estimate of the 
%   autocorrelation (calculated via the FFT):
%              1   N-k
%      r(k) = ---  sum x(m) x(m+k)
%             N-k  m=0
%   Thus, the expected value of the estimated PSD is related to the 
%   true PSD by a single convolution operation with one window. The 
%   mean is removed from the signal prior to estimation to prevent an 
%   impulse at 0 Hz from dominating the smoothed estimate.
%
%   The estimated PSD is scaled such that the estimate is 
%   asymptotically unbiased and Parseval's relation is approximately 
%   satisfied:
%                          +pi
%      var(x) ~= inv(2*pi) int p(w) dw ~= sum(p)/length(p).
%                          -pi
%
%   If only the window length is specified, the blackman window is 
%   used. The specified window length should be odd. If the specified
%   window length is even, 1 is added to make it odd. If the window 
%   itself is specified with an even number of elements, a zero is 
%   appended to make the window odd.
%
%   Example: Estimate the PSD of an electrocardiogram signal and 
%   plot the results. Include at least 5000 points in the PSD 
%   estimate and use a window length of 10 seconds.
%
%      load NOISYECG.mat;
%      x  = decimate(noisyecg(1:50e3),5);
%      fs = fs/5;
%      BlackmanTukey(x,fs,10,5000);
%
%   M. Hayes, Statistical Digital Signal Processing and Modeling. 
%   New York: John Wiley & Sons, 1996, pp.420-423.
%
%   J. G. Proakis, C. M. Rader, F. Ling, C. L. Nikias, M. Moonen, 
%   and I. K. Proudler, Algorithms for Statistical Signal Processing.
%   Saddle River, NJ: Prentice Hall, 2002, pp. 449-452.
%
%   Version 1.01 JM
%
%   See also SPECTRUM, WINDOW, SpectralAnalysis, and AutoCorrelate.
