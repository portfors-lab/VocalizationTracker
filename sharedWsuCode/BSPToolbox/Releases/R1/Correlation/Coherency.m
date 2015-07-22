function [c,f] = Coherency(x1,x2,fsa,wla,sna,ola,nfa,pfa)
%Coherency: Computes the estimated coherence of two signals.
%
%   [c,f] = Coherency(x1,x2,fs,wl,sn,ol,nf,pf)
%
%   x1   Input Signal
%   x2   Input Signal
%   fs   Sample rate (Hz). Default = 1 Hz.
%   wl   Length of window to use (sec). Default = (duration of x)/10.
%        If a vector, specifies entire window.
%   sn   Signal to noise ratio used to bias the coherence. 
%        Default = inf.
%   ol   Percent window overlap.  Default=50.
%   nf   Minimum no. of frequencies to estimate. Default = 200.
%   pf   Plot format: 0=none (default), 1=screen. 
%
%   c    Estimated coherency spectrum
%   f    Frequencies at which coh is estimated
%
%   The coherency spectrum is a function that gives an output between 
%   0 and 1 and is a measure of the correlation of frequency 
%   components of random signals x1 and x2. The coherency spectrum
%   is defined as
%
%        c = abs(Sxy)./sqrt(Sxx.*Syy)).
%
%   where Sxy is the cross-power spectral density and Sxx is the 
%   power spectral density of the signal x. If x and y are zero-mean 
%   white noise processes with a coefficient of correlation that is 
%   equal to p,
%
%        p = E[x*y]/sqrt(Var[x]*Var[y]),
%
%   the coherency is related to p by c = p^2.
%
%   This function removes the mean prior to spectral estimation. The
%   spectral components are estimated using Welch's method. If only 
%   the window length is specified, the Blackman window is used. The 
%   overlap argument must be a number between 0 and 100. The 
%   magnitude squared length FFT is calculated for each window and 
%   averaged to estimate Sxy, Sxx, and Syy. 
%   
%   Example: Calculate the coherence of an ABP data segment and an ICP
%   data segment, each with a sample rate of 125 Hz and 10 s data 
%   windows, using Hanning windows.  Plot the results to the screen.
%
%      load ABPICP.mat
%      x1  = abp(1:10000);
%      x2  = icp(1:10000);
%      [c,f] = Coherency(x1,x2,125,hanning(10*125),[],[],[],1);
%
%   Challis, R. E., and Kitney, R. I., "Biomedical Signal Processing
%   (in four parts), Part 3, The Power Spectrum and Coherence 
%   Function", Med. & Biol. Eng. & Comput., 29, pp.225-241, 1991.
%
%   Priestley, "Spectral Analysis and Time Series," Academic Press,
%   pp.556-557, 1981.
%
%   Version 1.00 JM
%
%   See Also AutoSpectra, CrossSpectra, and Cohereogram.
