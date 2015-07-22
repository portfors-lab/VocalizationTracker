function [cs,f] = CrossSpectra(x1,x2,fsa,wla,ola,fra,pfa)
%CrossSpectra: Computes the cross-spectra of two signals.
%
%   [cs,f] = CrossSpectra(x1,x2,fs,wl,ol,fr,pf)
%
%   x1   Input Signal
%   x2   Input Signal
%   fs   Sample Rate, hertz.
%   wl   Length of each of the windows (seconds).  Default=1001
%   ol   Percent window overlap.  Default=50.
%   fr   Minimum no. of frequencies to estimate. Default=200.
%   pf   Plot flag:  0=none (default), 1=screen.
%
%   cs   Estimated Crossspectra
%   f    Frequencies
%
%   CrossSpectra estimates the cross power spectral density of 
%   of two equal length input signals, x1 and x2.  The cross power 
%   spectral density is calculated using Welch's method of the 
%   averaged periodogram defined as:
%
%         cs = abs(X1 .* conj(X2))
%
%   where X1 and X2 are zero padded Fourier transforms of x1 and x2.  
%
%   This function removes the mean of each signal prior to spectral 
%   estimation.  If only the window length is specified, the 
%   blackman window is used. The overlap argument must be a number
%   between 0 and 100.  If no output arguments are specified, the 
%   output will plot to the screen.
%
%   Example: Calculate the cross-spectra of an ABP data segment and 
%   an ICP data segment with a sample rate of 125 Hz, window length
%   of 8 s, 50% overlap, hamming window, and plot the results to the 
%   screen.
%
%      load ABPICP.mat
%      x1 = abp(1:4000);
%      x2 = icp(1:4000);
%      CrossSpectra(x1,x2,125,hamming(8*125));
%
%   Challis, R. E., and Kitney, R. I., "Biomedical Signal Processing
%   (in four parts), Part 3, The Power Spectrum and Coherence 
%   Function", Med. & Biol. Eng. & Comput., 1991, 29, pp. 225-241.
%
%   Version 1.00 LJ
%
%   See Also AutoSpectra and Coherence.
