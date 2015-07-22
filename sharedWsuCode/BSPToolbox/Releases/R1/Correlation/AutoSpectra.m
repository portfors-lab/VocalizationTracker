function [as,f] = AutoSpectra(x,fsa,wla,ola,fra,pfa)
%AutoSpectra: Estimates the autospectra of a signal.
%
%   [as,f] = AutoSpectra(x,fs,wl,ol,fr,pf)
%
%   x    Input Signal
%   fs   Sample Rate (Hz).  Default=1 Hz.
%   wl   Length of each of the windows (seconds).  Default=1001
%   ol   Percent window overlap.  Default=50.
%   fr   Minimum no. of frequencies to estimate. Default=200.
%   pf   Plot flag:  0=none (default), 1=screen.
%
%   as   Estimated Autospectra
%   f    Frequencies (Hz) at which autospectra is estimated 
%
%   AutoSpectra estimates the power spectral density of a stationary
%   input signal 'x'.  The power spectral density is calculated 
%   using Welch's method of the averaged periodogram defined as:
%
%         as = sqrt(Sxx)
%
%   where Sxx is the power spectral density of the signal x.  
%
%   This function removes the mean prior to spectral estimation.  If
%   only the window length is specified, the blackman window is used.
%   The overlap argument must be a number between 0 and 100.  If no 
%   output arguments are specified, the output will plot to the 
%   screen.
%
%   Example: Calculate the autospectra of an ABP signal with a 
%   a sample rate of 125 Hz, window length of 8 s, 50% overlap,
%   hanning window, and plot the results to the screen.
%
%      load ABPICP.mat
%      [as,f] = AutoSpectra(abp,125,hanning(8*125),[],[],1);
%
%   Challis, R. E., and Kitney, R. I., "Biomedical Signal Processing
%   (in four parts), Part 3, The Power Spectrum and Coherence 
%   Function", Med. & Biol. Eng. & Comput., 29, pp.225-241, 1991.
%
%   Version 1.00 LJ
%
%   See Also CrossSpectra, Coherence, Spectrogram, and 
%   CrossSpectrogram.
