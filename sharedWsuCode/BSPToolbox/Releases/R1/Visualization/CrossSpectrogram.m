function [CS,t,f] = CrossSpectrogram(x1,x2,fsa,wla,fra,nfa,nsa,pfa);
% CrossSpectrogram:  Estimate and plot the cross-spectrogram of two
% signals.
%
%   [CS,t,f] = CrossSpectrogram(x1,x2,fs,wl,fr,nf,ns,pf);
%
%   x1   Input signal.
%   x2   Input signal.
%   fs   Sample rate, hertz.  Default = 1 Hz.
%   wl   Length of window to use (sec). Default = 1000 samples.
%        If a vector, specifies entire window.
%   fr   Array of minimum and maximum frequencies.  
%        Default = [0 fs/2].
%   nf   Number of frequencies to evaluate, default = wl/2.
%   ns   Requested number of times (horizontal pixels) to evaluate. 
%        Default = min(2^10,nx).
%   pf   Plot flag:  0=none (default), 1=screen.
%
%   CS   CrossSpectrogram output matrix.
%   t    Times at which the CrossSpectrogram was evaluated (s).
%   f    Frequencies at which the CrossSpectrogram was evaluated (Hz).
%
%   Plots the estimated cross-spectrogram of two non-stationary 
%   signals, x1 and x2.  The cross-spectrum is calculated using
%   modified periodogram.  The mean of each signal is removed before 
%   taking the FFT.  If only the window length is specified, the 
%   blackman window is used.
%
%   The x-axis on the contour plot represents the time in seconds 
%   or minutes, depending on the length of the signal.  The y-axis
%   represents the frequency, in hertz.  The colorbar represents 
%   the magnitude of the cross-spectral density.  The two original 
%   signals are plotted beneath the spectrogram.  
%
%   Example:  Plot the Cross-Spectrogram of an ABP and ICP signal, 
%   which are decimated to 12.5 Hz, using a Hanning window and a
%   window length of 100 s.
%
%      load ABPICP.mat
%      x1 = decimate(abp,10);
%      x2 = decimate(icp,10);
%      CrossSpectrogram(x1,x2,12.5,hanning(100*12.5));  
%
%   Hayes, M., "Statistical Digital Signal Processing and Modeling," 
%   John Wiley & Sons, 1996.
%
%   Version 1.00 LJ
%
%   See also Spectrogram and Cohereogram.
