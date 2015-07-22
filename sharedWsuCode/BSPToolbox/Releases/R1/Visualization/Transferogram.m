function [TF,PH,t,f] = Transferogram(x,y,fsa,wla,fra,nfa,nsa,pfa);
%Transferogram: Nonstationary visualization of transfer function.
%
%   [TF,PH,t,f] = Transferogram(x,y,fs,wl,fr,nf,ns,pf);
%
%   x    Input signal.
%   y    Input signal.
%   fs   Sample rate, hertz.  Default = 1 Hz.
%   wl   Length of window to use (sec). Default = 1000 samples.
%        If a vector, specifies entire window.
%   fr   Array of minimum and maximum frequencies.  
%        Default=[0 fs/2].
%   nf   Number of frequencies to evaluate, default=wl/2.
%   ns   Requested number of times (horizontal pixels) to evaluate. 
%        Default=min(2^10,nx).
%   pf   Plot flag:  0=none (default), 1=screen.
%
%   TF   Transfer function magnitude output matrix.
%   PH   Transfer function phase output matrix.
%   t    Times at which the transfer function were evaluated (s).
%   f    Frequencies at which the transfer function were evaluated 
%        (Hz).
%
%   Plots the estimated transfer function of two non-stationary 
%   signals, x and y, x being the input signal and y being the 
%   output signal.  The transfer function is calculated using
%   modified periodogram and is defined as:
%
%         tf = Sxy./abs(Sxx)
%
%   where Sxy is the cross-power spectral density and Sxx is the 
%   power spectral density of the signal x.
%
%   The mean of each signal is removed before taking the FFT.  If only
%   the window length is specified, the blackman window is used.  The
%   output, when plotted, will display two surface plots, one for the
%   magnitude of the transfer function and one for the phase.  If no
%   output arguments are specified, the default will print to the 
%   screen.
%
%   Example:  Plot the transfer function of an ABP and ICP signal, 
%   which are decimated to 25 Hz, where icp is the input signal and
%   abp is the output signal.  Use a Hanning window with a window 
%   length of 50 s and plot the output to the screen.
%
%      load ABPICP.mat
%      y = decimate(abp,5);
%      x = decimate(icp,5);
%      Transferogram(x,y,25,hanning(50*25));  
%
%   Bendat, J., and Piersol, A., "Engineering Applications of
%   Correlation and Spectral Analysis", John Wiley & Sons, pp.97-114,
%   1980.
%
%   Version 1.00 LJ
%
%   See also Spectrogram and Cohereogram.
