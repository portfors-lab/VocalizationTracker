function [tf,f] = TransferFunction(x,y,fsa,wla,ola,fra,pfa)
%TransferFunction:  Estimates transfer function of two signals.
%
%   [tf,f] = TransferFunction(x,y,fs,wl,ol,fr,pf)
%
%   x    Input Signal
%   y    Output Signal
%   fs   Sample rate (Hz). Default = 1 Hz.
%   wl   Length of window to use (sec). Default = 1001 samples.
%        If a vector, specifies entire window.
%   ol   Percent window overlap. Default=50.
%   fr   Minimum no. of frequencies to estimate. Default=275.
%   pf   Plot format: 0=none (default), 1=screen.
%
%   tf   Estimated transfer function.
%   f    Frequencies at which tf is estimated.
%   
%   TransferFunction estimates the transfer function of a system
%   where x is the input signal and y is the output signal.  This
%   function uses Welch's method of the averaged periodogram and is
%   defined as
%         
%         tf = Sxy./abs(Sxx)
%
%   where Sxy is the cross-power spectral density and Sxx is the 
%   power spectral density of the signal x.
%   
%   This function removes the mean and detrends the signal windows
%   prior to spectral estimation.  If only the window length is
%   specified, the blackman window is used.  The overlap argument 
%   must be a number between 0 and 100.  If no output arguments are
%   specified, the output will plot the magnitude transfer function
%   and the phase to the screen.
%   
%   Example: Calculate the transfer function using ICP as the input
%   signal and ABP as the output signal, decimated to 25 Hz, using
%   a hanning window with a window length of 30 s.  Plot the
%   results to the screen.
%
%      load ABPICP.mat
%      x = decimate(icp,5);
%      y = decimate(abp,5);
%      [tf,f] = TransferFunction(x,y,25,hanning(30*25),[],[],1);
%
%   Bendat, J., and Piersol, A., "Engineering Applications of
%   Correlation and Spectral Analysis", John Wiley & Sons, pp.97-114,
%   1980.
%
%   Version 1.00 LJ
%
%   See Also AutoSpectra and CrossSpectra.
