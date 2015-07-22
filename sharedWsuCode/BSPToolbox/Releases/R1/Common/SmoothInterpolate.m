function [y] = SmoothInterpolate(t,x,ti,wla,pfa);
%SmoothInterpolate: Kernel smoothing for non-uniform sampled signals.
%
%   [y] = SmoothInterpolate(t,x,ti,wl,pf);
%
%   t    Times of signal observations (sec).
%   x    Values of signal observations.
%   ti   Times to generate estimate of smoothed signal values (sec).
%   wl   Length of kernel window to use (sec). Default = 5 sec.
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   y    Estimated filtered signal values at times specified by ti.
%
%   Smoothes (lowpass filters) the signal specified by the (t,y) 
%   coordinates and evaluates the estimate at the times specified by
%   the vector ti. This function is faster than other kernel 
%   smoothing routines because it takes advantage of the fact that t 
%   and ti are sorted in increasing order (required). This routine 
%   uses a truncated guassian kernel with standard deviation 
%   specified by wl. Points more than 5 standard deviations away from 
%   the evaluation time, ti, are ignored. Although the units 
%   specified above are seconds, the signal times (t), estimation 
%   times (ti), and window length (wl) can actually be in any measure 
%   (e.g. samples), as long as they are consistent with one another.
%
%   Example: Do a smooth interpolation of interbeat intervals 
%   estimated from an electrocardiogram R detector at a constant 
%   sample rate of 5 Hz. Use a kernel width of 2 seconds.
%
%      load ECG.mat;
%      np = length(ecg);
%      ri = ECGDetectRInterbeat(ecg,fs,fs);
%      nr = length(ri);
%      ibi = diff(ri)/fs;
%      t   = (ri(1:nr-1) + ri(2:nr))/(2*fs);
%      fsi = 5;
%      ti  = 0:(1/fsi):((np-1)/fs); 
%      x   = SmoothInterpolate(t,ibi,ti,2,1);
%
%   M.P. Wand and M.C. Jones, Kernel Smoothing. New York: Chapman & 
%   Hall, 1995.
%
%   Version 1.00 JM
%
%   See also Detectors, Lowpass, and Smooth.
