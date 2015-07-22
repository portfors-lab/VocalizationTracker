function [pi,y] = PowerPeaks(x,fsa,fca,pa,xma);
%PowerPeaks: Finds peaks in the lowpass filtered signal power 
%
%   [pi,x] = PowerPeaks(x,fs,fc,pw,xm);
%
%   x     Input signal.
%   fs    Signal sample rate (Hz). Default = 125 Hz.      
%   fc    Cutoff frequency (Hz). Default = fs/4 Hz.
%   pw    Energy power. Default = 2.
%   xmn   Minimum value necessary to qualify as an energy peak. 
%         Default = 0.
%   pf    Plot flag: 0=none (default), 1=screen.
% 
%   pi    Vector of signal indices that contain the peaks.
%   y     Smoothed signal energy.
%
%   Returns the peaks in the local smoothed signal energy. Signal energy 
%   is defined here as the absolute value of the signal to the power of p. 
%   sorted vector x that are at least as big as xm (optional). Uses a 
%   Lowpass filter type 4 to prevent ringing after impulses.
%
%   Example: Find the energy peaks in an electrocardiogram signal.
%
%      load ECG;
%      x  = decimate(ecg,10);
%      [pi,y] = PowerPeaks(x,fs/10,2,4);
%      k = 1:length(y);
%      t = (k-1)/fs;
%      plot(t,x,'k',t,y,'b',t(pi),y(pi),'r.','MarkerSize',20);
%      xlim([10 10.5])
%      legend('Signal','Signal Power','Power Peaks');
%
%   Desmond J. Higham and Nicholas J. Higham, "MATLAB Guide," 2000.
%
%   Version 1.00 JM
%
%   See also DetectMaxima, Lowpass, and ECGDetectREnergy.
