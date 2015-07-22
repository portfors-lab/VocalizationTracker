function [y] = PeakRestore(x,fsa,fca,tha,ppa,gfa,pfa);
%PeakRestore: Nonlinear lowpass filter without peak attenuation
%
%   [y] = PeakRestore(x,fs,fc,th,pp,gc,pf);
%
%   x    Input signal.
%   fs   Sample rate (Hz). Default=2 Hz.
%   fc   Cutoff frequency (Hz). Default=fs/4.
%   th   Vector containing the two peak thresholds or indices. 
%   pp   Peak polarity. 1=Positive peaks (default), 2=Negative peaks,
%        3 = Both.
%   gf   Gain filter coefficient. Default=0.9.
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   y    Filter output.
%  
%   Linear-phase lowpass filters are often used as a preprocessing 
%   step in signal processing applications to eliminate high 
%   frequency components of noise that do not overlap with the signal 
%   spectrum. Signals with sharp peaks confound this type of noise 
%   reduction because sharp peaks have a broad spectrum and are 
%   significantly attenuated by lowpass filters. This implements 
%   a simple method that restores the signal peaks to their full 
%   amplitude by adding a filtered estimate of the peak residuals to 
%   the original lowpass filter output.
%
%   The fourth input argument is used differently depending on how
%   long the th vector is.
%
%   Scalar:           Any peak that is above th is labeled a peak.
%   2-element vector: Any peak that is between the two elements is
%                     labeled a peak.
%   n-element vector: vector containing the indices of the peaks.
%                     This enables the function to be used with any
%                     peak detector.
%
%   Example: Apply a linear and nonlinear filter to an 
%   electrocardiogram signal.
%
%      load('ECG.mat');   
%      x   = ecg(1:10000);
%      x   = Highpass(x,fs,5); % Drift filter
%      xlp = Lowpass(x,fs,20,4,2);
%      xpr = PeakRestore(x,fs,20);
%      k   = 1:length(x);
%      t   = (k-1)/fs;
%      plot(t,x,'g',t,xlp,'r',t,xpr,'b');
%      xlim([9.4 10.1]);
%      legend('No drift','Linear Lowpassed','Peak Restored',4);
%
%   J. McNames and B. Goldstien, "A Nonlinear Lowpass Filter that 
%   Eliminates Peak Attenuation," ICASSP 2002, in press.
%
%   Version 1.00 JM
% 
%   See also Lowpass.
