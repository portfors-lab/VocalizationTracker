function [y,n] = HighPass(x,fsa,fca,cfa,pfa)
%HighPass: Highpass filter
%
%   [y,n] = HighPass(x,fs,fc,cf,pf)
%
%   x    Input signal       
%   fs   Signal sample rate (Hz). Default=125 Hz      
%   fc   Cutoff frequency (Hz). Default=fs/4 Hz 
%   cf   Causality flag: 1 = causal, 2 = noncausal (default) for tp=1
%   pf   Plot flag:hig 0=none (default), 1=screen
%
%   y    Filtered Signal
%   n    Order of the filter
%
%   Filters the input signal x with a cutoff frequency fc using an
%   elliptical filter. The highpass filter can be causal or noncausal.
%   The causal implementation uses only the present and previous values
%   to determine the filter's output y, and therefore it is physically
%   realizable for realtime processing. The noncausal implementation 
%   filters the data in the forward direction, and the filtered sequence
%   is then reversed and run back through the filter; Y is the time 
%   reverse of the output of the second filtering operation.  The result
%   has precisely zero phase distortion and magnitude modified by the 
%   square of the filter's magnitude response.     
%
%   Example: Filter the raw intracranial pressure signal using a 
%   highpass filter with zero phase (noncausal) and with cutoff 
%   frequency 0.5 Hz. This will filter out the low frequency
%   components (frequencies below 0.5 Hz) and detrend the data:
%
%      load ICP; 
%      [y,n] = Highpass(icp,fs,0.5,2,1);
%
%   Version 1.00 MA
%
%   See also HighPass, Lowpass, filter, filtfilt, ellip, and butter.
