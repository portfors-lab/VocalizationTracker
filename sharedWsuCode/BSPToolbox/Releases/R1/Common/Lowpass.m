function [y,n] = LowPass(x,fsa,fca,fta,cfa,pfa);
%LowPass: Lowpass filter
%
%   [y,n] = LowPass(x,fs,fc,ft,cf,pf)
%
%   x    Input signal       
%   fs   Signal sample rate (Hz). Default=125 Hz      
%   fc   Cutoff frequency (Hz). Default=fs/4 Hz 
%   ft   Type: 1=Elliptic (default), 2=Butterworth, 
%        3=FIR based on Blackman Window, 4=Minimun Ringing
%   cf   Causality flag: 1 = causal, 2 = noncausal (default) for tp=1,2
%   pf   Plot flag: 0=none (default), 1=screen
%
%   y    Filtered Signal
%   n    Order of the filter
%
%   Filters the input signal x with a cutoff frequency fc using an 
%   elliptical butterworth filter. The lowpass filter can be causal 
%   or noncausal. The causal implementation uses only the present and
%   previous values to determine the filter's output y, and therefore
%   it is physically realizable for realtime processing. The noncausal 
%   implementation filters the data in the forward direction, and the 
%   filtered sequence is then reversed and run back through the filter.
%   The result has precisely zero phase distortion and magnitude 
%   modifiedby the square of the filter's magnitude response.     
%
%   Example: Filter the raw intracranial pressure signal using an 
%   elliptic lowpass filter with zero phase (noncausal) and with 
%   cutoff frequency fs/4 Hz.This will filter out the high frequency 
%   components (frequencies above fs/4 Hz) and smooth the data.
%
%      load ICP; 
%      [y,n] = LowPass(icp,fs,fs/4,1,2,1);
%
%   Version 1.00 MA
%
%   See also Smooth, HighPass, filter, filtfilt, ellip, and butter.
