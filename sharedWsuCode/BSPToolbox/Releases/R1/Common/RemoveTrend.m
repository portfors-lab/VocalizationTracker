function [y,n] = RemoveTrend(x,whpa,cfa,pfa)
%RemoveTrend: Highpass filter used to detrend data
%
%   [y,n] = RemoveTrend(x,whpa,tpa,cfa,pfa)
%
%   x     Input signal       
%   whp   Normalized cutoff frecuency 0.0 < whp < 1.0, with 1.0 
%         correponding to half the sampe rate. Default = 0.5     
%   cf    Causality flag: 1 = causal, 2 = noncausal (default) for tp=1
%   pf    Plot flag: 0=none (default), 1=screen
%
%   y    Filtered Signal
%   n    Order of the filter
%
%   Filters the input signal x with a cutoff frequency whp using an
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
%   Example: Detred the ICP signal
%
%      load ICP; 
%      [y,n] = RemoveTrend(icp,0.01,2,1);
%
%   Version 1.00 MA
%
%   See also HighPass, Smooth, filter, filtfilt, ellip, and butter.
