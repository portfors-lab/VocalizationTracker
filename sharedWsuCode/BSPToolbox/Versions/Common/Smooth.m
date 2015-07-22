function [y,n] = Smooth(x,wlpa,tpa,cfa,pfa)
%Smooth: Smoothing filter (lowpass)
%
%   [y,n] = Smooth(x,wlp,tp,cf,pf)
%
%   x     Input signal          
%   wlp   Normalized cutoff frecuency 0.0 < wlp < 1.0, with 1.0 
%         correponding to half the sampe rate. Default = 0.5 
%   tp    Type: 1=Elliptic (default), 2=Butterworth, 
%         3=FIR based on Blackman Window, 4=Minimun Ringing
%   cf    Causality flag: 1 = causal, 2 = noncausal (default) for tp=1,2
%   pf    Plot format. 0=none (default), 1=screen
%
%   y     Filtered Signal
%   n     Order of the filter
%
%   Filters the input signal x with a normalized cutoff frequency 
%   wlp. The lowpass filter can be causal or noncausal. The causal 
%   implementation usesonly the present and previous values to determine
%   the filter's output y, and therefore it is physically realizable for
%   real time processing. The noncausal implementation filters the data 
%   in the forward direction, and the filtered sequence is then reversed
%   and run back through the filter; Y is the time reverse of the output
%   of the second filtering operation.  The result has precisely zero 
%   phase distortion and magnitude modified by the square of the 
%   filter's magnitude response.     
%
%   Example: Smooth the ICP waveform:
%
%      load ICP; 
%      [y,n] = Smooth(icp,0.25,1,2,1);
%
%   Version 1.00 MA
%
%   See also Lowpass, HighPass, filter, filtfilt, ellip, and butter.

%=====================================================================
% Process function arguments
%=====================================================================
if nargin<1 | nargin>5,
    help Smooth;
    return;
    end;

wlp = 0.5;                              % Default sampling rate
if exist('wlpa') & ~isempty(wlpa),
    wlp = wlpa;
    end;
   
tp = 1;                                % Default type
if exist('tpa') & ~isempty(tpa),
    tp = tpa;
    end;
    
cf = 2;                                % Default flag
if exist('cfa') & ~isempty(cfa),
    cf = cfa;
    end;

pf = 0;                                % Default - no plotting
if nargout==0,                         % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;

%=====================================================================
% Process Inputs
%=====================================================================
x  = x(:);
LD = length(x);
k  = 1:LD;

%=====================================================================
% LowPass Filtering
%=====================================================================  
[y, n] = LowPass(x, 2, wlp, tp, cf, pf);