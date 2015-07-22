function [p,f] = EstimatePSD(x,fsa,wla,wta,pfa)
%EstimatePSD: Modified periodogram power spectral density estimation 
%
%   [p,f] = EstimatePSD(x,fs,w,pf)
%
%   x    Input signal.       
%   fs   Signal sample rate (Hz). Default=1 Hz.      
%   wt   Window type: 1=Blackman (default), 2=Bartlett, 
%        3=Hanning, 4 = Hamming, 5 = rectangular
%   wl   Window length (Samples). Default=1024 Samples
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   p    Estimated power spectral density 
%   f    Frequencies (Hz) at which p is estimated 
%
%   Estimates the power spectral density (PSD) of an input signal 
%   using the modified periodogram method. The PSD of the input 
%   signal is calculated by taking the absolute value of the FFT of 
%   the windowed signal squared. The user can specify the type of 
%   window (blackman, bartlett, or hanning), the sampling frequency, 
%   and whether or not to plot the results. 
%
%   Example: Estimate the PSD of an arterial blood pressure signal 
%   sampled at 125 Hz using a Blackman window and plot the results:
%
%      load bspdata.m; 
%      [p,f] = EstimatePSD(ABP, 125, 2);
%
%   Hayes, M., "Statistical Digital Signal Processing and Modeling," 
%   John Wiley & Sons, pp.408-412, 1996.
%
%   Version 1.00 MA
%
%   See also EstimateHR, Periodogram, boxcar, chebwin, hamming, and 
%   kaiser.

%====================================================================
% Process function arguments
%====================================================================
if nargin<1 | nargin>4,
    help EstimatePSD;
    return;
    end;

fs = 125;                               % Default sampling rate, Hz
if exist('fsa') & ~isempty(fs),
    fs = fsa;
    end;
    
wt = 1;                                 % Default window type
if exist('wta') & ~isempty(wta),
    wt = wta;
    end;
    
wl = 1;                                 % Default window length
if exist('wa') & ~isempty(wa),
    wl = wla;
    end;
    
pf = 0;                                 % Default - no plotting
if nargout==0,                          % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;

x = x(:);


