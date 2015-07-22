function [y,n] = Highpass(x,fsa,fca,cfa,pfa)
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

%=======================================================================
% Process function arguments
%=======================================================================

if nargin<1 | nargin>6,
    help Highpass;
    return;
    end;

fs = 125;                              % Default sampling rate, Hz
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
    end;
 
fc = fs/4;                             % Default cutoff frequency, Hz
if exist('fca') & ~isempty(fca),
    fc = fca;
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
    
%=======================================================================
% Process Inputs
%=======================================================================
x  = x(:);
LD = length(x);
k  = 1:LD;

%=======================================================================
% HighPass Filtering
%=======================================================================
nf  = fs/2;
whp = fc/nf; 
Wp  = whp*1.2;
Ws  = whp*0.8;
Rp  = 0.5;
Rs  = 40;
                                    
    if cf==1,                                    % causal
     	[n,Wn] = ellipord(Wp, Ws, Rp, Rs);
    	[B,A ] = ellip(n, Rp, Rs, Wn, 'high');
        y      = filter(B,A,x);
    else                                         % non-causal
    	[n,Wn] = ellipord(Wp, Ws, Rp, Rs);
    	[B,A ] = ellip(n, Rp, Rs, Wn, 'high');
        y      = filtfilt(B,A,x);
    end;
 


    
%=======================================================================
% Plotting
%=======================================================================
if pf == 1
    figure(1)
    figureset(1)
    h = plot(k./fs, x, 'b', k./fs, y, 'r');
    title('Raw Signal and Highpass Filtered Signal');
    xlabel('Time,s');
    ylabel('Amplitude');
    legend('Raw Signal', 'Highpass Filtered');
    box off; axisset
    
   % figure(2)
   % figureset(2)
   % freqz(B,A, 512, fs)
    

end;

%=======================================================================
% Take care of outputs
%=======================================================================
if nargout==0,
clear('y','n');
end;   