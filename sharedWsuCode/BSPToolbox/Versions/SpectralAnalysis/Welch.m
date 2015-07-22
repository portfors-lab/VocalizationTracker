function [psd,f] = Welch(x,fsa,wa,laga,overa,pa,pfa)
%Welch: Estimate the PSD using the welch method
%
%   [psd,f] = Welch(x,fs,w,lag,over,tp,pf)        
%
%   x      Input signal    
%   fs     Sampling frequency (Hz). Default=125      
%   w      Window type: 1=Blackman, 2=Bartlett, 3=Hanning. Default=1 
%   lag    Length of each overlaping signals. Default=50
%   over   Overlapping percentage. Default=10%
%   tp     Portion of serie to be tapered. Default=50%
%   pf     Plot format: 0=none, 1=screen. Default=0
%
%   psd    Power spectral density 
%   f      Frequencies (Hz) at which p is estimated 
% 
%   Estimates the power spectral density (PSD) of an input signal 
%   using the Welch method. The PSD of the input signal is 
%   estimated by averaging the FFT of the lag-point data window.
%   The user can specify the type of window (blackman, bartlett, 
%   or hanning), the sampling frequency, and whether or not 
%   to plot the results. 
%
%   If w is a scalar, then it specifies the window type as shown 
%   above. If w is a vector, then w is treated as a user-specified 
%   window. If it is a different length than x, zero padding is used 
%   to make the windowing possible.
%
%   Example: Estimate the PSD of an intracranial blood pressure 
%   signal sampled at 125 Hz using a Blackman window and plot the 
%   results
%
%      load ICP.mat; 
%      [psd,f] = Welch(icp,125,1);
%
%   Hayes M., "Statistical Digital Signal Processing and Modeling," 
%   John Wiley & Sons, pp.415-419, 1996.
%
%   Version 1.00 JB
%
%   See also BlackmanTukey, EigenVectorPSD, HeartRate,
%   MaximumEntropyPSD, MinimumNormPSD, ModifiedPeriodogram, Music,
%   and Wigner.



% =================================================================
% Default parameters
% =================================================================
if ( nargin < 2 | nargin >7)
    help Welch;
    return;
end;


fs = 125;                       % Default sampling rate, Hz
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
    end;

lag  = 50;
if exist('laga') & ~isempty(laga),
    lag = laga;
    end;

over  = 0.1;
if exist('overa') & ~isempty(overa),
    over = overa;
    end;

w = 1;
if exist('wa') & ~isempty(wa),
    w = wa;
    end;

tp  = 0.5;
if exist('tpa') & ~isempty(tpa),
    tp = tpa;
    end;

pf  = 0;                         % Default - no plotting
if exist('pfa'),
    pf = pfa;   
else
    if  ~isempty('pfa') 
    pf = 1;
end
end
   
 


% =================================================================
% Estimate PSD
% =================================================================
  n = length(x);
  n            = max(n, size(w,1));
  xx           = zeros(1,n);
  xx(1:length(x)) = x;
  x               = xx;

  n1 = 1;
  n0 = (1-over)*lag;
  nsect = 1+floor((n-lag)/(n0));
  Px    = 0;
  
  for i =1:nsect
      [per, fqs] = ModifiedPeriodogram(x(n1:n1+lag-1), fs, w, [], 0);
      Px = Px + per/nsect;
      n1 = n1 + n0;
  end
  psd = Px;

% =================================================================
% Plotting (TBC by James)
% =================================================================
if (pf ==1 | nargout == 0)
    figure(1);
    FigureSet(1);
   stem(fqs, Px);grid
   title('Welch Periodogram');
   xlabel('Hz');
   ylabel('PSD');
end