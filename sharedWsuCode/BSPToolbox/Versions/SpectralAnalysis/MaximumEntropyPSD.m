function [psd, f] = MaximumEntropyPSD(x,fsa,ora,pfa)
%MaximumEntropyPSD: Estimate the PSD using the Maximum Entropy Method
%
%   [psd,f] = MaximumEntropyPSD(x,fs,or,pf)        
%
%   x     Input signal    
%   fs    Sampling frequency (Hz). Default=125      
%   or    Order of variance estimate. Default=5
%   pf    Plot format: 0=none, 1=screen. Default=0
%
%   psd   Power spectral density 
%   f     Frequencies (Hz) at which p is estimated 
%
%   Estimates the power spectral density (PSD) of an input signal 
%   using the Maximum Entropy method. The PSD of the input 
%   signal is estimated by extrapolating the autocorrelation sequence
%   recursively and evaluating the spectrum of the given 
%   autocorrelation sequence. Order can vary between 2 to 5 without 
%   adding too much computational time. The user can specify whether
%   or not to plot the results. 
%
%   Example: Estimate the PSD of an intracranial pressure signal 
%   sampled at 125 Hz using an order = 2 and plot the results using
%   the maximum entropy method.
%
%      load ICP.mat; 
%      MaximumEntropyPSD(icp, 125, 2);
%
%   Hayes M., "Statistical Digital Signal Processing and Modeling," 
%   John Wiley & Sons, pp.437, 1996
%
%   Version 1.00 JB
%
%   See also EigenVectorPSD, HeartRate,
%   BlackmanTukey, MinimumNormPSD, ModifiedPeriodogram
%   Music, Welch, and Wigner.


% =======================================
% Default parameters
% =====================================
if ( nargin < 1 | nargin > 4)
    help MaximumEntropyPSD;
    break;
end;


fs = 125;
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
    end;    

or = 2;
if exist('ora') & ~isempty(ora),
    or = ora;
    end;    
    
pf  = 0;                                % Default - no plotting
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end; 

    
% ============================================
% Estimate PSD
% ============================================
n       = length(x);
n2      = ceil(log2(n));
n3      = 0:n-1;
N       = 2^n2;
f       = n3*(1/(N*(1/fs)));
fi      = f(1:N/2);
[a,e]   = EstimateAllPoleModel(x,or);
Px      = 20 * log10(e)-10*log10(abs(fft(a,N)));
Px      = Px(1:N/2);
f = fi;
psd = Px;
f = f';

if (pf ==1 | nargout == 0)
   figure(1);
   FigureSet(1);
   stem(fi, abs(Px));grid
    title('Maximum Entropy Method based PSD')
    xlabel('Hz');
   ylabel('PSD');
end

