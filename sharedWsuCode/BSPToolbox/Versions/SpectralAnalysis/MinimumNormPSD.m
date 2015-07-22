function [psd] = MinimumNormPSD(x,fsa,pa,ma,pfa)
%MinimumNorm: Estimate PSD using the Minimum Norm Algorithm
%
%   [psd,f] = MinimumNormPSD(x,fs,p,m,pf)        
%
%   x     Input signal    
%   fs    Sampling frequency (Hz). Default=125      
%   p     Order of variance estimate. Default=2
%   m     Order of eigen decomposition, must be > p+1. Default=5
%   pf    Plot format: 0=none, 1=screen. Default=0
%
%   psd   Power spectral density 
%   f     Frequencies (Hz) at which p is estimated 
%
%   Estimates the power spectral density (PSD) of an input signal 
%   using the Minimum Norm Decomposition method. The PSD of the 
%   input signal is by decomposition of the autocorrelation 
%   sequence in eigen vectors. The minimum norm of the vectors 
%   are used to estimate the PSD.
%
%   The user can specify the type of window (blackman, bartlett, 
%   or hanning or others), the sampling frequency, 
%   and whether or not to plot the results. 
%
%   Example: Estimate the PSD of an intracranial pressure signal 
%   sampled at 125 Hz and plot the results.
%
%      load ICP.mat; 
%      MinimumNormPSD(icp);
%
%   Hayes M., "Statistical Digital Signal Processing and Modeling," 
%   John Wiley & Sons, 1996, pp.466
%
%   Version 1.00 JB
%
%   See also EigenVectorPSD, HeartRate,
%   BlackmanTukey, MaximumEntropyPSD, ModifiedPeriodogram
%   Music, Welch, and Wigner.

% ================================================
% Default parameters
% ==================================================
if ( nargin < 1 | nargin > 5)
    help MinimumNormPSD
    break;
end;

fs = 125;
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
    end;    

p = 2;
if exist('pa') & ~isempty(pa),
    p = pa;
    end;    
   
m = 5;
if exist('ma') & ~isempty(ma),
    m = ma;
    end;    

pf  = 0;                                % Default - no plotting
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end; 


    
% ==================================================
% Default parameters
% ================================================
x       = x(:);
n       = length(x);
n2      = ceil(log2(n));
n3      = 0:n-1;
N       = 2^n2;
f       = n3*(1/(N*(1/fs)));
fi      = f(1:N/2);

if  m < p+1 | length(x) < m, error('Size of R is inappropriate'), end
R       = SetCovarMatrix(x,m);
[v,d]   = eig(R);
[y, i]  = sort(diag(d));
V       = [];
for   j = 1:m-p
    V   = [V, v(:,i(j))];
end;
a       = V*V(1,:)';
Px      = -20*log10(abs(fft(a, N)));
Px      = Px(1:N/2);
psd     = Px;
f       = fi;

if (pf ==1 | nargout == 0)
    figure(1);
    FigureSet(1);
   stem(fi, abs(Px));grid
   title('Minimum norm based PSD')
   xlabel('Hz');
   ylabel('PSD');
end


