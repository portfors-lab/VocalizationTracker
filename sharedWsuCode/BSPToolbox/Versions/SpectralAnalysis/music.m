function [psd,f] = music(x,fsa,pa,ma,pfa)
% Music: Estimate the PSD with the MUSIC algorithm
%
%   [psd,f] = music(x,fs,p,m,pf)
%
%   x     Input signal           
%   fsa   Sampling frequency (Hz). Default=125    
%   pa    Order of variance estimate. Default=2
%   ma    Order of eigen decomposition, M must be > p+1. Default=5
%   pf    Plot format: 0=none, 1=screen. Default=0
%  
%   psd   Power spectral density 
%   f     Frequencies (Hz) at which p is estimated 
%
%   This function estimates the power spectral density of an input 
%   signal using the music algothim. The user must specify the 
%   orders of the variance estimate and the eigen decompostion. 
%   Note that the order of eigen decomposition must be greater
%   than the order of variance estimate by at least one unit
%
%   Example: Estimate the PSD of an intracranial blood pressure 
%   signal sampled at 125 Hz using a Blackman window and plot the 
%   results
%
%      load ICP.mat; 
%      music(icp,fs);
%
%   Hayes M., "Statistical Digital Signal Processing and Modeling," 
%   John Wiley & Sons, pp.466, 1996.
%
%   Version 1.00 JB
%
%   See also BlackmanTukey, EigenVectorPSD, HeartRate, 
%   ModifiedPeriodogram, MaximumEntropyPSD, MinimumNormPSD, 
%   Welch, and Wigner.



% ===================================================================
% Default parameters
% ===================================================================
if ( nargin < 1 | nargin > 5)
    help music;
    return;
end;

fs = 125;
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
    end;    

p = 2;
if exist('pa') & ~isempty(pa),
    p = pArg;
    end;    
   
m = 5;
if exist('ma') & ~isempty(ma),
    m = ma;
    end;    

    
pf  = 0;                                % Default - no plotting
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end; 
    
% ===================================================================
% Estimate PSD
% ===================================================================
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
Px      = 0;
for   j = 1:m-p
    Px  = Px + abs(fft(v(:,i(j)),N));
end;
Px      = -20*log10(Px);

Px      = Px(1:N/2);
psd     = Px;
f       = fi;

% ==================================================================
% Plotting (TBC by James)
% ==================================================================
if (pf ==1 | nargout == 0)
    figure(1);
    FigureSet(1);
   stem(fi, abs(Px));grid
    title('Music based PSD')
   xlabel('Hz');
   ylabel('PSD');
end
