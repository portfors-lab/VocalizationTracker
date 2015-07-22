function [p,f,N] = EstimatePSD(x, fsa, wa, pfa)
%EstimatePSD Modified periodogram power spectral density estimation 
%
%   [p,f,N] = EstimatePSD(x, fs, w, pf)
%
%   x    Input signal       
%   fs   Signal sample rate (Hz), default=125 Hz      
%   w    Window type: 1=Blackman (default), 2=Bartlett, 3=Hanning 
%   pf   Plot format: 0=none (default), 1=screen
%
%   p    Estimated power spectral density 
%   f    Frequencies (Hz) at which p is estimated 
%   N    Window length used to calculate the FFT after zero-padding
%
%   Estimates the power spectral density (PSD) of an input signal 
%   using the modified periodogram method. The PSD of the input 
%   signal is calculated by taking the absolute value of the FFT of 
%   the windowed signal squared. The user can specify the type of 
%   window (blackman, bartlett, or hanning), the sampling frequency, 
%   and whether or not to plot the results. 
%
%   If w is a scalar, then it specifies the window type as shown 
%   above. If w is a vector, then w is treated as a user-specified 
%   window. If it is a different length than x, zero padding is used 
%   to make the windowing possible.
%
%   Example: Estimate the PSD of an intracranial pressure signal 
%   sampled at 125 Hz using a Blackman window and plot the results:
%      load ICP; 
%      [p,f] = EstimatePSD(icp, fs, 1, 1);
%
%   Hayes, M., "Statistical Digital Signal Processing and Modeling," 
%   John Wiley & Sons, pp.408-412, 1996.
%
%   Version 1.00
%
%   See also EstimateHR, Periodogram, boxcar, chebwin, hamming, and 
%   kaiser, blackman.

%--------------------------------------------------------------------
% Process function arguments
%--------------------------------------------------------------------
if nargin<1 | nargin>4,
    help EstimatePSD;
    return;
    end;

fs = 125;                               % Default sampling rate, Hz
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
    end;
 
w = 1;                                  % Default window
if exist('wa') & ~isempty(wa),
    w = wa;
    end;
    
pf = 0;                                 % Default - no plotting
if nargout==0,                          % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;

x  = x(:);
x  = x - mean(x);
n1 = 1;
n2 = length(x);
n3 = 0:n2-1;
n  = ceil(log2(n2));
N  = 2^n; 
w  = ones(N,1);

% ---------------------------------------------------------------------
% Select Window
% ---------------------------------------------------------------------
if w==1,
    w = blackman(n2); 
elseif w==2,
    w = bartlett(n2);
elseif w==3,
    w = hanning(n2);
    end;

% ---------------------------------------------------------------------
% Estimate PSD
% ---------------------------------------------------------------------
xw    = x.*w/norm(w);
PX    = abs(fft(xw, N)).^2/N;
hr    = n3*(1/(N*(1/fs)));
p     = PX(1:N/2);
f     = hr(1:N/2);

% ---------------------------------------------------------------------
% Plotting (TBC by James)
% ---------------------------------------------------------------------
if pf==1, 
    h = stem(f,p);
    set(h(1),'Marker','.');
    title('PSD Estimation');
    xlabel('Hz');
    ylabel('PSD');
 end
