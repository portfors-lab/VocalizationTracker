function [p,f] = BlackmanTukey(x,fsa,wla,nfa,pfa)
%BlackmanTukey: Estimate PSD using the Blackman-Tukey method
%
%   [p,f] = BlackmanTukey(x,fs,wl,nf,pf);        
%
%   x    Input signal.
%   fs   Sample rate (Hz). Default = 1 Hz.
%   wl   Length of window to use (sec). Default = (max. possible).
%        If a vector, specifies entire window.
%   nf   Number of frequencies to evaluate.
%        Default = max(128,round(wl/2)).
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   p    Power spectral density. 
%   f    Frequencies at which p is estimated (Hz).
%
%   Estimates the power spectral density (PSD) of an input signal 
%   using the Blackman-Tukey method. This method estimates the
%   autocorrelation sequence, multiplies the autocorrelation sequence
%   by a window, and then calculates the Fourier transform of the
%   windowed autocorrelation. The window multiplication in the time
%   domain is equivalent to convolution (a form of smoothing) in the
%   frequency domain and trades variance of the PSD estimate for
%   increased bias. 
%
%   This implementation uses an unbiased estimate of the 
%   autocorrelation (calculated via the FFT):
%              1   N-k
%      r(k) = ---  sum x(m) x(m+k)
%             N-k  m=0
%   Thus, the expected value of the estimated PSD is related to the 
%   true PSD by a single convolution operation with one window. The 
%   mean is removed from the signal prior to estimation to prevent an 
%   impulse at 0 Hz from dominating the smoothed estimate.
%
%   The estimated PSD is scaled such that Parseval's relation is 
%   approximately satisfied: var(x) = df*sum(p)*2, where df is the
%   frequency step size and the factor of 2 uses symmetry to 
%   estimate a two-sided integral of the PSD from -fs/2 to fs/2. 
%
%   If only the window length is specified, the blackman window is 
%   used. The specified window length should be odd. If the specified
%   window length is even, 1 is added to make it odd. If the window 
%   itself is specified with an even number of elements, a zero is 
%   appended to make the window odd.
%
%   Example: Estimate the PSD of an electrocardiogram signal and 
%   plot the results. Include at least 5000 points in the PSD 
%   estimate and use a window length of 10 seconds.
%
%      load NoisyECG.mat;
%      x  = decimate(ecg(1:50e3),5);
%      fs = fs/5;
%      BlackmanTukey(x,fs,10,5000);
%
%   M. Hayes, Statistical Digital Signal Processing and Modeling. 
%   New York: John Wiley & Sons, 1996, pp.420-423.
%
%   J. G. Proakis, C. M. Rader, F. Ling, C. L. Nikias, M. Moonen, 
%   and I. K. Proudler, Algorithms for Statistical Signal Processing.
%   Saddle River, NJ: Prentice Hall, 2002, pp. 449-452.
%
%   Version 1.01.01.23 JM
%
%   See also SPECTRUM and AutoCorrelate.

%====================================================================
% Error Checking
%====================================================================    
if ~IsScript(mfilename) & nargin<1,
    help BlackmanTukey;
    return;
    end;

nx = length(x);
if nx==0,
    error('Signal is empty.\n');
    end;
    
%====================================================================
% Process function arguments
%====================================================================
fs = 1; % Default sampling rate, Hz
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
    end;
    
m  = 32;
if exist('ma') & ~isempty(ma),
    m = ma;
    end;

wl = nx*2-1;                 % Default window length
if ~rem(wl,2),
    wl = wl + 1;             % Make odd
    end;
wn = blackman(wl)';          % Default window shape
if exist('wla') & ~isempty(wla),
    if length(wla)==1,       % User specified window length only
        wl = round(wla*fs);  % Convert from seconds to samples
        wl = max(wl,3);      % Must be at least 3 points
        wl = min(wl,nx*2-1); % Must not be more than nx*2-1
        if ~rem(wl,2),
            wl = wl + 1;     % Make odd
            end;
        wn = blackman(wl)';  % Default window shape
    else % User specified the actual window, not just the length
        wn = wla(:)';
        wl = length(wn);
        if ~rem(wl,2),
            wn = [wn 0]; % Append a zero to make odd
            wl = wl + 1;
            end;
        end;
    end;

nf = max(128,round(wl/2)); % Default No. of frequencies 
if exist('nfa') & ~isempty(nfa),
    nf = nfa;
    end;   
        
pf = 0;                                 % Default - no plotting
if IsScript(mfilename) | nargout==0,    % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;

%====================================================================
% Preprocessing
%====================================================================    
wn = wn(:);      % Make into a column vector
x  = x(:);       % Make into a column vector
mx = mean(x);    % Estimate the signal mean
vx = var(x);     % Estimate the signal variance
x  = x - mx;     % Remove mean

%====================================================================
% Calculate Autocorrelation Quickly using FFT
%====================================================================
nz = 2*nx-1;               % Minimum no. of points to use FFT    
nz = 2^(ceil(log2(nz)));   % Convert to power of 2 for FFT
xf = fft(x,nz);            % Calculate FFT of x
r  = ifft(xf.*conj(xf));   % Calculate biased autocorrelation
r  = real(r(1:nx));        % Eliminate superfluous zeros and nearly zero imaginary part
k  = (0:nx-1).';           % Generate scaling indices 
a  = 1./(nx-k);            % Calculate scaling coefficient
r  = r.*a;                 % Scale 
no = (wl+1)/2;             % Number of one-sided points to included in estimate 
r  = [r(no:-1:2);r(1:no)]; % Make two-sided

if 0, % Error checking routine to ensure autocorrelation was calculated correctly
    rb = zeros(nx,1); % Brute force technique
    for c1 = 1:nx,
        rb(c1) = sum(x(1:nx-(c1-1)).*x(c1:nx))/(nx-(c1-1));
        end;
    rb = [rb(no:-1:2);rb(1:no)]; % Make two-sided
    fprintf('Estimation Error: %f\n',max(abs(diff(rb-r))));
    end;
    
%====================================================================
% Variable Allocations, Indexing, Windowing, & Estimation
%====================================================================
nz   = 2*nf;                  % No. window points needed in FFT
nz   = 2^(ceil(log2(nz)));    % Convert to power of 2 for FFT
fi   = 1:(nz/2+1);            % Indices of frequencies to return
nf   = length(fi);            % No. of frequencies that PSD is evaluated at 
f    = (fi-1)*fs/nz;          % Frequencies in units of Hz
fr   = fs/nz;                 % Frequency resolution (frequency difference between adjacent estimates)
p    = abs(fft(r.*wn,nz)).^2; % Window autocorrelation and calculate the estimate
p    = p(fi);

%====================================================================
% Postprocessing
%====================================================================  
f = f(:);               % Ensure is column vector
p = p(:);               % Ensure is column vector
p = p*vx/(2*fr*sum(p)); % Scale to satisfy Parseval's relation

%====================================================================
% Plot PSD
%====================================================================
if pf==1,
    figure;
    FigureSet;
    plot(f,p);
    title('Blackman-Tukey Estimated PSD');
    xlabel('Frequency (Hz)');
    ylabel('PSD (scaled)');
    xlim([0 fs/2]);
    ylim([0 1.02*max(p)]);
    box off;
    zoom on;
    AxisSet;
    end

%====================================================================
% Process Return Arguments
%====================================================================
if ~IsScript(mfilename) & nargout==0,
    clear('p','f');
    end;

    