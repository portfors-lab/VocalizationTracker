function [p,f] = Welch(x,fsa,wla,ola,nfa,pfa)
%Welch: Estimate PSD using the Welch's method.
%
%   [p,f] = Welch(x,fs,wl,ol,nf,pf)
%
%   x    Input signal.
%   fs   Sample rate (Hz). Default = 1 Hz.
%   wl   Length of window to use (sec). Default = (max. possible).
%        If a vector, specifies actual window.
%   ol   Percent window overlap.  Default = 50%.
%   nf   Number of frequencies to evaluate.
%        Default = max(128,round(wl/2)).
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   p    Power spectral density. 
%   f    Frequencies at which p is estimated (Hz).
%
%   Estimates the power spectral density (PSD) of an input signal 
%   using Welch's method. This method estimates the PSD by 
%   calculating the FFT of overlapping windowed segments of the 
%   signal. The window multiplication in the time domain is 
%   equivalent to convolution (a form of smoothing) in the frequency 
%   domain and trades variance of the PSD estimate for increased 
%   bias. 
%
%   Since the PSD is an estimate of power, multiplication by the 
%   window in the time domain is equivalent to convolution by the
%   square of the window in the frequency domain. Thus, it may be
%   more intuitive to multiply the signal by the square root of
%   the popular windows.
%
%   The mean is removed from the signal prior to estimation to 
%   prevent an impulse at 0 Hz from dominating the smoothed estimate.
%   If only the window length is specified, the square root of the 
%   Blackman window is used. 
%
%   The estimated PSD is scaled such that the estimate is 
%   asymptotically unbiased and Parseval's relation is approximately 
%   satisfied:
%                          +pi
%      var(x) ~= inv(2*pi) int p(w) dw ~= sum(p)/length(p).
%                          -pi
%
%   Example: Estimate the PSD of an electrocardiogram signal and 
%   plot the results. Include at least 5000 points in the PSD 
%   estimate and use a window length of 10 seconds.
%
%      load NoisyECG.mat;
%      x  = decimate(ecg(1:50e3),5);
%      fs = fs/5;
%      Welch(x,fs,10,[],5000);
%
%   M. Hayes, Statistical Digital Signal Processing and Modeling. 
%   New York: John Wiley & Sons, 1996, pp. 415-420.
%
%   J. G. Proakis, C. M. Rader, F. Ling, C. L. Nikias, M. Moonen, 
%   and I. K. Proudler, Algorithms for Statistical Signal Processing.
%   Saddle River, NJ: Prentice Hall, 2002, pp. 447-449.
%
%   Version 1.00 JM
%
%   See also SPECTRUM, WINDOW, and SpectralAnalysis.


%====================================================================
% Error Checking
%====================================================================    
if nargin<1,
    help Welch;
    return;
    end;

nx = length(x);
if nx==0,
    error('Signal is empty.');
    end;
    
%====================================================================
% Process function arguments
%====================================================================
fs = 1;                              % Default sampling rate
if exist('fsa') & ~isempty(fsa)
    fs = fsa;
    end

wl = round(nx/10);                   % Default window length (samples) - guaranteed to be odd
wn = blackman(wl+2);                 % Default window shape
wn = wn(2:wl+1);                     % Truncate leading and ending zeros
if exist('wla') & ~isempty(wla),
    if length(wla)==1,               % User specified window length only
        wl = round(wla*fs);          % Convert seconds to samples
        wn = sqrt(blackman(wl));     % Default window shape
    else                             % User specified the actual window, not just the length
        wn = wla(:);
        wl = length(wn);
        end;
    if wl>nx,
        error('User-specified window length exceeds signal duration.');
        end;
    end;
        
ol = 0.5;                            % Default: 50% overlap
if exist('ola') & ~isempty(ola),
    ol = min(max(0,ola/100),1);
    end;

nf = max(128,round(wl/2));           % Default No. of frequencies 
if exist('nfa') & ~isempty(nfa)
    nf = nfa;
    end   
    
pf = 0;                              % Default - no plotting
if nargout==0, % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;

%====================================================================
% Preprocessing and Memory Allocation
%====================================================================
x   = real(x (:).');             % Convert to a row vector
wn  = real(wn(:).');             % Convert to a row vector
wn  = wn/sqrt(mean(wn.^2));      % Normalize window to make asymptotically unbiased 
mx  = mean(x);                   % Signal mean
vx  = var(x);                    % Signal variance
x   = x  - mean(x);              % Remove the signal mean

%====================================================================
% Initialize Variables
%====================================================================
ss   = wl - ol*wl;               % Step size (fractional samples)
ss   = max(1,ss);                % Must be at least 1 sample 
ns   = floor((nx-(wl-ss))./ss);  % No. of steps to take
ns   = max(1,ns);                % Take at least 2 steps
nz   = 2^(ceil(log2(nf))+1);     %  Zero vector
fi   = 1:(nz/2+1);               % Index of frequencies to evaluate at
f    = fs*((fi-1)/nz)';          % Frequencies
nf   = length(f);                % No. of frequencies
fr   = fs/nz;                    % Frequency resolution (frequency difference between adjacent estimates)

%====================================================================
% Estimate Coherence
%====================================================================
p = zeros(1,nf);
for cnt = 1:ns,
    i0  = (cnt-1)*round(ss) + 1; % Initial index
    i1  = (i0 + (wl-1));         % Final index
    i1  = min(i1,nx);            % Keep from hitting the limit
    i0  = i1-(wl-1);             % Adjust initial, if necessary    
    k   = (i0:i1);               % Specify signal indices   
    ps  = fft(x(k).*wn,nz);      % Calculate PSD estimate from segment
    p   = p + abs(ps(fi)).^2;    % Add to running sum
    end;    
    
p = p/ns;                        % Generate average 
p = p/wl;                        % Scale
    
%====================================================================
% Postprocessing
%====================================================================  
f = f(:);               % Ensure is column vector
p = p(:);               % Ensure is column vector

%====================================================================
% Plot PSD
%====================================================================
if pf==1,
    figure;
    FigureSet;
    plot(f,p);
    title(sprintf('Welch%cs Method Estimated PSD',char(39)));
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
if nargout==0,
    clear('p','f');
    end;

    

