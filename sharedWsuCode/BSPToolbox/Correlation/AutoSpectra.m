function [as,f] = AutoSpectra(x,fsa,wla,ola,fra,pfa)
%AutoSpectra: Estimates the autospectra of a signal.
%
%   [as,f] = AutoSpectra(x,fs,wl,ol,fr,pf)
%
%   x    Input Signal
%   fs   Sample Rate (Hz).  Default=1 Hz.
%   wl   Length of each of the windows (seconds).  Default=1001
%   ol   Percent window overlap.  Default=50.
%   fr   Minimum no. of frequencies to estimate. Default=200.
%   pf   Plot flag:  0=none (default), 1=screen.
%
%   as   Estimated Autospectra
%   f    Frequencies (Hz) at which autospectra is estimated 
%
%   AutoSpectra estimates the power spectral density of a stationary
%   input signal 'x'.  The power spectral density is calculated 
%   using Welch's method of the averaged periodogram defined as:
%
%         as = sqrt(Sxx)
%
%   where Sxx is the power spectral density of the signal x.  
%
%   This function removes the mean prior to spectral estimation.  If
%   only the window length is specified, the blackman window is used.
%   The overlap argument must be a number between 0 and 100.  If no 
%   output arguments are specified, the output will plot to the 
%   screen.
%
%   Example: Calculate the autospectra of an ABP signal with a 
%   a sample rate of 125 Hz, window length of 8 s, 50% overlap,
%   hanning window, and plot the results to the screen.
%
%      load ABPICP.mat
%      [as,f] = AutoSpectra(abp,125,hanning(8*125),[],[],1);
%
%   Challis, R. E., and Kitney, R. I., "Biomedical Signal Processing
%   (in four parts), Part 3, The Power Spectrum and Coherence 
%   Function", Med. & Biol. Eng. & Comput., 29, pp.225-241, 1991.
%
%   Version 1.00 LJ
%
%   See Also CrossSpectra, Coherence, Spectrogram, and 
%   CrossSpectrogram.

%====================================================================
% Error Check
%====================================================================
if nargin == 0 
    help AutoSpectra;
    return;
    end

%====================================================================
% Process Input Signal
%====================================================================
xr = x;
sx = std(x);
if sx==0,
    error('Signal ''1'' is constant.');
    end

mx = mean(x);
if abs(mx)/sx>1e-5,
    x = x - mx;
end

nx = length(x);
k = (1:nx)';
n = (0:nx-1)';
x = x(:);  %  Make sure x is a column vector

%====================================================================
% Process function arguments
%====================================================================
fs = 1; % Default: 1 Hz
if exist('fsa') & ~isempty(fsa),
     fs = abs(fsa);
 end

wl = min(1001,round(nx/2)); % Default window length
wn = blackman(wl); % Default window shape
if exist('wla') & ~isempty(wla),
    if length(wla)==1, % User specified window length only
        wl = round(wla*fs);
        if ~rem(wl,2),
            wl = wl + 1; % Make odd
        end
        wn = blackman(wl); % Default window shape
    else % User specified the actual window, not just the length
        wn = wla(:);
        wl = length(wn);
        if ~rem(wl,2), % Make odd
            wn = [wn;0];
            wl = wl + 1; 
        end
    end
end

ol = 0.5; % Default overlap
if exist('ola') & ~isempty(ola),
    ol = min(max(0,ola/100),1);
end

fr = 200; % Default no. frequencies to evaluate at
if exist('fra') & ~isempty(fra)
    fr = fra;
    end   
nz = 2^(ceil(log2(fr))+1);       %  Zero vector

if nargout == 0 & ~exist('pfa')
    pf = 1; % Plot to screen if no output arguments
else
    pf = 0;
    end

if exist('pfa') & ~isempty(pfa)
    pf = pfa;
    end

%====================================================================
% Variable Allocations & Initialization
%====================================================================
ss   = wl - ol*wl;                 % Step size (fractional samples)
ss   = max(1,ss);                  % Must be at least 1 sample 
ns   = floor((nx-(wl-ss))./ss);    % No. of steps to take
ns   = max(2,ns);                  % Take at least 2 steps
fi   = 1:(nz/2+1);                 % Index of frequencies to evaluate at
f    = fs*((fi-1)/nz)';            % Frequencies
nf   = length(f);                  % No. of frequencies

%====================================================================
% Estimate Autospectra
%====================================================================
S = zeros(ns,nf);

for cnt      = 1:ns,
    i0       = (cnt-1)*round(ss) + 1; % Initial index
    i1       = (i0 + (wl-1));         % Final index
    i1       = min(i1,nx);            % Keep from hitting the limit
    i0       = i1-(wl-1);             % Adjust initial, if necessary
    k        = (i0:i1);               % Signal window
    x1       = x(k).*wn;
    psd1     = ((abs(fft(x1)).^2)./length(x1))';
    S(cnt,:) = psd1(fi);
end

as = sqrt(sum(S)/ns)';

%====================================================================
% Plotting
%====================================================================
if pf==1
    figure;
    FigureSet;
    h = plot(f,as);
    set(h,'Marker','.');
    xlim([0 max(f)]);
    ylim([0 1.05*max(as)]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    title('Autospectra','FontWeight','bold');
    AxisSet;
    box off;
    end;

if nargout == 0,
     clear as;
     clear f;
end
