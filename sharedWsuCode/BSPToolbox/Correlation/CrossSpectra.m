function [cs,f] = CrossSpectra(x1,x2,fsa,wla,ola,fra,pfa)
%CrossSpectra: Computes the cross-spectra of two signals.
%
%   [cs,f] = CrossSpectra(x1,x2,fs,wl,ol,fr,pf)
%
%   x1   Input Signal
%   x2   Input Signal
%   fs   Sample Rate, hertz.
%   wl   Length of each of the windows (seconds).  Default=1001
%   ol   Percent window overlap.  Default=50.
%   fr   Minimum no. of frequencies to estimate. Default=200.
%   pf   Plot flag:  0=none (default), 1=screen.
%
%   cs   Estimated Crossspectra
%   f    Frequencies
%
%   CrossSpectra estimates the cross power spectral density of 
%   of two equal length input signals, x1 and x2.  The cross power 
%   spectral density is calculated using Welch's method of the 
%   averaged periodogram defined as:
%
%         cs = abs(X1 .* conj(X2))
%
%   where X1 and X2 are zero padded Fourier transforms of x1 and x2.  
%
%   This function removes the mean of each signal prior to spectral 
%   estimation.  If only the window length is specified, the 
%   blackman window is used. The overlap argument must be a number
%   between 0 and 100.  If no output arguments are specified, the 
%   output will plot to the screen.
%
%   Example: Calculate the cross-spectra of an ABP data segment and 
%   an ICP data segment with a sample rate of 125 Hz, window length
%   of 8 s, 50% overlap, hamming window, and plot the results to the 
%   screen.
%
%      load ABPICP.mat
%      x1 = abp(1:4000);
%      x2 = icp(1:4000);
%      CrossSpectra(x1,x2,125,hamming(8*125));
%
%   Challis, R. E., and Kitney, R. I., "Biomedical Signal Processing
%   (in four parts), Part 3, The Power Spectrum and Coherence 
%   Function", Med. & Biol. Eng. & Comput., 1991, 29, pp. 225-241.
%
%   Version 1.00 LJ
%
%   See Also AutoSpectra and Coherence.

%====================================================================
% Error Check
%====================================================================
if nargin == 0 
    help CrossSpectra;
    return;
end

%====================================================================
% Process Input Signals
%====================================================================
nx1 = length(x1);
xr1 = x1;
sx1 = std(x1);
if sx1==0,
    fprintf('ERROR: Signal ''1'' is constant.\n');
    return;
end

nx2 = length(x2);
xr2 = x2;
sx2 = std(x2);
if sx2==0,
    fprintf('ERROR: Signal ''2'' is constant.\n');
    return;
end

if nx1~=nx2,
    fprintf('ERROR: Signals must be of the same length.\n');
    return;
end

mx1 = mean(x1);
if abs(mx1)/sx1>1e-5,
    x1 = x1 - mx1;  % Remove mean
end

mx2 = mean(x2);
if abs(mx2)/sx2>1e-5,
    x2 = x2 - mx2;  % Remove mean
end

nx = nx1;
k  = (1:nx)';
n  = (0:nx-1)';
x1 = x1(:);  %  Make sure x1 is a column vector
x2 = x2(:);  %  Make sure x2 is a column vector

%====================================================================
% Process function arguments
%====================================================================
fs = 1; % Default sampling rate
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

if nargout == 0 & ~exist('pf')
    pf = 1; % Plot to screen if no output arguments
end

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
% Estimate Cross-spectra
%====================================================================
S1 = zeros(ns,nf);
S2 = zeros(ns,nf);

for cnt       = 1:ns,
    i0        = (cnt-1)*round(ss) + 1; % Initial index
    i1        = (i0 + (wl-1));         % Final index
    i1        = min(i1,nx);            % Keep from hitting the limit
    i0        = i1-(wl-1);             % Adjust initial, if necessary
    k         = (i0:i1);               % Signal window
    seg1      = x1(k).*wn;             % Windowing factor
    seg2      = x2(k).*wn;             % Windowing factor
    psd1      = fft(x1,nz);
    psd2      = fft(x2,nz);
    S1(cnt,:) = psd1(fi)';
    S2(cnt,:) = psd2(fi)';
    end

cs1 = (sum(S1)/ns)';
cs2 = (sum(S2)/ns)';
cs  = (abs(cs1.*conj(cs2)));

%====================================================================
% Plotting
%====================================================================
if pf==1,
    figure
    FigureSet;
    h = plot(f,cs);
    set(h,'Marker','.');
    box off;
    xlim([0 max(f)]);
    ylim([0 1.05*max(cs)]);
    xlabel('Frequency, Hz');
    ylabel('Magnitude');
    title('Cross-spectra','FontWeight','bold');
    AxisSet;
    zoom on;
    end

if nargout == 0,
    clear cs;
    clear f;
end
