function [tf,f] = TransferFunction(x,y,fsa,wla,ola,fra,pfa)
%TransferFunction:  Estimates transfer function of two signals.
%
%   [tf,f] = TransferFunction(x,y,fs,wl,ol,fr,pf)
%
%   x    Input Signal
%   y    Output Signal
%   fs   Sample rate (Hz). Default = 1 Hz.
%   wl   Length of window to use (sec). Default = 1001 samples.
%        If a vector, specifies entire window.
%   ol   Percent window overlap. Default=50.
%   fr   Minimum no. of frequencies to estimate. Default=275.
%   pf   Plot format: 0=none (default), 1=screen.
%
%   tf   Estimated transfer function.
%   f    Frequencies at which tf is estimated.
%   
%   TransferFunction estimates the transfer function of a system
%   where x is the input signal and y is the output signal.  This
%   function uses Welch's method of the averaged periodogram and is
%   defined as
%         
%         tf = Sxy./abs(Sxx)
%
%   where Sxy is the cross-power spectral density and Sxx is the 
%   power spectral density of the signal x.
%   
%   This function removes the mean and detrends the signal windows
%   prior to spectral estimation.  If only the window length is
%   specified, the blackman window is used.  The overlap argument 
%   must be a number between 0 and 100.  If no output arguments are
%   specified, the output will plot the magnitude transfer function
%   and the phase to the screen.
%   
%   Example: Calculate the transfer function using ICP as the input
%   signal and ABP as the output signal, decimated to 25 Hz, using
%   a hanning window with a window length of 30 s.  Plot the
%   results to the screen.
%
%      load ABPICP.mat
%      x = decimate(icp,5);
%      y = decimate(abp,5);
%      [tf,f] = TransferFunction(x,y,25,hanning(30*25),[],[],1);
%
%   Bendat, J., and Piersol, A., "Engineering Applications of
%   Correlation and Spectral Analysis", John Wiley & Sons, pp.97-114,
%   1980.
%
%   Version 1.00 LJ
%
%   See Also AutoSpectra and CrossSpectra.

%====================================================================
% Error Check
%====================================================================
if nargin < 2 
    help TransferFunction;
    return;
end

%====================================================================
% Process Input Signals
%====================================================================
x = x(:);
y = y(:);

nx = length(x);
xr = x;
sx = std(x);
if sx==0,
    fprintf('ERROR: Signal ''x'' is constant.\n');
    return;
end

ny = length(y);
yr = y;
sy = std(y);
if sy==0,
    fprintf('ERROR: Signal ''y'' is constant.\n');
    return;
end

if nx~=ny,
    fprintf('ERROR: Signals must be of the same length.\n');
    return;
end

mx = mean(x);
if abs(mx)/sx>1e-5,
    x = x - mx;  % Remove mean
end

my = mean(y);
if abs(my)/sy>1e-5,
    y = y - my;  % Remove mean
end

k  = (1:nx)';
n  = (0:nx-1)';
x  = x(:);  %  Make sure x1 is a column vector
y  = y(:);  %  Make sure x2 is a column vector

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

fr = 275; % Default no. frequencies to evaluate at
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
% Estimate Transfer Function
%====================================================================
for cnt = 1:ns,
    i0  = (cnt-1)*round(ss) + 1; % Initial index
    i1  = (i0 + (wl-1));         % Final index
    i1  = min(i1,nx);            % Keep from hitting the limit
    i0  = i1-(wl-1);             % Adjust initial, if necessary
    k   = (i0:i1);               % Signal window
    X   = fft(detrend(x(k)).*wn,nz);
    Y   = fft(detrend(y(k)).*wn,nz);
 
    Sxx(cnt,:) = (abs((X(fi)).^2))';
    Sxy(cnt,:) = (conj(X(fi)).*Y(fi))';
end

den = (sum(Sxx))';
num = (sum(Sxy))';
tf  = num./den;
H   = (abs(tf));  % Magnitude transfer function
ph  = angle(tf).*180/pi; % Adjust from radians to degrees

%====================================================================
% Plotting
%====================================================================
if pf==1,
    figure;
    FigureSet;
    subplot(2,1,1);
        h = plot(f,H,'r');
        set(h,'LineWidth',1.2);
        grid on;
        box off;
        xlim([0 max(f)]);
        ylim([0 1.05*max(H)]);
        ylabel('Magnitude (scaled)');
        title('Transfer Function','FontWeight','bold');
    subplot(2,1,2);
        h = plot(f,ph,'r');
        set(h,'LineWidth',1.2);
        grid on;
        box off;
        xlim([0 max(f)]);
        ylim([-185 185]);
        xlabel('Frequency (Hz)');
        ylabel('Phase (degrees)');
    AxisSet;
    end;

if nargout == 0,
    clear tf;
    clear f;
end