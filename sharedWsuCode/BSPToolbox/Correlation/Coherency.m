function [c,f] = Coherency(x1,x2,fsa,wla,sna,ola,nfa,pfa)
%Coherency: Computes the estimated coherence of two signals.
%
%   [c,f] = Coherency(x1,x2,fs,wl,sn,ol,nf,pf)
%
%   x1   Input Signal
%   x2   Input Signal
%   fs   Sample rate (Hz). Default = 1 Hz.
%   wl   Length of window to use (sec). Default = (duration of x)/10.
%        If a vector, specifies entire window.
%   sn   Signal to noise ratio used to bias the coherence. 
%        Default = inf.
%   ol   Percent window overlap.  Default=50.
%   nf   Minimum no. of frequencies to estimate. Default = 200.
%   pf   Plot format: 0=none (default), 1=screen. 
%
%   c    Estimated coherency spectrum
%   f    Frequencies at which coh is estimated
%
%   The coherency spectrum is a function that gives an output between 
%   0 and 1 and is a measure of the correlation of frequency 
%   components of random signals x1 and x2. The coherency spectrum
%   is defined as
%
%        c = abs(Sxy)./sqrt(Sxx.*Syy)).
%
%   where Sxy is the cross-power spectral density and Sxx is the 
%   power spectral density of the signal x. If x and y are zero-mean 
%   white noise processes with a coefficient of correlation that is 
%   equal to p,
%
%        p = E[x*y]/sqrt(Var[x]*Var[y]),
%
%   the coherency is related to p by c = p^2.
%
%   This function removes the mean prior to spectral estimation. The
%   spectral components are estimated using Welch's method. If only 
%   the window length is specified, the Blackman window is used. The 
%   overlap argument must be a number between 0 and 100. The 
%   magnitude squared length FFT is calculated for each window and 
%   averaged to estimate Sxy, Sxx, and Syy. 
%   
%   Example: Calculate the coherence of an ABP data segment and an ICP
%   data segment, each with a sample rate of 125 Hz and 10 s data 
%   windows, using Hanning windows.  Plot the results to the screen.
%
%      load ABPICP.mat
%      x1  = abp(1:10000);
%      x2  = icp(1:10000);
%      [c,f] = Coherency(x1,x2,125,hanning(10*125),[],[],[],1);
%
%   Challis, R. E., and Kitney, R. I., "Biomedical Signal Processing
%   (in four parts), Part 3, The Power Spectrum and Coherence 
%   Function", Med. & Biol. Eng. & Comput., 29, pp.225-241, 1991.
%
%   Priestley, "Spectral Analysis and Time Series," Academic Press,
%   pp.556-557, 1981.
%
%   Version 1.00 JM
%
%   See Also AutoSpectra, CrossSpectra, and Cohereogram.

%====================================================================
% Error Check
%====================================================================
if nargin < 2 
    help Coherency;
    return;
    end

%====================================================================
% Error Checking
%====================================================================
nx1 = length(x1);
nx2 = length(x2);

if nx1~=nx2,
    error('Signals must be of the same length.');
    end
nx = nx1;

sx1 = std(x1);
if sx1==0,
    error('Signal ''1'' is constant.');
    end

sx2 = std(x2);
if sx2==0,
    error('Signal ''2'' is constant.');
    end
    
%====================================================================
% Process function arguments
%====================================================================
fs = 1;                              % Default sampling rate
if exist('fsa') & ~isempty(fsa)
    fs = fsa;
    end

wl = round(nx/10);                   % Default window length (samples) - guaranteed to be odd
wn = blackman(wl);                   % Default window shape
if exist('wla') & ~isempty(wla),
    if length(wla)==1,               % User specified window length only
        wl = round(wla*fs);
        wn = blackman(wl);           % Default window shape
    else                             % User specified the actual window, not just the length
        wn = wla(:);
        wl = length(wn);
        end;
    end;
    
sn = inf;                   % Default: infinite signal to noise ration (SNR)
if exist('sna') & ~isempty(sna),
    sn = sna;                         % User-specified (must be a positive number)
    end;    
    
ol = 0.5;                            % Default: 50% overlap
if exist('ola') & ~isempty(ola),
    ol = min(max(0,ola/100),1);
    end;

nf = 200;                            % Default no. frequencies to evaluate at
if exist('nfa') & ~isempty(nfa)
    nf = nfa;
    end   
nz = 2^(ceil(log2(nf))+1);           %  Zero vector
    
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
x1  = x1(:).'; % Convert to a row vector
x2  = x2(:).'; % Convert to a row vector    
wn  = wn(:).'; % Convert to a row vector

wn  = wn*sqrt(wl*inv(nz*sum(wn.^2))); % Ensure window preserves energy 

m1  = mean(x1);
m2  = mean(x2);

x1  = x1 - mean(x1);
x2  = x2 - mean(x2);

%====================================================================
% Initialize Variables
%====================================================================
ss   = wl - ol*wl;               % Step size (fractional samples)
ss   = max(1,ss);                % Must be at least 1 sample 
ns   = floor((nx-(wl-ss))./ss);  % No. of steps to take
ns   = max(2,ns);                % Take at least 2 steps
fi   = 1:(nz/2+1);               % Index of frequencies to evaluate at
f    = fs*((fi-1)/nz)';          % Frequencies
nf   = length(f);                % No. of frequencies

%====================================================================
% Estimate Coherence
%====================================================================
S12 = zeros(ns,nf);
S11 = zeros(ns,nf);
S22 = zeros(ns,nf);

for cnt = 1:ns,
    i0  = (cnt-1)*round(ss) + 1; % Initial index
    i1  = (i0 + (wl-1));         % Final index
    i1  = min(i1,nx);            % Keep from hitting the limit
    i0  = i1-(wl-1);             % Adjust initial, if necessary    
    k   = (i0:i1);
    
    X1  = fft(x1(k).*wn,nz);
    X2  = fft(x2(k).*wn,nz);
            
    S12(cnt,:) = X1(fi).*conj(X2(fi));
    S11(cnt,:) = abs(X1(fi)).^2;
    S22(cnt,:) = abs(X2(fi)).^2;    
    end;    
    
num   = abs(sum(S12)).^2;
d1    = sum(S11);
d2    = sum(S22);
d1    = d1+(1/sn)*sum(d1)/var(x1); % Add equivalent signal noise
d2    = d2+(1/sn)*sum(d2)/var(x2); % Add equivalent signal noise
den   = d1.*d2;
%den   = sum(S11).*sum(S22);



id    = find(den~=0);
c     = zeros(nf,1);
c(id) = num(id)./den(id);
c     = sqrt(c);                 % Convert coherence to coherency
c     = c(:);                    % Convert to column vector

%====================================================================
% Plotting
%====================================================================
if pf==1,
    figure;
    FigureSet;
    plot(f,c);
    ylabel('Coherency');
    xlabel('Frequency, Hz');
    title('Coherency Spectrum','FontWeight','bold');
    box off;
    AxisSet;
    ylim([0 1]);
    xlim([0 max(f)]);
    end

%====================================================================
% Process Return Arguments
%====================================================================
if nargout == 0,
    clear c;
    clear f;
    end

