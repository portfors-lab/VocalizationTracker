function [c,ci,f] = ModulatedCoherence(x1,x2,fsa,wla,sna,aa,nfa,pfa)
%Coherencymod: Computes the estimated coherence of two signals vs frequency.
%
%   [c,f] = ModulatedCoherence(x1,x2,fs,wl,sn,nf,pf)
%
%   x1   Input Signal
%   x2   Input Signal
%   fs   Sample rate (Hz). Default = 1 Hz.
%   wl   Length of window to use (sec). Default = (duration of x)/10.
%        If a vector, specifies entire window.
%   sn   Signal to noise ratio used to bias the coherence. 
%        Default = inf.
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
    help ModulatedCoherence;
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
    
ol = 0;                              % Default: 0% overlap

nf = 200;                            % Default no. frequencies to evaluate at
if exist('nfa') & ~isempty(nfa)
    nf = nfa;
    end   
nz = 2^(ceil(log2(nf))+1);           % Zero vector
    
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
nd   = 2*length(fi)-1;           % No. of diagonals
di   = -(nd-1)/2:(nd-1)/2;       % Index of diagonals
ci   = zeros(nf,nf,2);           % Confidence intervals
nw   = floor(nx/wl);             % No. of independent windows
%====================================================================
% Estimate Coherence
%====================================================================
S12 = zeros(nf,nf);
S11 = zeros(nf,nf);
S22 = zeros(nf,nf);

for cnt = 1:ns,
    dx1 = zeros(nf,nf);
    
    i0  = (cnt-1)*round(ss) + 1; % Initial index
    i1  = (i0 + (wl-1));         % Final index
    i1  = min(i1,nx);            % Keep from hitting the limit
    i0  = i1-(wl-1);             % Adjust initial, if necessary
    k   = (i0:i1);

    X1  = fft(x1(k).*wn,nz);
    X2  = fft(x2(k).*wn,nz);
    X2f  = repmat(X2(fi).',[1 nf]);
    X1f  = repmat(X1(fi),[nf 1]);
    
    S12(:,:) = S12 + X2f.*conj(X1f).*(exp(-j*(fi-1)*i0*2*pi/nz)*exp(-j*(fi-1)*i0*2*pi/nz).');
    S11(:,:) = S11 + abs(X1f.*repmat(exp(-j*(fi-1)*(i0)*2*pi/nz),[nf 1])).^2;
    S22(:,:) = S22 + abs(X2f.*repmat(exp(-j*(fi-1)*(i0)*2*pi/nz).',[1 nf])).^2;
    end;

num   = squeeze(abs(S12).^2);
d1    = squeeze(abs(S11));
d2    = squeeze(abs(S22));
d1    = d1+(1/sn)*repmat(sum(d1,1),[nf 1])/var(x1); % Add equivalent signal noise
d2    = d2+(1/sn)*repmat(sum(d2,2),[1 nf])/var(x2); % Add equivalent signal noise
den   = d1.*d2;


id    = find(den~=0);
c     = zeros(nf);
c(id) = num(id)./den(id);

% Approximation of unbiased estimate
%c = c - (1/nw)*(1-c).^2.*(1+2*c/nw);
%c = min(1,max(0,c)); %0 <= c <= 1


% Approximation of confidence intervals
%ci = ones(size(c,1),size(c,2),2);
% z = acosh(c);
% ci(:,:,1) = norminv(0.025,z,sqrt(1/(2*nw-2)));
% ci(:,:,2) = norminv(0.975,z,sqrt(1/(2*nw-2)));
% ci = cosh(ci);
%e = zeros(size(c));
%k = find(c ~= 0);
%e(k) = (1-c(k)).*sqrt(2./(nw*c(k)));
%ci(:,:,1) = (1-2*e).*c;
%ci(:,:,2) = (1+2*e).*c;

% coherence must be 0 <= c <= 1
%ci = min(1,max(0,ci)); 
%figure;plot(c(:),reshape(ci(:,:,1),[],1),'.r',c(:),reshape(ci(:,:,2),[],1),'.b');
%figure;imagesc(abs(ci(:,:,1)-ci(:,:,2)),[0 1]);

%a = abs(ci(:,:,2) - ci(:,:,1));
%a = 1-a; %large ranges become transparent

%====================================================================
% Plotting
%====================================================================
if pf==1,
    figure;
    FigureSet;
    %image([f],[f],ones(nf,nf,3));
    %hold on;
    i = imagesc([f],[f],c,[0 1]);
    %ciImage(c,ci,f,f);
    set(gca, 'ydir', 'normal');
    colormap(colorspiral(256,1));
    %set(i, 'alphadata', a);
    ylabel('x2 Frequency (Hz)');
    xlabel('x1 Frequency (Hz)');
    title('Coherency Spectrum','FontWeight','bold');
    box off;
    AxisSet;
    %figure;
    %imagesc(a(:,:,1));
    colorbar;
    %ciBar;
    
    
    end

%====================================================================
% Process Return Arguments
%====================================================================
if nargout == 0,
    clear c;
    clear f;
end
end

function ih = ciImage(data,ci,x,y,pa)
%
%   Creates an image representing the data and confidence intervals
%
%   data = 2 dimensional data
%   ci   = confidence intervals for each data point
%   sd   = Flag to determine whether data should be scaled. If this is 
%          not set data values must be between 0 and 1(default=1)
%   pa   = Plot axis handle

%parameters
th = 0.1;    %lower intensity threshold - scale everything to be brighter than this

%error checking
if(ndims(data) ~= 2)
    error('Requires 2 dimentional data')
    return
end
if(ndims(ci) ~= 3 || size(ci,3) ~= 2)
    error('Requires upper and lower confidence bounds for each data point');
    return
end
if(max(data(:)) > 1 || min(data(:)) < 0)
    k = [find(data > 1) find(data < 0)];
    length(k)
    k(1:min(end,5))
    data(k)
    error('Data must be between 0 and 1');
    return
end
if(any(ci(:) > 1))
    error('Max confidence interval is 1');
    return
end
if(any(ci(:) < 0))
    error('Min confidence interval is 0');
    return
end


ra = squeeze(abs(ci(:,:,1)-ci(:,:,2)));  %range of confidence interval
mra = zeros(size(ra)); %maximum range allowed by the color space
%scale data to be above specified threshold
data = th+(data*(1-th));

%set the luminance
r = data;
g = data;
b = data;

%narrow range for small data values because there is not enough room for 
%the whole range
k = find(data <= 1/3);
ra(k) = 4*data(k).*ra(k);
mra(k) = 4*data(k);

k = find(data > 1/3);
ra(k) = 2*(1-data(k)).*ra(k);
mra(k) = 2*(1-data(k));

%shift to red for high ci range 
k = find(ra > 0.5*mra);
r(k) = r(k) + (ra(k)-0.5*mra(k));
g(k) = g(k) - 0.5*(ra(k)-0.5*mra(k));
b(k) = b(k) - 0.5*(ra(k)-0.5*mra(k));

% %debugging code
% if(any(b(:) < 0))
%     k = find(b < 0);
%     disp(['r:Subtracted too much from blue']);
%     b(k(1:min(length(k),5)))
% end

%shift to green for low ci range
k = find(ra < 0.5*mra);
g(k) = g(k) + (0.5*mra(k)-ra(k));
r(k) = r(k) - 0.5*(0.5*mra(k)-ra(k));
b(k) = b(k) - 0.5*(0.5*mra(k)-ra(k));

% %debugging code
% if(any(b(:) < 0))
%     k = find(b < 0);
%     disp(['g:Subtracted too much from blue']);
%     b(k(1:min(length(k),5)))
%     k(1:min(length(k),5))
% end


im = zeros([size(data,1) size(data,2) 3]);
im(:,:,1) = r;
im(:,:,2) = g;
im(:,:,3) = b;

%figure;
if exist('x') & exist('y') & ~isempty(x) & ~isempty(y)
    if exist('pa') & ~isempty(pa)
        set(pa,'cdata', im);
        ih = pa;
    else
        ih = image([x],[y],im);
    end
else
    if exist('pa') & ~isempty(pa)
        set(pa,'cdata', im);
        ih = pa;
    else
        ih = image(im);
    end
end
end

function h = ciBar

h = colorbar;
hc = get(h, 'children');

%make the ci/data map
x = linspace(0,1,64);
y = linspace(0,1,20);
[Y,X] = meshgrid(y,x);
Y(:,:,2) = zeros(size(Y));

ih = ciImage(X,Y,[],[],hc);
set(ih,'ydata',[0 1]);
set(get(ih,'parent'),'ylim',[0 1]);
end