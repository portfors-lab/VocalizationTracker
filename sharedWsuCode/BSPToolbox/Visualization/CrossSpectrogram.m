function [CS,t,f] = CrossSpectrogram(x1,x2,fsa,wla,fra,nfa,nsa,pfa);
% CrossSpectrogram:  Estimate and plot the cross-spectrogram of two
% signals.
%
%   [CS,t,f] = CrossSpectrogram(x1,x2,fs,wl,fr,nf,ns,pf);
%
%   x1   Input signal.
%   x2   Input signal.
%   fs   Sample rate, hertz.  Default = 1 Hz.
%   wl   Length of window to use (sec). Default = 1000 samples.
%        If a vector, specifies entire window.
%   fr   Array of minimum and maximum frequencies.  
%        Default = [0 fs/2].
%   nf   Number of frequencies to evaluate, default = wl/2.
%   ns   Requested number of times (horizontal pixels) to evaluate. 
%        Default = min(2^10,nx).
%   pf   Plot flag:  0=none (default), 1=screen.
%
%   CS   CrossSpectrogram output matrix.
%   t    Times at which the CrossSpectrogram was evaluated (s).
%   f    Frequencies at which the CrossSpectrogram was evaluated (Hz).
%
%   Plots the estimated cross-spectrogram of two non-stationary 
%   signals, x1 and x2.  The cross-spectrum is calculated using
%   modified periodogram.  The mean of each signal is removed before 
%   taking the FFT.  If only the window length is specified, the 
%   blackman window is used.
%
%   The x-axis on the contour plot represents the time in seconds 
%   or minutes, depending on the length of the signal.  The y-axis
%   represents the frequency, in hertz.  The colorbar represents 
%   the magnitude of the cross-spectral density.  The two original 
%   signals are plotted beneath the spectrogram.  
%
%   Example:  Plot the Cross-Spectrogram of an ABP and ICP signal, 
%   which are decimated to 12.5 Hz, using a Hanning window and a
%   window length of 100 s.
%
%      load ABPICP.mat
%      x1 = decimate(abp,10);
%      x2 = decimate(icp,10);
%      CrossSpectrogram(x1,x2,12.5,hanning(100*12.5));  
%
%   Hayes, M., "Statistical Digital Signal Processing and Modeling," 
%   John Wiley & Sons, 1996.
%
%   Version 1.00 LJ
%
%   See also Spectrogram and Cohereogram.

%====================================================================
% Error checking
%====================================================================
if  nargin < 2,
     help CrossSpectrogram;
     return;
end

nx1 = length(x1);
nx2 = length(x2);

if nx1~=nx2,
    fprintf('ERROR: Signals must be of the same length.\n');
    return;
end
nx = nx1;

if nx==0,
    fprintf('ERROR: Signal is empty.\n');
    return;
end

xr1 = x1;
mx1 = mean(x1);
sx1 = std(x1);
if sx1==0,
    fprintf('ERROR: Signal ''1'' is constant.\n');
    return;
end

xr2 = x2;
mx2 = mean(x2);
sx2 = std(x2);
if sx2==0,
    fprintf('ERROR: Signal ''2'' is constant.\n');
    return;
end

%====================================================================
% Preprocessing
%====================================================================    
x1  = x1(:)';  %  Make sure x1 is a row vector
x2  = x2(:)';  %  Make sure x2 is a row vector
 
if abs(mx1)/sx1>1e-5,
    x1 = x1 - mx1;  % Remove mean
end

if abs(mx2)/sx2>1e-5,
    x2 = x2 - mx2;  % Remove mean
end
    
%====================================================================
% Process function arguments, fill in defaults
%====================================================================
fs = 1; % Default sampling rate
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
end    

wl = 1000;          % Default window length
wn = blackman(wl)'; % Default window shape
if exist('wla') & ~isempty(wla),
    if length(wla)==1, % User specified window length only
        wl = max(3,round(wla*fs));
        if rem(wl,2),
            wl = wl - 1; % Make even
        end
        wn = blackman(wl)'; % Default window shape
    else % User specified the actual window, not just the length
        wn = wla(:)';
        if rem(length(wn),2),
            wn = wn(1:length(wn)-1); % Make even
        end
        wl = length(wn);
    end
end
    
Fmin = 0;
Fmax = fs/2;
if exist('fra') & ~isempty(fra),
    Fmin = max(fra(1),0);
    Fmax = min(fra(2),fs/2);
end   

nf = round(wl/2);
if exist('nfa') & ~isempty(nfa),
    nf = nfa;
end   
    
ns = min(2^10,nx);
if exist('nsa') & ~isempty(nsa),
    ns = min(nsa,nx);
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
st   = floor((nx-1)/(ns-1));                 % Step size
B    = st:st:nx;                             % Evaluation times
nb   = length(B);                            % No. points in B
wp   = ((fs/2)*(nf/(Fmax-Fmin))-1)*2;        % No. window points
wp   = 2^(ceil(log2(wp)));                   % Convert to power of 2 
b0   = floor(wp*(Fmin/fs))+1;                % Lower frequency index 
b1   = ceil (wp*(Fmax/fs))+1;                % Upper frequency index
fi   = b0:b1;                                % Frequency indices
f    = (fs*(fi-1)/wp)';                      % Frequencies
nf   = length(fi);

%====================================================================
% Main loop
%====================================================================
CS  = zeros(nf,nb);  % Power Spectral Density Matrix

for cnt       = 1:nb,
    b         = B(cnt);
    xi        = max(b-(wl/2-1),1);  % Signal initial index
    pi        = (wl/2-1)-(b-xi);    % Beginning zero index
    xe        = min(b+(wl/2)  ,nx); % Signal final index
    pe        = (wl/2  )-(xe-b);    % End zero index
    seg1      = [x1(xi)*zeros(1,pi) x1(xi:xe) x1(xe)*zeros(1,pe)];
    seg2      = [x2(xi)*zeros(1,pi) x2(xi:xe) x2(xe)*zeros(1,pe)];
    seg1      = seg1.*wn;           % Multiply by window type
    seg2      = seg2.*wn;           % Multiply by window type
	psd1      = fft(seg1,wp);  
    psd2      = fft(seg2,wp);
    csd       = abs(psd1 .* conj(psd2));
    CS(:,cnt) = csd(fi)';
end

%====================================================================
% Timing
%====================================================================
Fmax = max(f);
t    = ((B-1)/fs)';
tp   = t;
Tmax = (nx-1)/fs;
td   = 1; % Time Divider
        
if max(t)>2000,
    tp = tp/60;
    td = 60; % Use minutes
end
    
%====================================================================
% Plot Results
%==================================================================== 
if pf > 0
    figure;
    FigureSet(1);
    clf;

    %----------------------------------------------------------------
    % Cross-Spectrogram
    %----------------------------------------------------------------
    ha1 = axes('Position',[0.10 0.21 0.75 0.69]);
    w0 = pi*sqrt(2/log(2));
    s = reshape(CS,nf*nb,1);
    p = percentile(s,[0.001 0.975]);
    imagesc(tp,f,CS,p); %,[Smin Smax]);
    set(ha1,'XLim',[0 Tmax]./td);
    set(ha1,'XAxisLocation','Top');
    set(ha1,'YDir','normal');
    hold on;
        h = plot(((wl/2-1)/fs*[1 1])./td,[0 fs],'k');
        set(h,'LineWidth',2.0);
        h = plot(((wl/2-1)/fs*[1 1])./td,[0 fs],'w');
        set(h,'LineWidth',1.0);
        h = plot(((Tmax-(wl/2-1)/fs)*[1 1])./td,[0 fs],'k'); 
        set(h,'LineWidth',2.0);
        h = plot(((Tmax-(wl/2-1)/fs)*[1 1])./td,[0 fs],'w'); 
        set(h,'LineWidth',1.0);
        hold off;
    ylabel('Frequency (Hz)');
    title('Cross-Spectrogram');
    
    %----------------------------------------------------------------
    % Colorbar
    %----------------------------------------------------------------    
    ha2 = colorbar;    
    set(ha2,'Box','Off')
    set(ha2,'YTick',[]);
    set(ha2,'Position',[0.97 0.15 0.02 0.75]);
    set(ha1,'Position',[0.08 0.37 0.88 0.53]);
    
    %----------------------------------------------------------------
    % Signal 1
    %----------------------------------------------------------------      
    ha3 = axes('Position',[0.08 0.26 0.88 0.10]);
    nx  = length(x1);
    k   = 1:nx;
    tx  = (k-1)/fs;
    h   = plot(tx/td,xr1);
    set(h,'LineWidth',1.5);
    set(ha3,'xTickLabel',[]);
    ymin = min(xr1);
    ymax = max(xr1);
    yrng = ymax-ymin;
    ymid = (ymax+ymin)/2;
    ymin = ymid - 0.30*yrng;
    ymax = ymid + 0.30*yrng;
    rd   = 10^(3-floor(log10(ymax))); % 3 Significant figures
    ymin = (ymin*rd)/rd;
    ymax = (ymax*rd)/rd;
    ymid = (ymax+ymin)/2;
    yrng = ymax-ymin;
    if yrng==0,
        yrng = 1;
    end;
    set(ha3,'YTick',[ymin ymid ymax]);    
    ymin = min(xr1);
    ymax = max(xr1);
    yrng = ymax-ymin;
    ymin = ymin - 0.1*yrng;
    ymax = ymax + 0.1*yrng;
    set(ha3,'XLim',[0 Tmax]./td);
    set(ha3,'YLim',[ymin ymax]);
    ylabel('x1');
    
    %----------------------------------------------------------------
    % Signal 2
    %----------------------------------------------------------------      
    ha3 = axes('Position',[0.08 0.15 0.88 0.10]);
    nx  = length(x2);
    k   = 1:nx;
    tx  = (k-1)/fs;
    h   = plot(tx/td,xr2);
    set(h,'LineWidth',1.5);
    ymin = min(xr2);
    ymax = max(xr2);
    yrng = ymax-ymin;
    ymid = (ymax+ymin)/2;
    ymin = ymid - 0.30*yrng;
    ymax = ymid + 0.30*yrng;
    rd   = 10^(3-floor(log10(ymax))); % 3 Significant figures
    ymin = (ymin*rd)/rd;
    ymax = (ymax*rd)/rd;
    ymid = (ymax+ymin)/2;
    yrng = ymax-ymin;
    if yrng==0,
        yrng = 1;
    end;
    set(ha3,'YTick',[ymin ymid ymax]);    
    ymin = min(xr2);
    ymax = max(xr2);
    yrng = ymax-ymin;
    ymin = ymin - 0.1*yrng;
    ymax = ymax + 0.1*yrng;
    set(ha3,'XLim',[0 Tmax]./td);
    set(ha3,'YLim',[ymin ymax]);
    ylabel('x2');
    
    if td==1,
        xlabel('Time (seconds)');    
    else
        xlabel('Time (minutes)');
    end
    
    axes(ha1);
    AxisSet(8);
    set(gcf,'PaperPosition',[0.25 0.25 10.5 8])
end
    
if nargout == 0
    clear CS;
    clear t;
    clear f;
end
    