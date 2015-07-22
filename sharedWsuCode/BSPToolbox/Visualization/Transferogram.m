function [TF,PH,t,f] = Transferogram(x,y,fsa,wla,fra,nfa,nsa,pfa);
%Transferogram: Nonstationary visualization of transfer function.
%
%   [TF,PH,t,f] = Transferogram(x,y,fs,wl,fr,nf,ns,pf);
%
%   x    Input signal.
%   y    Output signal.
%   fs   Sample rate, hertz.  Default = 1 Hz.
%   wl   Length of window to use (sec). Default = 1000 samples.
%        If a vector, specifies entire window.
%   fr   Array of minimum and maximum frequencies.  
%        Default=[0 fs/2].
%   nf   Number of frequencies to evaluate, default=wl/2.
%   ns   Requested number of times (horizontal pixels) to evaluate. 
%        Default=min(2^10,nx).
%   pf   Plot flag:  0=none (default), 1=screen.
%
%   TF   Transfer function magnitude output matrix.
%   PH   Transfer function phase output matrix.
%   t    Times at which the transfer function were evaluated (s).
%   f    Frequencies at which the transfer function were evaluated 
%        (Hz).
%
%   Plots the estimated transfer function of two non-stationary 
%   signals, x and y, x being the input signal and y being the 
%   output signal.  The transfer function is calculated using
%   modified periodogram and is defined as:
%
%         tf = Sxy./abs(Sxx)
%
%   where Sxy is the cross-power spectral density and Sxx is the 
%   power spectral density of the signal x.
%
%   The mean of each signal is removed before taking the FFT.  If only
%   the window length is specified, the blackman window is used.  The
%   output, when plotted, will display two surface plots, one for the
%   magnitude of the transfer function and one for the phase.  If no
%   output arguments are specified, the default will print to the 
%   screen.
%
%   Example:  Plot the transfer function of an ABP and ICP signal, 
%   which are decimated to 25 Hz, where icp is the input signal and
%   abp is the output signal.  Use a Hanning window with a window 
%   length of 50 s and plot the output to the screen.
%
%      load ABPICP.mat
%      y = decimate(abp,5);
%      x = decimate(icp,5);
%      Transferogram(x,y,25,hanning(50*25));  
%
%   Bendat, J., and Piersol, A., "Engineering Applications of
%   Correlation and Spectral Analysis", John Wiley & Sons, pp.97-114,
%   1980.
%
%   Version 1.00 LJ
%
%   See also Spectrogram and Cohereogram.

%====================================================================
% Error Check
%====================================================================
if  nargin < 2,
     help Transferogram;
     return;
end

nx = length(x);
ny = length(y);

if nx~=ny,
    fprintf('ERROR: Signals must be of the same length.\n');
    return;
end;

if nx==0,
    fprintf('ERROR: Signal is empty.\n');
    return;
end;

%====================================================================
% Preprocessing
%====================================================================    
% Make x into a row vector
if size(x,2)==1,
    x = x';
    end;

if size(y,2)==1,
    y = y';
    end;
    
xr = x;
mx = mean(x);
sx = std(x);
if sx==0,
    fprintf('ERROR: Signal ''1'' is constant.\n');
    return;
end;

yr = y;
my = mean(y);
sy = std(y);
if sy==0,
    fprintf('ERROR: Signal ''2'' is constant.\n');
    return;
end;

if abs(mx)/sx>1e-5,
    x = x - mx;
end

if abs(my)/sy>1e-5,
    y = y - my;
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
b0   = floor(wp*(Fmin/fs))+1;                % Lower frequency bin index 
b1   = ceil (wp*(Fmax/fs))+1;                % Upper frequency bin index
fi   = b0:b1;                                % Frequency indices
f    = (fs*(fi-1)/wp)';                         % Frequencies
nf   = length(fi);

%====================================================================
% Main loop
%====================================================================
TF = zeros(nf,nb);  
PH = zeros(nf,nb);

for cnt       = 1:nb,
    b         = B(cnt);
    xb        = max(b-(wl/2-1),1);  % Signal initial index
    pb        = (wl/2-1)-(b-xb);    % Beginning zero index
    xe        = min(b+(wl/2)  ,nx); % Signal final index
    pe        = (wl/2  )-(xe-b);    % End zero index
    x1        = [x(xb)*zeros(1,pb) x(xb:xe) x(xe)*zeros(1,pe)];
    y1        = [y(xb)*zeros(1,pb) y(xb:xe) y(xe)*zeros(1,pe)];
    a1        = x1.*wn;             % Window type
    y1        = y1.*wn;             % Window type
    psdx      = fft(a1,wp);
    psdy      = fft(y1,wp);
    Sxx       = abs(psdx.^2);
    Sxy       = conj(psdx).*psdy;
    H         = Sxy./Sxx;
    TF(:,cnt) = abs(H(fi)');
    PH(:,cnt) = angle(H(fi)').*180/pi;
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
    % Magnitude Transferogram
    %----------------------------------------------------------------
    ha1 = axes('Position',[0.10 0.21 0.75 0.69]);
    w0 = pi*sqrt(2/log(2));
    s = reshape(TF,nf*nb,1);
    p = prctile(s,[2.5 97.5]); %p = prctile(s,[0.025 0.975]);
    imagesc(tp,f,TF,p); 
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
    title('Transfer Function - Magnitude');
    
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
    nx  = length(x);
    k   = 1:nx;
    tx  = (k-1)/fs;
    h   = plot(tx/td,xr);
    set(h,'LineWidth',1.5);
    set(ha3,'xTickLabel',[]);
    ymin = min(xr);
    ymax = max(xr);
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
    ymin = min(xr);
    ymax = max(xr);
    yrng = ymax-ymin;
    ymin = ymin - 0.1*yrng;
    ymax = ymax + 0.1*yrng;
    set(ha3,'XLim',[0 Tmax]./td);
    set(ha3,'YLim',[ymin ymax]);
    ylabel('x');
    
    %----------------------------------------------------------------
    % Signal 2
    %----------------------------------------------------------------      
    ha3 = axes('Position',[0.08 0.15 0.88 0.10]);
    nx  = length(y);
    k   = 1:nx;
    tx  = (k-1)/fs;
    h   = plot(tx/td,yr);
    set(h,'LineWidth',1.5);
    ymin = min(yr);
    ymax = max(yr);
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
    ymin = min(yr);
    ymax = max(yr);
    yrng = ymax-ymin;
    ymin = ymin - 0.1*yrng;
    ymax = ymax + 0.1*yrng;
    set(ha3,'XLim',[0 Tmax]./td);
    set(ha3,'YLim',[ymin ymax]);
    ylabel('y');
    
    if td==1,
        xlabel('Time (seconds)');    
    else
        xlabel('Time (minutes)');
    end
    
    axes(ha1);
    AxisSet(8);
    set(gcf,'PaperPosition',[0.25 0.25 10.5 8])
    figure;
    FigureSet(2);
    clf;
    
    %----------------------------------------------------------------
    % Phase Transferogram
    %----------------------------------------------------------------
    ha4 = axes('Position',[0.10 0.21 0.75 0.69]);
    w0 = pi*sqrt(2/log(2));
    s = reshape(PH,nf*nb,1);
    p = prctile(s,[2.5 97.5]); %p = prctile(s,[0.025 0.975]);
    imagesc(tp,f,PH,p); 
    set(ha4,'XLim',[0 Tmax]./td);
    set(ha4,'XAxisLocation','Top');
    set(ha4,'YDir','normal');
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
    ylabel('Phase, degrees');
    title('Transfer Function - Phase');
    
    %----------------------------------------------------------------
    % Colorbar
    %----------------------------------------------------------------    
    ha5 = colorbar;    
    set(ha5,'Box','Off')
    set(ha5,'YTick',[]);
    set(ha5,'Position',[0.97 0.15 0.02 0.75]);
    set(ha4,'Position',[0.08 0.37 0.88 0.53]);
    
    %----------------------------------------------------------------
    % Signal 1
    %----------------------------------------------------------------      
    ha6 = axes('Position',[0.08 0.26 0.88 0.10]);
    nx  = length(x);
    k   = 1:nx;
    tx  = (k-1)/fs;
    h   = plot(tx/td,xr);
    set(h,'LineWidth',1.5);
    set(ha6,'xTickLabel',[]);
    ymin = min(xr);
    ymax = max(xr);
    yrng = ymax-ymin;
    ymid = (ymax+ymin)/2;
    ymin = ymid - 0.30*yrng;
    ymax = ymid + 0.30*yrng;
    rd   = 10^(3-floor(log10(ymax))); % 3 Significant figures
    ymin = round((ymin*rd)/rd);
    ymax = round((ymax*rd)/rd);
    ymid = round((ymax+ymin)/2);
    yrng = ymax-ymin;
    if yrng==0,
        yrng = 1;
    end;
    set(ha6,'YTick',[ymin ymid ymax]);    
    ymin = min(xr);
    ymax = max(xr);
    yrng = ymax-ymin;
    ymin = ymin - 0.1*yrng;
    ymax = ymax + 0.1*yrng;
    set(ha6,'XLim',[0 Tmax]./td);
    set(ha6,'YLim',[ymin ymax]);
    ylabel('x');
    
    %----------------------------------------------------------------
    % Signal 2
    %----------------------------------------------------------------      
    ha6 = axes('Position',[0.08 0.15 0.88 0.10]);
    nx  = length(y);
    k   = 1:nx;
    tx  = (k-1)/fs;
    h   = plot(tx/td,yr);
    set(h,'LineWidth',1.5);
    ymin = min(yr);
    ymax = max(yr);
    yrng = ymax-ymin;
    ymid = (ymax+ymin)/2;
    ymin = ymid - 0.30*yrng;
    ymax = ymid + 0.30*yrng;
    rd   = 10^(3-floor(log10(ymax))); % 3 Significant figures
    ymin = round((ymin*rd)/rd);
    ymax = round((ymax*rd)/rd);
    ymid = round((ymax+ymin)/2);
    yrng = ymax-ymin;
    if yrng==0,
        yrng = 1;
    end;
    set(ha6,'YTick',[ymin ymid ymax]);    
    ymin = min(yr);
    ymax = max(yr);
    yrng = ymax-ymin;
    ymin = ymin - 0.1*yrng;
    ymax = ymax + 0.1*yrng;
    set(ha6,'XLim',[0 Tmax]./td);
    set(ha6,'YLim',[ymin ymax]);
    ylabel('y');
    
    if td==1,
        xlabel('Time (seconds)');    
    else
        xlabel('Time (minutes)');
    end
    
    axes(ha4);
    AxisSet(8);
    set(gcf,'PaperPosition',[0.25 0.25 10.5 8])
end

if nargout == 0
    clear TF;
    clear PH;
    clear f;
    clear t;
end
    