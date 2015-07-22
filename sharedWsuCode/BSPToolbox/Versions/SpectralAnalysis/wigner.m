function wigner(x,sra,fra,wla,nfa,nsa,tsa,pfa);
%Wigner: Estimate the smoothed wigner-ville distribution of a 
%signal
%
%   wigner(x,sr,fr,wl,nf,ns,ts,pf)
%
%   x    Input signal
%   sr   Sample rate. Default=125
%   fr   Array of minimum and maximum frequencies. 
%        Default=5*SR/NX to SR/2
%   wl   Length of window to use (default 128). Need not be power
%        of 2
%   nf   Number of frequencies to evaluate. Default=WL/2
%   ns   Requested number of times (horizontal pixels) to evaluate. 
%        Default=128)
%   ts   Time (in seconds) of the first element of the input signal 
%        Default=0
%   pf   Plot format: 0=none, 1=screen. Default=0.
%
%   This function estimates the wigner distribution of an input 
%   signal.The user can specify the horizontal pixels of the 
%   image to be generated. 
%
%   The Wigner Ville distribution is obtained by adding pieces made
%   up of the product of the signal at a past time multiplied by the 
%   signal at a future time, the time in the past being equal to the 
%   time in the future. 
%
%   Example: Estimate the wigner ville distribution of an intracranial
%   blood pressure at 125 Hz
%
%      load ICP.mat; 
%      wigner(icp, fs);
%
%   Leon Cohen, "Time-Frequency Analysis," Prentice Hall Signal 
%   Processing Series, pp.113-135, 1995.
%
%   Version 1.00 JB
%
%   See also BlackmanTukey, EigenVectorPSD, HeartRate,
%   MaximumEntropyPSD, MinimumNormPSD, ModifiedPeriodogram, Music,
%   and Welch.


% =================================================================
% Default parameters
% =================================================================
if ( nargin < 1 | nargin > 8)
    help wigner;
    return;
end;


%==================================================================
% Process function arguments, fill in defaults
%==================================================================
NX = length(x);
if NX==0,
    fprintf('ERROR: Signal is empty.\n');
    return;
    end;

sr = 1;
if exist('sra') & ~isempty(sra),
    sr = sra;
    end;    

Fmin = 0;
Fmax = sr/2;
if exist('fra') & ~isempty(fra),
    Fmin = max(fra(1),0);
    Fmax = min(fra(2),sr/2);
    end;    
    
%===================================================================
% Preprocessing
%===================================================================    
% Make x into a row vector
if size(x,2)==1,
    x = x';
    end;

xr = x;
mx = mean(x);
sx = std(x);
if sx==0,
    fprintf('ERROR: Signal is constant.\n');
    return;
    end;


if abs(mx)/sx>1e-5,
    %fprintf('WARNING: Removing mean from signal.\n');
    x = x - mx;
    end;

%==================================================================
% Finish processing function arguments, fill in defaults
%==================================================================
wl = 512;
if exist('wla') & ~isempty(wla),
    wl = wla; 
    end;
    
nf = round(wl/2);
if exist('nfa') & ~isempty(nfa),
    nf = nfa;
    end;   
    
ns = min(2^10,NX);
if exist('nsa') & ~isempty(nsa),
    ns = min(nsa,NX);
    end;     
    
ts = 0;
if exist('tsa') & ~isempty(tsa),
    ts = tsa;
    end;     
pf  = 0;                                % Default - no plotting
if exist('pfa') & ~isempty('pfa'),
    pf = pfa;   
end    

%===================================================================
% Variable Allocations & Initialization
%===================================================================
NX   = length(x);                         % No. points in X
ST   = floor((NX-1)/(ns-1));              % Step size
B    = ST+wl:ST:NX-wl;                    % Evaluation times
B1   = ST:ST:NX;
NB1  = length(B1);
NB   = length(B);                         % No. points in B
WP   = ((sr/2)*(nf/(Fmax-Fmin))-1)*2;     % No. window points
WP   = 2^(ceil(log2(WP)));                % Convert to power of 2 
b0   = floor(Fmin/(sr/2)*(1+WP)/2)+1;     % Lower frequency bin index 
b1   = ceil (Fmax/(sr/2)*(1+WP)/2)+1;     % Upper frequency bin index
fi   = b0:b1;                             % Frequency indices
f    = (fi-1)*(sr/2)/((1+WP));            % Frequencies
nf   = length(fi);

%====================================================================
% Main loop
%====================================================================
S  = zeros(nf,NB);  % Power Spectral Density
wn = blackman(wl/2)'; % Window
wn = wn;            % Ensure it passes DC
for m = 1:NB
    up   = (B(m) + (wl)/2)-1;
    down = (B(m) - (wl)/2)+1;
    
    
    s1   =      x(B(m) :-1 :down);
    s2   = conj(x(B(m) : 1 :up));
    ss   = s1.*s2;
    ss   = sqrt(ss);
   
    s    = ss.*wn;
    shat = abs(fft(s,WP));
    shat  = shat(fi);
    shat  = shat /length(fi);
    wvd(m,:) = shat;
end
upcorr = [];
downcorr = [];
diff1 = ST+wl:-ST:ST;
for   t = 1:length(diff1)
    upcorr = [wvd(1,:); upcorr];
end
diff2 = NX-wl:ST:NX;
for   t = 1:length(diff2)
    downcorr = [wvd(NB,:); downcorr];
end
   
wvd = [ upcorr ; wvd; downcorr];
wvd = wvd(1:NB1, :);

NB = NB1;
B = B1;
S  = wvd';
t = (B-1)/sr;

TD = 1; % Time Divider
if ts+max(t)>2000,
    TD = 60; % Use minutes
    end;

if 0, % Generate image of ranks - good technique
    t1    = reshape(S,nf*NB,1);
    [Y,I] = sort(t1);
    z     = zeros(size(I));
    z(I)  = (1:length(I))';
    S     = reshape(z,nf,NB);
    fprintf('WARNING: Using spectrogram ranks.\n');
    end;

%====================================================================
% Plot Results
%====================================================================  
if (pf ==1 | nargout == 0)
    Tmax = (NX-1)/sr;
    f1=figure;
    FigureSet(1);
    clf;

    %----------------------------------------------------------------
    % Spectrogram
    %----------------------------------------------------------------
    ha1 = axes('Position',[0.10 0.21 0.75 0.69]);
    w0 = pi*sqrt(2/log(2));
    s = reshape(S,nf*NB,1);
    p = percentile(s,[0.001 0.975]);
    imagesc((t+ts)/TD,f,abs(S),p); %,[Smin Smax]);
    set(ha1,'XLim',[ts ts+Tmax]./TD);
    set(ha1,'XAxisLocation','Top');
    set(ha1,'YDir','normal');
    hold on;
        h = plot((ts+(wl/2-1)/sr*[1 1])./TD,[0 sr],'k');
        set(h,'LineWidth',2.0);
        h = plot((ts+(wl/2-1)/sr*[1 1])./TD,[0 sr],'w');
        set(h,'LineWidth',1.0);
        h = plot((ts+(Tmax-(wl/2-1)/sr)*[1 1])./TD,[0 sr],'k'); 
        set(h,'LineWidth',2.0);
        h = plot((ts+(Tmax-(wl/2-1)/sr)*[1 1])./TD,[0 sr],'w'); 
        set(h,'LineWidth',1.0);
        hold off;
    ylabel('Frequency (Hz)');
    
    %----------------------------------------------------------------
    % Colorbar
    %----------------------------------------------------------------    
    ha2 = colorbar;    
    set(ha2,'Box','Off')
    set(ha2,'YTick',[]);
    set(ha2,'Position',[0.97 0.10 0.02 0.80]);
    set(ha1,'Position',[0.08 0.21 0.88 0.69]);
    
    %----------------------------------------------------------------
    % Signal
    %----------------------------------------------------------------      
    ha3 = axes('Position',[0.08 0.10 0.88 0.10]);
    NX  = length(x);
    k   = 1:NX;
    tx  = ts+(k-1)/sr;
    h   = plot(tx/TD,xr);
    set(h,'LineWidth',1.5);
    %set(h,'Marker','.');
    %set(h,'MarkerSize',6);
    ymin = min(xr);
    ymax = max(xr);
    yrng = ymax-ymin;
    ymid = (ymax+ymin)/2;
    ymin = ymid - 0.30*yrng;
    ymax = ymid + 0.30*yrng;
    rd   = 10^(3-floor(log10(ymax))); % 3 Significant figures
    ymin = round(ymin*rd)/rd;
    ymax = round(ymax*rd)/rd;
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
    set(ha3,'XLim',[ts ts+Tmax]./TD);
    set(ha3,'YLim',[ymin ymax]);
    if TD==1,
        xlabel('Time (seconds)');    
    else
        xlabel('Time (minutes)');
        end;
    
    axes(ha1);
    AxisSet(8);
    set(gcf,'PaperPosition',[0.25 0.25 10.5 8])
    title('BSP Laboratory -- Wigner-Ville Distribution');

    end;