function [A,t,d] = Autocorrelogram(x,SRarg,DRarg,WLarg,NDarg,NSarg,TSarg);
% [A,t,d] = Autocorrelogram(x,fs,DR,WL,ND,NS,TS);
%   Calculates estimates of the spectral content at the specified times
%   using a blackman window of the specified length.
%
%   x : Input signal.
%   fs: Sample rate (default = 1).
%   DR: Array of minimum and maximum delays in seconds (default 0 to (WL-1)/fs).
%   WL: Length of window to use (default 512). Must be even.
%   ND: Number of frequencies to evaluate (default WL)
%   NS: Requested number of times (horizontal pixels) to evaluate 
%       (default 512).
%   TS: Time (in seconds) of the first element of the input signal 
%       (default 0)

%clear all;
SCRIPT = 0;

if SCRIPT==1,
    SRarg  = 125;
    WLarg  = 256;
    NSarg  = 100;
    NFarg  = 256;
    TSarg  = 100;
    
    NK = 2^12;
    k = 1:NK;
    t = (k-1)/SRarg;
    x = 1 + Tweet(t,10,0.1);
    %x = 6*(k==1) + 6*(k==round(NK/2)) + 1;
    %x = 5+5*(k-1)/NK + 0.1*randn(size(k));
    %x = (k-1)/NK;
    %x = ((k-1)/NK).^2;
    %x = sin(2*pi*9.5*t);
    %x = sawtooth(2*pi*1*t);
    %x = square(2*pi*1*t);
    %x = 5*(k-1)/NK + chirp(t,1,max(t),5,'logarithmic');
    %x = k.*sin(2*pi*10*k/NK);
    %x = (k-1)/NK + sin(2*pi*10*k/NK);
    %x = sin(2*pi*10*k/NK) + 2*randn(size(k));
    %x = exp(-((k-NK/2)./(2*NK/25)).^2.);
    %x = cos(pi*(k-round(NK/2))/NK);
    %x = ones(size(k));
    %x = randn(size(k));
    end;

%====================================================================
% Process function arguments, fill in defaults
%====================================================================
NX = length(x);

SR = 1;
if exist('SRarg') & ~isempty(SRarg),
    SR = SRarg;
    end;    
    
WL = 512;
if exist('WLarg') & ~isempty(WLarg),
    WL = WLarg; 
    end;
    
ND = WL;
if exist('NDarg') & ~isempty(NDarg),
    ND = NDarg;
    end;   
    
NS = min(2^10,NX);
if exist('NSarg') & ~isempty(NSarg),
    NS = min(NSarg,NX);
    end;     
    
TS = 0;
if exist('TSarg') & ~isempty(TSarg),
    TS = TSarg;
    end;    
    
Dmin = 0;
Dmax = WL-1;
if exist('DRarg') & ~isempty(DRarg),
    if length(DRarg)==1,
        Dmin = 0;
        Dmax = min(round(DRarg(2)*SR),Dmax);
    else
        Dmin = max(round(DRarg(1)*SR),Dmin);
        Dmax = min(round(DRarg(2)*SR),Dmax);
        end;    
    end;    
    
%====================================================================
% Preprocessing
%====================================================================    
% Make x into a row vector
if size(x,2)==1,
    x = x';
    end; 

%====================================================================
% Variable Allocations & Initialization
%====================================================================
NX   = length(x);                % No. points in X
st   = floor((NX-1)/(NS-1));     % Step size
B    = st:st:NX;                 % Evaluation times
NB   = length(B);                % No. points in B
st   = floor((Dmax-Dmin)/ND);    % Step size
st   = max(1,st);                         
di   = Dmin:st:Dmax;             % Autocorrelation delay indices
ND   = length(di);

%====================================================================
% Main loop
%====================================================================
A  = zeros(ND,NB);  % Autocorrelogram
for c1 = 1:NB,
    b   = B(c1);
    xi  = max(b-(WL/2-1),1); % Signal initial index
    pi  = (WL/2-1)-(b-xi);
    xf  = min(b+(WL/2)  ,NX); % Signal final index
    pf  = (WL/2  )-(xf-b);
    y   = [x(xi)*ones(1,pi) x(xi:xf) x(xf)*ones(1,pf)];
    
    a = Autocorrelation(y,di,'fast')';
    %a = a(di); 
    
    if 0,
        subplot(2,1,1);
        plot(y);
        subplot(2,1,2);
        plot(a);
        fprintf('Pausing...\n'); pause;
        end;
    A(:,c1) = a';
	end;
t = (B-1)/SR;

TD = 1; % Time Divider
if max(t)>60*60,
    TD = 60; % Use minutes
    end;


%====================================================================
% Plot Results
%====================================================================  
if SCRIPT | nargout==0,
    Tmax = (NX-1)/SR;
    figure;
    FigureSet(1);
    clf;

    %----------------------------------------------------------------
    % Autocorrelogram
    %----------------------------------------------------------------
    ha1 = axes('Position',[0.10 0.21 0.75 0.69]);
    w0 = pi*sqrt(2/log(2));
    s = reshape(A,ND*NB,1);
    p = percentile(s,[0.001 0.975]);
    p = [-1 1];
    d = (di-1)/SR;
    imagesc((t+TS)/TD,d,A,p); %,[Smin Smax]);
    set(ha1,'XLim',[TS TS+Tmax]./TD);
    set(ha1,'XAxisLocation','Top');
    set(ha1,'YDir','normal');
    hold on;
        h = plot((TS+(WL/2-1)/SR*[1 1])./TD,[0 SR],'k');
        set(h,'LineWidth',2.0);
        h = plot((TS+(WL/2-1)/SR*[1 1])./TD,[0 SR],'w');
        set(h,'LineWidth',1.0);
        h = plot((TS+(Tmax-(WL/2-1)/SR)*[1 1])./TD,[0 SR],'k'); 
        set(h,'LineWidth',2.0);
        h = plot((TS+(Tmax-(WL/2-1)/SR)*[1 1])./TD,[0 SR],'w'); 
        set(h,'LineWidth',1.0);
        hold off;
    ylabel('Delay (sec)');

    %----------------------------------------------------------------
    % Colorbar
    %----------------------------------------------------------------    
    ha2 = colorbar;    
    set(ha2,'Box','Off')
    set(ha2,'Position',[0.95 0.10 0.02 0.80]);
    set(ha1,'Position',[0.10 0.21 0.84 0.69]);
    
    %----------------------------------------------------------------
    % Signal
    %----------------------------------------------------------------      
    ha3 = axes('Position',[0.10 0.10 0.84 0.10]);
    
    NX  = length(x);
    k   = 1:NX;
    tx  = TS+(k-1)/SR;
    h   = plot(tx/TD,x);
    set(h,'LineWidth',1.5);
    ymin = min(x);
    ymax = max(x);
    rd   = 10^(3-floor(log10(ymax))); % 3 Significant figures
    ymin = round(ymin*rd)/rd;
    ymax = round(ymax*rd)/rd;
    ymid = (ymax+ymin)/2;
    yrng = ymax-ymin;
    if yrng==0,
        yrng = 1;
        end;
    set(ha3,'YTick',[ymin ymid ymax]);    
    ymin = ymin - 0.1*yrng;
    ymax = ymax + 0.1*yrng;
    set(ha3,'XLim',[TS TS+Tmax]./TD);
    set(ha3,'YLim',[ymin ymax]);
    if TD==1,
        xlabel('Time (seconds)');    
    else
        xlabel('Time (minutes)');
        end;
    
    axes(ha1);
    AxisSet(8);
    set(gcf,'PaperPosition',[0.25 0.25 10.5 8])
    end;