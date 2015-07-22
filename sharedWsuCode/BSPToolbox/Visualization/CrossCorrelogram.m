function [C,t,d] = CrossCorrelogram(x1,x2,fsa,dra,wla,nfa,nsa,tsa,pfa);
% CrossCorrelogram: Estimate and plot the cross-correlation of 
% two non-stationary signals.
%
%   [C,t,d] = CrossCorrelogram(x1,x2,fs,dr,wl,nf,ns,ts,pf);
%
%   x1   Input signal.
%   x2   Input signal.
%   fs   Sample rate.  Default=1.
%   dr   Maximum delays in seconds on verticle axis.  Default=wl/2.
%   wl   Length of window to use.  Default=512. 
%   nf   Number of frequencies to evaluate.  Default=wl.
%   ns   Requested number of times (horizontal pixels) to evaluate. 
%        default = 1024.
%   ts   Time (in seconds) of the first element of the input signal, 
%        Default=0.
%   pf   Plot flag:  0=none (default), 1=screen.
% 
%   C    CrossCorrelogram output matrix
%   t    Times at which the CrossCorrelogram was evaluated (s).
%   d    Delays at which the CrossCorrelogram was evaluated (s).
%
%   Plots the estimated cross-correlation of two non-stationary 
%   signals, x1 and x2.  The cross-correlation is calculated 
%   across overlapping windows of the signals and plotted using 
%   a surface plot.  The top and bottom 2.5 percentile of each 
%   cross-correlation is removed before plotting.  The x-axis 
%   represents the time in seconds or minutes, depending on the 
%   length of the signal.  The y-axis represents the lag in seconds.
%   The colorbar represents the magnitude of the cross-correlation, 
%   and the two original signals are plotted beneath the correlogram.  
%
%   Example:  Plot the Cross-Correlogram of an ABP and ICP signal 
%   with a sampling rate of 125 Hz and delay window of -1:1.
%
%      load ABPICP.mat
%      [C,t,d] = CrossCorrelogram(abp,icp,125,1,[],[],[],[],1);  
%
%   Challis, R. E., and Kitney, R. I., "Biomedical Signal Processing
%   (in four parts), Part 1, Time-domain Methods", Med. & Biol. Eng.
%   & Comput., 1990, 28, pp. 509-524.
%
%   Version 1.00 LJ
%
%   See also AutoCorrelogram, Autocorrelation, and CrossCorrelation.

%====================================================================
% Error Check
%====================================================================
if nargin < 2 
    help CrossCorrelogram;
    return;
end

%====================================================================
% Process Input Signals
%====================================================================
nx1 = length(x1);
nx2 = length(x2);

if nx1~=nx2,
    fprintf('ERROR: Signals must be of the same length.\n');
    return;
end;
nx = nx1;

if nx==0,
    fprintf('ERROR: Signal is empty.\n');
    return;
end

x1 = x1(:)';
x2 = x2(:)';

%====================================================================
% Process function arguments, fill in defaults
%====================================================================
fs = 1; % Default sampling rate
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
end   

wl = 512; % Default window length
if exist('wla') & ~isempty(wla),
    wl = wla; 
end

nf = wl; % Default number of frequencies
if exist('nfa') & ~isempty(nfa),
    nf = nfa;
end   
 
ns = min(2^10,nx); % Default horizontal pixel
if exist('nsa') & ~isempty(nsa),
    ns = min(nsa,nx);
end     

ts = 0; % Default time 
if exist('tsa') & ~isempty(tsa),
    ts = tsa;
end    

dmin = -wl./2;   % Default Neg Delay index
dmax =  wl./2-1; % Default Pos Delay index
if exist('dra') & ~isempty(dra),
    dmax = max(round(dra*fs/2),dmin);
    dmin = min(round(-dra*fs/2-1),dmax);
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
st   = floor((nx-1)/(ns-1));         % Step size
et   = st:st:nx;                     % Evaluation times
nb   = length(et);                    % No. points in B
st   = floor((dmax-dmin)/nf);        % Step size
st   = max(1,st);                         
di   = dmin:st:dmax;                 % CrossCorrelation delay indices
nf   = (length(di)*2)-1;             % Length of -delay:+delay
nl   = floor(nf./2)+1;               % Length of 1:delay to use for
                                     % CrossCorrelate function
                                     
%====================================================================
% Main loop
%====================================================================
C  = zeros(nf,nb);                   % CrossCorrelogram matrix
for c1   = 1:nb,
    b    = et(c1);
    xi   = max(b-(wl/2-1),1);        % Signal initial index
    vi   = (wl/2-1)-(b-xi);
    xf   = min(b+(wl/2)  ,nx);       % Signal final index
    vf   = (wl/2  )-(xf-b);
    y1   = [x1(xi)*zeros(1,vi) x1(xi:xf) x1(xf)*zeros(1,vf)];
    y2   = [x2(xi)*zeros(1,vi) x2(xi:xf) x2(xf)*zeros(1,vf)];
    c    = CrossCorrelation(y1,y2,nl,'fast')';
    C(:,c1) = c';
end

%====================================================================
% Timing
%====================================================================
t   = ((et-1)/fs)';
d   = (2.*(di-1)/fs)';
td  = 1; % Time Divider
if max(t)>60*60,
    td = 60; % Use minutes
end

%====================================================================
% Plot Results
%====================================================================  
if pf > 0
    tmax = (nx-1)/fs;
    tw   = [(wl/2-1)/fs,(tmax-(wl/2)/fs)];
    figure;
    FigureSet(1);

    %----------------------------------------------------------------
    % Correlogram
    %----------------------------------------------------------------
    ha1 = axes('Position',[0.10 0.21 0.75 0.69]);
    w0  = pi*sqrt(2/log(2));
    s   = reshape(C,nf*nb,1);
    p   = prctile(s,[2.5 97.5]);
    imagesc((t+ts)/td,d,C,p); 
    set(ha1,'XLim',[ts ts+tmax]./td);
    set(ha1,'XAxisLocation','Top');
    set(ha1,'YDir','normal');
    title('CrossCorrelogram');
    hold on;
        h = plot([1;1]*tw,[dmin dmax],'k');
        set(h,'LineWidth',2.0);
        h = plot([1;1]*tw,[dmin dmax],'w');
        set(h,'LineWidth',1.0);
        hold off;
    ylabel('Delay (sec)');

    %----------------------------------------------------------------
    % Colorbar
    %----------------------------------------------------------------    
    ha2 = colorbar;    
    set(ha2,'Box','Off')
    set(ha2,'Position',[0.90 0.15 0.02 0.75]);
    set(ha1,'Position',[0.08 0.37 0.81 0.53]);
    
    %----------------------------------------------------------------
    % Signal 1
    %----------------------------------------------------------------      
    ha3  = axes('Position',[0.08 0.26 0.81 0.10]);
    set(ha3,'Visible','Off');
    k    = 1:nx;
    tx   = ts+(k-1)/fs;
    h    = plot(tx/td,x1(k));
    set(h,'LineWidth',1.5);
    set(ha3,'xTickLabel',[]);
    ymin = min(x1);
    ymax = max(x1);
    rd   = 10^(3-floor(log10(ymax))); % 3 Significant figures
    ymin = (ymin*rd)/rd;
    ymax = (ymax*rd)/rd;
    ymid = (ymax+ymin)/2;
    yrng = ymax-ymin;
    if yrng==0,
        yrng = 1;
    end;
    set(ha3,'YTick',[ymin ymid ymax]);    
    ymin = ymin - 0.1*yrng;
    ymax = ymax + 0.1*yrng;
    set(ha3,'XLim',[ts ts+tmax]./td);
    set(ha3,'YLim',[ymin ymax]);
    ylabel('x1');
    
    %----------------------------------------------------------------
    % Signal 2
    %----------------------------------------------------------------      
    ha4  = axes('Position',[0.08 0.15 0.81 0.10]);
    k    = 1:nx;
    tx   = ts+(k-1)/fs;
    h    = plot(tx/td,x2(k));
    set(h,'LineWidth',1.5);
    ymin = min(x2);
    ymax = max(x2);
    rd   = 10^(3-floor(log10(ymax))); % 3 Significant figures
    ymin = (ymin*rd)/rd;
    ymax = (ymax*rd)/rd;
    ymid = (ymax+ymin)/2;
    yrng = ymax-ymin;
    if yrng==0,
        yrng = 1;
        end;
    set(ha4,'YTick',[ymin ymid ymax]);    
    ymin = ymin - 0.1*yrng;
    ymax = ymax + 0.1*yrng;
    set(ha4,'XLim',[ts ts+tmax]./td);
    set(ha4,'YLim',[ymin ymax]);
    ylabel('x2');
        
    if td==1,
        xlabel('Time (seconds)');    
    else
        xlabel('Time (minutes)');
        end;
    
    axes(ha1);
    AxisSet(8);
    end
    
if nargout==0
    clear C;
    clear t;
    clear d;
    end