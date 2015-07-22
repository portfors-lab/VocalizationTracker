function [di] = DetectionPlots(x,fsa,thia,eda,fcfa);
%DetectionPlots: Displays plots to inspect detection performance.
%   
%   [di,th] = DetectionPlots(x,fs,thi,ed,fcf);
%
%   x    Input signal.
%   fs   Sample rate (Hz). Default = 1 Hz.
%   thi  If scalar, the initial threshold. Default = Detected beat 
%        indices (samples). Default = 4*RSD (RSD is a robust estimate
%        of the signal standard deviation).
%        If vector, lists the detection indices.
%   ed   Event duration (s). Default = median inter-event interval.
%   fcf  Flag to create new figure? 1 = Yes (default), 0 = No.
%
%   di   Detected event indices.
%
%   A handly tool for manually selecting a threshold for 
%   threshold detection. Plots both the histogram of 
%   detected indices and an overlap plot of the detected
%   events.
%
%   Example: Edit detected beats of a noisy electrocardiogram. 
%
%      load MER.mat; 
%      DetectionPlots(x,fs,si);
%
%   Task Force of the European Society of Cardiology and the American 
%   Society of Pacing and Electrophysiology, "Heart rate variability:
%   Standards of measurement, physiological interpretation, and 
%   clinical use," Circulation, vol. 93, no. 5, pp. 1043-1065, Aug. 
%   1996.
%
%   Version  0.00.01.01 JM
%
%   See also EditAnnotations. 

%====================================================================
% Error Checking
%====================================================================    
if nargin<1,
    help DetectionPlots;
    return;
    end;
    
if exist('eda') & length(eda)>1,
    error('ed argument must be a scalar.');
    end;

%=====================================================================
% Process function arguments
%=====================================================================
fs = 1;
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
    end;        
        
thi = 0;
if exist('thia') & ~isempty(thia),
    thi = thia;
else
    rsd = 1.483*median(abs(x-median(x)));
    thi = 4*rsd;
    end;

ed = 0;
if exist('eda') & ~isempty(eda),
    ed = eda;
    end;      
       
fcf = 1;
if exist('fcfa') & ~isempty(fcfa),
    fcf = fcfa;
    end;

%====================================================================
% Author-specified Parameters
%====================================================================
nbph = 100;                                                % Number of bins for the peak histogram
nbop =  50;                                                % Number of bins for the overlap histogram
ncs  =  10;                                                % Number of close spikes to plot

%====================================================================
% Preprocessing
%====================================================================
x   = x(:);                                                % Make into a column vector
nx  = length(x);
pi  = DetectMaxima(x,mean(x),1);

if length(thi)==1,
    th = thi;
    di = pi(find(x(pi)> th));
    ni = pi(find(x(pi)<=th));
else
    di  = thi;                                             % Function was called with detection indices
    ni  = setdiff(pi,di);                                  % Not-detection indices are those peaks not in di
    
    xmin = min(x([pi]));
	xmax = max(x([pi]));
    bc   = linspace(xmin,xmax,nbph)';                          % Bin centers
    te   = zeros(nbph,1);
    for c1=1:length(bc),
        te(c1) = sum(x(di)<bc(c1)) + sum(x(ni)>bc(c1)); 
        end;
    [jnk,imin] = min(te);
    th  = bc(imin); 
    end;

di = di(:);
ni = ni(:);
    
if ed==0,
    if length(di)>1,
        ed = median(diff(sort(di))/fs);
    else
        ed = 50/fs; 
        end;
    end;

%====================================================================
% Figure Initialization & Axis Layout
%====================================================================
if fcf,
    figure;
    FigureSet(1);
    end;
clf;

le = 0.10;                                                 % Left edge (plots along left edge)
pw = 0.99-le;                                              % Plot width (top 2 plots)
ph = 0.25;                                                 % Plot height

ax1 = axes('Position',[le   0.74 pw   ph]);
ax2 = axes('Position',[le   0.41 pw   ph]);
ax3 = axes('Position',[le   0.08 0.25 ph]);
ax4 = axes('Position',[0.42 0.08 0.25 ph]);
ax5 = axes('Position',[0.74 0.08 0.25 ph]);
         
%====================================================================
% Time-Series Plot
%====================================================================
axes(ax1);
ss = max(1,ceil(nx/100e3));                               % Step size - don't plot more than 100k points
kd = 1:ss:nx;                                             % Index of samples to plot  
td = (kd-1)/fs;                                           % Times of points plotted
k  = 1:nx;
t  = (k-1)/fs;
h  = plot(td,x(kd),'b');
set(h(1),'Color',0.5*[1 1 1]);
hold on;
    h = plot(t(ni),x(ni),'b.');
    set(h,'Color','r');
    h = plot(t(di),x(di),'r.');
    set(h,'Color','g');
    h = plot([0 (nx-1)/fs],th*[1 1],'b-');
    set(h,'LineWidth',2);    
    hold off;
xlim([0 (nx-1)/fs]);    
xlabel('Time (s)');
ylabel('Signal');

%====================================================================
% Peak Histogram
%====================================================================
axes(ax2);
xmin = min(x([pi]));
xmax = max(x([pi]));
bc   = linspace(xmin,xmax,nbph)';                          % Bin centers
bw   = min(diff(sort(unique(bc))));                        % Bin width
bhn  = hist(x(ni),bc);                                     % Bin heights for non-peaks
bhp  = hist(x(di),bc);                                     % Bin heights for peaks
bhn  = bhn(:);                                             % Convert to column vector
bhp  = bhp(:);                                             % Convert to column vector
ymax = max(1,min([3*max(bhp),1.02*max([bhn;bhp])]));
h    = bar(bc,bhn,1.0);
set(h,'FaceColor',0.7*[1 0 0]);
set(h,'EdgeColor',0.7*[1 0 0]);
hold on;
    h = bar(bc,bhp,1.0);
    set(h,'FaceColor',0.7*[0 1 0]);
    set(h,'EdgeColor',0.7*[0 1 0]);
    h = plot(th*[1 1],[0 ymax],'b-');
    set(h,'LineWidth',2);
    hold off;
box off;
xlabel('Signal');
ylabel('Frequency');
xlim([xmin-bw xmax+bw]);
ylim([0 ymax]);
box off;
AxisSet;

%====================================================================
% Overlap Plot
%====================================================================
axes(ax3);
idi = ceil(ed*fs/2);                                       % Duration between events
id = (-idi:idi)';
t   = id/fs;
hold on;
    if length(di)>1,
        I   = id*ones(1,length(di)) + ones(length(id),1)*di';
        I(I<=0) = 1;
        I(I>nx) = nx;
        h   = plot(t,x(I),'k');             
        set(h,'LineWidth',0.1);
        set(h,'Color',0.5*[1 1 1]);
        end;
    [jnk,is] = sort(x(di));
    nid = length(is);
    dic = di(is(1:min(10,nid)));                               % Indices of closest spikes

    [jnk,is] = sort(-x(ni));
    nid = length(is);
    nic = ni(is(1:min(10,nid)));                               % Indices of closest spikes

    if length(dic)>1,
        Id  = id*ones(1,length(dic)) + ones(length(id),1)*dic';
        Id(Id<=0) = 1;
        Id(Id>nx) = nx;
        h   = plot(t,x(Id),'g');             
        set(h,'LineWidth',0.1);
        end;
    
    if length(nic)>1,
        In  = id*ones(1,length(nic)) + ones(length(id),1)*nic';
        In(In<=0) = 1;
        In(In>nx) = nx;
        h   = plot(t,x(In),'r');             
        set(h,'LineWidth',0.1);
        end;
    
    h = plot([min(t) max(t)],th*[1 1],'b-');
    set(h,'LineWidth',2);
    
    hold off;
xlim([min(t) max(t)]);
ylim([min(x) max(x)]);
box off;
xlabel('Time (s)');
ylabel('Signal');
AxisSet;

%====================================================================
% Overlap Plot
%====================================================================
axes(ax4);
if length(di)~=0,
    idi = ceil(ed*fs/2);                                       % Duration between events
    did = (-idi:idi)';
    t   = did/fs;
    I   = did*ones(1,length(di)) + ones(length(did),1)*di';
    I(I<=0) = 1;
    I(I>nx) = nx;
    xmin = min(x);
    xmax = max(x);
    bc   = linspace(xmin,xmax,nbop)';                          % Bin centers
    bw   = min(diff(sort(unique(bc))));                        % Bin width
    N    = hist(x(I)',bc);
    nmax = max(max(N));
    cm   = jet(nmax);
    cm   = [1 1 1;cm];
    colormap(cm);
    imagesc(t,bc,N);
    set(gca,'YDir','Normal');
    hold on;
        h = plot([min(t) max(t)],th*[1 1],'w-');
        set(h,'LineWidth',3);
        h = plot([min(t) max(t)],th*[1 1],'b-');
        set(h,'LineWidth',2);    
        hold off;
    xlim([min(t) max(t)]);
    ylim([min(bc) max(bc)]);
    box off;
    xlabel('Time (s)');
    AxisSet;
    end;
    
%====================================================================
% Closest Plot
%====================================================================
axes(ax5);
idi = ceil(ed*fs/2);                                       % Duration between events
id = (-idi:idi)';
t   = id/fs;
hold on;
    if length(di)>1,
        I   = id*ones(1,length(di)) + ones(length(id),1)*di';
        I(I<=0) = 1;
        I(I>nx) = nx;
        h   = plot(t,x(I),'k');             
        set(h,'LineWidth',0.1);
        set(h,'Color',0.5*[1 1 1]);
        end;
    [jnk,is] = sort(x(di));
    nid = length(is);
    dic = di(is(1:min(10,nid)));                               % Indices of closest spikes

    [jnk,is] = sort(-x(ni));
    nid = length(is);
    nic = ni(is(1:min(10,nid)));                               % Indices of closest spikes

    if length(dic)>1,
        Id  = id*ones(1,length(dic)) + ones(length(id),1)*dic';
        Id(Id<=0) = 1;
        Id(Id>nx) = nx;
        h   = plot(t,x(Id),'g');             
        set(h,'LineWidth',0.1);
        end;
    
    if length(nic)>1,
        In  = id*ones(1,length(nic)) + ones(length(id),1)*nic';
        In(In<=0) = 1;
        In(In>nx) = nx;
        h   = plot(t,x(In),'r');             
        set(h,'LineWidth',0.1);
        end;
    
    h = plot([min(t) max(t)],th*[1 1],'b-');
    set(h,'LineWidth',2);
    
    hold off;
xlim([min(t) max(t)]);
yrg = max(x([dic;nic])) - min(x([dic;nic]));
ylim([th-1*yrg th+1*yrg]);
%ylim([min(x) max(x)]);
box off;
xlabel('Time (s)');
ylabel('Signal');
AxisSet;

%====================================================================
% Post-Processing
%====================================================================
if nargout==0,
    clear('di');
    end;
