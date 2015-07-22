function [] = DetectionPlot(x,fs,di,pta,eda,fcf);
%DetectionPlots: Displays plots to inspect detection performance.
%   
%   DetectionPlot(x,fs,di,pt,ed,fcf);
%
%   x    Input signal.
%   fs   Sample rate (Hz). 
%   di   Detection indices (samples).
%   pt   Plot type. 1=Overlap plot (default), 2=Overlap histogram.
%   ed   Event duration (s). Default = median inter-event interval.
%   fcf  Flag to create new figure? 1 = Yes (default), 0 = No.
%
%   A handly tool for plotting the morphology of detected events.
%
%   Example: Edit detected beats of a noisy electrocardiogram. 
%
%      load MER.mat; 
%      DetectionPlot(x,fs,si);
%
%   Task Force of the European Society of Cardiology and the American 
%   Society of Pacing and Electrophysiology, "Heart rate variability:
%   Standards of measurement, physiological interpretation, and 
%   clinical use," Circulation, vol. 93, no. 5, pp. 1043-1065, Aug. 
%   1996.
%
%   Version  0.00.01.01 JM
%
%   See also ThresholdDetectionPlots and EditAnnotations. 

%====================================================================
% Error Checking
%====================================================================    
if nargin<1,
    help DetectionPlot;
    return;
    end;
    
if exist('eda') & length(eda)>1,
    error('ed argument must be a scalar.');
    end;

%=====================================================================
% Process function arguments
%=====================================================================
pt = 1;
if exist('pta') & ~isempty(pta),
    pt = pta;
    end;         

ed = 0;
if exist('eda') & ~isempty(eda),
    ed = eda;
else    
    if length(di)>1,
        ed = median(diff(sort(di))/fs);
    else
        ed = 50/fs; 
        end;
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
x   = x(:)-mean(x);                                        % Make into a column vector
nx  = length(x);

%====================================================================
% Figure Initialization & Axis Layout
%====================================================================
if fcf,
    figure;
    FigureSet(1);
    end;
clf;

switch pt,
    %----------------------------------------------------------------
    % Time Series Plot
    %----------------------------------------------------------------
    case 1,         
    mi = DetectMaxima(x,0,1);
    ss = max(1,ceil(nx/100e3));                               % Step size - don't plot more than 100k points
    kd = 1:ss:nx;                                             % Index of samples to plot  
    td = (kd-1)/fs;                                           % Times of points plotted
    k  = 1:nx;
    t  = (k-0.5)/fs;
    h  = plot(td,x(kd),'b');
    set(h(1),'Color',0.5*[1 1 1]);
    hold on;
        h = plot(t(mi),x(mi),'r.');
        set(h,'Color','r');
        set(h,'MarkerSize',2);
        h = plot(t(di),x(di),'g.');
        set(h,'Color','g');
        set(h,'MarkerSize',10);        
        hold off;
    xlim([0 nx/fs]);
    AxisSet;
    xlabel('Time (s)');
    ylabel('Signal (units?)');
    box off;
    
    %----------------------------------------------------------------
    % Overlap Histogram Image
    %----------------------------------------------------------------
    case 2,
    idi = ceil(ed*fs/2);                                       % Duration before and after each event
    did = (-idi:idi)';                                         % Segment indices (samples)
    t   = did/fs;                                              % Segment indices (seconds)
    I   = did*ones(1,length(di)) + ones(length(did),1)*di';    % Matrix of segment indices
    I(I<=0) = 1;                                               % Handle left edge condition
    I(I>nx) = nx;                                              % Handle right edge condition 
    xmin = min(x);                                             % Select lower range of bins
    xmax = max(x);                                             % Select upper range of binx 
    bc   = linspace(xmin,xmax,nbop)';                          % Bin centers
    bw   = min(diff(sort(unique(bc))));                        % Bin width
    N    = hist(x(I)',bc);                                     % Get histogram of values
    nmax = max(max(N));                                        % Maximum value
    N    = 256*N/nmax;
    cm   = jet(255);
    cm   = [1 1 1;cm];
    colormap(cm);
    imagesc(t,bc,N);
    set(gca,'YDir','Normal');
    hold on;
        xt = median(x(I'));                                % Template
        h = plot(t,xt,'k',t,xt,'w');
        set(h(1),'LineWidth',3.0);
        set(h(2),'LineWidth',1.5);        
        hold off;
    xlim([min(t) max(t)]);
    ylim([min(bc) max(bc)]);
    box off;
    xlabel('Time (s)');
    AxisSet;
    
    %----------------------------------------------------------------
    % Overlap Histogram Image
    %----------------------------------------------------------------
    case 3,    
    figure;
    FigureSet(1);
    tl = tdr(1)+tdr(2)+1;                                  % Template length (samples)
    ns = [-tdr(1)-round(0.25*tl):tdr(2)+round(0.25*tl)];   % Number of samples prior and posterior to a peak of a spike
    iv = find(ipi+min(ns)>0);                              % Valid (plottable) indices
    ti = ipi(iv)*ones(1,length(ns))+ones(length(iv),1)*ns;% Template indices matric
    h = plot(ns/fs,x(ti));
    set(h,'Color',0.8*[1 1 1]);                            % Light gray
    set(h,'LineWidth',0.1);                                % Very thin lines
    xlim([ns(1) ns(end)]/fs);
    ymx = max(tmp);
    ymn = min(tmp);
    yrg = ymx-ymn;
    ylim([ymn-0.05*yrg ymx+0.05*yrg]);   
    hold on;
        h=plot(xlim,0*ylim,'k:',0*xlim,ylim,'k:');
        h=plot([-tdr(1):tdr(2)]/fs,tmp,'g');
        set(h,'LineWidth',5);
        set(h,'Color','g');
        hold off;
    box off;
    xlabel('Time (s)');
    ylabel('Microelectrode Amplitude (mV)');
    AxisSet;    
    end;

%====================================================================
% Post-Processing
%====================================================================
if nargout==0,
    end;
