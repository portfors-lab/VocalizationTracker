function [di,mi]=SpikeDetectTemplate(x,fs,frla,kwa,tda,mfa,pfa);
% SpikeDetectTemplate: Automatic template matching spike detector
%
%   [di] = SpikeDetectTemplate(x,fs,frl,sm,kw,mf,pf);
%
%   x    Input signal.
%   fs   Sample rate (Hz).
%   frl  Firing rate limits. Default = [5 250] Hz.
%   kw   Kernel width. Default = 0.7 (x IQR). 
%   td   Template durations. Default = [0.45 0.55]*1e-3 s. 
%   mf   Manual flag: 0=automatic (default), 1=manual thresholds.
%   pf   Plot flag: 0=none (default), 1=screen, 2=current figure.
%
%   id   Indices of detected spikes
%   mi   Cell array of the sorted similarity maxima indices detected 
%        during each iteration
%
%   Uses automatic template matching to detect spikes in neural
%   signals. Primarily designed for extracellular neural signals such
%   as microelectrode recordings.
%
%   Example: Detect spikes in a MER signal.
%      load MER.mat;
%      SpikeDetectTemplate(x,fs,[5 250]);
%
%   Version 1.10 JM
%
%   See also DetectMaxima, PowerPeaks, Threshold,
%   DetectionErrorIndices, SpikeDetectThreshold.


%==========================================================
% Process function arguments
%==========================================================
if nargin<2,
    help SpikeDetectTemplate;
    return;
    end;

frmax  = 250;                                              % Maximum firing rate (Hz)
frmin  = 5;                                               % Minimum firing rate (Hz)
if exist('frla') & ~isempty(frla),
    frmin = frla(1);
    frmax = frla(2);
    end
    
kwp    = 0.7;                                              % Kernel width: 0.6 of IQR
if exist('kwa') & ~isempty(kwa),
    kwp = kwa;
    end;
    
tdur = [0.45 0.55]*1e-3;                                   % Template duration before & after peaks, respectively (seconds)
if exist('tda') & ~isempty(tda),
    tdur = tda;
    end;
    
mf = 0;                                                    % Manual flag
if exist('mfa') & ~isempty(mfa),
    mf = mfa;
    end;
    
pf = 0;                                                    % Default - no plotting
if nargout==0,                                             % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;

%====================================================================
% Author-Specified Parameters
%====================================================================
sm     = 2;                                                % Cross-correlation similarity measure
va     = 0.50;                                             % Valley depth argument
ni     = 3;                                                % Number of iterations
power  = 2;                                                % Signal power
detect = 1;                                                % Type of detection 0=all spikes, 1=most prominent 

%====================================================================
% Preprocessing
%====================================================================
x  = x(:);                                                 % Convert signal to column vector, if not already
x  = x-mean(x);                                            % Make the signal's mean zero
x  = x/std(x);                                             % Normalized input signal  
nx = length(x);
t  = [0:nx-1]/fs;                                          % time scale
fcp   = frmax;                                             % Cutoff frequency for lowpass filter applied to power signal
nsmin = ceil(frmin*(nx-1)/fs);                             % Minimum number of spikes
nsmax = ceil(frmax*(nx-1)/fs);                             % Maximum number of spikes
tdr   = [floor(tdur(1)*fs) ceil(tdur(2)*fs)];              % Convert template duration to units of samples
mi    = cell(ni+1,1);                                      % Cell array of sorted similarity maxima indices

%====================================================================
% Coarse Spike Detection
%====================================================================
pxid      = find(x>=0);
px        = zeros(size(x));
px(pxid)  = x(pxid);
[pi,xsp]  = PowerPeaks(px,fs,fcp,power,[],0);              % Locate peaks in smoothed signal power
[xsps,id] = sort(-xsp(pi));                                % Get IDs of peaks in descending order

pis  = pi(id);                                             % Convert to actual peak indices in descending order
np   = length(pi);                                         % Number of peaks detected 
i1   = min(np,nsmax);                                      % Maximum index of peak array to use 
ipi  = pis(1:nsmin);                                       % Select the largest 'nsmin' number of spikes to create the first spike template
mi{1} = pis;

%====================================================================
% Algorithm 3~5: Matched Filter Iterations
%====================================================================
for mit=1:ni,     
    [pis,fom,tmp] = TemplateDetect(x,fs,ipi,tdr,sm,pf);     % Use template matching to detect peak indices (sorted) : fom is the crosscorrelation      
    np  = length(pis);                                     % Number of peaks detected 
    
    mi{mit+1} = pis;
    
    tpf = 0; % Threshold plot flag
    if pf & mit==ni,
        tpf = 1;
        end;
    
    if ~mf,
        [th,thmax]= ThresholdPickAuto(fom(pis),nsmin,nsmax,sm,kwp,va,detect,tpf); % Pick up an optimal threshold automatically for the matched filter method
    else
        [th,thmax]= ThresholdPickAuto(fom(pis),nsmin,nsmax,sm,kwp,va,detect,1);  % Pick up an optimal threshold automatically for the matched filter method
        gi = ginput(1);
        th = gi(1);
        end;        
    if mit==ni & th==thmax,                                % If threshold equals thmax,
        ipi = [];                                          % no spikes declared
    else
        if sm==2,
            ipi = pis(find(fom(pis)>th));                      % The finalized indices of spikes' peaks
        else
            ipi = pis(find(fom(pis)<th));                      % The finalized indices of spikes' peaks
            end;
        ipi = ipi(:);
        end    
    
    if pf==1 | pf==2,
        k   = -tdr(1):tdr(2);                                      % Row vector of offset indices
        ntp = length(k);                                           % Number of template points
        nip = length(ipi);                                         % Number of initial points
        SI  = ipi*ones(1,ntp) + ones(nip,1)*k;                     % Matrix of signal indices : one row for one spike
        vri = find(all(SI>0 & SI<=nx,2));                          % Identify the valid row indices (those than don't exceed the signal limits)
        SI  = SI(vri,:);                                           % Strip out invalid rows
        tmpn = median(x(SI));                                       % Estimated template vector(row)    

        figure;
        FigureSet;
        h = plot(k/fs,x(SI),'k');
        set(h,'Color',0.7*[1 1 1]);
        hold on;
            h2 = plot(k/fs,tmp,'g',k/fs,tmpn,'b');
            set(h,'LineWidth',3);
            set(h,'Marker','.');
            set(h,'MarkerSize',15);
            hold off;
        xlim([k(1) k(end)]/fs);
        box off;
        xlabel('Time (s)');
        title(sprintf('Iteration %d',mit));    
        AxisSet;
        legend(h2,'Prior-Detection Template','Post-Detection Template');
        end;
    end
  
di = ipi;                                                  % Final outputs of this function

%==========================================================================
% Display the detected spike peaks by the iterated matched filter method
%==========================================================================
if (pf==1 || pf==2) & ~isempty(di),
    figure;    
    FigureSet(2);
    h = plot(t,x,t(di),x(di),'r.');
    box off;
    xlabel('Time (s)');
    ylabel('Microelectrode Amplitude (mV)');
    AxisSet;
    
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
elseif (pf==1 || pf==2) & isempty(di),                     % If there is no spike based on the estimated PDF decsion logic
    figure;                                                % Just plot the original MER
    FigureSet(2);
    h = plot(t,x);
    box off;
    xlabel('Time (s)');
    ylabel('Microelectrode Amplitude (mV)');
    legend('No spikes found');
    AxisSet(12);
    end
    
%=========================================================================
% Process the Output Arguments
%=========================================================================
if nargout==0,
    clear('di','mi');
    end;    
