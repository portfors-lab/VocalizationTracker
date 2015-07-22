function [di,th,te,mi] = SpikeDetectOptimalThreshold(x,fs,ti,aia,pfa);
% SpikeDetectOptimalThreshold: Detects spikes based on the threshold
%
%   [di,th,te,mi] = SpikeDetectOptimalThreshold(x,fs,ti,ai,pf);
%
%   x    Input signal.
%   fs   Sample rate (Hz). 
%   ti   True indices of events.
%   ai   Acceptance interval. Default = 0.001 s.
%   pf   Plot flag: 0=none (default), 1=screen, 2=current figure.
%
%   di   Indices of detected events.
%   th   Optimal threshold.
%   te   Total error of optimal threshold.
%   mi   Indices of maxima, sorted in decreasing order.
%
%   Uses the "true" indices of events to identify the optimal 
%   threshold that could be used for threshold detection. This serves
%   as an upper bound on the best performance that could be achieved
%   with this type of detector. This is especially commonly used for
%   detecting spikes in extracellular neural signals.
%
%   The loop through all possible thresholds stops when the total 
%   error exceeds 10 times the minimum error seen so far.
%
%   The total error is only accounted for with an approximate 
%   acceptance interval. When the detected indices are finally
%   calculated, the first peak above the "optimal" threshold is 
%   selected and all subsequent peaks above the threshold that are
%   within ai of these peaks are ignored. 
%
%   Example: Detect spikes in a microelectrode recording.
%
%      load MER.mat;
%      SpikeDetectOptimalThreshold(x,fs,si);
%
%   Version 1.10 JM
%
%   See also DetectMaxima, PowerPeaks, PickThreshold,
%   DetectionErrorIndices, SpikeDetectTemplate, DetectionPlots.


%=========================================================================
% Process function arguments
%=========================================================================
if nargin<3,
    help SpikeDetectOptimalThreshold;
    return;
    end;
  
ai = 1e-3;                                                 % Default acceptance interval of 1 ms
if exist('aia','var') && ~isempty(aia),
    ai = aia;
    end;         
    
pf = 0;                                                    % Default - no plotting
if nargout==0,                                             % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;

%=========================================================================
% Preprocessing
%=========================================================================
x   = x(:);
ti  = ti(:);
nt  = length(ti);                                          % No. true events

%=========================================================================
% Main Loop
%=========================================================================
[mi] = DetectMaxima(x,[],1);                               % Detect all local maxima as candidates of spike peaks
mx   = x(mi);
[mx,si] = sort(mx,1,'descend');                            % Sort maxima in decreasing order
mi   = mi(si);                                             % Sort indices in the same manner
nmx  = length(mx);

tp     = zeros(nt,1);                                      % Initialize vector with total number of events        
fp     = 0;
tea    = zeros(nmx+1,1);                                   % Total error array
tha    = zeros(nmx+1,1);                                   % Threshold array
tmin   = inf;
tea(1) = nt;                                              
tha(1) = mx(1)+0.01*std(x);
for c1 = 2:nmx,
    tha(c1) = mean(mx(c1-1:c1));
    [ds,id] = min(abs(mi(c1-1)-ti));                       % Index distance between one spike candidate and all true spikes    
    if ds<=ceil(ai*fs/2),
        tp(id) = 1;            
    else 
        fp = fp + 1;            
        end;    
    fn = nt-sum(tp);
    tea(c1) = (fn+fp)/nt;
    if tea(c1)<tmin,
        tmin = tea(c1);
        end;
    if tea(c1)>10*tmin,
        tha = tha(1:c1);
        tea = tea(1:c1);
        break;
        end;
    end;

%=========================================================================
% Postprocessing
%=========================================================================
tea = tea*100;                                             % Convert to percentage
[te,imin] = min(tea);
th = tha(imin);
di = sort(mi(1:imin-1));

nd = length(di);
ei = zeros(nd,1);                                          % Indices to eliminate
id = 1;
while id<=nd
    in = id + 1;
    while in<=nd & di(in)-di(id)<ai*fs,
        ei(in) = 1;
        in = in + 1;
        end;
    id = in;
    end;
di = di(find(~ei));       

%=========================================================================
% Display The Results
%=========================================================================
if pf,
    if pf==1,
        figure;
        FigureSet(1);
        end;
    subplot(2,1,1);    
    nx = length(x);
    k = 1:nx;    
    t = (k-0.5)/fs;
    h1 = plot(t,x);
    set(h1,'LineWidth',0.2);
    set(h1,'Color',0.8*[1 1 1]);
    hold on;
        h2 = plot(t(mi),x(mi),'r.');
        set(h2,'MarkerSize',5);
        
        h3 = plot(t(di),x(di),'g.');
        set(h3,'MarkerSize',6);
        
        h4 = plot([t(1) t(end)],th*[1 1],'k--');
        set(h4,'LineWidth',1.5);
        hold off;
    xlim([t(1) t(end)]);
    xlabel('Time (s)');
    ylabel('Signal (?)');
    box off;
    AxisSet;
    legend([h1;h3;h4],'Signal','Detected Events','Optimal Threshold','Location','SouthEast');
    
    subplot(2,1,2);
    h = plot(tha,tea);
    set(h,'LineWidth',1.5);
    hold on;
        h = plot(th,te,'g.');
        set(h,'MarkerSize',8);
        hold off;
    box off;
    xlim([ tha(end) tha(1)]);
    ylim([0 100]);
    ylabel('Total Error (%)');
    xlabel('Threshold');
    AxisSet;
    end;
    
%=========================================================================
% Process the Output Arguments
%=========================================================================
if nargout==0,
    clear('di','th','te','mi');
    end;
    