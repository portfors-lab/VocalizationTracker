function [brs,bei,npr] = BRSSequence(abp,fsa,pi,ri,fse,tha,pfa)
%SpectralHRV: Standard spectral heart rate variability metrics.
%
%   [brs,bei,npr] = BRSSequence(abp,fsa,pi,ri,fse,tha,pfa)
%
%   abp   Arterial blood pressure signal.
%   fsa   Arterial blood pressure sample rate (Hz).
%   pi    Indices of systolic peaks in pressure (samples).
%   ri    Indices of R peaks in electrocardiogram signal (samples).
%   fse   Sample rate of electrocardiogram signal.
%   th    Thresholds for ramps (mmHg,samples). Default = [1 1].
%   pf    Plot flag: 0=none (default), 1=screen.
%
%   brs   Baroreflex sensitivity estimated by sequence technique.
%   bei   Baroreflex effectiveness index.
%   npr   Number of blood pressure ramps.
%
%   This function calculates the baroreflex sensitivity using the
%   sequence technique.
%     
%   Example: To be completed.
%
%   Marco Di Rienzo, Gianfranco Parati, Paolo Castiglioni, Roberto 
%   Tordi, Giuseppe Mancia, and Antonio Pedotti, "Baroreflex 
%   effectiveness index: an additional measure of baroreflex control 
%   of heart rate in daily life" Am J Physiol Regulatory Integrative 
%   Comp Physiol, vol. 280: pp. R744-R751, 2001.
%
%   Version 0.00.00.23 JM
%
%   See also Detectors and Metrics.

%====================================================================
% Error Checking
%====================================================================    
if nargin<5,
    help BRSSequence;
    return;
    end;
    
%====================================================================
% Process function arguments
%====================================================================    
th = [1 1]; % Default sampling rate, Hz
if exist('tha') & ~isempty(tha) & length(tha)==2,
    th = tha
    end;

pf = 0;                                 % Default - no plotting
if nargout==0,    % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
    
%====================================================================
% Preprocessing
%====================================================================  
abp = abp(:); % Convert to column vector

ar = zeros(length(pi),1); % Allocate memory for pressure indices
si = zeros(length(pi),1); % Allocate memory for sequence indices
ss = zeros(length(pi),1); % Allocate memory for sequence slopes (BRS)

% Identify canidates
%pd = diff(abp(pi));
%si = find(abs(pd)>th(1)); % Only consider those that are above the threshold 

%====================================================================
% Calculate Metrics
%====================================================================  
ac = 0; % Count of ABP ramps
sc = 0; % Count of RRI ramps & sequences

for c1 = 1:length(pi)-2, % Consider each peak - is it the start of a ramp
    %ia = pi(si(c1));
    %ia = pi(c1);
    x  = abp(pi(c1 + (0:2)));
    
    if all(diff(x)>th(1)), % Criteria for ABP ramp
        ac     = ac + 1;
        ar(ac) = pi(c1);
        
        % Check to determine if there is a corresponding RRI ramp
        tp  = (pi(c1)-1)/fsa;        % Time of ABP peak
        eri = tp*fse+1;          % Expected r index (not an integer)       
        ie  = min(find(ri>eri)); % Index of closest r that occurs after ABP peak 
       
        if ie<=length(ri)-4, 
            rri = diff(ri(ie + (0:3))); % RR intervals            
            if all(diff(rri)>th(2)), % Criteria for RRI ramp
                sc     = sc + 1;
                si(sc) = pi(c1);
                y      = rri;
                ss(sc) = sum((x-mean(x)).*(y-mean(y)))/sum((x-mean(x)).^2);
                end;            
            end;
        end;
    if all(-diff(x)>th(1)), % Criteria for ABP ramp
        ac     = ac + 1;
        ar(ac) = pi(c1);
        
        % Check to determine if there is a corresponding RRI ramp
        tp  = (pi(c1)-1)/fsa;        % Time of ABP peak
        eri = tp*fse+1;          % Expected r index (not an integer)       
        ie  = min(find(ri>eri)); % Index of closest r that occurs after ABP peak 
       
        if ie<=length(ri)-4, 
            rri = diff(ri(ie + (0:3))); % RR intervals            
            if all(-diff(rri)>th(2)), % Criteria for RRI ramp
                sc     = sc + 1;
                si(sc) = pi(c1);
                y      = rri;
                ss(sc) = sum((x-mean(x)).*(y-mean(y)))/sum((x-mean(x)).^2);
                end;            
            end;
        end;
        
        
    end;
ar = ar(1:ac);
si = si(1:sc);
ss = ss(1:sc);
    
%====================================================================
% Postprocessing
%====================================================================    
brs = mean(abs(ss));
bei = sc/ac; 
npr = ac;

%====================================================================
% Plot 
%====================================================================
if pf==1,
    figure;
    FigureSet;
    k = 1:length(abp);
    t = (k-1)/fsa;
    rrt = (ri(1:length(ri)-1) + ri(2:length(ri)))/2;
    rrt = (rrt-1)/fse;
    rri = diff(ri)/fse;
    rri = 0.2*std(abp(pi))*(rri-mean(rri))/std(rri) + mean(abp(pi));
    
    rt = (ri-1)/fse;
    
    h = plot(t,abp,'r',rrt,rri,'g',t(si),abp(si),'b.',rt,mean(abp(pi))*ones(length(rt),1),'k.');
    set(h(2),'Marker','.');
    set(h(3),'MarkerSize',15)    
    xlim([0 max(t)]);
    xlabel('Time (sec)');
    ylabel('ABP (mmHg)');
    AxisSet;
    box off;
    end
    
%====================================================================
% Process Return Arguments
%====================================================================
if nargout==0,
    clear('brs');
    clear('bei');
    end;
