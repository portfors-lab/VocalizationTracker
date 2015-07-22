function [m] = TimeDomainHRV(di,fs,pfa)
%TimeDomainHRV: Standard heart rate variability metrics.
%
%   [m] = TimeDomainHRV(di,fs,pf);
%
%   di   Detected indices of normal contractions (samples).
%   fs   Signal sample rate (Hz).
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   m    Vector of the four HRV metrics recommended by the task 
%        force: [SDNN; HRV-TI; SDANN; RMSDD; SDSD; MHR].
%
%   This function calculates the four heart rate variability (HRV)
%   metrics recommended by the task force:
%      
%      SDNN   Standard deviation (sec).
%      HRV-TI HRV Triangular index (unitless ratio).
%      SDANN  Std. dev. of the 5 minute average intervals (sec).
%      RMSDD  Square root of mean squared interval differences (sec).
%      SDSD   Standard deviation of suc. interval differences (sec).
%      MHR    Mean heart rate (Hz).
%  
%   The standard does not completely specify how the HRV-TI index
%   should be calculated. This implementation uses a bin width
%   of 1/fs, as suggested by the task force. The bins are aligned 
%   such that the center of the first bin occurs at 0 (this was not 
%   addressed by the task force). The task force recommended that 
%   this only be used with recordings 20 minutes
%   or longer.
%
%   Similarly, the standard does not completely specify how the
%   window of segments should be calculated for the SDANN metric.
%   This implementation uses all of the beats within 2.5 minutes
%   before and after each individual beat. The metric is then
%   calculated as the standard deviation of these averages. 
%   Hence, this metric is most appropriate for segments much
%   longer than 5 minutes (e.g. 24 hours).
%
%   See reference for details. Other metrics (not implemented) are
%   also described in the reference, but not implemented in the 
%   toolbox.
%
%   Example: Calculate the HRV metrics of the first five minutes
%   of an ECG.
%
%      load ECG;
%      x  = ecg(1:fs*5*60);
%      di = ECGDetectREnergy(x,fs);
%      m  = TimeDomainHRV(di,fs,1);
%
%   Task Force of the European Society of Cardiology and the American 
%   Society of Pacing and Electrophysiology, "Heart rate variability:
%   Standards of measurement, physiological interpretation, and 
%   clinical use," Circulation, vol. 93, no. 5, pp. 1043-1065, Aug. 
%   1996.
%
%   Version 0.00.00.23 JM
%
%   See also Detectors.

%====================================================================
% Error Checking
%====================================================================    
if nargin<2,
    help TimeDomainHRV;
    return;
    end;

nb = length(di); % No. beats
if nb==0,
    error('Signal is empty.\n');
    end;
    
%====================================================================
% Process function arguments
%====================================================================       
pf = 0;                                                    % Default - no plotting
if nargout==0,                       % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
    
%====================================================================
% Author-Specified Parameters
%====================================================================    
wl = 5*60*fs;                                              % Surrounding window length (samples) for SDANN

%====================================================================
% Preprocessing
%==================================================================== 
nnis = diff(di);                                           % NN intervals (samples)
nni  = nnis/fs;                                            % NN intervals (seconds)

%====================================================================
% Calculate Metrics
%====================================================================    
sdnn = std(nni);                                           % Standard deviation of the NN intervals

a  = zeros(nb,1);
ib = 1;                                                    % Beginning segment index
ie = 1;                                                    % Ending segment index
for c1 = 1:nb,
    while         di(c1  )-di(ib+1)>wl/2,                  % Ensure lower index is no more than 2.5 min from index
        ib = ib + 1;
        end;
    while ie<nb & di(ie+1)-di(c1  )<wl/2,                  % Ensure upper index is no more than 2.5 min from center indexm but as far out as possible
        ie = ie + 1;
        end;
    a(c1) = mean(diff(di(ib:ie))/fs);
    end;
sdann = std(a);

rmssd = sqrt(mean((diff(nni).^2)));                        % Square root of the mean squared differences of successive NN intervals

sdsd = std(diff(nni));                                     % Standard deviation of successive interval differences

bc = (min(nnis):max(nnis))/fs;                             % Bin centers
ic = hist(nni,bc);                                         % Count of intervals in histogram
ti = sum(ic)/max(ic);                                      % Traingular index

mhr = mean(1./nni);

%====================================================================
% Postprocessing
%====================================================================    
m = [sdnn;ti;sdann;rmssd;sdsd;mhr];

%====================================================================
% Plot PSD
%====================================================================
if pf==1,
    figure;
    FigureSet;
    y = 100*bc/sum(bc);
    h = bar(be,y,'histc');
    set(h,'FaceColor',0.5*[1 1 1]);
    title('Histogram of NN Intervals');
    xlabel('Time (sec)');
    ylabel('Percent');
    xlim([min(be) max(be)]);
    ylim([0 1.02*max(y)]);
    box off;
    zoom on;
    AxisSet;
    end
    
%====================================================================
% Process Return Arguments
%====================================================================
if nargout==0,
    clear('m');
    end;
