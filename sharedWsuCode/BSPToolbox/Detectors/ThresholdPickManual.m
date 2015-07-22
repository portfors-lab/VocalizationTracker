function [di,th,ef] = SelectThreshold(x,fsa,thia,eda);
%SelectThreshold: Select a threshold based on histograms and overlap.
%   
%   [di,th,ef] = SelectThreshold(x,fs,thi,ed);
%
%   x    Input signal.
%   fs   Sample rate (Hz). Default = 1 Hz.
%   thi  Initial threshold. Default = Detected beat indices (samples). 
%        Default = 4*RSD (RSD is a robust estimate of the
%        signal standard deviation).
%   ed   Event duration (s). Default = median inter-event interval.
%
%   di   Detected event indices.
%   th   Detection threshold.
%   ef   Edit flag. 1=manual annotation recommended, 0=annotation was 
%        good.
%   
%   A handly tool for manually selecting a threshold for threshold 
%   detection. Plots both the histogram of detected indices and an 
%   overlap plot of the detected events.
%
%   Example: Edit detected beats of a microelectrode recording. 
%
%      load MER.mat; 
%      di = SelectThreshold(x,fs);
%
%   Task Force of the European Society of Cardiology and the American 
%   Society of Pacing and Electrophysiology, "Heart rate variability:
%   Standards of measurement, physiological interpretation, and 
%   clinical use," Circulation, vol. 93, no. 5, pp. 1043-1065, Aug. 
%   1996.
%
%   Version  0.00.01.23 JM
%
%   See also DetectionPlots and EditAnnotations. 

%=====================================================================
% Process function arguments
%=====================================================================
if nargin<1,
    help SelectThreshold;
    return;
    end;

fs = 1;
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
    end;        
    
ed = 0;
if exist('eda') & ~isempty(eda),
    ed = eda;
    end;      
    
thi = 0;
if exist('thia') & ~isempty(thia),
    thi = thia;
else
    rsd = 1.483*median(abs(x-median(x)));
    thi = 4*rsd;
    end;
       
%====================================================================
% Preprocessing
%====================================================================
th  = thi;

figure;
FigureSet(1,'Wide');

%====================================================================
% Main Loop
%====================================================================
df = 0;
ef = 0;
while 1,      
    di = DetectionPlots(x,fs,th,ed,0);
    if ~df,
        fprintf('Select a new threshold on the histogram plot.\n');
        fprintf('Hit e on the keyboard if you want to edit the detected spikes.\n');
        %fprintf('If no spikes are in this signal, hit n on the keyboard.\n');
        fprintf('Hit right mouse button or x when done.\n');
        df = 1;
        end;
    drawnow;
    
    [ut, mag, ui] = ginput(1);
    
    if isempty(ui),
        fprintf('Empty!\n');
        continue;
        end; 

    switch ui,
    case {1,2}, % Left mouse button
        th = ut;
    case {3,'x','X'}
        break;
    case {'e','E'}
        ef = 1;        
    case {'n','N'} 
        di = [];
        break;
        end;
    end;

%====================================================================
% Post Processing
%====================================================================
if nargout==0,
    clear('di','th','ef');
    end;


