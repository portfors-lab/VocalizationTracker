function [te,tp,fp,fn] = DetectionErrorIndices(fs,di,ti,aia);
%DetectionErrorIndices: Replaced by DetectionMetrics.

fprintf('This function has been replaced by DetectionMetrics.\n');
return;

%DetectionErrorIndices: Calculates error indices for detectors
%
%   [te,tp,fp,fn] = DetectionErrorIndices(fs,di,ti,ai);
%
%   fs    Sampling frequency (Hz)
%   did   Indices of detected spike peaks
%   tid   Indices of true spike peaks
%   ai    Acceptance interval (ms). Default = 1 ms.
%
%   te    Total error (0-100%)
%   tp    No. true positives
%   fp    No. false positives
%   fn    No. false negatives
%
%   Correctly calculates five detection error indices. If multiple 
%   events are detected within the acceptance interval, only counts
%   as one true positive and no false positives.
%
%   The total error is defined as te = 100*(fp+fn)/length(ti).
%
%   Example: Calculate the error statistics on a spike detection
%   algorithm.
%
%      load MER.mat;
%      di = SpikeDetectTemplate(x,fs);
%      [te,tp,fp,fn] = DetectionErrorIndices(fs,di,si)
%
%   Version 1.10 JM

%=========================================================================
% Process function arguments
%=========================================================================
if nargin<3,
    help DetectionErrorIndices;
    return;
    end;

ai = 1e-3;                                                 % Default acceptance interval, 1 s
if exist('aia') & ~isempty(aia),
    ai = aia;
    end;

%=========================================================================
% Preprocessing
%=========================================================================
di  = di(:);                                               % Make di into a column vector
nd  = length(di);                                          % No. detected events

ti  = ti(:);                                               % Make ti into a column vector
nt  = length(ti);                                          % No. events

fp  = 0;                                                   % Initialize count of false positives
tp  = zeros(nt,1);                                         % Initialize vector with total number of events        

if nt==0,
    tp = 0;
    fp = nd;
    fn = 0;
    te = 100;
    end;

%=========================================================================
% Main Loop
%========================================================================= 
for c1=1:nd,
    [ds,id] = min(abs(di(c1)-ti));                         % Index distance between one spike candidate and all true spikes
    if ds <= ceil(ai*fs/2),
        tp(id) = 1;            
    else 
        fp = fp + 1;            
        end;
    end;
tp = sum(tp);                                          % No. of true positives
fp = fp;                                               % No. of false positives
fn = nt-tp;                                            % No. of false negatives
te = 100*(fp+fn)/nt;                                   % Total Error (%)
