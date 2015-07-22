function [te,tp,fp,fn] = DetectionMetrics(ti,di,ai);
%DetectionMetrics: Calculates performance metrics for detection
%
%   [te,tp,fp,fn] = DetectionMetrics(ti,di,ai);
%
%   ti   True event indices (gold standard).
%   di   Detection indices.
%   ai   Acceptance interval (samples). 
%
%   tp   Number of true positives.
%   fp   Number of false positives.
%
%   Correctly calculates five detection error indices. If multiple 
%   events are detected within the acceptance interval, only counts
%   as one true positive and all others are counted as false 
%   positives.
%
%   The acceptance interval is symmetric about each true event. Thus
%   if the acceptance interval is 10 samples a detected index that is
%   5 samples before or 5 samples after the true index would be 
%   counted as a true positive. 
%
%   This is an O(nt) + O(nd) algorithm that does not do a find for
%   every index, where nt=length(ti) and nd=length(di). Thus, it is 
%   much more efficient than straightforward algorithms for 
%   calculating these indices.
%
%   The total error is defined as te = 100*(fp+fn)/length(ti).
%
%   Example: Compare a threshold detector to true spikes in a 
%   microelectrode recording signal with an acceptance interval
%   of 1 ms.
%
%      load MER.mat; 
%      im = DetectMaxima(x);
%      [jnk,is] = sort(x(im),'descend');
%      im = im(is(1:1000));
%      [te,tp,fp,fn] = DetectionMetrics(si,im,1e-3*fs) 
%
%   Hayes, M., "Statistical Digital Signal Processing and Modeling," 
%   John Wiley & Sons, 1996.
%
%   Version 1.00 JM
%
%   See also DetectMaxima, EditAnnotations, and Detectors.

%====================================================================
% Error Checking
%====================================================================    
if nargin<3,
    help DetectionMetrics;
    return;
    end;
    
%====================================================================
% Preprocessing
%====================================================================      
nt = length(ti);
nd = length(di);

if nd==0,
    tp = 0;
    fp = 0;
    fn = nt;
    te = 100*(fn+fp)/nt;
    return;
    end;

ti = sort(ti);
di = sort(di);

%====================================================================
% Calculate True Positives
%====================================================================  
i0 = 1;
i1 = 1;
tp = 0;
for c1=1:nt,
    while i0+1<nd & di(i0+1)<ti(c1),                       % Increment i0 until di(i0) is just to the left of ti(c1)
        i0 = i0 + 1;
        end;
    while di(i1)<ti(c1) & i1<nd,                           % Increment i1 until di(i1) is just to the right of ti(c1)
        i1 = i1 + 1;
        end;
    if min(abs(ti(c1)-di(i0:i1)))<=ai/2,
        tp = tp + 1;
        end;
    end;    

% %====================================================================
% % Calculate False Positives
% %====================================================================    
% i0 = 1;
% i1 = 1;    
% fp = 0;
% for c1=1:length(di),
%     while ti(i0+1)<di(c1) & i0+1<nt,                       % Increment i0 until di(i0) is just to the left of ti(c1)
%         i0 = i0 + 1;
%         end;
%     while ti(i1)<di(c1) & i1<nt,                           % Increment i1 until di(i1) is just to the right of ti(c1)
%         i1 = i1 + 1;
%         end;
%     if min(abs(di(c1)-ti(i0:i1)))>ai,
%         fp = fp + 1;
%         end;
%     end;  

%====================================================================
% Postprocessing
%==================================================================== 
fp = nd-tp;
fn = nt-tp;
te = 100*(fn+fp)/nt;
