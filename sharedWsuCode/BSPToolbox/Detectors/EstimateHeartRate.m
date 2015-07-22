
function [HR, ind] = EstimateHeartRate(x, fsa, wa, pfa)
%EstimateHeartRate Estimate the heart rate (hr)
%
%   [hr] = EstimateHR(x, fs, w, pf)
%
%   x    Input signal       
%   fs   Signal sample rate (Hz), default=125 Hz      
%   w    Window length (s), default = 10 s, w must be greater than 10 
%   pf   Plot format. 0:none (default), 1:screen
%
%   hr   heart rate estimate 
%
%   Estimates the heart rate from an input time seres. The input time series is 
%   expected to be a pressure waveform such as ABP or ICP. The output hr is column
%   vector containing the heart rate time series. The function estimates the heart 
%   rate within each window using the PSD.
%
%   Example: Estimate the HR of an intracranial pressure signal 
%   sampled at 125 Hz and plot the results:
%      load ICP; 
%      [hr] = EstimateHeartRate(icp, fs, 10, 1);
%
%   Version  1.0
%
%   See also EstimatePSD, Blackman, Periodogram


%--------------------------------------------------------------------
% Process function arguments
%--------------------------------------------------------------------

if nargin<1 | nargin>4,
    help EstimateHeartRate;
    return;
    end;

fs = 125;                                 % Default sampling rate, Hz
if exist('fsa') & ~isempty(fsa), 
    fs = fsa;
    end;
 
w = 10*fs;                               % Default window
if exist('wa') & ~isempty(wa),
    w = wa*fs;
    end;
    
pf = 0;                                  % Default - no plotting
if nargout==0,                           % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;

x  = x(:);
n1 = 1;
n2 = length(x);
n3 = 0:n2-1;
n  = ceil(log2(w));
k1 = 1;

% --------------------------------------------------------------------
% Create Windows
% --------------------------------------------------------------------
for k=1:w:n2-w
    if ((k+w-1)<n2-w)
    [PSD, hertz, N] = EstimatePSD(x(k:k+w-1),fs,1,0);
    n4 = ceil(0.5*N*(1/fs));
    n5 = ceil(4.0*N*(1/fs));
    [M, I] = max(PSD(n4:n5));
    HR(k1) = (n4+I)*(1/(N*(1/fs)));
    index(k1) = k;   
    k1 = k1+1;
    else
    [PSD, hertz, N] = EstimatePSD(x(k:n2),fs,1,0);
    n4 = ceil(0.5*N*(1/fs));
    n5 = ceil(4.0*N*(1/fs));
    [M, I] = max(PSD(n4:n5));
    HR(k1) = (n4+I)*(1/(N*(1/fs)));
    index(k1) = k;      
        
    end
end

% ---------------------------------------------------------------------
% Plotting 
% ---------------------------------------------------------------------
if pf==1
    t1    = n3./(fs*60);
    t2    = (index+round(w/2))./(fs*60);
    xmin1 = min(t1);
    xmax1 = max(t1);
    ymin1 = min(x);
    ymax1 = max(x);
    xmin2 = min(t2);
    xmax2 = max(t2);
    ymin2 = min(HR);
    ymax2 = max(HR);


    subplot(2,1,1)
    plot(t1, x)
    title('Input Time Series');
    ylabel('Amplitude');
    axis([xmin1 xmax1 ymin1 ymax1]);
    grid on;

    subplot(2,1,2)
    plot(t2,HR)
    title('Heart Rate Estimation');
    xlabel('Time, minutes');
    ylabel('Heart Rate, Hz');
    axis([xmin2 xmax2 ymin2 ymax2]);
    grid on;
end;
