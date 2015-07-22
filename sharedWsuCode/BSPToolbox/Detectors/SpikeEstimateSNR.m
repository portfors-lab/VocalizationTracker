function [sh,sp,np] = SpikeEstimateSNR(x,fs,si,xra,pfa);
% SpikeEstimateSNR: Estimates the signal-to-noise ratio of a MER.
%
%   [sh,np] = SpikeEstimateSNR(x,fs,si,xr,pf);
%
%   x    Input signal.
%   fs   Sample rate (Hz). 
%   si   Spike indices.
%   xr   Exclusion range [before after] (s). Default = [0.001 0.001].
%   pf   Plot flag: 0=none (default), 1=screen, 2=current figure.
%
%   sh   Estimated Signal to Noise ratio
%   sp   Estimated signal power (variance).
%   np   Estimated noise power (variance).
%
%   The signal x is assumed to contain WSS background noise and 
%   discrete, isolated events indicated by vector of the spike indices 
%   si. We assume segments that contain spikes include both the signal
%   and the additive background noise and that these are uncorrelated.
%   The background noise power is estimated from all of the portions
%   of the signal that are outside the exclusion range around each
%   spike.
%   
%   Example: Estimate the SNR of a MER signal.
%
%      load MER.mat;
%      [sh,sp,np] = SpikeEstimateSNR(x,fs,si);
%
%   Version 1.10 JM   

%=========================================================================
% Process function arguments
%=========================================================================
if nargin<3,
    help SpikeEstimateSNR;
    return;
    end;

xr = [0.001 0.001];                                        % Default exclusion range (s)
if exist('xra') & ~isempty(xra),
    xr = xra;
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
x = x(:) - mean(x);
nx = length(x);
ns = length(si);

%=========================================================================
% Estimate SNR
%=========================================================================
iv  = zeros(nx,1);                                         % Indicator variables for exclusion range
bid = si - ceil(fs*xr(1));                                 % Indices 0.5 ms prior to spikes
aid = si + ceil(fs*xr(2));                                 % Indices 2 ms posterior to spikes

bid(bid<1)  = 1;                                           % Handle edge conditions
aid(aid>nx) = nx;

for c1=1:ns,
    iv(bid(c1):aid(c1)) = 1;
    end;
tp = mean(x(find( iv)).^2);                                % Estimated total power
np = mean(x(find(~iv)).^2);                                % Estimated noise power
sp = tp - np;                                              % Estimated signal power
sh = sp/np;                                                % Estimated SNR

%====================================================================
% Plot Results
%====================================================================  
if pf,
    if pf==1,
        figure;
        FigureSet(1);
        end;
    k = 1:nx;
    t = (k-0.5)/fs;
    xe = x;
    xe(find(~iv)) = nan;
    h = plot(t,xe);
    set(h,'LineWidth',0.2);
    set(h,'Color','g');
    hold on;
        xn = x;
        xn(find(iv)) = nan;
        h = plot(t,xn,'b');
        set(h,'LineWidth',0.2);
        set(h,'Color','b');
        
        h = plot(t(si),x(si),'r.');
        set(h,'MarkerSize',8);
        hold off;
    xlabel('Time (s)');
    ylabel('Signal (?)');
    box off;
    AxisSet;
    legend('Exclusion Regions','Noise Regions','Detected Spikes');
    end;
