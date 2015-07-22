function [y,ni] = ImpulseReject(x,xf,fsa,tha,pfa);
%ImpulseReject: Rejects impulses in the signal
%
%   [y] = ImpulseReject(x,xf,fs,th,pf);
%
%   x    Unfiltered input signal.
%   fs   Sample rate (Hz). Default = 1 Hz.
%   xf   Filtered signal.
%   th   Threshold for impulse rejection. Default = 4.
%   df   Divide by zero flag. 0 = filter (default), 1 = don't filter.
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   y    Output signal.
%   ni   Number of detected impulses.
%
%   Uses normalized difference,
%   
%   d = |x-xf|/s
%
%   where s is a robust estimate of the standard deviation given by
%
%   s = 1.483 median(|x-median(x)|)
%   
%   to determine the location of outliers. Specifically, at a given
%   sample index n
%
%   x(n) = {x(n)    if d(n)<th
%          {xf(n)   Otherwise
%
%   Example: Apply filter to the interbeat interval series for the
%   noisy ECG signal. 
%
%      load ecg.mat;
%      x = diff(ei);
%      xf = MedianFilter(x,10);
%      ImpulseReject(x,xf);
%
%   Astola, J. and Kuosmanen, P., "Fundamentals of Nonlinear Digital 
%   Filtering," CRC Press, pp.194-201, 1997.
%
%   Version 1.00 JM
%
%   See also MedianFilter.

%====================================================================
% Error Checking
%====================================================================    
if nargin<1,
    help ImpulseReject;
    return;
    end;
    
if length(x)==0 | length(xf)==0,
    error('Signal is empty.');
    end;

if length(x)~=length(xf),
    error('Input signals are different lengths.');
    end;    
    
if var(x)==0,
    error('Signal is constant.');
    end;
    
%====================================================================
% Calculate Basic Signal Statistics
%====================================================================    
nx  = length(x);
xr  = x;

%====================================================================
% Process Function Arguments
%====================================================================
fs = 1;
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
    end;    
        
th = 4;
if exist('tha') & ~isempty(tha),
    th = tha;
    end;  
    
pf = 0;                                 % Default - no plotting
if nargout==0,    % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
    
%====================================================================
% Calculate Output
%====================================================================
mx    = median(x);
mad   = median(abs(x-mx));
if length(th)==1,
	s     = 1.483*mad;
	if s==0,
        y = x; 
	else
		y     = x;
		sr    = abs(x-xf)/s; % Scaled Residual
		id    = find(sr>th);
        ni    = length(id);
		y(id) = xf(id);
        end;
    else
        y   = x;
        id  = find(abs(x-xf)>th);
        y(id) = xf(id);
        end;

%=========================================================================
% Plotting
%=========================================================================
if pf==1
    figure;
    FigureSet(1)
    subplot(2,1,1);
        k = 1:nx;
        t = k/fs;
        h = plot(t,x,'g',t,xf,'r',t,y,'b');
        set(h(1),'LineWidth',3);
        set(h(2:3),'LineWidth',1.2);
        xlabel('Time (s)');
        ylabel('Unscaled');
        title('Signals');
        xlim([0 (nx-1)/fs]);
        yrg = max(y)-min(y);
        ylim([min(y)-0.05*yrg max(y)+0.05*yrg]);
        box off; 
        AxisSet;
        legend(h,'Unfiltered','Filtered','Impulses Removed',2);
    subplot(2,1,2);
        k = 1:nx;
        t = k/fs;
        h = plot(t,sr,'r',t,abs(x-y)/s,'b',[0 (nx-1)/fs],th*[1 1],'k:');
        xlabel('Time (s)');
        ylabel('Unscaled');
        title('Absolute Residuals');
        xlim([0 (nx-1)/fs]);
        ylim([0 1.1*th]);
        box off; 
        AxisSet;
        legend('Filtered','Impulses Removed','Threshold',2);
    end;

%=========================================================================
% Take care of outputs
%=========================================================================
if nargout==0,
    clear('y','ni');
    end;   