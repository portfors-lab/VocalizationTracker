function [y] = SmoothSeries(t,x,ti,wla,pfa);
%SmoothSeries: Kernel smoothing for non-uniform sampled signals.
%
%   [y] = SmoothSeries(t,x,ti,wl,pf);
%
%   t    Times of signal observations (sec).
%   x    Values of signal observations.
%   ti   Times to generate estimate of smoothed signal values (sec).
%   wl   Length of kernel window to use (sec). Default = 5 sec.
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   y    Estimated filtered signal values at times specified by ti.
%
%   Smoothes (lowpass filters) the signal specified by the (t,y) 
%   coordinates and evaluates the estimate at the times specified by
%   the vector ti. This function is faster than other kernel 
%   smoothing routines because it takes advantage of the fact that t 
%   and ti are sorted in increasing order (required). 
%
%   This implementation uses a truncated guassian kernel with 
%   standard deviation specified by wl. Points more than 5 standard 
%   deviations away from the evaluation time, ti, are ignored. 
%   Although the units specified above are seconds, the signal times 
%   (t), estimation times (ti), and window length (wl) can actually 
%   be in any measure (e.g. samples), as long as they are consistent 
%   with one another.
%
%   Example: Do a smooth interpolation of interbeat intervals 
%   estimated from an electrocardiogram R detector at a constant 
%   sample rate of 5 Hz. Use a kernel width of 2 seconds.
%
%      load ECG.mat;
%      np = length(ecg);
%      ri = ECGDetectRInterbeat(ecg,fs,fs);
%      nr = length(ri);
%      ibi = diff(ri)/fs;
%      t   = (ri(1:nr-1) + ri(2:nr))/(2*fs);
%      fsi = 5;
%      ti  = 0:(1/fsi):((np-1)/fs); 
%      x   = SmoothSeries(t,ibi,ti,2,1);
%
%   M.P. Wand and M.C. Jones, Kernel Smoothing. New York: Chapman & 
%   Hall, 1995.
%
%   Version 0.00.00.23 JM
%
%   See also Detectors, Lowpass, and Smooth.

%====================================================================
% Error Checking
%====================================================================    
if nargin<3,
    help SmoothSeries;
    return;
    end;
    
if length(x)~=length(t),
    error('Signal values and times must be of the same length.');
    end;
    
%====================================================================
% Process function arguments
%====================================================================     
wl = 5;                                                    % Default window length
if exist('wla') & ~isempty(wla),
    wl = max(0,wla);
    if wla<=0,
        error('Window length argument was zero or negative.');
        end;
    end;

pf = 0;                                                    % Default - no plotting
if nargout==0,                                             % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
    
%====================================================================
% Preprocessing and Memory Allocation
%====================================================================  
t  = t(:);                                                 % Make into a column vector
x  = x(:);                                                 % Make into a column vector
ti = ti(:);                                                % Make into a column vector

nx = length(x);
ni = length(ti);
y  = zeros(ni,1);

i1 = 1;
i2 = 1;

oct = 0;

for c = 1:ni,
	while (i1<nx) & (ti(c)-t(i1+1)>5*wl) & (i1<i2),
		i1 = i1 + 1;
		end;
	while (i2<nx) & (i2<=i1 | t(i2)-ti(c)<5*wl),
		i2 = i2 + 1;
		end;

	u = (t(i1:i2)-ti(c))/wl;
	w = exp(-u.^2/2);                                      % Gausian kernel weights 
	if sum(w)==0,                                          % No samples within the kernel width, use uniform weights
		oct = oct + 1;
		w = ones(size(w));
		end;
	y(c) = sum(x(i1:i2).*w)/sum(w);
	end;
   
%====================================================================
% Plot Results
%====================================================================  
if pf==1,
    figure;
    FigureSet(1);
    h = plot(t,x,'r',ti,y,'b');
    set(h,'Marker','.');
    set(h(1),'MarkerSize',15);
    set(h(1),'LineWidth',1.5);
    set(h(2),'MarkerSize',10);
    set(h(2),'LineWidth',1.2);
    xlabel('Time (sec)');
    ylabel('Signal (scaled)');
    box off;
    AxisSet;
    legend('Original','Smoothed & Interpolated');
    end;

%====================================================================
% Process Return Arguments
%====================================================================
if nargout==0,
    clear('y');
    end;

