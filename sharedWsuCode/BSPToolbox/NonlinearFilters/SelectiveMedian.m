function [y] = SelectiveMedian(x, wla, pfa)
%SelectiveMedian: Selective Median Filter 
%
%   [y] = SelectiveMedian(x,wl,pf)
%
%   x       Input signal
%   wl      Width of the sliding window in samples (odd integer).
%           Default=11.
%   pf      Plot format: 0=none (default), 1=screen.
%
%   y       Filtered signal
%
%   Filters the signal using a Selective Median Filter. A window
%   of width wl is placed at the beginning of the input vector. The
%   input values within the window are divided into three groups:
%   the central value, the values before the central value (lower 
%   values), and the values after the central value (upper values). 
%   The median of the lower values and the median of the upper values 
%   are then computed and compared to the central value. The one whose
%   value is closest to the central value is stored in the output 
%   vector y. If both of them are at the same distance from the central 
%   value, the median of the whole window is stored in the output vector
%   y. The same procedure is repeated as the window slides through the
%   data, advancing one sample at a time.
%
%   Example: Filter the nonlinear filters test signal using a Selective 
%   Median Filter with wl = 31, and plot the results.
%
%      load NFSignal.mat;
%      [y] = SelectiveMedian(x, 31, 1);
%
%   Astola, J. and Kuosmanen, P., "Fundamentals of Nonlinear Digital 
%   Filtering," CRC Press, pp.108-109, 1997.
%
%   Version 1.02 CC
%
%   See also SelectiveAverage and ComparisonSelection.

%--------------------------------------------------------------------
% Process function arguments
%--------------------------------------------------------------------
if nargin<1 | nargin>3,
    help SelectiveMedian;
    return;
    end;

wl = 11; % Default window size = 30 samples
if exist('wla') & ~isempty(wla),
    wl = wla;
end;

    
pf = 0; % Default - no plotting
if nargout==0, % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;

%--------------------------------------------------------------------
% Define function variables
%--------------------------------------------------------------------
s    = size(x);

if s(2) ~= 1,
    x = x';
end

if rem(wl,2) == 0,
    wl = wl+1;
end

N    = length(x);
xap  = [ x(1)*ones((wl-1)/2,1) ; x ; x(end)*ones((wl-1)/2,1)];


%--------------------------------------------------------------------
% Filter signal
%--------------------------------------------------------------------
for i = 1:N
    xr = xap(i:i+(wl-1));
    medianl = median(xr(1:(wl-1)/2));
    medianu = median(xr((wl-1)/2 +2:wl));
    if abs(medianl - xr((wl-1)/2 +1)) < abs(medianu - xr((wl-1)/2+1))
       y(i,1) = medianl;
   elseif abs(medianl - xr((wl-1)/2 +1)) > abs(medianu - xr((wl-1)/2+1))
       y(i,1) = medianu;
   else
       y(i,1) = median(xr);
   end
end
    
    

%--------------------------------------------------------------------
% Plot Results
%--------------------------------------------------------------------
if exist('pf') & pf == 1,
    figure
    subplot(2,1,1)
    plot((0:N-1), x, 'r', (0:N-1), y, 'b');
    axis([0 N min(x)-4 max(x)+4]);
    xlabel('Samples');
    ylabel('Amplitude');
    grid on;
    title('Selective Median Filter');
    subplot(2,1,2)
    plot((0:N-1), (x-y), 'k');
    axis([0 N min(x)-4 max(x)+4]);
    xlabel('Samples');
    ylabel('Error');
    grid on;
    subplot(2,1,1)
    AxisSet(10); 
    FigureSet(1);
end
