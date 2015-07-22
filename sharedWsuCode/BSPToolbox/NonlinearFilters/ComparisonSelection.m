function [y] = ComparisonSelection(x, wla, da, pfa)
%ComparisonSelection: Comparison and Selection Filter 
%
%   [y] = ComparisonSelection(x,w,d,pf)
%
%   x       Input signal
%   wl      Width of the sliding window in samples (odd integer). 
%           Default = ceiling of length(x)/100.
%   d       Distance in rank from the median that the output value is
%           to be taken from. Default = w/10.
%   pf      Plot format: 0=none (default), 1=screen.
%
%   y       Filtered signal
%
%   Filters the signal using a Comparison and Selection Filter. A 
%   window of width w is placed at the beginning of the input vector. 
%   The input values within the window are sorted in ascending order, 
%   and the mean and median of the values within the window are 
%   calculated. If the mean is greater than the median, the value 
%   which is d places before the central value (the median value)is
%   stored in the output vector y. Conversely, if the median is 
%   greater then the mean, the value which is d places after the 
%   central value is stored in the output vector y. If both mean and 
%   median are equal, the central value is stored into the output 
%   vector y. The same procedure is repeated as the window slides 
%   through the data, advancing one sample at a time.
%
%   Example: Filter the nonlinear filters test signal using a  
%   Comparison Selection Filter with window length w = 31 and distance
%   value d = 3 samples, and plot the results.
%
%      load NFSignal.mat
%      [y] = ComparisonSelection(x, 31, 3, 1);
%
%   Astola, J. and Kuosmanen, P., "Fundamentals of Nonlinear Digital 
%   Filtering," CRC Press, pp.107, 1997.
%
%   Version 1.02 CC
%
%   See also SelectiveAverage and SelectiveMedian.

%--------------------------------------------------------------------
% Process function arguments
%--------------------------------------------------------------------
if nargin<1 | nargin>4,
    help ComparisonSelection;
    return;
    end;

wl = ceil(length(x)/100);                  % Default window size 
if exist('wla') & ~isempty(wla),
    wl = wla;
end;

d = ceil(wl/10);                     %Default value for d
if exist('da') & ~isempty(da),
    d=da;
end;
    
pf = 0;                             % Default - no plotting
if nargout==0,                      % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;

%--------------------------------------------------------------------
% Define function variables
%--------------------------------------------------------------------
s    = size(x);
if s(2) ~= 1,                   % Convert input to a Nx1 vector
    x = x';
end

if rem(wl,2) == 0,               % If window width is even, make it odd
    wl = wl+1;
end

N    = length(x);
xap  = [ x(1)*ones((wl-1)/2,1) ; x ; x(end)*ones((wl-1)/2,1)];


%--------------------------------------------------------------------
% Filter signal
%--------------------------------------------------------------------
for i = 1:N
    xr = xap(i:i+(wl-1));
    xs = sort(xr);
    xmean = mean(xr);
    xmedian = median(xr);
    if xmean > xmedian
       y(i,1) = xs((wl-1)/2-d);
   elseif xmean < xmedian
       y(i,1) = xs((wl-1)/2+d);
   else
       y(i,1) = xs((wl-1)/2+1);
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
    title('Comparison and Selection Filter');
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




