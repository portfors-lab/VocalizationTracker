function [y] = WeightedMedian(x, wla, aa, pfa)
%WeightedMedian: Weighted Median Filter 
%
%   [y] = WeightedMedian(x,wl,a,pf)
%
%   x       Input signal
%   wl      Width of the sliding window in samples (odd integer).
%           Default=11.
%   a       Weight vector (must have the same length as wl). 
%           Default=vector of ones.
%   pf      Plot format: 0=none (default), 1=screen.
%
%   y       Filtered signal
%
%   Filters the signal using a Weighted Median Filter. A window
%   of width wl is placed at the beginning of the input vector. The
%   values within the window are weighted by computing a(i)<>x(i)
%   for values of i from 1 to N, where a(i)<>x(i) is a vector which
%   contains the value x(i) a(i) times (for instance, 3<>2 = [2 2 2])
%   The vectors obtained from the different i values are concatenated
%   and the median of the new vector is computed and stored in the 
%   output vector y. The same procedure is repeated as the window 
%   slides through the data, advancing one sample at a time.
%
%   Example: Filter the nonlinear filters test signal using a Weighted 
%   Median Filter with wl = 7 and a = [0 1 1 2 5 2 1 1 0], and plot the 
%   results.
%
%      load NFSignal.mat;
%      [y] = WeightedMedian(x, 31, [0 1 1 2 5 2 1 1 0], 1);
%
%   Astola, J. and Kuosmanen, P., "Fundamentals of Nonlinear Digital 
%   Filtering," CRC Press, pp.73-77, 1997.
%
%   Version 1.00 CC
%
%   See also MedianFilter, WeightedOrderStatistics, and RankOrder.

%--------------------------------------------------------------------
% Process function arguments
%--------------------------------------------------------------------
if nargin<1 | nargin>4,
    help WeightedMedian;
    return;
    end;

wl = 11; % Default window size = 11 samples
if exist('wla') & ~isempty(wla),
    wl = wla;
    end;
 
a = ones(1,wl); % Default value of a = all ones
if exist('aa') & ~isempty(aa),
    a = aa;
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

if s(2) ~= 1,               % Convert input to a Nx1 vector
    x = x';
end

if rem(wl,2) == 0,           % If window width is even, make it odd
    wl = wl+1;
end

N    = length(x);
xap  = [ x(1)*ones((wl-1)/2,1) ; x ; x(end)*ones((wl-1)/2,1)];


%--------------------------------------------------------------------
% Filter signal
%--------------------------------------------------------------------
for i = 1:N
    xr = xap(i:i+(wl-1));
    xwl = [];
    for j = 1:wl
        xwl = [xwl ; xr(j)*ones(a(j),1)];
    end
    y(i,1) = median(xwl);
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
    title('Weighted Median Filter');
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