function [y] = RankOrder(x, wla, ra, pfa)
%RankOrder: Rank Order Filter 
%
%   [y] = RankOrder(x,wl,r,pf)
%
%   x       Input signal
%   wl      Width of the sliding window in samples (odd integer)
%   r       Rank value (number between 0 and 100)
%   pf      Plot format: 0=none (default), 1=screen.
%
%   y       Filtered signal
%
%   Filters the signal using a Rank Order Filter. A window of width
%   wl is placed at the beginning of the input vector. The value 
%   associated with the rth percentile is computed and stored in the 
%   output vector y. The same procedure is repeated as the window 
%   slides through the data, advancing one sample at a time.
%
%   Example: Filter the nonlinear filters test signal using a Rank
%   Order Filter with r = 75 and wl = 31, and plot the results.
%
%      load NFSignal.m
%      [y] = RankOrder(x, 31, 75, 1);
% 
%   Astola, J. and Kuosmanen, P., "Fundamentals of Nonlinear Digital 
%   Filtering," CRC Press, pp.78-81, 1997.
%
%   Version 1.00 CC
%
%   See also WeightedOrderStatistics and MedianFilter. 

%--------------------------------------------------------------------
% Process function arguments
%--------------------------------------------------------------------
if nargin<1 | nargin>4,
    help RankOrder;
    return;
    end;

wl = 31; % Default window size = 31 samples
if exist('wla') & ~isempty(wla),
    wl = wla;
    end;
 
r = 50; % Default value of a = 50th percentile
if exist('ra') & ~isempty(ra),
    r = ra;
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
    y(i,1) = prctile(xr,r);
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
    title('Rank Order Filter');
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
