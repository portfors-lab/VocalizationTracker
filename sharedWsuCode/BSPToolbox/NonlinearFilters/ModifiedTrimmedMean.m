function [y] = ModifiedTrimmedMean(x, wla, qa, pfa)
%ModifiedTrimmedMean: Modified Trimmed Mean Filter
%
%   [y] = ModifiedTrimmedMean(x,wl,q,pf)
%
%   x       Input signal
%   wl      Width of the sliding window in samples (odd integer).
%           Default=31.
%   q       Maximum distance allowed between the median and each of
%           the data points within a window so that the data point
%           will not be trimmed when computing the output. Default=0.2.
%   pf      Plot format: 0=none (default), 1=screen.
%
%   y       Filtered signal
%
%   Filters the signal using a Modified Trimmed Mean Filter. A window
%   of width wl is placed at the beginning of the input vector. The
%   input values within the window are sorted in ascending order, and
%   all the values whose distance from the median is greater than q are 
%   trimmed. The remaining values are averaged, and their mean is stored
%   in the output vector y. The same procedure is repeated as the window 
%   slides through the data, advancing one sample at a time.
%
%   Example: Filter the nonlinear filters test signal using a Modified
%   Trimmed Mean Filter with q = 0.2 and wl = 21, and plot the results.
%
%      load NFSignal.mat;
%      [y] = ModifiedTrimmedMean(x, 21, 0.2, 1);
%
%   Astola, J. and Kuosmanen, P., "Fundamentals of Nonlinear Digital 
%   Filtering," CRC Press, p.56, 1997.
%
%   Version 1.00 CC
%
%   See also TrimmedMean, WinsorizedTrimmedMean, and DWModifiedTrimmedMean.

%--------------------------------------------------------------------
% Process function arguments
%--------------------------------------------------------------------
if nargin<1 | nargin>4,
    help ModifiedTrimmedMean;
    return;
    end;

wl = 31; % Default window size = 31 samples
if exist('wla') & ~isempty(wla),
    wl = wla;
    end;
 
q = 0.2; % Default value of q = 0.2
if exist('qa') & ~isempty(qa),
    q = qa;
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
    xs    = sort(xap(i:i+(wl-1)));
    xmed  = median(xs);
    sum1  = 0;
    sum2  = 0;
    for m = 1:wl
        if abs(xs(m)-xmed) < q
        sum1 = sum1 + xs(m);
        sum2 = sum2 + 1;
        end
    end
    y(i,1)  = (sum1/sum2);
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
    title('Modified Trimmed Mean Filter');
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


