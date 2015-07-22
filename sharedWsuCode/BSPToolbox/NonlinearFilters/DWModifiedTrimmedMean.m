function [y] = DWModifiedTrimmedMean(x, w1a, w2a, qa, pfa)
%DWModifiedTrimmedMean: Double-Window Modified Trimmed Mean Filter
%
%   [y] = DWModifiedTrimmedMean(x,w1,w2,q,pf)
%
%   x       Input signal
%   w1      Width of the median window in samples (odd integer).
%           Default=11.
%   w2      Width of the mean window in samples (odd integer)
%           (w2 must be >= w1). Default=31.
%   q       Maximum distance allowed (positive integer). Default=0.2
%   pf      Plot format: 0=none (default), 1=screen.
%
%   y       Filtered signal
%
%   Filters the signal using a Double-Window Modified Trimmed Mean Filter. 
%   A window of width w2 is placed at the beginning of the input vector. 
%   The input values within the window are sorted in ascending order, and
%   a window of width w1 (w1 <= w2) is placed in the center of the sorted 
%   vector. The median of these w1 values is calculated. All the values in   
%   the initial window w2 whose distance from the median is greater than q 
%   are trimmed. The remaining values are averaged, and their mean is stored
%   in the output vector y. The same procedure is repeated as the window 
%   slides through the data, advancing one sample at a time.
%
%   Example: Filter the nonlinear filters test signal using a Double-Window
%   Modified Trimmed Mean Filter with q = 0.2, w1 = 11 and w2 = 21, and plot
%   the results.
%
%      load NFSignal.mat;
%      [y] = DWModifiedTrimmedMean(x, 11, 21, 0.2, 1);
%
%   Astola, J. and Kuosmanen, P., "Fundamentals of Nonlinear Digital 
%   Filtering," CRC Press, p.56, 1997.
%
%   Version 1.02 CC
%
%   See also TrimmedMean, WinsorizedTrimmedMean, and ModifiedTrimmedMean.

%--------------------------------------------------------------------
% Process function arguments
%--------------------------------------------------------------------
if nargin<1 | nargin>5,
    help DWModifiedTrimmedMean;
    return;
    end;

w1 = 11; % Default window1 size = 11 samples
if exist('w1a') & ~isempty(w1a),
    w1 = w1a;
end;
    
w2 = 31; % Default window2 size = 31 samples
if exist('w2a') & ~isempty(w2a),
    w2 = w2a;
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

if s(2) ~= 1,                  % Convert input to a Nx1 vector
    x = x';
end

if rem(w1,2) == 0,             % If window1 width is even, make it odd
    w1 = w1+1;
end
if rem(w2,2) == 0,             % If window2 width is even, make it odd
    w2 = w2+1;
end

N    = length(x);
xap  = [ x(1)*ones((w2-1)/2,1) ; x ; x(end)*ones((w2-1)/2,1)];


%--------------------------------------------------------------------
% Filter signal
%--------------------------------------------------------------------
for i = 1:N
    xs    = sort(xap(i:i+(w2-1)));
    xm    = xs(1+(w2-w1)/2:end-(w2-w1)/2);
    xmed  = median(xm);
    sum1  = 0;
    sum2  = 0;
    for m = 1:w2
        if abs(xs(m)-xmed) < q
        sum1 = sum1 + xs(m);
        sum2 = sum2 + 1;
        end
    end
    if sum2 == 0,
        y(i,1) = 0;
    else
        y(i,1)  = (sum1/sum2);
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
    title('Double Window Modified Trimmed Mean Filter');
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

