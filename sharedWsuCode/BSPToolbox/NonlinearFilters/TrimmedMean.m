function [y] = TrimmedMean(x, wla, alphaa, pfa)
%TrimmedMean: Trimmed Mean Filter
%
%   [y] = TrimmedMean(x,wl,alpha,pf)
%
%   x       Input signal
%   wl      Width of the sliding window in samples (odd integer).
%           Default=31.
%   alpha   Proportion of values to be trimmed at both ends of the 
%           window (0 <= alpha <=0.48). Default=0.20.
%   pf      Plot format: 0=none (default), 1=screen.
%
%   y       Filtered signal
%
%   Filters the signal using an alpha-Trimmed Mean Filter. A window
%   of width wl is placed at the beginning of the input vector. The
%   input values within the window are sorted in ascending order, and
%   a number of alpha*wl values are trimmed at each end. The remaining
%   values are averaged, and their mean is stored in the output vector 
%   y. The same procedure is repeated as the window slides through the
%   data, advancing one sample at a time.
%
%   Example: Filter the nonlinear filters test signal using a Trimmed 
%   Mean Filter with alpha = 0.15 and wl = 31, and plot the results.
%
%     load NFSignal.mat;
%     [y] = TrimmedMean(x, 31, 0.15, 1);
%
%   Astola, J. and Kuosmanen, P., "Fundamentals of Nonlinear Digital 
%   Filtering," CRC Press, pp.52-55, 1997.
%
%   Version 1.00 CC
%
%   See also WinsorizedTrimmedMean, ModifiedTrimmedMean, and 
%   DWModifiedTrimmedMean.

%--------------------------------------------------------------------
% Process function arguments
%--------------------------------------------------------------------
script = 0;
if ~script
if nargin<1 | nargin>4,
    help TrimmedMean;
    return;
    end;

wl = 31; % Default window size = 31 samples
if exist('wla') & ~isempty(wla),
    wl = wla;
    end;
 
alpha = 0.20; % Default value of alpha = 0.20
if exist('alphaa') & ~isempty(alphaa),
    alpha = alphaa;
    end;
    
pf = 0; % Default - no plotting
if nargout==0, % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
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
    xs    = sort(xap(i:i+(wl-1)));
    y(i,1)  = mean(xs(1+floor(alpha*wl)):xs(floor(end-alpha*wl)));
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
    title('Trimmed Mean Filter');
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

