function [y] = MedianFilter(x,wla,pfa)
%MedianFilter:  Median Filter 
%
%   [y] = MedianFilter(x,w,pf)
%
%   x       Input signal
%   wl      Width of the sliding window in samples (odd integer).
%           Default=21.
%   pf      Plot format: 0=none (default), 1=screen.
%
%   y       Filtered signal
%
%   Filters the signal using a Median Filter. A window of width w
%   is placed at the beginning of the input vector. The median of
%   the values within the window is calculated and stored in the 
%   output vector y. The same procedure is repeated as the window 
%   slides through the data, advancing one sample at a time.
%
%   Example: Filter the nonlinear filters test signal using a Median 
%   Filter with  wl = 31, and plot the results.
%
%      load NFSignal.mat;
%      [y] = MedianFilter(x, 31, 1);
%
%   Astola, J. and Kuosmanen, P., "Fundamentals of Nonlinear Digital 
%   Filtering," CRC Press, pp.52-55, 1997.
%
%   Version 1.00 CC
%
%   See also WeightedMedian, RankOrder, and WeightedOrderStatistics. 

%====================================================================
% Error Checking
%====================================================================    
if nargin<1,
    help MedianFilter;
    return;
    end;

wl = 21; % Default window size = 21 samples
if exist('wla') & ~isempty(wla),
    wl = wla;
    end;
    
pf = 0;        % Default - no plotting
if nargout==0, % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;

%====================================================================
% Define function variables
%====================================================================
x = x(:);                                                  % Convert into a column vector
nx = length(x);                                            % Length of x

if rem(wl,2) == 0,                                         % If window width is even, make it odd
    wl = wl+1;
    end

%====================================================================
% Proprocessing
%====================================================================
hw = (wl-1)/2;                                             % Half window length
wi = -hw:hw;                                               % Window offset Indices

%====================================================================
% Filter signal
%====================================================================
y = zeros(nx,1);
for c1 = 1:nx,
    i0 = max( 1,c1-hw);
    i1 = min(nx,c1+hw);
    y(c1) = median(x(i0:i1));
    end;
    
%--------------------------------------------------------------------
% Plot Results
%--------------------------------------------------------------------
if exist('pf') & pf == 1,
    figure;
    FigureSet(1);   
    subplot(2,1,1)
        h=plot((0:nx-1), x, 'c', (0:nx-1), y, 'b');
        axis tight;
        title('Median Filter');
        xlabel('Samples');
        ylabel('Amplitude');
        grid on; 
        box off; 
        AxisSet;
    subplot(2,1,2)
        h=plot((0:nx-1), (x-y), 'k');
        axis tight;
        xlabel('Samples');
        ylabel('Error');
        grid on; 
        box off; 
        AxisSet;
        set(gcf, 'PaperOrientation', 'landscape');
    end





