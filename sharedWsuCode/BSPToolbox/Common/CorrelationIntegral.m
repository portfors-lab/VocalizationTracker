function [p,Cr,r] = CorrelationIntegral(X,pfa)
%CorrelationIntegral: Estimate the correaltion dimension
%
%   [Cr,r] = CorrelationIntegral(X,pf);        
%
%   X    Delay embedding of the attractor.
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   Cr   Correlation function. 
%   r    Distance.   
%
%   Calculates the correlation dimension of the attractor
%   embedded in X. The matrix X contains the n-dimensional 
%   embedding of the attractor, where every row of Y is a 
%   vector of delay coordinates:
%  
%           |x(t1), x(t1-T), x(t1-2T), ... , x(t1-(n-1)T)|
%           |x(t2), x(t2-T), x(t2-2T), ... , x(t2-(n-1)T)|
%       X = |x(t3), x(t3-T), x(t3-2T), ... , x(t3-(n-1)T)|
%           | ...    ...       ...     ...      ...      |
%           |x(tN), x(tN-T), x(tN-2T), ... , x(tN-(n-1)T |
%
%   The correlation dimension is defined as a ratio:
%                                  log(C(r))
%               cordim(X) = lim ---------------
%                           r->0    log(r)
%   where C(r), the correlation function, is defined to be
%   the proportion of pairs of orbit points within r units
%   of one another.
%
%   Example: Estimate the correlation dimension of the 
%   attractor reconstructed from the ICP time series. 
%       
%       load ICP;
%       icplp = LowPass(icp,fs,10,1,2,0);
%       icps  = HighPass(icplp(10^4:10^4+10*fs),fs,0.6,2,0);
%       X     = AttractorReconstruction(icps,5,10,1);
%       [Cr,r]= CorrelationIntegral(X,1);
%
%   K. T. Alligood, T. D. Sauer, J. A. Yorker, CHAOS: An
%   Introduction to Dynamical Systems. 
%   Springer-Verlag, 1996, pp.180-183.
%
%   Version 0.00.00.21 MA
%
%   See also AttractorReconstruction.

%=========================================================
% Process function arguments
%=========================================================

if nargin<1 | nargin>2,
    help CorrelationDimension;
    return;
    end;
    
pf = 0;                     % Default - no plotting
if nargout==0,              % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;

    
%=========================================================
% Calculate Euclidian Distances
%=========================================================
L  = size(X);
D  = zeros(L);

for j = 1:L(1,1)-1
    for i = (j+1):L(1,1),
        D(j,i) = sum((X(j,:) - X(i,:)).^2);
    end
end;
D = sqrt(D);


%=========================================================
% Calculate Correlation Function
%=========================================================
% The correlation function C(r) is defined to be the 
% proportion of pairs of orbit points with r units of one
% another.

%e     = -1:-0.1:-4;
ed = D(:);
e     = linspace(0,5);
r     = 2.^e;
Crn   = zeros(1,length(e));

for i=1:length(e);
    Crni    = find(ed < r(i));
    Crn(i)  = length(Crni);
    if Crn(i) == 0;
        Crn(i) = 1;
    end;
end;

Cr    = Crn./length(ed);
p = linefit(log(r),log(Cr));

%=========================================================
% Plot
%=========================================================
if pf ==1,
   % figure; 
   % FigureSet(1);
    plot(log(r), log(Cr), 'b.');
    title('Correlation Dimension');
    xlabel('ln(r)');
    ylabel('ln(C(r))');
    AxisSet;
    box off;
end;

%====================================================================
% Process Return Arguments
%====================================================================
if nargout==0,
    clear('Cr', 'r');
    end;










