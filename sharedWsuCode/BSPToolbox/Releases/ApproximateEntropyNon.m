function [ varargout ] = ApproximateEntropyNon(varargin)
%ApproximateEntropyNon: Calculates approximate entropy of time-series.
%
%   [y] = AproximateEntropyNon(x,fs,m,r,N,off,pf);
%
%   x     Input signal.
%   fs    Sampling frequency (Hz).
%   m     The length of the vectors compared in the calculation, 
%         default = 2.
%   r     The threshold value for determining whether the scalar 
%         components of a vector are close, deault = 0.2*std(x).
%   N     Length of the intervals for which the individual apEn 
%         values will be calculated, for relatively good performance
%         and accuracy N should be about 200-300, default = 250.
%   off   Offset between the beginning point of succesive intervals,
%         default = 50.
%   pf    Plot flag: 0=none (default), 1=screen.
%
%   y     A real valued signal that represents the approximate 
%         entropy of the signal at various consecutive intevals
%
%   This function estimates the approximate entropy of a signal at 
%   specified intervals of the signal given. The calculation of the
%   approximate entropy statistical measure is a quantification of 
%   the degree of regularity and complexity present in a given signal.
%     
%   The approximate entropy is a real valued signal that is calculated
%   using the methods decribed in [PINCUS]. For the result of the 
%   calculation to be good estimate of the actual entropy the 
%   following conditions should be satisfied.
%
%   - m should be chosen such that it is within the range 
%     log(n) <= m <= 3.4*log(n) where n is length of x
%   - r = about 15% of the signals standard deviation
%   - r should be at least 3 times greater than the estimated maximum 
%     of the noise component of the signal
%
%   Pincus, Steven M and Goldberger, Ary L. "Physiological time-
%   series analysis: what does regularity quantify?" American 
%   Physiological Society, vol 94. pp H1643-H1656, 1994.
%
%   Version 0.02.00.21 DT
%
%   See also : ApproximateEntropy.


%=====================================================================
% Proccess Arguments
%=====================================================================


    if nargin < 2 | nargin > 7 | nargout > 1
        
        help('ApproximateEntropyNon')
        return
        
    end
    
    % extract input arguments
    x   = varargin{1};
    fs  = varargin{2};
   
    if nargin >= 3
        m = varargin{3};
    else
        m = 2;
    end
    
    if nargin >= 4
        r   = varargin{4};
    else
        r = 0.2*std(x);
    end
    
    if nargin >= 5
        N   = varargin{5};
    else
       if length(x) >= 1000
           N = 250;
       else
           N = floor(length(x)/10);
       end
    end
    
    if nargin >= 6
        off = varargin{6};
    else
       if length(x) >= 1000
           off = 50;
       else
           off = floor(length(x)/10);
       end
   end
    
    if nargin == 7
        pf = varargin{7};
    else
        if nargout == 0
            pf = 1;
        else
            pf = 0;
        end
    end

%=====================================================================
% Error Checking
%=====================================================================


    if length(N)~=1
        error('N must be scalar.');
        return
    elseif isa(N,'numeric') ~= 1
        error('N must be numeric');
        return
    elseif N<=0
        error('N must be positive.');
        return
    elseif floor(N) ~= N
        error('N must be an integer.');
        return
    elseif length(off)~=1
        error('off must be scalar.');
        return
    elseif isa(off,'numeric') ~= 1
        error('off must be numeric');
        return
    elseif off<1
        error('off must be positive.');
        return
    elseif floor(off) ~= off
        error('off must be an integer.');
        return
    end
    
%=====================================================================
% Initialize Variables
%=====================================================================

    
    lenX     = length(x);            % lemgth of signal 
    numVals  = floor((lenX-N)/off);  % number of values to find  
    apen     = zeros(numVals,1);     % time series to output 
    intInd   = 1:N;                  % interval indices
    
%=====================================================================
% Calculate apEn time series
%=====================================================================

    k = 0; % initial interval starting point

    for i = 1:numVals
        
        apen(i) = ApproximateEntropy(x(k+intInd),m,r);
        k = k + off;
        
    end
        
%=====================================================================
% Plotting
%=====================================================================

    if pf

        % lots of cryptic plot formatting commands
        figH  = figure;
        plotH = plot( (0:off:numVals*off-1)/fs, apen, '.r' );
        axH   = get(figH,'CurrentAxes');
        titH  = get(axH,'Title');
        ylabH = get(axH,'YLabel');
        xlabH = get(axH,'XLabel');
        set(titH,'String','Approximate Entropy');
        set(ylabH,'String','Entropy');
        set(xlabH,'String','Time (s)');
        set(axH,'YLim',[0 max(apen)*1.25]);
        
        FigureSet(1)
        
    end
    
%=====================================================================
% Output
%=====================================================================
    
    if nargout == 1
        
        varargout{1} = apen;
        
    end