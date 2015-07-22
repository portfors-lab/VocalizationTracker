function [varargout] = ApproximateEntropy(varargin)
%ApproximateEntropy: Calculates the approximate entropy of a signal.
%
%   [y] = AproximateEntropy(x,fs,m,r,pf);
%
%   x    Input signal.
%   fs   The sample frequency (Hz).
%   m    The length of the vectors compared in the calculation,
%        default = 2.
%   r    The threshold value for determining whether the scalar 
%        components of a vector are close, default = 0.2*std(x).
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   y    A real value that represents the approximate entropy of the
%        signal, a measure of regularity.
%
%   This function estimates the approximate entropy of a signal. The
%   calculation of the approximate entropy statistical measure is a 
%   quantification of the degree of regularity and complexity present 
%   in a given signal.
%     
%   The approximate entropy is a real value that is calculated using 
%   the methods decribed in [PINCUS]. For the result of the 
%   calculation to be good estimate, the following conditions should 
%   be satisfied.
%
%   - m should be chosen such that it is within the range 
%     log(n) <= m <= 3.4*log(n) where n is length of x
%   - r = about 15% of the signals standard deviation
%   - r should be at least 3 times greater than the estimated maximum 
%     of the noise component of the signal
%
%   Example: Calculate approximate entropy part of ECG signal, plot
%   condititional probability.
%
%      load ECG
%      R  = ECGDetectRInterbeat(ecg,fs,fs);
%      rr = [ R(1)-1 ; R(2:length(R))-R(1:length(R)-1)-1 ]/fs;
%      x  = SmoothSeries(R/fs,rr,(0:200:length(ecg)-1)/fs,2);
%      r  = 0.2*std(x);
%      apEn = ApproximateEntropy(x(1:5000),fs/250,3,r,1);
%
%   Pincus, Steven M and Goldberger, Ary L. "Physiological time-
%   series analysis: what does regularity quantify?," American 
%   Physiological Society, vol 94. pp H1643-H1656, 1994.
%
%   Version 0.03.00.21 DT
%
%   See also ApproximateEntropyNon.


%=====================================================================
% Proccess Arguments
%=====================================================================
if nargin < 2 | nargin > 5 | nargout > 1
    help('ApproximateEntropy')
    return;
    end;
    
%=====================================================================
% Extract input arguments
%=====================================================================
x  = varargin{1};
fs = varargin{2};

if nargin >= 3
    m  = varargin{3};
else
    m = 2;
    end

if nargin >= 4
    r  = varargin{4};
else
    r = 0.2*std(x);
    end

if nargin == 5 
    pf = varargin{5};
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
sigSize = size(x);

if sigSize(1) == 0
    error('x is empty.');
    return
elseif sigSize(2) ~= 1
    error('x is not a column vector.');
    return
elseif isa(x,'numeric') ~= 1   
    error('x must be numeric');
    return
elseif length(fs)~=1
    error('fs must be scalar.');
    return
elseif isa(fs,'numeric') ~= 1 
    error('fs must be numeric');
    return
elseif fs <= 0
    error('fs must be positive.');
    return
elseif length(m)~=1
    error('m must be scalar.');
    return
elseif isa(m,'numeric') ~= 1 
    error('m must be numeric');
    return
elseif m<1 | m>sigSize(1)
    error('m is out of range.');
    return
elseif floor(m) ~= m
    error('m must be an integer.');
    return
elseif length(r)~=1
    error('r must be scalar.');
    return
elseif isa(r,'numeric') ~= 1  
    error('r must be numeric');
    return
elseif r<=0
    error('r must be positive.');
    return 
elseif length(pf) ~= 1  
    error('pf must be scalar.');
    return
elseif isa(pf, 'numeric') ~= 1
    error('pf must be numeric');
    return
elseif pf ~= 0 & pf ~= 1
    error('pf must be a Boolean Value.');
    return
    end

%=====================================================================
% Variables
%=====================================================================
xLen       = length(x);        % Signal length
numVect    = xLen-m;           % Number of vectors
condProb   = zeros(numVect,1); % Conditional probabilities
vectInd    = 0:m-1;            % Vector indices
closeInd   = [];               % indices of close vectors
maxEn      = 1/12;             % small number to use for no match entropy
    
%=====================================================================
% Find componentwise close sets of vectors and use
% to calculate respective conditional probabilities 
%=====================================================================
if m==1,
    xp = x(1:numVect);
    for i = 1:numVect,
        ci = abs(x(i)-xp)<=r;
		closeInd = find(ci);	
        totalClose = length(closeInd);
        if totalClose==1,                                      % Is there a single close unique close vector (the point itself)
            totalNextClose = maxEn;
        else                                                   % number of vectors that are close and whose next values are also close
            totalNextClose  = sum(abs(x(i+m)-x(closeInd+m))<=r);
            end
        condProb(i) = totalNextClose/totalClose;               % calculate conditional probability   
        end;
else
    for i = 1:numVect
        for j = 1:numVect,                                     % Find component wise close vectors
            close = all(abs(x(vectInd+i)-x(vectInd+j))<=r);
            if close,
                closeInd = [closeInd  j];
                end
            end
        totalClose = length(closeInd);
        if totalClose == 1,                                    % Is there a close unique close vector
            totalNextClose = maxEn;
        else                                                   % number of vectors that are close and whose next values are also close
            totalNextClose  = sum( abs( x(i+m) - x(closeInd+m) ) <= r );
            end
        closeInd = [];                                         % reset close indices
        condProb(i) = totalNextClose/totalClose;               % calculate conditional probability
        end
    end;
%=====================================================================
% Approximate Entropy is the mean of the logarithm of the 
% conditional probabilities, plot conditional probability
%=====================================================================
apEn = -sum(log(condProb))/(numVect);
    
%=====================================================================
% Plotting
%=====================================================================
if pf      
    % lots of cryptic plot formatting commands
    figH  = figure;
    plotH = plot( (0:numVect-1)/fs, condProb, '.r' );
    axH   = get(figH,'CurrentAxes');
    titH  = get(axH,'Title');
    ylabH = get(axH,'YLabel');
    xlabH = get(axH,'XLabel');
    set(titH,'String','Conditional Probability');
    set(ylabH,'String','Probability');
    set(xlabH,'String','Time (s)');
    set(axH,'YLim',[0 1]);
    FigureSet(1);
    end

%=====================================================================
% Output
%=====================================================================
    
    if nargout == 1
        
        varargout{1} = apEn;
        
    end
    
%=====================================================================