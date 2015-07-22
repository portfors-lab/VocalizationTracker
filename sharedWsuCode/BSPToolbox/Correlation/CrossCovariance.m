function [cc,lg] = CrossCovariance(x1,x2,lag,method,pfa)
%CrossCovariance: Estimate the cross-covariance of two signals.
%
%   [cc,lg] = CrossCovariance(x1,x2,lag,method,pf)
%
%   x1       Input signal 1
%   x2       Input signal 2
%   lag      Length of lag (samples).  Default=length(x).
%   method   'biased', 'unbiased' (default), 'fast'.
%   pf       Plot format:  0=none (default), 1=screen.
%
%   cc       Estimated cross-covariance
%   lg       X-axis lags
%
%   CrossCovariance is a function that gives an unnormalized 
%   covariance of two equal length random stationary signals, x1 
%   and x2.  The cross-covariance is defined as
%
%      cc(m) = sum((x1(m) - mean(x1)).*(x2(m+nx) - mean(x2)))
%
%   where m is the covariance index and nx is the length of
%   signal.  The range of the output is specified by the number
%   of lags chosen, and will be from -lag:+lag.
%
%   The method chosen determines how the output will be normalized.
%   The 'unbiased' method is the most accurate, but slowest method,
%   and uses a normalization factor of 1/(nx-abs(m)).  The 'biased'
%   method is also slow, but will give a positive definite output and
%   uses a normalization factor of 1./nx.  The 'fast' method is the
%   least accurate method.  It calculates the cross-covariance by 
%   taking the inverse fft of (X1 .* conj(X2)).  If no output 
%   argument is specified, the default will plot to the screen.  
%
%   Example:  Calculate the cross-covariance of an ABP signal and an
%   ICP signal with a lag of 1000 and a biased output printed to the 
%   screen.
%
%      load ABPICP.mat
%      x1 = abp(1:2000);
%      x2 = icp(1:2000);
%      [cc,lg] = CrossCovariance(x1,x2,1000,'biased',1);
%
%   Challis, R. E., and Kitney, R. I., "Biomedical Signal Processing
%   (in four parts), Part 1, Time-domain Methods", Med. & Biol. Eng.
%   & Comput., 1990, 28, pp. 509-524.
%
%   Version 1.00 LJ
%
%   See Also AutoCov, Autocorrelate, and CrossCorrelate.

%====================================================================
% Error Check
%====================================================================
if nargin < 2 
    help CrossCovariance;
    return;
end

%====================================================================
% Process Input Signals
%====================================================================
nx1 = length(x1);
xr1 = x1;
mx1 = mean(x1);
sx1 = std(x1);
if sx1==0,
    fprintf('ERROR: Signal ''1'' is constant.\n');
    return;
end

nx2 = length(x2);
xr2 = x2;
mx2 = mean(x2);
sx2 = std(x2);
if sx2==0,
    fprintf('ERROR: Signal ''2'' is constant.\n');
    return;
end

if nx1~=nx2,
    fprintf('ERROR: Signals must be of the same length.\n');
    return;
end

nx   = length(x1);
mlag = nx; % Default lag
n    = (-nx+1:nx-1)';
lg   = n;
x1   = x1(:);  % Convert to a column vector
x2   = x2(:);  % Convert to a column vector

%====================================================================
% Process function arguments
%====================================================================
if exist('lag') & ~isempty(lag)
   mlag = min(lag,mlag);
   lg   = (-mlag+1:mlag-1)';
end

if nargin < 4
   method = 'Unbiased';
end

if nargout == 0 & ~exist('pfa')
    pf = 1; % Plot to screen if no output arguments
else
    pf = 0;
end

if exist('pfa') & ~isempty(pfa)
    pf = pfa;
end

%====================================================================
% Biased Method
%====================================================================
if strcmpi(method, 'Biased')
    a = 0;
    x1mean = mean(x1);
    x2mean = mean(x2);
    for c = 1:length(n)
        m = n(c);
        if c < nx
           cov1 = x1(1:c);
           cov2 = x2(end+1-c:end);
        else
           cov1 = x1(c-nx+1:end);
           cov2 = x2(1:c-2*a);
           a    = a + 1;
        end
        cc(c) = 1./nx.*sum((cov1-x1mean).*(cov2-x2mean));
     end
     cc = cc(nx-mlag+1:nx+mlag-1);
     cc = cc(:);
     
%====================================================================
% Unbiased Method
%====================================================================
elseif strcmpi(method, 'Unbiased')
    a = 0;
    for c = 1:length(n)
        m = n(c);
        if c < nx
           cov1 = x1(1:c);
           cov2 = x2(end+1-c:end);
        else
           cov1 = x1(c-nx+1:end);
           cov2 = x2(1:c-2*a);
           a    = a + 1;
        end
     x1mean = mean(cov1);
     x2mean = mean(cov2);
     cc(c)  = 1./(nx-abs(m)).*sum((cov1-x1mean).*(cov2-x2mean));
     end
     cc = cc(nx-mlag+1:nx+mlag-1);
     cc = cc(:);
     
%====================================================================
% Fast Method
%====================================================================
elseif strcmpi(method, 'Fast')
    x1    = x1 - mean(x1);
    x2    = x2 - mean(x2);
    zpad  = 2^(ceil(log2(nx))+1);   %  Zero pad  
    xpad1 = [x1;zeros(zpad-nx,1)];  %  Pad 'x1' with zeros 
    xpad2 = [x2;zeros(zpad-nx,1)];  %  Pad 'x2' with zeros 
    Xpad1 = fft(xpad1);
    Xpad2 = fft(xpad2);
    cc    = real(ifft(Xpad1.*conj(Xpad2)));
    cc    = [cc(end-mlag+1:end);cc(1:mlag-1)];
    cc    = cc*nx./((nx-abs(lg)).*(nx-1));
    cc    = cc(:);
   
%====================================================================
% Error Message, Unknown Method Specified
%====================================================================
else
   error('Method must be ''Biased'', ''Unbiased'', or ''Fast''');
end

%====================================================================
% Plotting
%====================================================================
if pf==1,
    figure;
    FigureSet;    
    h = plot(lg,cc);
    set(h,'LineWidth',1.2);
    box off;
    xlabel('Lag');
    ylabel('Cross-Covariance');
    title('Cross-Covariance','FontWeight','bold');
    AxisSet;
    end

if nargout == 0
    clear cc;
    clear lg;
    end