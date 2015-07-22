function [cc,lg] = CrossCorrelate(x1,x2,lag,method,pfa)
%CrossCorrelate: Estimate the cross-correlation of two signals.
%
%   [cc,lg] = CrossCorrelate(x1,x2,lag,method,pf)
%
%   x1       Input Signal 1
%   x2       Input Signal 2
%   lag      Length of lag (samples).  Default=length(x).
%   method   'biased', 'unbiased' (default), 'fast'.
%   pf       Plot format:  0=none (default), 1=screen
%
%   cc       Estimated cross-correlation
%   lg      X-axis lags
%
%   CrossCorrelate is a function that gives a normalized correlation
%   of two equal length random stationary signals, x1 and x2.  The 
%   cross-correlation is defined as:
%
%      cc(m) = sum((x1(m) - mean(x1)).*(x2(m+nx) - mean(x2)))
%
%   where m is the correlation index and nx is the length of the
%   signal.  This figure is divided by the var(x1).*var(x2) so that
%   -1<=cc<=+1.  The range of the output is specified by the number
%   of lags chosen, and will be from -lag:+lag.
%
%   The method chosen determines how the output will be normalized.
%   The 'unbiased' method is the most accurate, but slowest method,
%   and uses a normalization factor of 1/(nx-abs(m)).  The 'biased'
%   method is also slow, but will give a positive definite output and
%   uses a normalization factor of 1./nx.  The 'fast' method is the
%   least accurate method.  It calculates the cross-correlation by 
%   taking the inverse fft of (X1 .* conj(X2)).
%
%   If no output argument is specified, the default will plot to the 
%   screen.  
%
%   Example:  Calculate the cross-correlation of an ABP data segment 
%   and an ICP data segment with a lag of 500, a biased output, and 
%   output plotted to the screen.
%
%      load ABPICP.mat
%      x1 = abp(1:1000);
%      x2 = icp(1:1000);
%      [cc lg] = CrossCorrelate(x1,x2,500,'Biased',1);  
%
%   Priestley, M., "Spectral Analysis and Time Series," Academic 
%   Press Limited, pp. 330-333, 1981.
%
%   Version 1.00 LJ
%
%   See Also Autocorrelate, CrossCovariance, AutoCovariance, and 
%   CrossCorrelogram.

%====================================================================
% Error Check
%====================================================================
if nargin < 2 
    help CrossCorrelate;
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
    error('Signal ''1'' is constant.');
    end

nx2 = length(x2);
xr2 = x2;
mx2 = mean(x2);
sx2 = std(x2);
if sx2==0,
    error('Signal ''2'' is constant.');
    end

if nx1~=nx2,
    fprintf('Signals must be the same length.');
    end

nx   = length(x1);
mlag = nx; % Default lag
n    = (-nx+1:nx-1)';
lg  = n;
x1   = x1(:);  %  Make sure x1 is a column vector
x2   = x2(:);  %  Make sure x2 is a column vector

%====================================================================
% Process function arguments
%====================================================================
if exist('lag') & ~isempty(lag)
    mlag = min(lag,mlag);
    lg  = (-mlag+1:mlag-1)';
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
if strcmpi(method,'Biased')  
    a = 0;
    x1mean    = mean(x1);
    x2mean    = mean(x2);
    variance1 = var(x1);
    variance2 = var(x2);
    for c = 1:length(n)
        m = n(c);
        if c < nx
            cor1 = x1(1:c);
            cor2 = x2(end-c+1:end);
        else
            cor1 = x1(c-nx+1:end);
            cor2 = x2(1:c-2*a);
            a    = a + 1;
        end
        num(c) = 1./nx.*sum((cor1-x1mean).*(cor2-x2mean));
    end
    cc = num/(sqrt(variance1*variance2));
    cc = cc(nx-mlag+1:nx+mlag-1);
    cc = cc(:);
    
%====================================================================
% Fast Method
%====================================================================
elseif strcmpi(method,'Fast')  
    cor1  = x1 - mean(x1);
    cor2  = x2 - mean(x2);
    zpad  = 2^(ceil(log2(nx))+1);     %  Zero vector 
    xpad1 = [cor1;zeros(zpad-nx,1)];  %  Pad 'x1' with zeros 
    xpad2 = [cor2;zeros(zpad-nx,1)];  %  Pad 'x2' with zeros 
    Xpad1 = fft(xpad1);
    Xpad2 = fft(xpad2);
    cc    = real(ifft(Xpad1.*conj(Xpad2)));
    cc    = cc./((nx-1)*sqrt(var(x1)*var(x2)));
    cc    = [cc(end-mlag+1:end);cc(1:mlag-1)];
    cc    = cc*nx./(nx-abs(lg));
    cc    = cc(:);
    
%====================================================================
% Unbiased Method
%====================================================================
elseif strcmpi(method,'Unbiased')  
    a = 0;
    for c = 1:length(n)
        m = n(c);
        if c < nx
            cor1 = x1(1:c);
            cor2 = x2(end-c+1:end);
        else
            cor1 = x1(c-nx+1:end);
            cor2 = x2(1:c-2*a);
            a    = a + 1;
        end
        x1mean  = mean(cor1);
        x2mean  = mean(cor2);
        num(c)  = 1./(nx-abs(m)).*sum((cor1-x1mean).*(cor2-x2mean));
        den1(c) = 1./(nx-abs(m)).*sum((cor1-x1mean).^2);
        den2(c) = 1./(nx-abs(m)).*sum((cor2-x2mean).^2);
        if den1(c)>0 & den2(c)>0
        cc(c) = num(c)./sqrt(den1(c).*den2(c));
        end
    end
    cc = [cc';zeros((length(n)-length(cc)),1)];  
    cc = cc(nx-mlag+1:nx+mlag-1);
    cc = cc(:);
   
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
    axis tight;
    %axis([-mlag mlag -1 1]);
    ylabel('Cross-Correlation');
    xlabel('Lag');
    title('Cross-Correlation','FontWeight','bold');
    AxisSet;
    end

if nargout == 0,
     clear cc;
     clear lg;
end
