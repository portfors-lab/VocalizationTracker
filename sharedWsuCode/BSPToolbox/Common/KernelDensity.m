function [f,xt,bh] = KernelDensity(x,xta,wla,pfa);
%KernelDensity: Kernel density estimation
%
%   [f,xt] = KernelDensity(x,xt,wl,pf);
%
%   x    Values of observations.
%   xt   Values of x where the PDF is to be estimated. Default =
%        100 points uniformly spaced from min(x) to max(x).
%   wl   Length of kernel window to use (sec). Default = 0.1*std(x).
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   xt   Values of x where PDF is estimated.
%   f    Estimated PDF of x at values xt.
%
%   Estimates the probability density function (PDF) using a 
%   traditional kernel smoother.
%
%   The current version uses a truncated guassian kernel with 
%   standard deviation specified by wl. Points more than 5 standard 
%   deviations away from the evaluation point are ignored. 
%
%   Example: Estimate the PDF of spike indices.
%
%      load MER.mat;
%      x   = x(1:round(fs*5));
%      sp  = (Lowpass(x.^2,fs,200,4)).^(1/2);
%      pi  = DetectMaxima(sp);
%      f   = KernelDensity(sp(pi),[],[],1);
%
%   M.P. Wand and M.C. Jones, Kernel Smoothing. New York: Chapman & 
%   Hall, 1995.
%
%   Version 0.00.00.23 JM
%
%   See also Detectors, Lowpass, and Smooth.

%====================================================================
% Error Checking
%====================================================================    
if nargin<1,
    help KernelDensity;
    return;
    end;
    
%====================================================================
% Process function arguments
%====================================================================   
xmin = min(x);
xmax = max(x);
xrng = xmax-xmin;
xt = linspace(xmin-0.025*xrng,xmax+0.025*xrng,100);        % Default estimation points
if exist('xta') & ~isempty(xta),
    xt = xta;
    end;

wl = 0.1*std(x);                                           % Default window length
if exist('wla') & ~isempty(wla),
    wl = max(eps,wla);
    end;

pf = 0;                                                    % Default - no plotting
if nargout==0,                                             % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
    
%====================================================================
% Preprocessing and Memory Allocation
%====================================================================  
xt = xt(:);                                                % Make into a column vector
xs = sort(x(:));                                           % Sort x in increasing order

nx = length(x);
ni = length(xt);
f  = zeros(ni,1);

i1 = 1;
i2 = 1;

oct = 0;

for c = 1:ni,
	while (i1<nx) & (xt(c)-xs(i1+1)>5*wl) & (i1<i2),
		i1 = i1 + 1;
		end;
	while (i2<nx) & (i2<=i1 | xs(i2)-xt(c)<5*wl),
		i2 = i2 + 1;
        end;
	u = (xs(i1:i2)-xt(c))/wl;                              % Arguments for the Gaussian kernels
	w = exp(-u.^2/2);                                      % Gausian kernel weights 
	f(c) = sum(w)/(nx*sqrt(2*pi*wl^2));
	end;

%====================================================================
% Plot Results
%====================================================================  
if pf,
    figure;
    FigureSet(1);
    bc = linspace(xmin-0.025*xrng,xmax+0.025*xrng,120)';    % Bin centers
    bh = hist(x,bc);                                       % Bin height
    bs = min(diff(bc));                                    % Bin spacing
    bh = bh/(sum(bh)*bs);                                  % Scale to have unit area
    h = bar(bc,bh,1);
    set(h,'FaceColor',0.7*[1 1 1]);
    set(h,'EdgeColor',0.5*[1 1 1]);
    hold on;
        h = plot(xt,f);
        set(h,'Marker','.');
        set(h,'LineWidth',2);
        hold off;
    xlabel('Data');
    ylabel('Density');
    xlim([min(xt) max(xt)]);
    ylim([0 max([bh(:);f])*1.025]);
    box off;
    AxisSet;
    end;

%====================================================================
% Process Return Arguments
%====================================================================
if nargout==0,
    clear('f','xt');
    end;

