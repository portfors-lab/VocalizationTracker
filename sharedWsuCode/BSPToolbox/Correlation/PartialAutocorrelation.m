function [pc,t] = PartialAutocorrelation(x,fsa,mla,pfa)
%PartialAutocorrelation: Estimate of signal partial autocorrelation.
%
%   [ac,pt] = PartialAutocorrelation(x,fs,ml,pf)
%
%   x    Input signal.
%   fs   Sample rate (Hz). Default = 1 Hz.
%   ml   Maximum lag (seconds). Default = (length(x)-1)/(5*fs).
%   pf   Plot flag:  0=none (default), 1=screen.
%
%   pc   Estimated autocorrelation.
%   t    Lag of estimates (sec).
%
%   Uses the Durbin-Levinson algorithm to estimate the partial
%   autocorrelation. The partial autocorrelation is defined as 
%
%      p(h) = Corr(r(h+1),r(1))
%
%   where r(h) is the residual of
% 
%      r(h) = x(h) - projection(x(1),...,x(h-1))
%
%   the signal x after it has been regressed on the h-1 previous
%   values of the signal. Since it is a measure of correlation, it
%   is bounded by 1 in magnitude.
%
%   The standard estimator is based on building a series of 
%   autoregressive models of order 1 through ml. The models are 
%   constructed from the estimated autocovariance of the signal.
%   This implementation allows either the biased on unbiased 
%   estimate of the autocovariance to be used (see the Autocovariance
%   description).
%
%   The bias and variance of both estimators is given in the 
%   references. These have not been implemented because they require
%   knowledge of the true autocorrelation. However, the plot displays
%   the large sample 95% estimated confidence intervals for a white 
%   Gaussian noise process with the same variance as the sample 
%   variance of the signal.
%
%   Example: Plot the estimated partial autocorrelation of an ABP 
%   signal with a maximum lag of 2 seconds.
%
%      load ABPICP.mat
%      x = abp(1:2000);
%      [pc,t] = PartialAutocorrelation(x,125,2,1);
%
%   M. B. Priestley, Spectral Analysis and Time Series, 
%   San Francisco, CA: Academic Press, 1981, pp. 371-372.
%
%   P. J. Brockwell, R. A. Davis, Time Series: Theory and Methods, 
%   2nd ed., New York, NY: Springer, 1991, pp. 242-245.
%
%   Version 2.00.22 JM
%
%   See also LEVINSON, SIGNAL, and Correlation. 

%====================================================================
% Error Check
%====================================================================
if nargin<1,
    help PartialAutocorrelation;
    return;
    end
        
if var(x)==0,
    error('Signal is constant.');
    end;
    
%====================================================================
% Process Function Arguments
%====================================================================
fs = 1;                                 % Default sample rate
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
    end; 
    
ml = round(length(x)/5);                % Default sample rate
if exist('mla') & ~isempty(mla),
    ml = round(mla*fs);
    ml = min(ml,length(x));
    end; 
        
pf = 0;                                 % Default - no plotting
if nargout==0,    % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end; 
    
%====================================================================
% Preprocessing
%====================================================================    
nx = length(x);
mx = mean(x);
vx = var(x);
x  = x(:);                 % Convert to a column vector
x  = x-mx;                 % Remove the mean
pc = zeros(ml,1);          % Partial autocorrelation
pt = zeros(ml,1);          % Partial autocorrelation times
pv = zeros(ml,1);          % Predictive variance
mc = zeros(ml,1);          % Model coefficients

%====================================================================
% Main Loop
%====================================================================
ac = Autocovariance(x,fs,ml);

pc(1) = 1;
mc(1) = 1;
pv(1) = vx;

pc(2) = ac(2)/ac(1);
mc(2) = pc(2);
pv(2) = ac(1)*(1-pc(2).^2);

for c1 = 3:ml,
    pc(c1)     = (ac(c1) - mc(2:c1-1).'*ac((c1-1):-1:2))/pv(c1-1);
    mc(2:c1-1) = mc(2:c1-1) - pc(c1)*mc(c1-1:-1:2); 
    mc(c1)     = pc(c1);
    %lc = levinson(ac(1:c1))';
    %[-mc(1:c1) lc(1:c1)]
    pv(c1)     = pv(c1-1)*(1-pc(c1).^2);
    end;

%====================================================================
% Postprocessing
%====================================================================    
t = (0:ml-1).'/fs;          % Lags in units of samples

%====================================================================
% Plotting
%====================================================================
if pf==1,
    cl = norminv(0.025,0,ones(ml,1)*sqrt(1/nx));
    cu = norminv(0.975,0,ones(ml,1)*sqrt(1/nx));
    
    %fprintf('Percentage Outside: %5.3f\n',sum(pc>cu | pc<cl)/length(pc));
    
    figure;
    FigureSet;
    h = plot(t,pc);
    set(h,'LineWidth',1.2);
    if ml<100,
        set(h,'Marker','.');
        set(h,'MarkerSize',12);
        end;
    hold on;
        h = plot([0 max(t)],[0 0],'k:');
        h = plot(t,cl,'r:'); set(h,'LineWidth',2.5);
        h = plot(t,cu,'r:'); set(h,'LineWidth',2.5);
        hold off;
    xlim([0 max(t)]);
    box off;
    zoom on;
    xlabel('Lag (s)');
    ylabel('$\rho$');
    title('Partial Autocorrelation');
    AxisSet;
    end

if nargout==0,
    clear pc;
    clear t;
    end



