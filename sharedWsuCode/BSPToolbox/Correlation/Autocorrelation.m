function [ac,t] = Autocorrelation(x,fsa,mla,ata,pfa)
%Autocorrelation: Esimation of signal autocorrelation.
%
%   [ac,t] = Autocorrelation(x,fs,ml,at,pf)
%
%   x        Input signal.
%   fs       Sample rate (Hz). Default = 1 Hz.
%   ml       Maximum lag (seconds). Default = (length(x)-1)/(5*fs).
%   at       Autocorrelation type. 1=Biased (default), 2=Unbiased.
%   pf       Plot flag:  0=none (default), 1=screen.
%
%   ac       Estimated autocorrelation.
%   t        Lag of estimates (sec).
%
%   This is sometimes called the autocorrelation function, the 
%   sample autocorrelation function, the normalized autocorrelation,
%   and the normalized autocovariance.
%
%   The biased and unbiased estimators of autocorrelation are 
%   related to those of the autocovariance by the simple relation
%
%      p(h) = ac(h)/ac(0)
%
%   where ac is the estimated autocovariance,
%
%      ac(h) = inv(n)*sum_{k=1}^{n-|h|} (x(k) - mx)(x(k+h) - mx)
%
%   for the biased estimator. The unbiased estimator of ac is
%
%      ac(h) = inv(n-|h|)*sum_{k=1}^{n-|h|} (x(k) - mx)(x(k+h) - mx)
%
%   The bias and variance of both estimators is given in the 
%   references. These have not been implemented because they require
%   knowledge of the true autocorrelation. However, the plot displays
%   the large sample 95% estimated confidence intervals for a white 
%   Gaussian noise process with the same variance as the sample 
%   variance of the signal.
%
%   Both methods use the fast Fourier transform (FFT) to generate
%   the estimates quickly.
%
%   Example: Plot the estimated autocorrelation of an ABP signal with
%   a maximum lag of 2 seconds.
%
%      load ABPICP.mat
%      x = abp(1:2000);
%      [ac,t] = Autocorrelation(x,125,2,[],1);
%
%   M. B. Priestley, Spectral Analysis and Time Series, 
%   San Francisco, CA: Academic Press, 1981, pp. 321-330.
%
%   P. J. Brockwell, R. A. Davis, Time Series: Theory and Methods, 
%   2nd ed., New York, NY: Springer, 1991, pp. 218-236.
%
%   Version 2.00.22 JM
%
%   See also XCOV, CONV, and Correlation. 

%====================================================================
% Error Check
%====================================================================
if nargin<1,
    help Autocorrelation;
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
    
at = 1;                                 % Default estimator type: biased
if exist('ata') & ~isempty(ata) & ata==2,
    at = 2;
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
ac = zeros(ml,1); 

%====================================================================
% Main Loop
%====================================================================
np = 2^nextpow2(2*nx-1);   % Figure out how much to pad the signal
X  = fft(x,np);
xc = ifft(abs(X).^2);
ac = real(xc(1:ml));

if at==1,                  % Biased estimate
    ac = ac./nx;   
elseif at==2,              % Unbiased estimate
    k  = (0:ml-1).';
    ac = ac./(nx-k);
    end;    
    
ac = ac/ac(1); % This is the only real difference between autocorrelation and autocovariance

%====================================================================
% Postprocessing
%====================================================================    
t = (0:ml-1).'/fs;          % Lags in units of samples

%====================================================================
% Plotting
%====================================================================
if pf==1,
    if at==1,
        cl    = norminv(0.025,0,ones(ml,1)*sqrt(1/nx));
        cu    = norminv(0.975,0,ones(ml,1)*sqrt(1/nx));
    elseif at==2,
        k     = 0:(ml-1);
        cl    = norminv(0.025,0,1./sqrt(nx-k)).';
        cu    = norminv(0.975,0,1./sqrt(nx-k)).';     
        end;
    %fprintf('Percentage Outside: %5.3f\n',sum(ac>cu | ac<cl)/length(ac));
    
    figure;
    FigureSet;
    h = plot(t,ac);
    set(h,'LineWidth',1.2);
    if ml<100,
        set(h,'Marker','.');
        set(h,'MarkerSize',12);
        end;    
    hold on;
        h = plot([0 max(t)],[0 0],'k');
        h = plot(t,cl,'r');
        h = plot(t,cu,'r');
        if at==1,
            h = plot(t, (nx-t*fs)/nx,'g');
            h = plot(t,-(nx-t*fs)/nx,'g');
            end;
        hold off;
    xlim([0 max(t)]);
    ylim([-1 1]);
    box off;
    zoom on;
    xlabel('Lag (s)');
    ylabel('$\rho$');
    title('Autocorrelation');
    AxisSet;
    end

if nargout==0,
    clear ac;
    clear t;
    end



