function [ac,t] = Autocovariance(x,fsa,mla,ata,pfa)
%Autocovariance: Esimation of signal autocovariance.
%
%   [ac,t] = Autocovariance(x,fs,ml,at,pf)
%
%   x        Input signal.
%   fs       Sample rate (Hz). Default = 1 Hz.
%   ml       Maximum lag (seconds). Default = (length(x)-1)/fs.
%   at       Autocovariance type. 1=Biased (default), 2=Unbiased.
%   pf       Plot flag:  0=none (default), 1=screen.
%
%   ac       Estimated autocovariance.
%   t        Lag of estimates (sec).
%
%   The biased (default) estimator is given by
%
%      ac(h) = inv(n)*sum_{k=1}^{n-|h|} (x(k) - mx)(x(k+h) - mx)
%
%   where mx = mean(x) and x is assumed to be n samples long.
%
%   This estimator is biased, but its asymptotic distribution has a 
%   mean that equals the true autocorrelation under mild assumptions.
%
%   The unbiased estimator is given by
%
%      ac(h) = inv(n-|h|)*sum_{k=1}^{n-|h|} (x(k) - mx)(x(k+h) - mx)
%
%   This estimator is not unbiased since the signal mean is estimated
%   from the data, but it is asymptotically unbiased and the bias is 
%   a smaller order of magnitude than the variance.
%
%   Many researchers prefer the biased estimate because it has lower
%   variance and is guaranteed to be semi-postive definite. This 
%   later property results in a power spectral density estimate that 
%   is always positive and autoregressive models that are always 
%   stable. It is also known that this property is shared with the 
%   true autocovariance of stationary processes. It has also been 
%   claimed that the biased estimate has a lower mean squared error 
%   due to the decreased variance of the estimate for large lags.
%
%   The bias and variance of both estimators is given in the 
%   references. These have not been implemented because they require
%   knowledge of the true autocovariance. However, the plot displays
%   the large sample 95% estimated confidence intervals for a white 
%   Gaussian noise process with the same variance as the sample 
%   variance of the signal.
%
%   Both methods use the fast Fourier transform (FFT) to generate
%   the estimates quickly. However, for small maximum lags relative
%   to the signal length this may actually be slower than the direct
%   implementation.
%
%   Example: Plot the estimated autocovariance of an ABP signal with
%   a maximum lag of 5 seconds.
%
%      load ABPICP.mat
%      x = abp(1:2000);
%      [ac,t] = Autocovariance(x,125,5,[],1);
%
%   M. B. Priestley, Spectral Analysis and Time Series, 
%   San Francisco, CA: Academic Press, 1981, pp. 321-330.
%
%   Version 2.00.22 JM
%
%   See also XCOV, CONV, and Correlation. 

%====================================================================
% Error Check
%====================================================================
if nargin<1,
    help Autocovariance;
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
    
ml = length(x);                       % Default sample rate
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
    
%====================================================================
% Postprocessing
%====================================================================    
t = (0:ml-1).'/fs;         % Lags in units of samples

%====================================================================
% Plotting
%====================================================================
if pf==1,
    if at==1,
        k     = 0:(ml-1);
        av    = (1/nx)*(1-k/nx).*vx.^2;
        av(1) = 2*inv(nx)*vx^2;
    elseif at==2,
        k     = 0:(ml-1);
        av    = (1./(nx-k)).*vx.^2;
        av(1) = 2*inv(nx)*vx^2;        
        end;
    cl    = norminv(0.025,0,sqrt(av)).';
    cu    = norminv(0.975,0,sqrt(av)).';
    cl(1) = cl(1)+ac(1);
    cu(1) = cu(1)+ac(1);
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
        h = plot([0 max(t)],[0 0],'k:');
        h = plot(t,cl,'r:');
        h = plot(t,cu,'r:');
        hold off;
    xlim([0 max(t)]);
    ymn = min(ac);
    ymx = max(ac);
    yrg = ymx - ymn;
    ylim([(ymn-.02*yrg) (ymx+.02*yrg)]);
    box off;
    zoom on;
    xlabel('Lag (s)');
    ylabel('Autocovariance');
    title('Autocovariance');
    AxisSet;
    end

if nargout==0,
    clear ac;
    clear t;
    end



