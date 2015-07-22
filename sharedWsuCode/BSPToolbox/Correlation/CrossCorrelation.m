function [cc,t] = CrossCorrelation(x1,x2,fsa,mla,ata,pfa)
%CrossCorrelation: Esimation of signal cross correlation.
%
%   [ac,t] = CrossCorrelation(x1,x2,fs,ml,at,pf)
%
%   x1       Input signal 1.
%   x2       Input signal 2.
%   fs       Sample rate (Hz). Default = 1 Hz.
%   ml       Maximum lag (seconds). Default = (length(x)-1)/(5*fs).
%   at       Crosscorrelation type. 1=Biased (default), 2=Unbiased.
%   pf       Plot flag:  0=none (default), 1=screen.
%
%   ac       Estimated crosscorrelation.
%   t        Lag of estimates (sec).
%
%   This is sometimes called the crosscorrelation function, the 
%   sample crosscorrelation function, the normalized 
%   crosscorrelation.
%
%   The biased and unbiased estimators of crosscorrelation are 
%   related to those of the autocovariance by the simple relation
%
%      p(h) = ac(h)/(s1 s2)
%
%   where s1 and s2 are the sample standard deviations of x1 and x2.
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
%      x1 = abp(1:2000);
%      x2 = icp(1:2000); 
%      [ac,t] = CrossCorrelation(x1,x2,125,2,[],1);
%
%   M. B. Priestley, Spectral Analysis and Time Series, 
%   San Francisco, CA: Academic Press, 1981, pp. 692-693.
%
%   Version 2.00.00 JM
%
%   See also XCOV, CONV, and Correlation. 

%====================================================================
% Error Check
%====================================================================
if nargin<1,
    help CrossCorrelation;
    return;
    end
        
if var(x1)==0 | var(x2)==0,
    error('Signal is constant.');
    end;
    
%====================================================================
% Process Function Arguments
%====================================================================
fs = 1;                                 % Default sample rate
if exist('fsa') & ~isempty(fsa),
    fs = fsa;
    end; 
    
ml = round(length(x1)/5);                % Default sample rate
if exist('mla') & ~isempty(mla),
    ml = round(mla*fs);
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
nx = min(length(x1),length(x2));
ml = min(ml,nx-1);
mx1 = mean(x1);
mx2 = mean(x2);
vx1 = var(x1);
vx2 = var(x2);
x1  = x1(:);                 % Convert to a column vector
x2  = x2(:);                 % Convert to a column vector
x1  = x1-mx1;                % Remove the mean
x2  = x2-mx2;                % Remove the mean
cc = zeros(2*ml+1,1);        % Allocate memory
nc = length(cc);

%====================================================================
% Main Loop
%====================================================================
np = 2^nextpow2(2*nx-1);   % Figure out how much to pad the signal
X1 = fft(x1,np);
X2 = fft(x2,np);
cp = ifft(X1.*conj(X2));
cc = real([cp(np-(ml-1:-1:0));cp(1:ml+1)]);

if at==1,                  % Biased estimate
    cc = cc./nx;   
elseif at==2,              % Unbiased estimate
    k  = (-ml:ml).';
    cc = cc./(nx-abs(ml));
    end;    
    
cc = cc/sqrt(vx1*vx2); % This is the only real difference between crosscorrelation and crosscovariance

if 0, % Brute force confirmation that it's doing the right thing
    cc2 = zeros(2*ml+1,1);
    for c1 = -ml:ml,
        if c1<0,
            cc2(c1+ml+1) = sum(x1(1:nx-abs(c1)).*x2(1+abs(c1):nx));
        else
            cc2(c1+ml+1) = sum(x1(1+abs(c1):nx).*x2(1:nx-c1));
            end;
        end;
    cc2 = cc2/nx;
    cc2 = cc2/sqrt(vx1*vx2);
    end;

%====================================================================
% Postprocessing
%====================================================================    
t = (-ml:ml).'/fs;          % Lags in units of samples

%====================================================================
% Plotting
%====================================================================
if pf==1,   
    figure;
    FigureSet;
    h = plot(t,cc);
    set(h,'LineWidth',1.2);
    if ml<100,
        set(h,'Marker','.');
        set(h,'MarkerSize',12);
        end;    
    hold on;
        if at==1,
            h = plot(t, (nx-abs(t)*fs)/nx,'g');
            h = plot(t,-(nx-abs(t)*fs)/nx,'g');
            end;
        hold off;
    xlim([min(t) max(t)]);
    ylim([-1 1]);
    AxisLines;
    box off;
    zoom on;
    xlabel('Lag (s)');
    ylabel('\rho');
    title('CrossCorrelation');
    AxisSet(8);
    end

if nargout==0,
    clear cc;
    clear t;
    end



