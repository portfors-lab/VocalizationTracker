function [b,e] = MAInnovations(x,nb,pfa)
%MAInnovations: Estimate the moving average process parameters.
%
%   [b,e] = MAInnovations(x,nb,pf)
%
%   x    Input signal.
%   nb   Model order (number of zeros). 
%   pf   Plot flag:  0=none (default), 1=screen.
%
%   b    Vector of model coefficients.
%   e    Estimated variance of white noise driving the process.
%
%   Uses the innovations algorithm to do a preliminary estimation
%   of a moving average model,
%
%      x(t) = z(t) + b(1) z(t-1) + ... + b(nb) z(t-nb)
%
%   where z(t) is a white noise process with zero mean and an 
%   estimated variance of e. 
% 
%   If pf=1, an image is shown of the model coefficients and 
%   estimated noise variance for model orders ranging from 1 to nb. 
%   This is a useful visualization technique for picking the model 
%   order. A good choice is to pick the smallest value of nb such 
%   that the model coefficients are not much different for larger 
%   values.
%
%   Example: Plot the estimated spectrum of a 25th order MA model 
%   of an ABP signal.
%
%      load ABPICP.mat
%      x = abp(1:2000);
%      [b,e] = MAInnovations(x,25);
%      freqz(b,1,1024,fs);
%      FigureSet;
%
%   M. B. Priestley, Spectral Analysis and Time Series, 
%   San Francisco, CA: Academic Press, 1981, pp. 245-250.
%
%   Version 2.00.22 JM
%
%   See also SIGNAL, IDENT, and Models. 

%====================================================================
% Error Check
%====================================================================
if nargin<1,
    help MAInnovations;
    return;
    end
        
if var(x)==0,
    error('Signal is constant.');
    end;
    
%====================================================================
% Process Function Arguments
%====================================================================        
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
B  = zeros(nb,nb);         % Model coefficients at various estimation stages
mv = zeros(nb,1);          % Estimated white noise variance

%====================================================================
% Main Loop
%====================================================================
ac     = Autocovariance(x,1,nb+1);
mv(1)  = ac(1);            % First estimate of noise power

for c1 = 1:nb,                        % This represents m in the reference 
    B(c1,c1) = inv(mv(1))*ac(c1+1);   % k = 0
    for c2 = 1:c1-1,                  % This represents k in the reference
        j           = 0:c2-1;
        B(c1,c1-c2) = inv(mv(c2+1))*(ac(c1-c2+1) - sum(B(c1,c1-j).*B(c2,c2-j).*mv(j+1).'));
        end;
    j        = 0:(c1-1);
    mv(c1+1) = ac(1) - sum(B(c1,c1-j).^2.*mv(j+1).');
    end;

%====================================================================
% Post Processing
%====================================================================
b = [1;B(nb,1:nb).'];    
e = mv(nb+1);

%====================================================================
% Plotting
%====================================================================
if pf==1,    
    figure;
    FigureSet;
    ha1 = axes('Position',[0.08 0.08 0.76 0.84]);
    imagesc(B);
    xlabel('Coefficient');
    ylabel('Iteration (Model Order)');
    title('Moving Average Coefficients and Noise Power versus Model Order');
    box off;
    
    ha2 = axes('Position',[0.86 0.08 0.05 0.84]);
    imagesc(mv(2:nb+1));
    title('s^2');
    set(ha2,'Visible','Off');
    box off;
    
    ha3 = axes('Position',[0.93 0.08 0.05 0.84]);
    colorbar(ha3,'peer',ha1);
    set(ha3,'Visible','Off');
    box off;
    
    AxisSet;
    end

if nargout==0,
    clear b;
    clear e;
    end



