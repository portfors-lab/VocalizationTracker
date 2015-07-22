function [a,b,e] = ARMABrockwell(x,na,nb,pfa)
%ARMABrockwell: Preliminary estimation of an ARMA process.
%
%   [a,b,e] = ARMABrockwell(x,na,nb,pf)
%
%   x    Input signal.
%   na   Denominator order (no. of poles)
%   nb   Numerator order (no. of zeros)
%   pf   Plot flag:  0=none (default), 1=screen.
%
%   a    Vector of model denominator coefficients.
%   b    Vector of model numerator coefficients.
%   e    Estimated variance of white noise driving the process.
%
%   Uses the method described in the reference to build a preliminary 
%   model of an ARMA process. The model is not guaranteed to be 
%   stable and is recommended only for low order models (e.g. na =
%   nb = 1). The model is consistent in probability. 
%
%   Example: Plot the estimated spectrum of a 3rd order ARMA model 
%   of an ABP signal.
%
%      load ABPICP.mat
%      x = abp(1:2000);
%      [a,b,e] = ARMABrockwell(x,3,3);
%      freqz(b,a,1024,fs);
%      FigureSet;
%
%   M. B. Priestley, Spectral Analysis and Time Series, 
%   San Francisco, CA: Academic Press, 1981, pp. 245-250.
%
%   Version 0.00.00.23 JM
%
%   See also SIGNAL, IDENT, and Models. 

%====================================================================
% Error Check
%====================================================================
if nargin<1,
    help ARMABrockwell;
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
[h,e] = MAInnovations(x,nb+na+15); % Estimated impulse response of ARMA process

A = zeros(na,na);
for c1 = 1:na,
    id = (nb+c1) + (0:-1:-(na-1));
    nz = find(id>0);
    id = id(nz);
    A(c1,nz) = h(id).';
    end;   
p = h(nb+1+(1:na));

ap = pinv(A)*p;            % Calculate partial denominator coefficients

bp = zeros(nb,1);          % Allocate memory for partial numerator coefficients
for c2 = 1:nb,
    id = 1:min(c2,na);
    bp(c2) = h(c2+1) - sum(ap(id).*h(c2-id+1));
    end;
    
%====================================================================
% Postprocessing
%====================================================================      
a = [1;-ap];              % Append 1 to final coefficients   
b = [1; bp];              % Append the 1

%====================================================================
% Plotting
%====================================================================
if pf==1,    
    figure;
    FigureSet;
    z = roots(b);
    p = roots(a);
    h = plot(real(z),imag(z),'bo',real(p),imag(p),'rx');
    set(h,'LineWidth',1.5);
    set(h,'MarkerSize',10);
    Circle;
    AxisLines;
    amx = max(abs([1.02;real(z);imag(z);real(p);imag(p)]));
    xlim([-amx amx]);
    ylim([-amx amx]);
    axis square;
    xlabel('Real');
    ylabel('Imaginary');
    title('Pole-Zero Plot of Model Coefficients');
    AxisSet;
    end

if nargout==0,
    clear a;
    clear b;
    clear e;
    end



