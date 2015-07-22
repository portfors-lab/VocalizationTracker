function [icc,ci,s2b,s2w] = IntraclassCorrelation(x,losa,pfa);
%IntraclassCorrelation: Calculates intraclass correlation coefficient
%
%   [icc,ci,s2b,s2w] = IntraclassCorrelation(x,los,pf);
%
%   x    Input data (cell array).
%   los  Level of significance. Default = 0.05.
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   icc  Scalar containing the intraclass correlation coefficient.
%   ci   Confidence interval (lower,upper).
%   s2b  Variance between classes. 
%   s2w  Variance within classes. 
%
%   This is function calculates the intraclass correlation 
%   coefficient using a one-way random effects model as described
%   by McGraw and Wong. This method is more general than that 
%   described in the paper, however, because it does not require
%   the same number of observations for each class. 
%
%   This implementation ignores elements that contain only one
%   measurement. 
%
%   The argument matrix must be a cell array. Each element of the
%   cell must contain a vector of measurements for each class.
%
%   Example: Generate the ICC on random data.
%
%      x = cell(3,1);
%      x{1} = randn(5,1);
%      x{2} = randn(3,1)+1;
%      x{3} = randn(4,1)-1;
%      ICC = IntraclassCorrelation(x,1);
%
%   McGraw, Kenneth O., Wong, S.P., "Forming Inferences About Some 
%   Intraclass Correlation Coefficients," Psychological Methods, 
%   Vol. 1, No. 1, 30-46, 1996.
%
%   Version 1.00 JM
%
%   See also anova1.

%====================================================================
% Error Checking
%====================================================================    
if nargin<1,
    help IntraclassCorrelation;
    return;
    end;
    
if size(x,1)<2,
    error('Cell array must contain at least 2 elements (classes).');
    end;
    
%====================================================================
% Calculate Basic Signal Statistics
%====================================================================    
nx = size(x,1);
    
%====================================================================
% Process Function Arguments
%====================================================================
los = 0.05;
if exist('losa') & ~isempty(losa) & losa>0 & losa<1,
    los = losa;
    end;    

pf = 0;                                 % Default - no plotting
if nargout==0,    % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
    
%====================================================================
% Build up cell array that only contains elements with repeated 
% measurements
%====================================================================   
y = {};                                                    % Declare new cell array for repeated observations
ny = 0;                                                    % Number of elements in cell array (r in ALSM)
for c1=1:nx,
    if iscell(x),
        z = x{c1};
    else
        z = x(c1,:);
        end;    
    if length(z)>1,
        ny     = ny + 1;
        z      = z(:);                                     % Convert to column vector, if not already  
        y{ny}  = z;
        end;
    end;

if ny<2,
    error('Cell array must contain at least 2 elements with repeated measures.');
    end;
    
%y = x;
%ny = nx;

%====================================================================
% Main Loop
%====================================================================
gmy = 0;                                                   % Grand mean of y 
no  = zeros(ny,1);                                         % Vector of the number of observations for each class (ni in ALSM)
my  = zeros(ny,1);                                         % Class means (Yi. in ALSM)
for c1=1:ny,
    no(c1) = length(y{c1});                                % Number of observations for class c1
    my(c1) = mean(y{c1});                                  % Mean of observations for class c1
    gmy    = gmy + sum(y{c1});                             % Intermediate step to calculate grand mean of y
    end;
nt  = sum(no);                                             % Total number of observations (nT in ALSM)
gmy = gmy/nt;                                              % Calculate the grand mean of all observations (Y.. in ALSM)

nyp = (sum(no)-sum(no.^2)/sum(no))/(ny-1);                 % The equivalent expression for number of observations when the number of observations is not equal for all classes (see 24.10a in ALSM, pg. 966)

SST = 0;                                                   % Sum of squares total (SSTO in ALSM)
SSR = 0;                                                   % Sum of squares between classes (rows) (SSTR in ALSM)
SSW = 0;                                                   % Sum of squares within classes or rows (SSE in ALSM)
for c1=1:ny,
    SST = SST + sum((y{c1}-gmy).^2);                       % Update SST
    SSR = SSR + no(c1)*(my(c1)-gmy).^2;                    % Update SSR
    SSW = SSW + sum((y{c1}-my(c1)).^2);                    % Update SSW 
    end;

if abs(SSR+SSW-SST)/SST>100*eps,                           % Error checking on partitioning of SST
    error(sprintf('abs(SSR+SSW-SST)/SST:%f too big. Check code for bug.',abs(SSR+SSW-SST)/SST));
    end;
    
msr = SSR/(ny-1);                                          % Mean sum of squares between rows (MSTR in ALSM)
msw = SSW/(nt-ny);                                         % Mean sum of squares within rows (MSE in ALSM)

s2b = max(0,(msr-msw)/nyp);                                % Estimate of variance between class (row) means (see 24.23 in ALSM, pg. 971)
s2w = msw;                                                 % Estimate of variance within class (row)

icc = s2b/(s2b + s2w);                                     % Calculate the ICC

msb = msr;

%====================================================================
% Calculate Confidence Interval
%====================================================================
Lr = (1/nyp)*((msr/msw)*(1/finv(1-los/2,ny-1,ny*(nyp-1)))-1); % Lower confidence interval on s2r/s2w
Ur = (1/nyp)*((msr/msw)*(1/finv(  los/2,ny-1,ny*(nyp-1)))-1); % Upper confidence interval on s2r/s2w

ci = zeros(2,1);
ci(1) = Lr/(1+Lr);
ci(2) = Ur/(1+Ur);
ci(1) = max(0,ci(1));                                      % Make lower confidence interval is no less than zero
ci(2) = min(1,ci(2));                                      % Make upper confidence interval is no less than zero

%====================================================================
% Plot Results
%====================================================================  
if pf>=1,
    figure;
    FigureSet(1);
    [jnk,is] = sort(my);
    ymax = 0;
    ymin = inf;
    for c1=1:ny,
        id = is(c1);
        hold on;    
        if no(id)>1,
            h = plot([c1 c1],[min(y{id}) max(y{id})],'k');
            end;
        h = plot(c1,my(id),'ks');
        set(h,'MarkerSize',5);
        set(h,'MarkerFaceColor','k');
        if no(id)>1,
            h = plot(c1,y{id},'ro');
            set(h,'MarkerSize',3);
            set(h,'MarkerFaceColor','r');        
            end;
        ymax = max([ymax max(y{id})]);
        ymin = min([ymin min(y{id})]);
        end;
    yrg = ymax-ymin;
    ylim([ymin-0.025*yrg ymax+0.025*yrg]);
    xlim([1 ny]);
    h = plot(xlim,gmy*[1 1],'b');
    hold off;
    xlabel('Ordered Class Index');
    ylabel('Data (units?)');
    title(sprintf('ICC:%5.3f (%5.3f-%5.3f)  gm:%4.2f  sb:%4.2f  sw:%4.2f  sw/gm:%4.2f',...
          icc,ci(1),ci(2),gmy,sqrt(s2b),sqrt(s2w),sqrt(s2w)/gmy));
    box off;
    AxisSet;
    end;

%====================================================================
% Process Return Arguments
%====================================================================
if nargout==0,
    clear('icc','ci','msw','msb');
    end;

%====================================================================
% Example to Validate Correct Calculation
%====================================================================
% d = randn(3,50);
% d(1,:) = d(1,:) + 3;
% d(2,:) = d(2,:);
% d(3,:) = d(3,:) - 3;
% 
% x = cell(3,1);
% x{1} = d(1,:);
% x{2} = d(2,:);
% x{3} = d(3,:);
% 
% ICC = IntraclassCorrelation(x,1);
% 
% k = size(d,2);
% n = size(d,1);
% dm = mean(d,2);
% gm = mean(mean(d));
% SSR = sum(k*(dm-gm).^2);
% MSR = SSR/(n-1);
% SSW = 0;
% for c1=1:n,
%     SSW = SSW + sum((d(c1,:)-dm(c1)).^2);
%     end;
% MSW = SSW/(n*(k-1));
% 
% ICC2 = (MSR-MSW)/(MSR+(k-1)*MSW);
% 
% disp([SSR SSW MSR MSW ICC2]);
% 
% fprintf('ICC:%5.3f ICC2:%5.3f\n',ICC,ICC2);

%return;






