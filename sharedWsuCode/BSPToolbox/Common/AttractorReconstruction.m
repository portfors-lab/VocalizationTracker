function [Y] = AttractorReconstruction(x,ta,na,pfa) 
%AttractorReconstruction: Attractor reconstruction from data
%
%   [Y] = AttractorReconstruction(x,t,n,pf);        
%
%   x    Input signal.
%   t    Delay (samples). Default = 5 samples.
%   n    Embedding dimension. Default = 18
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   Y    Attractor reconstruction / delay plot. 
%
%   Generates a trayectory in an n-dimensional space from 
%   a time series. The attractor reconstruction can be
%   visualized in 2 or 3 dimensions as a delay plot. The 
%   matrix Y contains the n-dimensional embedding of the 
%   attractor, where every row of Y is a vector of delay 
%   coordinates:
%  
%    Y(i,:) = [x(t), x(t-T), x(t-2T), ... , x(t-(n-1)T)]
%
%   Example: Generate the attractor reconstruction for a 
%   10 s. segment of ICP data. The data is first smoothed 
%   with a lowpass filter and detrended to ensure a proper
%   attractor reconstruction.
%       
%       load ICP;
%       icplp = LowPass(icp,fs,10,1,2,0);
%       icps  = HighPass(icplp(10^4:10^4+10*fs),fs,0.6,2,0);
%       Y     = AttractorReconstruction(icps,5,10,1);
%
%   K. T. Alligood, T. D. Sauer, J. A. Yorker, CHAOS: An
%   Introduction to Dynamical Systems. 
%   Springer-Verlag, 1996, pp.537-553.
%
%   Version 0.00.00.21 MA
%
%   See also CorrelationDimension.

%=========================================================
% Process function arguments
%=========================================================
if nargin<1 | nargin>4,
    help AttractorReconstruction;
    return;
    end;

t = 5;                 % Default delay = 5 samples
if exist('ta') & ~isempty(ta),
    t = ta;
    end;
    
n = 18;                % Default n = 12 dimensions
if exist('na') & ~isempty(na),
    n = na;
    end;
    
pf = 0;                % Default - no plotting
if nargout==0,         % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
    
%=========================================================
% Attractor Reconstruction
%=========================================================
ld = length(x);
Y  = zeros(ld-(n-1)*t,n);   

for i=1:n;
    Y(:,i) = x((n-i)*t+1:ld-(i-1)*t);
end;

%=========================================================
% Plotting
%=========================================================
if pf==1,
    %figure;
    %FigureSet(1) %'wide');
    %subplot(1,2,1)
    plot(Y(:,1),Y(:,2));
    title('Attractor Reconstruction (Delay Plot)');
    xlabel('x(t)');
    ylabel('x(t-T)');
    box off;
    zoom on;
    AxisSet;
    
   
    %subplot(1,2,2)
    %plot3(Y(:,1),Y(:,2), Y(:,3));
    %title('Attractor Reconstruction (Delay Plot)');
    %xlabel('x(t)');
    %ylabel('x(t-T)');
    %zlabel('x(t-2T)');
    %box on;
    %grid on;
    %zoom on;
    %AxisSet;
end

%====================================================================
% Process Return Arguments
%====================================================================
if nargout==0,
    clear('Y');
    end;