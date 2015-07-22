function [y] = Hyperbola(x,ymaxa,pfa);
%Hyperbola: Generates a hyperbolic window function
%
%   [y] = Hyperbola(x);
%
%   x    If scalar, window length. If vector, indices of window.
%   ymax Maximum value of the window amplitude. Default = 0.95.
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   y    Window.
%
%   This function returns a vector that represents a hyperbolic 
%   window function. Visually, this is very similar to a triangular
%   or tent window function. However, the hyperbola is analytic
%   (all it's derivatives exist at all points) and has a rounded
%   peak. The parameter ymax controls how rounded the peak is.
%   This function is used in the ColorSpiral colormap to prevent
%   a discontinuity at the midpoint of the colormap.
%
%   Example: Generate the spectrogram of an intracranial pressure
%   signal using a Hyperbola window that is 45 s in duration.
%
%      load ICP.mat; 
%      icpd = decimate(icp,15);
%      wl   = round(45*fs/15);
%      Spectrogram(icpd,fs/15,Hyperbola(wl));
%
%   C. H. Edwards, D. E. Penney, "Calculus and Analytic Geometry,"
%   2nd edition, Prentice-Hall, 1986.
%
%   Version 1.00 JM
%
%   See also triang, window, and ColorSpiral.

%   See http://mathworld.wolfram.com/Hyperbola.html for details. 

%====================================================================
% Error Checking
%====================================================================    
if nargin<1,
    help Hyperbola;
    return;
    end;
    
%====================================================================
% Process Function Arguments
%====================================================================     
if length(x)==1,                                           % If is an integer, make it into an array
    x = 1:x;
    end;

ymax = 0.95;
if exist('ymaxa') & ~isempty(ymaxa),
    ymax = ymaxa;
    end;    

pf = 0;                                                    % Default - no plotting
if nargout==0,                                             % Plot if no output arguments
    pf = 1;
    end;  
if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
    
%====================================================================
% Preprocessing
%====================================================================   
a    = sqrt((1-ymax).^2/(1-(1-ymax).^2));
xmin = min(x);
xmax = max(x);
xs   = 2*(x-xmin)/(xmax-xmin) - 1;                         % Scale so it ranges from -1 to 1
nx   = length(x);

%====================================================================
% Main Routine
%====================================================================   
y        = 1-sqrt(xs.^2+a^2)/sqrt(1+a^2);
y(y<0)   = 0;

ya       = ones(nx,1);
ya(xs<0) =  xs(xs<0)+1;
ya(xs>0) = -xs(xs>0)+1;

%====================================================================
% Postprocessing
%====================================================================   
y      = y(:);                                             % Convert into a column vector

%====================================================================
% Plot Default Figure
%====================================================================
if pf,
    figure;
    FigureSet(1);
    h = plot(x,y);
    set(h,'LineWidth',1.5);
    hold on;
        h = plot(x,ya,'r');
        set(h,'LineWidth',0.5);
        h = plot((xmax-xmin)/2 + xmin,ymax,'k.');
        set(h,'MarkerSize',15);
        hold off;
    xlim([xmin xmax]);
    ylim([min(y) 1.03]);
    box off;
    AxisSet;
    end;
    
%====================================================================
% Process Return Arguments
%====================================================================
if nargout==0,
    clear('y');
    end;    