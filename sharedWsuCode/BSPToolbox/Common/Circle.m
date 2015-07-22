function [h] = Circle(rda,sta)
%Circle: Draws a circle around the origin with specified radius
%
%   [h] = Circle(rd,st)
%
%   rd   Radius. Default = 1.
%   st   Style of circle. Default = 'k:' (Black dotted line).
%
%   h    Handle to line object.
%
%   Plots a circle on the current plot around the origin. This is
%   useful for pole-zero plots.
%
%   Example: Generate a scatter plot of a complex random variable 
%   and add the unit circle to it. Use a red dashed line for the
%   circle and set the linewidth to 2.0. Note that the aspect 
%   ratio must be one or the circle will appear to be an oval.
%
%      x   = 1/sqrt(2)*(randn(50,1) + j*randn(50,1));
%      figure;
%      FigureSet; 
%      plot(real(x),imag(x),'.');
%      h = Circle([],'r--');
%      set(h,'LineWidth',2.0);
%      axis(1.2*[-1 1 -1 1]);
%      axis('equal');
%
%   Desmond J. Higham and Nicholas J. Higham, "MATLAB Guide," 2000.
%
%   Version 1.00 JM
%
%   See also AxisLines, grid, polar, axis, and pzmap.

%====================================================================
% Process function arguments
%====================================================================
rd = 1;                               % Default sampling rate, Hz
if exist('rda') & ~isempty(rda),
    rd = abs(rda);
    end;
   
st = 'k:';    
if exist('sta') & ~isempty(sta),
    st = sta;
    end;
    
%====================================================================
% Plot
%====================================================================
hf = ishold; % Hold flag

if ~hf,
    hold on;
    end;

th = (0:0.01:1)*2*pi;
r  = rd*ones(size(th));
h  = polar(th,r,st);
set(h,'Clipping','Off');

if ~hf, 
	hold off;
	end;
    
%====================================================================
% Process Output Arguments
%====================================================================      
if nargout==0,
    clear h;
    return;
    end;    
