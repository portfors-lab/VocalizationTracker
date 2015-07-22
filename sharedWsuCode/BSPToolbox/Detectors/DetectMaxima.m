function [im] = DetectMaxima(x,thrsa,pltfa,pfa)
%DetectMaxima: Detect maxima (peaks) in the input signal x
%
%   im = DetectMaxima(x,thrs,pltf,pfa)
%
%   x      Input signal       
%   thrs   Threshold value (optional)
%   pltf   Plateau flag: 0=count plateaus (default),
%                        1=only find distinct maxima 
%   pf     Plot flag: 0=none (default), 1=screen
%  
%   im     Location (index) of all maxima in x
%
%   Finds the location (index) of all the maxima (peaks) in the 
%   input signal x and returns them in a sorted vector (y). 
%   If the thrs input is used, the funtion returns those peaks that
%   are greater than thrs. 
%
%   Example: Detect all the peaks in the ICP signal provided with the 
%   bsp toolbox data (ICP.m) and plot the results.
%
%      load ICP; 
%      DetectMaxima(icp);
%
%   Version 1.0 MA
%
%   See also DetectMinima, ManualDetector, and PressureDetect.

%=========================================================================
% Process function arguments
%=========================================================================
if nargin<1 | nargin>4,
    help DetectMaxima;
return;
end;

thrs = min(x)-1;                                           % Default threshold
if exist('thrsa') & ~isempty(thrsa),
    thrs = thrsa;
end;

pltf = 0;                             % Default threshold
if exist('pltfa') & ~isempty(pltfa),
    pltf = pltfa;
    end;
    
pf = 0;                                    % Default - no plotting
if nargout==0,                             % Plot if no output arguments
    pf = 1;
    end;  

if exist('pfa') & ~isempty(pfa),
    pf = pfa;
    end;
    
%=========================================================================
% Detect Maxima
%=========================================================================
x  = x(:);
nx = length(x);
k  = 1:nx;

if nx==1,
    id = [1];
elseif pltf==0,                                            % Count plateaus
    id = find((x(1:nx-2)<=x(2:nx-1))&(x(2:nx-1)>=x(3:nx)))+1;
    if x(1)>=x(2),
        id = [1;id];
        end;
    if x(nx)>=x(nx-1),
        id = [id;nx];
        end;
else                                                       % Only find distinct maxima (no plateaus)
    id = find(((x(1:nx-2)<=x(2:nx-1))&(x(2:nx-1)>x(3:nx))) | ((x(1:nx-2)<x(2:nx-1))&(x(2:nx-1)>=x(3:nx))))+1;
    if x(1)>x(2),
        id = [1;id];
        end;
    if x(nx)>x(nx-1),
        id = [id;nx];
        end;
    end;
    
im = id(x(id)>thrs(1));                                       % Only include maxima that are above the threshold 

%=========================================================================
% Plotting  
%=========================================================================
if pf ==1,
    figure;
    FigureSet(1);
    plot(k, x, 'b', k(im), x(im), 'r.');
    title('Detect Maxima');
    xlabel('Samples');
    ylabel('Magnitude');
    end;

%=========================================================================
% Take care of outputs
%=========================================================================
if nargout==0,
	clear('im');
    end;   