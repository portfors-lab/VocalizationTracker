function im = DetectMinima(x,thrsa,pltfa,pfa)
%DetectMinima: Detect minima in the input signal x
%
%   im = DetectMinima(x,thrs,pltf,pfa)
%
%   x      Input signal       
%   thrs   Threshold value (optional)
%   pltf   Plateau flag: 0=count plateaus (default),
%                        1=only find distinct minima (default)
%   pf     Plot flag: 0=none (default), 1=screen
%
%   y      Location (index) of the mimima  
%
%   Finds the location (index) of all the minima in the 
%   input signal x and returns them in a sorted vector (MIN). 
%   If the thrs input is used, the funtion returns those peaks that
%   are smaller than thrs. 
%
%   Example: Detect all the minima in the ICP signal provided with the 
%   bsp toolbox data (ICP.m) and plot the results.
%
%      load ICP; 
%      DetectMinima(icp);
%
%   Version 1.0 MA
%
%   See also DetectMaxima, ManualDetector, and PressureDetector.


%=========================================================================
% Process function arguments
%=========================================================================
if nargin<1 | nargin>3,
    help DetectMinima;
    return;
    end;

thrs = min(x)-1;                                           % Default threshold
if exist('thrsa') & ~isempty(thrsa),
    thrs = thrsa;
    end;

pltf = 1;                                                  % Default threshold
if exist('pltfa') & ~isempty(pltfa),
    pltf = pltfa;
    end;
    
pf = 0;                                                    % Default - no plotting
if nargout==0,                                             % Plot if no output arguments
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
    id = [];
elseif pltf==0,                                                % Count plateaus
    id = find((x(1:nx-2)>=x(2:nx-1))&(x(2:nx-1)<=x(3:nx)))+1;
    if nx>1 & x(1)<=x(2),
        id = [1;id];
        end;
    if nx>1 & x(nx)<=x(nx-1),
        id = [id;nx];
        end;
else                                                       % Only find distinct minima (no plateaus)
    id = find(((x(1:nx-2)>=x(2:nx-1))&(x(2:nx-1)<x(3:nx))) | ((x(1:nx-2)>x(2:nx-1))&(x(2:nx-1)<=x(3:nx))))+1;
    if nx>1 & x(1)<x(2),
        id = [1;id];
        end;
    if nx>1 & x(nx)<x(nx-1),
        id = [id;nx];
        end;
    end;
    
im = id(x(id)<thrs);                                       % Only include maxima that are above the threshold 

%=========================================================================
% Plotting  
%=========================================================================
if pf ==1,
    figure;
    FigureSet(1,'LTX');
    plot(k, x, 'b', k(im), x(im), 'r.')
    title('Detect Minima')
    xlabel('Samples')
    ylabel('Magnitude')
    end;

%=========================================================================
% Take care of outputs
%=========================================================================
if nargout==0,
	clear('im');
	end;   