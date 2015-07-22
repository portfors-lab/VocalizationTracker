function [h] = PixelPlot(x,y,I,lma,cma,nfa);
%PixelPlot: Generates image with nonuniform pixel size.
%
%   [h] = PixelPlot(x,y,I,lm,cm)
%
%   x    Vector of horizontal locations of pixels.
%   y    Vector of vertical locations of pixels.
%   I    Matrix containing the values of the image.
%   lm   Vector with 2 elements specifing the smallest and largest
%        values to map to the minimum and maximum color map 
%        indices. Default = [Imin Imax]. 
%   cm   Colormap to use. Default = figure's colormap. 
%   nf   New figure flag. 1=New figure (default), 0=current figure.
%
%   h    Handle to all of the patches that compose the pixelplot.
%
%   This is designed to be an alternative to MATLAB's image and 
%   imagesc functions which 1) require uniform horizontal and 
%   vertical spacing of pixels and 2) use the figure's colormap.
%   Surf can also be used for this purpose, but the location of
%   pixels with this command is difficult to control and understand.
%
%   For the element I(i,j), PixelPlot colors a patch of the image
%   with a horizontal range from the midpoint between x(j-1),x(j) to 
%   the midpoint between x(j),x(j+1). Similarly, the vertical range
%   is from the midpoint between y(i-1),y(i) to the midpoint between
%   y(i),y(i+1). The left edge of the left-most pixel is chosen such 
%   that that x(1) is located in the center of the pixel. The right 
%   edge is chosen in the same manner. The vertical pixel edges are 
%   chosen in the same manner as the horizontal pixel edges. Note 
%   that the rows of I correspond to the elements of the y vector
%   and the columns of I correspond to the elements of the x vector.
%   Thus, the matrix represents a vertically inverted version of the
%   image.
%
%   One consequence of chosing pixel locations in this way is that 
%   the pixels are not necessarily centered about the pixel 
%   locations. Another consequence is that pixels in regions of
%   high pixel density comprise a smaller portion of the image than
%   pixels in low pixel density. Thus, a larger portion of the image
%   area is dedicated to displaying I where the pixel density is the
%   lowest. For data visualization applications, this may be exactly
%   the opposite of what is desired. 
%
%   Example: Plot the morphology of QRS complexes over a segment of
%   an electrocardiogram.
%
%      load MER.mat;
%      si = si(1:50);
%      di = -round(2e-3*fs):round(2e-3*fs);
%      id = ones(length(di),1)*si' + (di)'*ones(1,length(si));
%      I  = x(id);
%      PixelPlot((si-1)/fs,di,I); 
%
%   R. L. Harris, "Information Graphics. A comprehensive illustrated
%   reference," Oxford University Press, 1999.
%
%   Version 1.00 JM
%
%   See also specgram, window, decimate, Beatogram, and Cohereogram.

%====================================================================
% Error Checking
%====================================================================    
if nargin<3,
    help PixelPlot;
    return;
    end;

%====================================================================
% Process Function Arguments
%====================================================================
lm = [min(min(I)) max(max(I))];
if exist('lma') & ~isempty(lma),
    lm = lma;
    end;    
    
nf = 1;
if exist('nfa') & ~isempty(nfa),
    nf = nfa;
    end;     
    
if nf,
    figure;
    end;

if exist('cma','var') & ~isempty(cma),
    cm = cma;
else
    cm = colormap;    
    end;     
    
%====================================================================
% Preprocessing
%====================================================================
ny = size(I,1);
nx = size(I,2);

cmax = lm(2);
cmin = lm(1);
nc   = size(cm,1);

X  = zeros(4,ny*nx);
Y  = zeros(4,ny*nx); 
%C  = zeros(1,ny*nx); 
C  = zeros(1,ny*nx,3);
cp = 0;

for c1=1:nx,
    if c1==1,
        xl   = x(c1) - (x(c1+1)-x(c1  ))/2;
        xr   = x(c1) + (x(c1+1)-x(c1  ))/2;
        xmin = xl;
    elseif c1==nx,
        xl   = x(c1) - (x(c1  )-x(c1-1))/2;
        xr   = x(c1) + (x(c1  )-x(c1-1))/2;
        xmax = xr;
    else
        xl   = x(c1) - (x(c1  )-x(c1-1))/2;
        xr   = x(c1) + (x(c1+1)-x(c1  ))/2;
        end;
    for c2=1:ny
        if c2==1,
            yb   = y(c2) - (y(c2+1)-y(c2  ))/2; 
            yt   = y(c2) + (y(c2+1)-y(c2  ))/2;
            ymin = yb;
        elseif c2==ny,
            yb   = y(c2) - (y(c2  )-y(c2-1))/2; 
            yt   = y(c2) + (y(c2  )-y(c2-1))/2;
            ymax = yt;
        else
            yb   = y(c2) - (y(c2  )-y(c2-1))/2; 
            yt   = y(c2) + (y(c2+1)-y(c2  ))/2;
            end;
        cid = round(nc*(I(c2,c1)-cmin)/(cmax - cmin));
        cid = min(max(1,cid),nc);
        
        cp = cp + 1;
        
        X(:,cp)   = [xl;xr;xr;xl];
        Y(:,cp)   = [yb;yb;yt;yt];
        %C(:,cp)   = cid;       
        C(1,cp,:) = cm(cid,:);       
        end;
    end;

FigureSet(1);
h = patch(X,Y,C);
set(h,'LineStyle','None');
drawnow;
xlim([xmin xmax]);
ylim([ymin ymax]);