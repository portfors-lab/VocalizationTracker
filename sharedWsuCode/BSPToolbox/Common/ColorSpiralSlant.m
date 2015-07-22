function [mp] = ColorSpiral(nca,npa,pfa);
%ColorSpiral: Generates a monotonic colormap with maximum color depth
%
%   [m] = ColorSpiral(n,np,pf);
%
%   nc   Number of colors (length of the colormap). Default = 64.
%   np   Number of sinusoidal periods. Default = 2.
%   pf   Plot flag: 0=none (default), 1=screen.
%
%   m    Color map.
%
%   This function returns an n x 3 matrix containing the RGB entries
%   used for colormaps in MATLAB figures. The colormap is designed
%   to have a monotonically increasing intensity, while maximizing
%   the color depth. This is achieved by generating a spiral through
%   the RGB cube that ranges from RBG = [0 0 0] to RGB = [1 1 1].
%
%   Example: Generate the spectrogram of an intracranial pressure
%   signal using a Blackman-Harris window that is 45 s in duration.
%
%      load ICP.mat; 
%      icpd = decimate(icp,15);
%      wl   = round(45*fs/15);
%      Spectrogram(icpd,fs/15,45);
%      colormap(ColorSpiral(256));
%
%   C. Ware, "Information Visualization: Percetpion for Design," 2nd edition,
%   Morgan Kaufmann, 2004.
%
%   Version 1.00 JM
%
%   See also DistinctColors, colormap, jet, and caxis.

%====================================================================
% Error Checking
%====================================================================    
if nargin>3,
    help ColorSpiral;
    return;
    end;
    
%====================================================================
% Process Function Arguments
%====================================================================    
nc = 64;                                                   % Number of colors in the colormap
if exist('nca') & ~isempty(nca),
    nc = nca;
    end;

np = 2;                                                    % Number of sinusoidal periods
if exist('npa') & ~isempty(npa),
    np = npa;
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
%wn = sqrt(3/8)*[0;triang(nc-2);0];
wn = sqrt(3/8)*Hyperbola(nc);
a12 = asin(1/sqrt(3));                                     % First  rotation angle
a23 = pi/4;                                                % Second rotation angle

T = inv([1.0 0.956 0.621; 1.0 -0.272 -0.647; 1.0 -1.106 1.703]);
luminanceHSV = T(1,:)';
luminanceSum = [1 1 1];

[thetaHSV,phiHSV,r] = cart2sph(luminanceHSV(1),luminanceHSV(2),luminanceHSV(3));
[thetaSum,phiSum,r] = cart2sph(luminanceSum(1),luminanceSum(2),luminanceSum(3));

thetaRotate = thetaHSV - thetaSum;
phiRotate   = phiHSV   - thetaSum;

%====================================================================
% Main Routine
%====================================================================    
r0 = linspace(sqrt(3),0,nc).';
g0 = cos(((r0-sqrt(3)/2)*np*2*pi/sqrt(3)));
b0 = sin(((r0-sqrt(3)/2)*np*2*pi/sqrt(3)));

theta  = pi/2*ones(nc,1);
radius = 1   *ones(nc,1);
phi    = ((r0-sqrt(3)/2)*np*2*pi/sqrt(3));

%[theta,phi,R] = cart2sph(zeros(nc,1),g0,b0);
thetaRotated = theta + thetaRotate;
phiRotated   = phi   + phiRotate;

[r0a,g0a,b0a] = sph2cart(thetaRotated,phiRotated,radius);
[thetaB,phiB,radiusB] = cart2sph(r0a,g0a,b0a);
thetaC = thetaB - thetaRotate;
phiC   = phiB   - phiRotate;
radiusC = radiusB;
[r0c,g0c,b0c] = sph2cart(thetaC,phiC,radiusC);

r0 = r0 + r0a;
g0 = wn.*g0a;
b0 = wn.*b0a;

[ag,rd] = cart2pol(r0,g0);
[r1,g1] = pol2cart(ag+a12,rd);            
b1      = b0;

[ag,rd] = cart2pol(r1,b1);
[r2,b2] = pol2cart(ag+a23,rd);            
g2      = g1;

%====================================================================
% Postprocessing
%====================================================================    
r  = max(min(r2,1),0);
g  = max(min(g2,1),0);
b  = max(min(b2,1),0);

mp = [r g b];

%====================================================================
% Plot Default Figure
%====================================================================
if pf,
    figure;
    FigureSet(1);
    h = plot([mp,sum(mp,2)/3,mean(rgb2gray(mp),2)]);
    set(h(1),'Color','r');
    set(h(2),'Color','g');
    set(h(3),'Color','b');
    set(h(4),'Color','k');
    set(h,'LineWidth',1.5);
    xlim([1 nc]);
    ylim([0 1]);
    box off;
    xlabel('Map Index');
    legend('Red','Green','Blue','Intensity','Gray');
    if nc<=256,                                            % MATLAB doesn't display colormaps with more than 256 colors correctly
        colormap(mp);
        colorbar;
        end;
    ylim([0 1.03]);
    AxisSet;
    end;
    
%====================================================================
% Process Return Arguments
%====================================================================
if nargout==0,
    clear('mp');
    end;    
    
% 
% figure;
% FigureSet(2);
% ha = axes;
% set(ha,'Visible','off');
% colormap(M);
% h = colorbar('Location','Manual');
% set(h,'Position',[0.1 0.1 0.8 0.8]);


