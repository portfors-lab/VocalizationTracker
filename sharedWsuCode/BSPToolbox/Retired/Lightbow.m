function [mp] = Lightbow(nca,npa,pfa);
%Lightbow: Generates a monotonic colormap with maximum color depth
%
%   [m] = Lightbow(n,np,pf);
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
%      colormap(Lightbow(256));
%
%   C. Ware, "Information Visualization: Percetpion for Design," 2nd edition,
%   Morgan Kaufmann, 2004.
%
%   Version 1.00 JM
%
%   See also colormap, jet, and caxis.

%====================================================================
% Error Checking
%====================================================================    
error('The function Lightbow has been replaced by ColorSpiral.');

if nargin<1,
    help Lightbow;
    return;
    end;
    
%====================================================================
% Process Function Arguments
%====================================================================    
nc = 64;                                                   % Number of colors in the colormap
if exist('nca') & ~isempty(nca),
    nc = nca;
    end;

np = 1.5;                                                    % Number of sinusoidal periods
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
% Main Loop
%====================================================================    
t  = linspace(sqrt(3),0,nc).';
wn = [0;sqrt(3/8)*triang(nc-2);0];
%ce = wn.*exp(j*(t*np*2*pi/sqrt(3)+pi/3));
ce = wn.*exp(j*((t-sqrt(3)/2)*np*2*pi/sqrt(3)));
x  = real(ce);
y  = imag(ce);

r  = t;
g  = x;
b  = y;

r0 = t;
g0 = x*0;
b0 = y*0;

% figure;
% FigureSet(1);
% plot3(r,g,b,'k');
% xlabel('r');
% ylabel('g');
% zlabel('b');
% hold on;
%     plot3(r0,g0,b0,'m');
%     plot3(r  ,g    ,b*0-1,'b');
%     plot3(r  ,g*0-1,b,'g');
%     plot3(r*0,g    ,b,'r');
%     hold off;
% %grid on;
% axis([0 2 -1 1 -1 1]);
% rotate3d on;
% AxisLines;

%==========================================================
% Rotate 45 degrees about the g axis
%==========================================================
drb  = (r.^2+b.^2).^(1/2);
arb  = angle(r+j*b);
r    = drb.*cos(asin(1/sqrt(3))+arb);
g    = g;
b    = drb.*sin(asin(1/sqrt(3))+arb);

drb0 = (r0.^2+b0.^2).^(1/2);
arb0 = angle(r0+j*b0);
r0   = drb0.*sin(asin(1/sqrt(3))+arb0);
g0   = g0;
b0   = drb0.*cos(asin(1/sqrt(3))+arb0);

% figure;
% FigureSet(2);
% plot3(r,g,b,'k');
% xlabel('r');
% ylabel('g');
% zlabel('b');
% hold on;
%     plot3(r0,g0,b0,'m');
%     plot3(r  ,g    ,b*0-1,'b');
%     plot3(r  ,g*0-1,b,'g');
%     plot3(r*0,g    ,b,'r');
%     hold off;
% axis([0 2 -1 1 -1 1]);
% %grid on;
% rotate3d on;
% AxisLines;

%return;

%==========================================================
% Rotate 45 degrees about the b axis
%==========================================================
drg  = (r.^2+g.^2).^(1/2);
arg  = angle(r+j*g);
r    = drg.*sin(pi/4+arg);
g    = drg.*cos(pi/4+arg);
b    = b;

drg0 = (r0.^2+g0.^2).^(1/2);
arg0 = angle(r0+j*g0);
r0   = drg0.*sin(pi/4+arg0);
g0   = drg0.*cos(pi/4+arg0);
b0   = b0;

% figure;
% FigureSet(4);
% plot3(r,g,b,'k');
% xlabel('r');
% ylabel('g');
% zlabel('b');
% hold on;
%     plot3(r0,g0,b0,'m');
%     plot3(r  ,g    ,b*0-1,'b');
%     plot3(r  ,g*0-1,b,'g');
%     plot3(r*0,g    ,b,'r');
%     hold off;
% %grid on;
% axis([0 2 -1 1 -1 1]);
% rotate3d on;
% AxisLines;

r  = max(min(r,1),0);
g  = max(min(g,1),0);
b  = max(min(b,1),0);

%====================================================================
% Swap red and green so that the colors more closely match jet
%====================================================================
u = r;
r = g;
g = u; 

mp = [r g b];

%====================================================================
% Plot Default Figure
%====================================================================
if pf,
    figure;
    FigureSet(1);
    h = plot([mp,sum(mp,2)/3]);
    set(h(1),'Color','r');
    set(h(2),'Color','g');
    set(h(3),'Color','b');
    set(h(4),'Color','k');
    set(h,'LineWidth',1.5);
    xlim([1 nc]);
    ylim([0 1]);
    box off;
    xlabel('Map Index');
    legend('Red','Green','Blue','Intensity');
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


