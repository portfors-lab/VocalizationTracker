function [h] = VArrow(ys,yt,xc,BothWays,AHWarg,AHLarg);
% [h] = VArrow(y-start, y-tip, x-location, BothWays,Width,Length);
% This function draws a horizontal arrow extending 
% from the point specified by (y-start,x) to 
% (y-tip, x). Uses the axes to determine the
% arrow-head length and height, so the limits on the
% x and y axes should be set before calling this 
% function.
%
% The last two arguments are optional.

yl   = get(gca,'YLim'); % Y-axis limits
xl   = get(gca,'XLim'); % X-axis limits
Yrng = yl(2) - yl(1);
Xrng = xl(2) - xl(1);

un   = get(gca,'Units');
set(gca,'Units','inches');
p    = get(gca,'Position');
set(gca,'Units',un);
ar   = p(3)/p(4);    % Aspect ratio

AHW = 0.005*Xrng;    % Arrow head width
AHL = 0.020*Yrng*ar; % Arrow head length

if exist('AHWarg') & ~isempty(AHWarg),
    AHW = AHWarg;
    end;
    
if exist('AHLarg') & ~isempty(AHLarg),
    AHL = AHLarg;
    end;
    

if exist('BothWays') & ~isempty(BothWays),% Arrow points both ways
    h1 = plot(xc*[1 1],[min([ys yt])+AHL max([ys yt])-AHL],'k');
    hs = ishold;
    hold on;
    h2 = patch(xc + [0 AHW -AHW],max([yt ys]) - [0 AHL AHL],[1 1 1],[0 0 0]);
    h3 = patch(xc + [0 AHW -AHW],min([yt ys]) + [0 AHL AHL],[1 1 1],[0 0 0]);
    h2 = [h2;h3];
    set(h2,'LineWidth',0.001);
elseif ys>yt,        % Arrow points down
    h1 = plot(xc*[1 1],[min([ys yt])+AHL max([ys yt])],'k');
    hs = ishold;
    hold on;
    h2 = patch(xc + [0 AHW -AHW],min([yt ys]) + [0 AHL AHL],[1 1 1],[0 0 0]);
    set(h2,'LineWidth',0.001);
else                 % Arrow points up 
    h1 = plot(xc*[1 1],[min([ys yt]) max([ys yt])-AHL],'k');
    hs = ishold;
    hold on;
    h2 = patch(xc + [0 AHW -AHW],max([yt ys]) - [0 AHL AHL],[1 1 1],[0 0 0]);
    set(h2,'LineWidth',0.001);
    end;
if ~hs,
    hold off;
    end;
    
h = [h1;h2];