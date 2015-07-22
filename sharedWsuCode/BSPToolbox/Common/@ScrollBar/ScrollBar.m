function sb = ScrollBar(varargin)
%ScrollBar: Creates a versitile scrollbar on the current figure.
%  
%   sb = ScrollBar('Property', Value, ...);        
%  
%   Property         |   Value
%   -----------------+------------------------------------------------------
%   'Orientation'    |   [{Horizontal} | Vertical]
%   'Position'       |   'auto' -or- [1x4] RECT -or- axis handle
%   'Callback'       |   function handle -or- axis handle
%   'Limit'          |   [1x2] scrollbar axis limits
%   'Range'          |   [1x2] scrollbar range
%   'ZoomButtons'    |   [{On} | Off]
%
%   ScrollBar creates an axis, three patches and optionally zoom buttons.
%   Moving the middle patch moves the entire bar, while moving either of
%   the end patches will stretch the bar in that direction. Whenever one
%   of the patches changes position, a user supplied callback is called
%   with a two element vector of the new scrollbar range as the only
%   argument. If an axis handle is provided instead of a callback
%   function, ScrollBar will change the axis limits of that axis to match
%   the new scrollbar range in whatever direction the scrollbar is
%   oriented.
%
%   If anything is plotted in the axis containing the scrollbar, it may
%   prevent the patch from moving when the plot is clicked. To avoid this
%   the HitTest property of the plot should be turned off.
%
%   The object returned is not a graphics handle. If it is displayed, the
%   handle of the axis containing the scrollbar is displayed instead. 
%  
%   Example: Create a plot of an MER signal and a vertical and horizontal
%   scrollbar on the same figure.
%  
%      load MER;
%      figure;
%      a = axes;
%      t = (1:length(x))/fs;
%      plot(t,x,'b',t(si),x(si),'.r');
%      hsb = ScrollBar;
%      vsb = ScrollBar('Orientation','Vertical', ...
%            'Callback',a, 'Limit', [min(x) max(x)]*1.1);
%      set(vsb, 'ZoomButtons', 'Off');
%      k = get(hsb, 'Axis');
%      axes(k);
%      hold on;
%      p = plot(t, (x-min(x))./(max(x)-min(x)));
%      set(p, 'HitTest', 'off');
%
%   See also Colorbar.

%allocate fields
sb.a            = [];   %scrollbar axis
sb.cb           = [];   %callback function handle
sb.cba          = [];   %axis to modify for auto callback
sb.placement    = [];   %auto or manual placement
sb.orientation  = [];   %scrollbar orientation
sb.patchWidth   = [];   %"arrow" patch width
sb.lim          = [];   %scrollbar axis limits
sb.range        = [];   %scrollbar range
sb.zb           = [];   %zoom button flag
sb.pp           = [];   %previous point
sb.bdcb         = [];   %button down callback 
sb.wbmf         = [];   %WindowButtonMotionFcn
sb.zbh          = [];   %zoom button handles
sb.ph           = [];   %patch height

f = gcf;                %scrollbar's figure

%process input arguments
if(mod(nargin,2) ~= 0)
    error('Unexpected number of arguments');
end
for arg = 1:2:nargin
    if(~ischar(varargin{arg}))
        error('Error parsing argument %f. Expected string', arg);
    end
    switch(lower(varargin{arg}))
        case 'orientation'
            if(strcmpi(varargin{arg+1}, 'Horizontal'))
                sb.orientation = 'Horizontal';
            elseif(strcmpi(varargin{arg+1}, 'Vertical'))
                sb.orientation = 'Vertical';
            else
                error('Unexpected value for %s : %s', ...
                    varargin{arg}, varargin{arg+1});
            end
        case 'position'
            if(prod(size(varargin{arg+1})) == 4)
                %The argument is the position rectangle for a new axis
                sb.a = axes('position', varargin{arg+1});
                sb.placement = 'Manual';
            elseif(ishandle(varargin{arg+1}))
                %The argument is an axis handle
                sb.a = varargin{arg+1};
                sb.placement = 'Manual';
            else
                error('Unexpected value for %s : %s', ...
                    varargin{arg}, varargin{arg+1});
            end
        case 'callback'
            if(isa(varargin{arg+1}, 'function_handle'))%function handle
                sb.cb = varargin{arg+1};
            elseif(ishandle(varargin{arg+1}))%axis handle
                sb.cba = varargin{arg+1};
                sb.cb = 'auto';
            else
                error('Unexpected value for %x : %s', ...
                    varargin{arg}, varargin{arg+1});
            end
        case 'limit'
            if(~isnumeric(varargin{arg+1}))
                error('Axis limits must be numeric');
            elseif(prod(size(varargin{arg+1})) == 2)
                sb.lim(1,1:2) = varargin{arg+1};
            else
                error('Limits expected to be numeric.');
            end
        case 'range'
            if(~isnumeric(varargin{arg+1}))
                error('Scrollbar range must be numberic');
            elseif(prod(size(varargin{arg+1})) == 2)
                sb.range(1,1:2) = varargin{arg+1};
            else
                error('Range expected to be numeric.');
            end
        case 'zoombuttons'
            if(strcmp(lower(varargin{arg+1}), 'on'))
                sb.zb = 'On';
            elseif(strcmp(lower(varargin{arg+1}), 'off'))
                sb.zb = 'Off';
            else
                error('Unexpected value for %x : %s', ...
                    varargin{arg}, varargin{arg+1});
            end
        otherwise
            error('Unexpected argument: %s', varargin{arg});
    end
end

%Guess good values for anything not provided
if isempty(sb.a)
    sb.a = [];
end
if isempty(sb.cb)
    sb.cb = 'auto';
    %find major axis
    k = findobj(f, 'Type', 'axes');
    if(isempty(k))
        warning('No axis on current figure');
        k = axes; %create an axis
    end
    ms = 0; %max size
    ma = 0; %major axis
    for(c1 = 1:length(k))
        pu = get(k(c1), 'units');
        set(k(c1), 'units', 'normalized');
        pos = get(k(c1), 'position');
        set(k(c1), 'units', pu);
        if pos(3)*pos(4) > ms
            ms = pos(3)*pos(4);
            ma = k(c1);
        end
    end
    sb.cba = ma;
end
if isempty(sb.cba)
    sb.cba = [];
end
if isempty(sb.placement)
    sb.placement = 'Auto';
end
if isempty(sb.orientation)
    if(strcmp(lower(sb.placement),'auto'))
        %No helpful info provided by user
        k = findobj('Tag', 'ScrollBar');
        if length(k) == 1%assume first scrollbar is horizontal
            sb.orientation = 'Vertical';
        else
            sb.orientation = 'Horizontal';
        end        
    else
        %use axis dimensions to determine orientation
        pos = get(sb.a, 'position');
        if(pos(3)>pos(4))                       %long axis
            sb.orientation = 'Horizontal';
        else                                    %tall axis
            sb.orientation = 'Vertical';
        end
    end
end
if isempty(sb.patchWidth)
    sb.patchWidth = .01;
end
if isempty(sb.lim)
    if(~isempty(sb.cba))%use major axis limit
        if(strcmp(lower(sb.orientation), 'horizontal'))
            sb.lim = get(sb.cba, 'xlim');
        else
            sb.lim = get(sb.cba, 'ylim');
        end
    else
        sb.lim = [0 1];
    end
end
if isempty(sb.range)
    sb.range = sb.lim;
end
if isempty(sb.zb)
    sb.zb = 'On';
end
sb.pp = [];%previous point
sb.bdcb = [];%button down callback 
sb.wbmf = [];%WindowButtonMotionFcn

if strcmp(lower(sb.placement), 'auto')
    %Try and guess a good spot for the scrollbar

    %find major axis
    k = findobj(f, 'Type', 'axes');
    if(length(k) == 0)
        warning('No axis on current figure');
        k = axes; %create an axis
    end
    ms = 0; %max size
    ma = 0; %major axis
    for(c1 = 1:length(k))
        pu = get(k(c1), 'units');
        set(k(c1), 'units', 'normalized');
        pos = get(k(c1), 'position');
        set(k(c1), 'units', pu);
        if pos(3)*pos(4) > ms
            ms = pos(3)*pos(4);
            ma = k(c1);
        end
    end
    if(strcmp(lower(sb.orientation), 'horizontal'))
        %push major axis up and place below
        pu = get(ma, 'units');
        set(ma, 'units', 'normalized');
        pos = get(ma, 'position');
        %should take up about 10% of the previous axis
        set(ma, 'position', pos+[0 0.1 0 -0.1]);
        sb.a = axes('position', [pos(1) pos(2) pos(3) 0.05]);
        set(ma, 'units', pu);
        
        %see if there are any other axes lined up with the major axis
        k = findobj('type', 'axes');
        k = k(find(k ~= sb.a));
        for(c1 = 1:length(k))
            apos = get(k(c1), 'position');
            if(apos(2) == pos(2))
                set(k(c1), 'position', apos+[0 0.1 0 -0.1]);
                %if it is a scrollbar, reset zoombuttons
                if(strcmp(get(k(c1), 'tag'), 'ScrollBar'))
                    asb = get(k(c1), 'userdata');
                    if(strcmp(lower(asb.zb), 'on'))
                        delete(asb.zbh);
                        asb.zbh = addZoomButtons;
                    end
                end
            end
        end
    else
        %push major axis right and place to the left
        pu = get(ma, 'units');
        set(ma, 'units', 'normalized');
        pos = get(ma, 'position');
        %should take up about 10% of the previous axis
        set(ma, 'position', pos+[0.1 0 -0.1 0]);
        sb.a = axes('position', [pos(1) pos(2) 0.05 pos(4)]);
        set(ma, 'units', pu);

        %see if there are any other axes lined up with the major axis
        k = findobj('type', 'axes');
        k = k(find(k ~= sb.a));
        for(c1 = 1:length(k))
            apos = get(k(c1), 'position');
            if(apos(1) == pos(1))
                set(k(c1), 'position', apos+[0.1 0 -0.1 0]);
                %if it is a scrollbar, reset zoombuttons
                if(strcmp(get(k(c1), 'tag'), 'ScrollBar'))
                    asb = get(k(c1), 'userdata');
                    if(strcmp(lower(asb.zb), 'on'))
                        delete(asb.zbh);
                        asb.zbh = addZoomButtons(asb);
                    end
                end
            end
        end
    end
end

if(ischar(sb.cb) && strcmp(lower(sb.cb), 'auto'))
    %use default callback function
    a = sb.cba;
    if(strcmp(lower(sb.orientation), 'horizontal'))
        sb.cb = @(lim) set(a,'xlim',lim);
    else
        sb.cb = @(lim) set(a,'ylim',lim);
    end
end

%%------------------------------------------------------------------------%
%% Set up the scrollbar                                                   %
%%------------------------------------------------------------------------%

%opengl is the only renderer that will do transparencies
set(f, 'Renderer', 'opengl');

%set axis properties
if(strcmp(lower(sb.orientation), 'vertical'))
    set(sb.a, 'ylim', sb.lim, 'xlimmode', 'manual');
    set(sb.a, 'xtick', [], 'xticklabel', []);
else
    set(sb.a, 'ylimmode', 'manual', 'xlim', sb.lim);
    set(sb.a, 'ytick', [], 'yticklabel', []);
end
set(sb.a, 'units', 'normalized', 'HitTest', 'on', 'box', 'on');
set(sb.a, 'Tag', 'ScrollBar');

%%set up patches
%Main bar
ra = sb.range;
pw = sb.patchWidth;
if strcmp(lower(sb.orientation), 'horizontal')
    sb.ph = get(sb.a, 'ylim');
    sb.p = patch([ra(1)+pw*(ra(2)-ra(1)) ra(2)-pw*(ra(2)-ra(1)) ...
        ra(2)-pw*(ra(2)-ra(1)) ra(1)+pw*(ra(2)-ra(1))],  ...
        [sb.ph(1) sb.ph(1) sb.ph(2) sb.ph(2)],'g');
else 
    sb.ph = get(sb.a, 'xlim');
    sb.p = patch([sb.ph(1) sb.ph(1) sb.ph(2) sb.ph(2)], ...
        [ra(1)+pw*(ra(2)-ra(1)) ra(2)-pw*(ra(2)-ra(1)) ...
        ra(2)-pw*(ra(2)-ra(1)) ra(1)+pw*(ra(2)-ra(1))],'g');
end
set(sb.p, 'ButtonDownFcn', {@ScrollBarClick, sb,sb.p});
set(sb.p, 'EdgeColor', 'g', 'FaceAlpha', 1, 'LineWidth', 2);
%set(sb.p, 'EdgeColor', 'g', 'LineWidth', 2);%strcat error on 'FaceAlpha' modification
set(sb.p, 'Tag', 'MainPatch');

%Left or down "arrow"
if strcmp(lower(sb.orientation), 'horizontal')
    xd = [ra(1), ra(1)+pw*(ra(2)-ra(1)), ra(1)+pw*(ra(2)-ra(1)), ra(1)];
    yd = [sb.ph(1) sb.ph(1) sb.ph(2) sb.ph(2)];
    sb.pl = patch(xd, yd, 'r');
else
    xd = [sb.ph(1) sb.ph(1) sb.ph(2) sb.ph(2)];
    yd = [ra(1), ra(1)+pw*(ra(2)-ra(1)), ra(1)+pw*(ra(2)-ra(1)), ra(1)];
    sb.pl = patch(xd, yd, 'r');
end
set(sb.pl, 'ButtonDownFcn', {@ScrollBarClick, sb,sb.pl});
set(sb.pl, 'EdgeColor', 'r', 'FaceAlpha', 1, 'LineWidth', 2);
%set(sb.pl, 'EdgeColor', 'r', 'LineWidth', 2);
set(sb.pl, 'Tag', 'LeftPatch');

%Right or up "arrow"
if strcmp(lower(sb.orientation), 'horizontal')
    xd = [ra(2)-pw*(ra(2)-ra(1)), ra(2), ra(2), ra(2)-pw*(ra(2)-ra(1))];
    yd = [sb.ph(1) sb.ph(1) sb.ph(2) sb.ph(2)];
    sb.pr = patch(xd,yd, 'r');
else
    xd = [sb.ph(1) sb.ph(1) sb.ph(2) sb.ph(2)];
    yd = [ra(2)-pw*(ra(2)-ra(1)), ra(2), ra(2), ra(2)-pw*(ra(2)-ra(1))];
    sb.pr = patch(xd, yd, 'r');
end
set(sb.pr, 'ButtonDownFcn', {@ScrollBarClick, sb, sb.pr});
set(sb.pr, 'EdgeColor', 'r', 'FaceAlpha', 1, 'LineWidth', 2);
%set(sb.pr, 'EdgeColor', 'r', 'LineWidth', 2);
set(sb.pr, 'Tag', 'RightPatch');

%move all of the patches to the back so they don't obscure anything
ao = get(sb.a, 'Children');
ao = ao(ao == sb.pr || ao == sb.pr || so == sb.p);
set(ao, 'hittest', 'off'); %scrollbar patches should be the only clickables
ao = [sb.p;sb.pl;sb.pr;ao];

if strcmp(lower(sb.zb), 'on')
    sb.zbh = addZoomButtons(sb);
end
fcn = get(f, 'WindowButtonMotionFcn');
set(f, 'WindowButtonMotionFcn', {@SetPointerType, sb,fcn});

axes(sb.a);%makes clicks on the edge of the axis trigger callback
           %MATLAB bug?
           
set(sb.a, 'userdata', sb);%store object in axis userdata so we can access
                          %it through a handle
%register class
sb = class(sb,'ScrollBar');

