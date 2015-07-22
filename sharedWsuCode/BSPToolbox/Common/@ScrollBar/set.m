function sb = set(sb, varargin)

if ~nargin
    error('Not enough input arguments');
end
sb = get(sb.a, 'userdata'); %make sure we're using the real version
if nargin == 1
    disp(sb);
%    disp(['CallbackFcn       ' sb.cb]);
%    disp(['Placement         ' sb.placement]);
    
    if(~nargout)
        clear sb;
    end
    return;
end

if(iscell(varargin{1}))
    args = varargin{1};
else
    args = varargin;
end
nargs = length(args);
for arg=1:2:nargs
    switch(lower(args{arg}))
        case 'callbackfcn'
            sb.cb = args{arg+1};
        case 'placement'
            sb.placement = args{arg+1};
        case 'axis'
            error('The axis cannot be modified.');
        case 'range'
            sb.range = args{arg+1};
            ChangeRange(sb);
        case 'limit'
            ChangeLimit(sb,args{arg+1});
        case 'zoombuttons'
            sb.zb = args{arg+1};
            if(strcmp(lower(sb.zb), 'on'))
                sb.zbh = addZoomButtons(sb);
            else
                delete(sb.zbh)
                sb.zbh = [];
            end
        otherwise
            error('There is no ''%s'' property in the ''ScrollBar'' class.', args{arg});
    end
end
set(sb.a,'userdata',sb);
if ~nargout
    clear sb;
end