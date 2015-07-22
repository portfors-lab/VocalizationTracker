function out = get(sb,prop)

if ~nargin
    error('Not enough input arguments');
end
sb = get(sb.a, 'userdata'); %make sure we're using the real version
if nargin == 1
    disp(sb);
    return;
end

if(~ischar(prop))
    error('Property name argument must be a string');
end
switch(lower(prop))
    case 'callbackfcn'
        if(ischar(sb.cb))
            out = sb.cb;
        else
            out = func2str(sb.cb);
        end
    case 'placement'
        out = sb.placement;
    case 'axis'
        out = sb.a;
    case 'range'
        out = sb.range;
    case 'limit'
        out = sb.lim;
    case 'orientation'
        out = sb.orientation;
    otherwise
        error('There is no ''%s'' property in the ''ScrollBar'' class.', prop);
end

