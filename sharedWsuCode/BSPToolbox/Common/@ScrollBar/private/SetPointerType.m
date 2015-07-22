function out = SetPointerType(cbo, foo,sb,fcn);
ptr = 0;
if(exist('fcn') && ~isempty(fcn))
    if(length(fcn)>1)
        ptr = feval(fcn{1},cbo,foo,fcn{2:end});
    else
        ptr = feval(fcn{1},cbo,foo);
    end
end
sb = get(sb.a, 'userdata');
out = 1; %1 means the pointer was over a patch 
cp = get(gcf, 'CurrentPoint');
units = get(sb.a, 'units');
set(sb.a, 'units', 'pixels');
pos = get(sb.a, 'position');
set(sb.a, 'units', units);
if(cp(1) > pos(1) && cp(1) < pos(1)+pos(3) && cp(2) > pos(2) && cp(2) < pos(2)+pos(4))
    cp = get(sb.a, 'CurrentPoint');
    xd = get(sb.pl, 'XData');
    yd = get(sb.pl, 'YData');
    ll = [min(xd) min(yd)];
    ur = [max(xd) max(yd)];
    if(cp(1,1) > ll(1) && cp(1,2) > ll(2) && cp(1,1) < ur(1) && cp(1,2) < ur(2))
        if(strcmp(lower(sb.orientation), 'horizontal'))
            set(gcf, 'Pointer', 'left');
            return;
        else
            set(gcf, 'Pointer', 'bottom');
            return;
        end
    end
    xd = get(sb.pr, 'XData');
    yd = get(sb.pr, 'YData');
    ll = [min(xd) min(yd)];
    ur = [max(xd) max(yd)];
    if(cp(1,1) > ll(1) && cp(1,2) > ll(2) && cp(1,1) < ur(1) && cp(1,2) < ur(2))
        if(strcmp(lower(sb.orientation), 'horizontal'))
            set(gcf, 'Pointer', 'right');
            return;
        else
            set(gcf, 'Pointer', 'top');
            return;
        end
    end
    xd = get(sb.p, 'XData');
    yd = get(sb.p, 'YData');
    ll = [min(xd) min(yd)];
    ur = [max(xd) max(yd)];
    if(cp(1,1) > ll(1) && cp(1,2) > ll(2) && cp(1,1) < ur(1) && cp(1,2) < ur(2))
        set(gcf, 'Pointer', 'fleur');
        return;
    end
end
if(~ptr)
    set(gcf, 'pointer', 'arrow');
    out = 0;
end

if ~nargout
    clear out;
end