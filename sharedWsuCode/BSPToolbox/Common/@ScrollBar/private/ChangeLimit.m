function ChangeLimit(sb,lim)

%modify the range of the axis
%sb = get(sb.a, 'userdata');
plim = sb.lim;
sb.lim = lim;
ra = lim(2) - lim(1);
prevra = plim(2) - plim(1);

%get patch data
if(strcmp(lower(sb.orientation), 'horizontal'))
    set(sb.a, 'xlim', lim);
    xdl = get(sb.pl, 'xdata');
    xdr = get(sb.pr, 'xdata');
    xd  = get(sb.p , 'xdata');
else
    set(sb.a, 'ylim', lim);
    xdl = get(sb.pl, 'ydata');
    xdr = get(sb.pr, 'ydata');
    xd  = get(sb.p , 'ydata');
end

%modify left patch width
pra = xdl(2) - xdl(1);
pra = pra*ra/prevra;
%xdl([1 4]) = xdl(2:3) - pra;
xdl = xdl - min(min(xdl-lim(1), 0));%keep in axis limits
xdl(2:3) = xdl([1 4]) + pra;


%modify right patch width
pra = xdr(2) - xdr(1);
pra = pra*ra/prevra;
%xdr([2:3]) = xdr([1 4]) + pra;
xdr = xdr - max(max(xdr-lim(2), 0));%keep in axis limits
xdr([1 4]) = xdr(2:3) - pra;


%modify main patch width
xd([1 4]) = max(xdl);
xd(2:3) = min(xdr);

%update patch position data
if(strcmp(lower(sb.orientation), 'horizontal'))
    set(sb.pl, 'xdata', xdl);
    set(sb.pr, 'xdata', xdr);
    set(sb.p , 'xdata', xd );
else
    set(sb.pl, 'ydata', xdl);
    set(sb.pr, 'ydata', xdr);
    set(sb.p , 'ydata', xd );
end
set(sb.a, 'userdata', sb);
try
    if(~isempty(sb.cb))%call user's callback fcn with new points as the arg
        feval(sb.cb,[xdl(1) xdr(2)]);
    end
catch
    disp('Error in user-supplied callback function');
    disp(lasterr);
end
