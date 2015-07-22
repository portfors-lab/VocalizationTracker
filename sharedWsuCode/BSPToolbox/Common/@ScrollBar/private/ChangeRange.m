function ChangeRange(sb)
%Change the scrollbar range
%DOES NOT UPDATE sb stored in sb.a

%modify left patch width
if(strcmp(lower(sb.orientation), 'horizontal')
    xdl = get(sb.pl, 'xdata');
    xdr = get(sb.pr, 'xdata');
    xd  = get(sb.p , 'xdata');
else
    xdl = get(sb.pl, 'ydata');
    xdr = get(sb.pr, 'ydata');
    xd  = get(sb.p , 'ydata');
end

%modify left patch position
pra = xdl(2) - xdl(1);
xdl([1 4]) = sb.range(1);
xdl(2:3) = xdl(1) + pra;

%modify right patch position
pra = xdr(2) - xdr(1);
xdr(2:3) = sb.range(2);
xdr([1 4]) = sb.range(2) - pra;

%modify main patch position
xd([1 4]) = xdl(2:3);
xd(2:3) = xdr([1 4]);

if(strcmp(lower(sb.orientation), 'horizontal'))
    set(sb.pl, 'xdata', xdl);
    set(sb.pr, 'xdata', xdr);
    set(sb.p , 'xdata', xd );
else
    set(sb.pl, 'ydata', xdl);
    set(sb.pr, 'ydata', xdr);
    set(sb.p , 'ydata', xd );
end
