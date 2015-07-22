function zocb(cbo,foo,sb,ah)
%zoom out callback
%sb     scrollbar object
%ah     axis handle object containing updated sb

sb = get(ah, 'userdata');%get updated object
ra = sb.range;
lim = sb.lim;
rng = (lim(2)-lim(1));
ctr = (ra(2)-ra(1))/2+ra(1);
lim(1) = ctr-rng;
lim(2) = ctr+rng;
ChangeLimit(sb,lim);
