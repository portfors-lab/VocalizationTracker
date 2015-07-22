function zbh = addZoomButtons(sb)

units = get(sb.a, 'units');
set(sb.a, 'units', 'pixels');
pos = get(sb.a, 'position');
set(sb.a, 'units', units);
if(strcmp(lower(sb.orientation), 'horizontal'))
    zobpos = [pos(1)-pos(4)/2 pos(2) pos(4)/2 pos(4)/2];
    zibpos = [pos(1)-pos(4)/2 pos(2)+pos(4)/2 pos(4)/2 pos(4)/2];
else
    zobpos = [pos(1) pos(2)-pos(3)/2 pos(3)/2 pos(3)/2];
    zibpos = [pos(1)+pos(3)/2 pos(2)-pos(3)/2 pos(3)/2 pos(3)/2];
end

zo = uicontrol('units', 'pixels', 'position', zobpos, ...
              'tag', 'zobutton', 'string', '-', 'callback', {@zocb,sb,sb.a});
zi = uicontrol('units', 'pixels', 'position', zibpos, ...
              'tag', 'zibutton', 'string', '+', 'callback', {@zicb,sb,sb.a});
set(zo, 'units', 'normalized');
set(zi, 'units', 'normalized');
zbh = [zo zi];