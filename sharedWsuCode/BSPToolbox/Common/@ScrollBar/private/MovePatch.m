function MovePatch(cbo, foo, sb, patch)
%move the patch

sb = get(sb.a, 'userdata');%get updated scrollbar

cp = get(sb.a, 'currentpoint'); %get current position of mouse
mvmt = cp-sb.pp;                %figure out how far it moved          
sb.pp = cp;
if(strcmp(lower(sb.orientation), 'horizontal'))
    xpos = get(sb.p, 'xdata');      %get current patch position
    xposr = get(sb.pr, 'xdata');
    xposl = get(sb.pl, 'xdata');
    mvmt = mvmt(1,1);
    lim = get(sb.a, 'xlim');
else
    xpos = get(sb.p, 'ydata');
    xposr = get(sb.pr, 'ydata');
    xposl = get(sb.pl, 'ydata');
    mvmt = mvmt(1,2);
    lim = get(sb.a, 'ylim');
end
sb.patchWidth = xposl(3)-xposl(1);



%figure out which patch moved
if(patch == sb.p)
    %don't allow to exceed axis limits
    if(xposr(2) + mvmt > lim(2))
        mvmt = lim(2)-xposr(2);
    elseif(xposl(1) + mvmt < lim(1))
        mvmt = lim(1)-xposl(1);
    end
    xpos = xpos + mvmt(1);              %move patches
    xposl = xposl + mvmt(1);
    xposr = xposr + mvmt(1);
elseif(patch == sb.pr)
    %don't allow to exceed axis limits
    if(xposr(2) + mvmt > lim(2))
        mvmt = lim(2)-xposr(2);
    end
    xpos(2:3) = ones(2,1)*max(xpos(2) + mvmt(1), xpos(1));
    xposr([1 4]) = ones(2,1)*max(xposr(1) + mvmt(1), xpos(2));
    xposr(2:3) = xposr([1 4])+xposl(2:3)-xposl([1 4]);
elseif(patch == sb.pl)
    %don't allow to exceed axis limits
    if(xposl(1) + mvmt < lim(1))
        mvmt = lim(1)-xposl(1);
    end
    xpos([1 4]) = ones(2,1)*min(xpos(1) + mvmt(1), xpos(2));
    xposl(2:3) = ones(2,1)*min(xposl(2) + mvmt(1), xpos(1));
    xposl([1 4]) = xposl(2:3)+xposr([1 4])-xposr(2:3);
end

if(strcmp(lower(sb.orientation), 'horizontal'))
    set(sb.p , 'xdata', xpos );
    set(sb.pr, 'xdata', xposr);
    set(sb.pl, 'xdata', xposl);
else
    set(sb.p , 'ydata', xpos );
    set(sb.pr, 'ydata', xposr);
    set(sb.pl, 'ydata', xposl);
end
newpts = [xposl(1) xposr(2)];
sb.range = newpts;
set(sb.a,'userdata',sb);
try
    if(~isempty(sb.cb))
        feval(sb.cb,newpts);%call user's callback fcn with newpts as the arg
    end
catch
    disp('Error in user-supplied callback function');
    disp(lasterr);
end
drawnow;                            %redraw