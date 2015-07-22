function ScrollBarClick(cbo, foo, sb, patch)
%Callback for when a patch is clicked
%Enables/disables callbacks used for moving the patch
sb = get(sb.a, 'userdata');%get updated scrollbar

if(~isempty(get(patch, 'ButtonDownFcn')))
    %start moving the patch
    sb.pp = get(sb.a, 'CurrentPoint');
    sb.bdcb = get(patch, 'ButtonDownFcn');
    set(patch, 'ButtonDownFcn', []);
    sb.wbmf = get(gcf, 'WindowButtonMotionFcn');
    set(gcf, 'WindowButtonMotionFcn', {@MovePatch,sb,patch});
    set(gcf, 'WindowButtonUpFcn', {@ScrollBarClick,sb,patch});
%    set(gcf, 'WindowButtonDownFcn', {@ScrollBarClick,sb,patch});
else
    %stop moving the patch
    set(patch, 'ButtonDownFcn', sb.bdcb);
    set(gcf, 'WindowButtonMotionFcn', sb.wbmf);
    set(gcf, 'WindowButtonUpFcn', []);
%    set(gcf, 'WindowButtonDownFcn', []);
end

set(sb.a, 'userdata', sb);%update scrollbar