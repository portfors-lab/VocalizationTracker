function SetColorBar(label,range)
%
%function GetColorScale(color_axis_range,range)
%
%   INPUT ARGUMENTS
%   label   The label used on the colorbar
%           Default : ''
%   range   The range of data for the colorbar to span.
%           Three modese are possible:
%               0 : Don't rescale (Default)
%               1 : Rescale to cover the full range of
%                   data but to center the range on 0
%               2 : Reslcale between [0 max]
%               3 : Rescale between [-max 0]
%               4 : Rescale to [0 1] for normalized data
%               5 : Rescale to [-1 1] for normalized data
%               [a b] : Rescale between a and b

if exist('range','var')
    if isempty(range)
        range = 0;
    end
else
    range = 0;
end

if exist('label','var')
    if isempty(label)
        label = [];
    end
else
    label = [];
end

%Set colorbar
if length(range) == 2
    caxis(range);
elseif range == 1
    max_val = max(abs(caxis));
    caxis([-max_val max_val]);
elseif range == 2
    max_val = max(max((caxis)),0);
    caxis([0 max_val]);
elseif range == 3
    min_val = min(min((caxis)),0);
    caxis([min_val 0]);
elseif range == 4
    caxis([0 1]);
elseif range == 5
    caxis([-1 1]);
end
c = colorbar;

%Set colorbar axis label
if ~isempty(label)
    y = ylabel(c,label);
    set(y,'Rotation',270);
    set(y,'VerticalAlignment','bottom');
end

