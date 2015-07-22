function h = PlotWithConfidenceIntervals(index,values,upper_bound,lower_bound,lineSpec)

%Make all vectors column vectors
index = index(:)';
values = values(:)';
upper_bound = upper_bound(:)';
lower_bound = lower_bound(:)';

if strfind(lineSpec,'r');
    fill_color = [1 0.7 0.7];
elseif strfind(lineSpec,'g');
    fill_color = [0.7 1 0.7];
elseif strfind(lineSpec,'b');
    fill_color = [0.7 0.7 1];
elseif strfind(lineSpec,'k');
    fill_color = [0.7 0.7 0.7];
elseif strfind(lineSpec,'m');
    fill_color = [1 0.8 1];
else
    fill_color = [0.6 0.6 1];
end

if strfind(lineSpec,'--')
    lineStyle = ':';
elseif strfind(lineSpec,'-.')
    lineStyle = '-.';
elseif strfind(lineSpec,':')
    lineStyle = ':';
else
    lineStyle = '-';
end

hold_figure = false;
if ~ishold
    hold on
    hold_figure = true;
end
f = fill([index fliplr(index)],[upper_bound fliplr(lower_bound)],fill_color);
set(f,'LineStyle',lineStyle);
h = plot(index,values,lineSpec);
set(h,'LineWidth',2);

if hold_figure
    hold off
end

