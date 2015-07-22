function PrintFigure(p_name,p_format,p_width,p_height,p_resolution,figures)
%   PrintFigure: Print a figure to a file
%
%   printFigure(name,format,width,height,resolution,figures)
%   
%   name        Name of file (string)
%               Default: 'figure'
%   format      Image format ('eps','epsc','eps2','epsc2','jpeg','tiff','png',...)
%               Default: 'jpeg'
%   width       Width of figure (inches)
%               Default: 8
%   height      Height of figure (inches)
%               Default: 6
%   resolution  Resolution of figure (dpi)
%               Default: 300
%   figures     A vector of pointers to figures to apply the function to
%               Default: current figure

name = 'figure';
if exist('p_name') && ~isempty(p_name)
   name = p_name;
end

format = '-djpeg';
if exist('p_format') && ~isempty(p_format)
   format = ['-d' p_format];
end

width = 8;
if exist('p_width') && ~isempty(p_width)
   width = p_width;
end

height = 6;
if exist('p_height') && ~isempty(p_height)
   height = p_height;
end

resolution = '-r300';
if exist('p_resolution') && ~isempty(p_resolution)
   resolution = ['-r' num2str(p_resolution)];
end

if exist('figures','var')
    if isempty(figures)
        figures = gcf;
    end
else
    figures = gcf;
end
if figures == 0
    figures = findall(0,'type','figure');
end

for figure_num = 1:length(figures)
     hFigure = figures(figure_num);
%     if strmatch(p_format, 'epsc2')
%        set(figures(figure_num), 'Renderer', 'painter');
%     end
%     set(hFigure,'Units','Inches');
%     set(hFigure,'PaperSize',[width height]);
%     set(hFigure,'PaperPosition',[0 0 width height]);
%     resolutionString = ['-r' int2str(resolution)];
     figureName = name;
     if length(figures) > 1
         figureName = [name '_' int2str(figure_num)];
     end
%     if strmatch(p_format,'pdf')
%         saveas(hFigure,figureName,'pdf');
%     else
%         print(hFigure,figureName,format,resolutionString);
    %     end
    if strmatch(p_format, 'fig')
        saveas(hFigure, figureName, 'fig');
    else
        print(hFigure, figureName, format, resolution);
    end
end


