function [] = FigureLatex(fha)
%FigureLatex: Uses the LaTeX interpreter for all text objects in the figure
%
%   [] = FigureLatex(fh)
%
%   fh   Figure handle. Default=Current figure. 
%
%   Causes all text objects in the figure to use the LaTeX interpreter. 
%   Handy for figures containing greek letters or mathematical formulas.
%
%   There is a known bug in MATLAB (R14SP3) that prevents figures using the
%   LaTeX interpreter on empty strings from rendering properly when using 
%   the OpenGL renderer. If you need to use the OpenGL renderer for some 
%   reason (like transparencies), make sure all of the text objects (like 
%   axis titles) in the figure are nonempty.
%
%   Example: Generate a standard plot with labels.
%
%      n = (1:31)';
%      y1 = 3*randn(1,31)+80;
%      y2 = repmat(mean(y1), [1 31]);
%      figure;
%      plot(n,y1,'b');
%      hold on;
%      p = plot(n,y2,'r');
%      xlabel('Day (n)');
%      ylabel('High Temperature ($^\circ$F)');
%      FigureSet;
%      AxisSet;
%      legend('$T$, daily temperature', ...
%             '$\bar{T} = \frac{1}{N}\sum^{N}_{n=1}T_n$');
%      AxisSet(12);
%      FigureLatex;
%
%   See also AxisSet, FormatTicks, and FigureSet.

%====================================================================
% Process function arguments
%====================================================================
fh = gcf; % Default figure
if exist('fha') & ~isempty(fha)
    fh = fha;
end

hv = get(fh, 'Children');% vector of graphics object handles

%====================================================================
% Main Loop
%====================================================================    
while(~isempty(hv))
    %change xlabels, ylabels, zlabels, and titles.
    if strcmp(get(hv(1), 'type'), 'axes')
        h = get(hv(1), 'XLabel');
        set(h, 'Interpreter', 'LaTeX');
        
        h = get(hv(1), 'YLabel');
        set(h, 'Interpreter', 'LaTeX');
        
        h = get(hv(1), 'ZLabel');
        set(h, 'Interpreter', 'LaTeX');
        
        h = get(hv(1), 'Title');
        set(h, 'Interpreter', 'LaTeX');
    else %this may be a graphics object with a string property
        try
            get(hv(1), 'String');
            set(hv(1), 'Interpreter', 'LaTeX');
        catch
            %Guess not. No problem though.
        end
    end
    hv = [get(hv(1), 'children'); hv(2:end)]; %remove current object
                                              %and add children to list
end
