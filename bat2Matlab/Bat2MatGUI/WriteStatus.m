function WriteStatus(text, color, fig)
    try
        if ~exist('fig','var')
            fig = findobj('name', 'Visualization Machine');
            if isempty(fig)
                if any(findobj('type', 'figure'))
                    fig = gcf;
                end
            end
        end

        if ~exist('color', 'var')
            color = 'black';
        end

        textPane = findjobj(fig, 'nomenu', 'class', 'javax.swing.JTextPane');
        if ~isempty(textPane) %status box present on figure
            statusText = textPane.getStyledDocument();
            statusText.putLine(text, color);
            textPane.setCaretPosition(statusText.getLength());
        else %status box not present, write to command line
            if strcmp(color, 'orange')
                color = 'systemcommands';
            end
            cprintf(color, [text '\n']);
        end
    catch
        %if all else fails
        disp(text)
    end
end

