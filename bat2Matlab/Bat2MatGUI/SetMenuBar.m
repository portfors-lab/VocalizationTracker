function SetMenuBar(figureHandle, stimPath)
%sets up menu bar for the visualization machine program.  Stimpath is the
%path of the folder which contains the stimulus call files.

fileMenu = uimenu(figureHandle, 'label', 'File');
compareMenu = uimenu(figureHandle, 'label', 'Comparison Analysis');

uimenu(fileMenu, 'label', 'Help', 'callback', @helpFun);
uimenu(fileMenu, 'label', 'Exit', 'callback', @exitFun);
uimenu(compareMenu, 'label', 'tuning curve difference', 'callback', 'TCDiff');
uimenu(compareMenu, 'label', 'Correlation Matrix', 'callback', 'CorrMetricGUI');
uimenu(compareMenu, 'label', 'Spike Predictions', 'callback', @predictFun);

    function helpFun(hObj, eventdata)
    %open up a help document or display to command window
    helpdlg('TODO: write help file')
    end

    function predictFun(hObj, evtdata)
        VisualizePredictions(stimPath)
    end

    function exitFun(hObj, eventdata)
        close(figureHandle);
    end

end