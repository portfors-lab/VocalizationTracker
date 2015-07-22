function TraceChooser(prefs, experimentData)

    fh = figure('position', [500 350 400 225], 'resize', 'off');

    bcolor = [0.8 0.8 0.8];


    uicontrol(fh, 'style', 'text', 'position', [50 185 250 30], 'string', 'Select test/trace to view', 'fontsize', 14, 'backgroundcolor', bcolor);

    uicontrol(fh, 'style', 'text', 'position', [15 150 50 25], 'string', 'test #:', 'fontsize', 12, 'backgroundcolor', bcolor);
    uicontrol(fh, 'style', 'text', 'position', [135 150 60 25], 'string', 'trace #:', 'fontsize', 12, 'backgroundcolor', bcolor);

    testBox = uicontrol(fh, 'style', 'edit', 'position', [75 150 40 25], 'backgroundcolor', bcolor);
    traceBox = uicontrol(fh, 'style', 'edit', 'position', [200 150 40 25], 'backgroundcolor', bcolor);

    bgh = uibuttongroup('parent', fh,  'units', 'pixels','position', [15 60 375 75], 'backgroundcolor', bcolor, 'title', 'Repitions');

    indButton = uicontrol(bgh, 'style', 'radiobutton', 'string', 'sweep #:', 'position', [10 15 75 25], 'backgroundcolor', bcolor);

    allButton = uicontrol(bgh, 'style', 'radiobutton', 'string', 'all', 'position', [150 15 50 25], 'backgroundcolor', bcolor, 'value', 1);

    superButton = uicontrol(bgh, 'style', 'radiobutton', 'string', 'all*', 'position', [200 15 45 25], 'backgroundcolor', bcolor);
    
    superduperButton = uicontrol(bgh, 'style', 'radiobutton', 'string', 'all traces', 'position', [260 15 100 25], 'backgroundcolor', bcolor);
    
    sweepBox = uicontrol(bgh, 'style', 'edit', 'position', [90 15 40 25], 'backgroundcolor', bcolor);

    uicontrol(fh, 'style', 'pushbutton', 'position', [150 10 90 30], 'string', 'Plot', 'callback', @okFun);

    uicontrol(fh, 'style', 'pushbutton', 'position', [250 10 90 30], 'string', 'Close', 'callback', 'close(gcf)' );


    function okFun(hObj, eventdata)
        testNum = str2num(get(testBox, 'string'));
        traceNum = str2num(get(traceBox, 'string'));
        if isempty(testNum)
           errordlg('must specify test number')
           return
        end
        if testNum > length(experimentData.test)
            errordlg('test number out of range')
            return
        end
        if any(get(bgh, 'SelectedObject') == [indButton, allButton, superButton])
            if isempty(traceNum)
                errordlg('must select trace number')
                return
            end
            rawData = ExtractRawData(experimentData, prefs, testNum, traceNum);
            rawData = rawData{testNum,traceNum};
        end
        figure('name', [prefs.cell_id num2str(testNum)])
        switch(get(bgh, 'SelectedObject'))
            case{indButton}
                sweepNum = str2num(get(sweepBox, 'string'));
                if isempty(sweepNum)
                    errordlg('must select sweep number')
                    return
                end
                plot(rawData(sweepNum,:));
            case{allButton}
                plot(1:length(rawData),rawData);                
            case{superButton}
                num_plot_rows = ceil(sqrt(size(rawData,1)));
                num_plot_cols = num_plot_rows;
                for sweep = 1:size(rawData,1)
                    subplot(num_plot_rows, num_plot_cols, sweep)
                    plot(rawData(sweep,:))
                    %ylim([-0.5 0.5])
                    
                end
            case{superduperButton}
                rawData = ExtractRawData(experimentData, prefs, testNum);
                rawData = rawData(testNum,:);
                num_plot_rows = ceil(sqrt(size(rawData,2)));
                num_plot_cols = num_plot_rows;
                for trace = 1:size(rawData,2)
                    subplot(num_plot_rows, num_plot_cols, trace)
                    data = rawData{trace};
                    plot(1:length(data),data)
                end
        end
%         if isempty(sweepNum)
%             rawData = ExtractRawData(experimentData, prefs, testNum, traceNum);
%             rawData = rawData{testNum,traceNum};
%             figure
%             num_plot_rows = ceil(sqrt(size(rawData,1)));
%             num_plot_cols = num_plot_rows;
%             for sweep = 1:size(rawData,1)
%                 subplot(num_plot_rows, num_plot_cols, sweep)
%                 plot(rawData(sweep,:))
%             end
%         end
%         if any([isempty(testNum), isempty(traceNum), isempty(sweepNum)])
%             errordlg('input field(s) empty')
%         else
%             rawData = ExtractRawData(experimentData, prefs, testNum, traceNum);
%             rawData = rawData{testNum,traceNum};
%             figure;
%             if sweepNum == 0
%                 plot(1:length(rawData),rawData);
%             else
%                 sweepNum = str2num(get(sweepBox, 'string'));
%                 plot(rawData(sweepNum,:));
%             end
%         end
    end
end