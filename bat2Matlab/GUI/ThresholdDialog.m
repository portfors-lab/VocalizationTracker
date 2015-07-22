function thresholds = ThresholdDialog(testNums, thresholds)
%presents a dialog box to get spike threshold value for each test

%Amy Boyle Jan 2011

    bcolor = [0.8 0.8 0.8];
    
    if isempty(thresholds)
        thresholds = ones(length(testNums),1)*0.2; %default value
    end
    figLength = 75 + length(testNums)*25;
    fh = figure('position', [300 250 325 figLength], 'resize', 'off', 'windowstyle', 'modal');
    titlePos = [50 figLength-40 150 30];
    uicontrol(fh, 'style', 'text', 'position',  titlePos, 'string', 'Thresholds', 'fontsize', 14, 'backgroundcolor', bcolor);
    yPos = figLength - 65;
    for idx = 1:length(testNums)
        uicontrol(fh, 'style', 'text', 'position', [145 yPos 75 20], 'string', ['Test ' num2str(testNums(idx))], 'horizontalalignment', 'left', 'backgroundcolor', bcolor);
        thresholdInputs(idx) = uicontrol(fh, 'style', 'edit', 'position', [225 yPos 50 20], 'string', num2str(thresholds(idx)));
        yPos = yPos - 25;
    end
    
    allInput = uicontrol(fh, 'style', 'edit', 'position', [25 figLength/2 85 20], 'string', '0.2');
    uicontrol(fh, 'style', 'pushbutton', 'position', [25 figLength/2-25 85 20], 'string', 'apply to all', 'callback', @applyFun);

    uicontrol(fh, 'style', 'pushbutton', 'position', [150 10 50 20], 'string', 'OK', 'callback', @okFun);
    
    uicontrol(fh, 'style', 'pushbutton', 'position', [220 10 50 20], 'string', 'Close', 'callback', 'close(gcf)');
     
    uiwait
    
    function applyFun(hObj, eventData)
        allVal = get(allInput, 'string');
        set(thresholdInputs, 'string', allVal);
    end
    
    function okFun(hObj, eventData)
       thresholds =  str2double(get(thresholdInputs, 'string'));
       close(fh);
    end
end