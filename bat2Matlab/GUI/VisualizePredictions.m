function VisualizePredictions(stimPath)

animalPath = [];
rootPath = pwd;
outPath = '';
%default saveing settings
saveLocation=pwd;
format = 'tiff';
resolution = 600;

% stimPath = '/home/amy/src/MATLAB/vocal_aux/2008Experiments/';
if exist('stimPath')
    if isempty(stimPath)
        stimPath = uigetdir(rootPath, 'Select folder containing audio stimulus files');
    end
else
    stimPath = uigetdir(pwd, 'Select folder containing audio stimulus files');
end

if ~isequal(stimPath(end), filesep)
   stimPath = [stimPath filesep];
end

tableIcon = imread('table.jpg');
traceIcon = imread('trace.jpg');

fh = figure('position', [300 250 425 450], 'resize', 'off', 'windowstyle', 'normal', 'menubar', 'none');
%fh = figure('position', [500 500 375 225], 'resize', 'off', 'windowstyle', 'modal');
bcolor = [0.8 0.8 0.8];

uicontrol(fh, 'style', 'text', 'position',  [100 400 250 30], 'string', 'Spike Predictions', 'fontsize', 14, 'backgroundcolor', bcolor);

uicontrol(fh, 'style', 'text', 'position', [15 350 200 20], 'string', 'Animal folder', 'horizontalalignment', 'left', 'backgroundcolor', bcolor);

pathBox = uicontrol(fh, 'style', 'edit', 'position', [15 325 225 25], 'string', animalPath, 'tooltipstring', 'folder where .pst and .raw files reside (same name as folder)');
uicontrol(fh, 'style', 'pushbutton', 'position', [240 325 100 25], 'string', 'Browse...', 'callback', @browseFun);

% uicontrol(fh, 'style', 'pushbutton', 'position', [350 325 25 25],'CData', tableIcon, 'tooltipstring', 'view file contents', 'callback', @exploreFun);
uicontrol(fh, 'style', 'pushbutton', 'position', [380 325 25 25],'CData', traceIcon, 'tooltipstring', 'view spike traces',  'callback', @previewFun);
uicontrol(fh, 'style', 'pushbutton', 'position', [115 275 25 25],'CData', tableIcon, 'tooltipstring', 'select test numbers', 'callback', @exploreFun);

uicontrol(fh, 'style', 'text', 'position', [15 275 100 20], 'string', 'Test numbers:', 'horizontalalignment', 'left', 'backgroundcolor', bcolor);
testNumBox =  uicontrol(fh, 'style', 'edit', 'position', [150 275 150 20], 'tooltipstring', 'Test numbers to run prediction on');
jhandle = findjobj(testNumBox);
set(jhandle, 'MouseClickedCallback', @testFieldClick);

% minBox = uicontrol(fh, 'style', 'edit', 'position', [115 275 100 20]);
% maxBox = uicontrol(fh, 'style', 'edit', 'position', [200 275 50 20]);

% uicontrol(fh, 'style', 'text', 'string', 'to', 'position', [170 275 20 20], 'backgroundcolor', bcolor);

bgh = uibuttongroup('parent', fh,  'units', 'pixels','position', [15 190 400 75], 'backgroundcolor', bcolor);% 'selectionchangefcn', @selectionFun);
sameRadio = uicontrol(bgh, 'style', 'radiobutton', 'string', 'Same Traces', 'position', [10 40 175 20], 'backgroundcolor', bcolor);
diffRadio = uicontrol(bgh, 'style', 'radiobutton', 'string', 'Different Traces', 'position', [125 40 150 20], 'backgroundcolor', bcolor);
allRadio = uicontrol(bgh, 'style', 'radiobutton', 'string', 'All Traces', 'position', [260 40 100 20], 'backgroundcolor', bcolor);
traceBox = uicontrol(bgh, 'style', 'edit', 'position', [25 10 150 20]);
jhandle = findjobj(traceBox);

set(jhandle, 'MouseClickedCallback', @traceFieldClick);

uicontrol(fh, 'style', 'text', 'position', [15 150 100 20], 'string', 'Threshold:', 'horizontalalignment', 'left', 'backgroundcolor', bcolor);
threshBox = uicontrol(fh, 'style', 'edit', 'position', [115 150 50 20], 'string', '0.2', 'tooltipstring', 'Spike peak threshold');

uicontrol(fh, 'style', 'text', 'position', [200 150 200 20], 'string', 'Model test nums:', 'horizontalalignment', 'left', 'backgroundcolor', bcolor);
modelBox = uicontrol(fh, 'style', 'edit', 'position', [315 150 100 20], 'tooltipstring', 'Test numbers used to generate prediction model (e.g. tuning curve)');
jhandleMB = findjobj(modelBox);
set(jhandleMB, 'MouseClickedCallback', @modelFieldClick);

txtCheckBox = uicontrol(fh, 'style', 'checkbox', 'position',[15 100 250 20], 'string', 'Save Figures', 'backgroundcolor', bcolor, 'callback', @checkFun);
% uicontrol(fh, 'style', 'text', 'position', [35 100 200 20], 'string', 'output location:', 'horizontalalignment', 'left', 'backgroundcolor', bcolor);
% txtPathBox = uicontrol(fh, 'style', 'edit', 'position', [35 75 225 25], 'string', outPath, 'enable', 'off');
% txtPathButton = uicontrol(fh, 'style', 'pushbutton', 'position', [260 75 100 25], 'string', 'Browse...', 'callback', @txtBrowseFun, 'enable', 'off');
saveOptButton = uicontrol(fh, 'style', 'pushbutton', 'position', [15 75 150 25], 'string', 'Save options', 'callback', @saveOptions, 'enable', 'off');

legendTickBox = uicontrol(fh, 'style', 'checkbox', 'position',[215 75 250 20], 'string', 'Include legend', 'backgroundcolor', bcolor);

okButton = uicontrol(fh, 'style', 'pushbutton', 'position', [200 10 90 30], 'string', 'OK', 'callback', @okFun);
uicontrol(fh, 'style', 'pushbutton', 'position', [300 10 90 30], 'string', 'Close', 'callback', 'close(gcf)');

experiment_data = [];
prefs = [];

    function browseFun(jObj,evtdata)
        animalPath = uigetdir(rootPath, 'folder containing data');
        if animalPath ~= 0
           set(pathBox, 'string', [animalPath filesep]);
           experiment_data = [];
        end
    end

    function txtBrowseFun(hObj, eventdata)
        thisPath = uigetdir(rootPath, 'Output text files location');
        set(txtPathBox, 'string', [thisPath filesep]);
    end

    function checkFun(hObj, eventdata)
        if(get(hObj, 'value'))
%             set(txtPathBox, 'enable', 'on');
%             set(txtPathButton, 'enable', 'on');
            set(saveOptButton, 'enable', 'on');
        else
%             set(txtPathBox, 'enable', 'off');
%             set(txtPathButton, 'enable', 'off');
            set(saveOptButton, 'enable', 'off');
        end
    end

    function exploreFun(hObj, eventData)
        %view file contents, get test numbers
        if isempty(animalPath)
            errordlg('! Must select animal file first, before able to select test numbers')
            return
        end
        if isempty(experiment_data)
            [prefs experiment_data] = GetExpData(animalPath);
            if isempty(prefs)
                return
            end
        end
        %ExploreTestData(experiment_data);
        testNums = SelectTestNumbers(experiment_data);
        if ~isempty(testNums)
            set(testNumBox, 'string', num2str(testNums'));
        end
    end

    function previewFun(jObj, eventData)
        if isempty(animalPath)
            errordlg('! Must select animal file first, before able to select test numbers')
            return
        end
        if isempty(experiment_data)
            [prefs experiment_data] = GetExpData(animalPath);
            if isempty(prefs)
                return
            end
        end
        if ~isempty(experiment_data)
            TraceChooser(prefs, experiment_data);
        end
    end

    function testFieldClick(hObj, evtdata)
        if get(evtdata, 'clickCount') == 2
            exploreFun
        end
    end

    function traceFieldClick(hObj, evtdata)
        if get(evtdata, 'clickCount') == 2
            selectedRadio = get(bgh, 'selectedobject');
            tnums = str2num(get(testNumBox, 'string'));
            switch selectedRadio
                case allRadio
                    return
                    %do nothing
                case sameRadio
                    if isempty(tnums)
                        WriteStatus('Select test numbers first');
                        return
                    end
                    traceCounts = [];
                    for t = tnums
                        traceCounts = [size(experiment_data.test(t).trace,2) traceCounts];
                    end
                    maxtrace = min(traceCounts);
                    liststr = {};
                    for tracenum = 1:maxtrace
                        liststr = [liststr {num2str(tracenum)}];
                    end
                    traceNums = listdlg('liststring', liststr, 'listsize', [150 150]);
                    if isempty(traceNums)
                        return
                    end
                    traceNums = num2str(traceNums);
                case diffRadio
                    if isempty(tnums)
                        WriteStatus('Select test numbers first', err);
                        return
                    end
                    traceNums = [];
                    for t = tnums
                        nums = num2str(SelectTraceNumbers(experiment_data, t)');
                        traceNums = [traceNums '; ' nums];
                        if isempty(nums)
                            return
                        end
                    end
                    traceNums = traceNums(2:end); %remove leading ;
                otherwise
                    return
            end
            set(traceBox, 'string', traceNums); 
        end
    end

    function modelFieldClick(hObj, evtdata)
        if get(evtdata, 'clickCount') == 2
            if isempty(animalPath)
                errordlg('! Must select animal file first, before able to select test numbers')
                return
            end
            if isempty(experiment_data)
                [prefs experiment_data] = GetExpData(animalPath);
                if isempty(prefs)
                    return
                end
            end
            %ExploreTestData(experiment_data);
            testNums = SelectTestNumbers(experiment_data);
            if ~isempty(testNums)
                set(modelBox, 'string', num2str(testNums'));
            end
        end
    end

    function saveOptions(hobject, eventdata)
        %allows user to select save parameters
        [saveLocation format resolution] = SaveDialog(saveLocation, format, resolution);
    end

    function okFun(hObj, evtdata)
    %load the experimental data
        if isempty(experiment_data)
            [prefs experiment_data] = GetExpData(animalPath);
            if isempty(prefs)
                return
            end
        end
        prefs.spike_time_peak_threshold = str2double(get(threshBox, 'string'));
        prefs.audio_directory = stimPath;
%         saveLocation = get(txtPathBox, 'string');
        saveFigs = get(txtCheckBox, 'value');

        testNums = str2num(get(testNumBox, 'string'));
        if isempty(testNums)
            errordlg('Select test numbers')
            return;
        end
                         
        %saveas(gcf,figname);
        
        selectedRadio = get(bgh, 'selectedobject');
        traceStr = get(traceBox, 'string');
        if selectedRadio == sameRadio
            traceNums = str2num(traceStr);
            if isempty(traceNums)
                errordlg('Select trace numbers')
                return;
            end
        elseif selectedRadio == diffRadio
            traceCell = textscan(traceStr, '%s', 'delimiter', ';');
            traceCell = cellfun(@str2num, traceCell{1}, 'uniformoutput', false);
            traceIdx = 1;
            if isempty(traceCell{1})
                errordlg('Select trace numbers')
                return;
            end
        end
        
        includeLegend = get(legendTickBox, 'value');
        
        modelNums = str2num(get(modelBox, 'string'));
        if isempty(modelNums)
            train_data = testNums;
        else
            train_data = modelNums;
        end
        
        %Create the model
        [constrModel model] = CreateModel(experiment_data,prefs,train_data);
    %   model = CreateModel(experiment_data,prefs,train_data);
        numStr = num2str(train_data);
%         figname = [saveLocation prefs.cell_id '_model1_' numStr '.pdf'];
        
        formatStr = ['-d' format];
        resolutionStr = ['-r' num2str(resolution)];
        for testNum = testNums
            if selectedRadio == diffRadio
                traceNums = traceCell{traceIdx};
                traceIdx = traceIdx +1;
                if isempty(traceNums)
                    WriteStatus('Invalid trace number input','orange');
                    continue
                end
            elseif selectedRadio == allRadio
                traceNums = 1:length(experiment_data.test(testNum).trace);
            end
                
            for traceNum = traceNums
                VisualizeTracePredictions(experiment_data, ...
                                          prefs, ...
                                          model, ...
                                          testNum, ...
                                          traceNum);
                if ~includeLegend
                    legend('off')
                end
                testStr = num2str(testNum);
                traceStr = num2str(traceNum);
                figname = [saveLocation prefs.cell_id '_vocal1_' testStr '_' traceStr];  
                
                if saveFigs 
                    if strmatch(format, 'fig')
                        saveas(gcf, figname, 'fig');
                    else
                        print(gcf, figname, formatStr, resolutionStr)
                    end
                    close gcf
                end
            end
        end
    end
end