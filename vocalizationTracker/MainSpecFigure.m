function MainSpecFigure(filePath, inputs, inHandles, plotdata)
% creates the main GUI figure for the vocalization tracker.
%
% MainSpecFigure prompts user for input file and creates a spectrogram.

% MainSpecFigure(filePath, inputs, inHandles, plotdata) takes in parameters
%   to load a previously saved session. All these parameters are saved in
%   the subfunction saveProj
%filePath       location of audio file
%inputs         structure containing the inputs for the uicontrols
%inHandles      structure containing handles and other plot info
%plotdata       struture containing modified signal and labels

fileTypes = {'*.call; *.call1; *.kanwal; *.wav', 'Audio Files (*.call, *.kanwal, *.wav)'; '*.*', 'All Files'};
if ~exist('filePath', 'var')
   [fileName, filePath, filter] = uigetfile(fileTypes);
   filePath = [filePath fileName];
    if fileName == 0
        return
    end    
end
callName = filePath;

%check to see if file is too large (only for wav)
maxfilesizeMB = 20;
fileinfo = dir(filePath);
if (fileinfo.bytes/1048576) > maxfilesizeMB
    [signal sampleRate realtime] = ParsePartialAudio(filePath);
else
    [signal sampleRate] = ParseAudioData(filePath);
    realtime = [0 length(signal)/sampleRate];
end
signalLength = length(signal);
signal = signal - mean(signal);
cropSignal = signal;
cropTime = realtime;
signalTime = signalLength/sampleRate;

%==========================================================================
%Display initial spectrogram
%==========================================================================
if ~exist('inputs', 'var')
    inputs = SetDefaults(1, sampleRate);
end
if exist('plotdata', 'var')
    cropSignal = plotdata.cropSignal;
end

%these two values match the ones in PlotHarmonicsStates ...
plotFreqSamples = 2^10;
plotTimeSamples = 2^10;

%window duration: time resolution & frequency resolution tradeoff. narrow
%window = fine time, but coarse frequency; narrow windows have short time
%duration, but a wide bandwidth. And the inverse.
windowDuration = inputs.windowDuration;
frequencyRange = inputs.values.frequencyRange; %use set frequency range for comparison between spectrograms
dbRange = inputs.values.dbRange;

inputs.sampleRate = sampleRate;

[spectrogram,t,f, handles] = NonparametricSpectrogram(cropSignal, ...
                                             sampleRate, ...
                                             'nTimes',plotTimeSamples, ...
                                             'nFrequencies',plotFreqSamples, ...
                                             'windowDuration',windowDuration, ...
                                             'plotType',1,...
                                             'spectrogramType',2, ...
                                             'dbRange', dbRange,...
                                             'frequencyRange', frequencyRange, ...
                                             'colormap', inputs.colormap, ...
                                             'visible', 'off');
fig = gcf;  
set(fig, 'name', callName);

if numel(get(0, 'MonitorPositions')) > 4
   set(fig, 'units', 'normalized','position', [1.1 0 1 1], 'menubar', 'none')
%    disp('two monitors')
else
   set(fig, 'units', 'normalized','position', [0 0 1 1], 'menubar', 'none') 
%    disp('one monitor')
end

%menu bar
% SetupMenubar(fig);
fileMenu = uimenu(fig, 'label', 'File');
uimenu(fileMenu, 'label', 'Save Project As...', 'callback', @saveProj);
uimenu(fileMenu, 'label', 'Exit', 'callback', 'close all');
fixMenu = uimenu(fig, 'label', 'Spectrogram Tools');
uimenu(fixMenu, 'label', 'Restore original spectrogram', 'callback', @resetSpecFun);
uimenu(fixMenu, 'label', 'Load different Audio file', 'callback', @loadNewSpec);
uimenu(fixMenu, 'label', 'write current to wav', 'callback', @outputWav);

%save original size/positions of spectrogram
specPosBig = get(handles.haSpectrogram, 'position');
signalPosBig = get(handles.haSignal, 'position');
colorbarPosBig = get(handles.haColorbar, 'position');
psdPosBig = get(handles.haPowerSpectralDensity, 'position');

%reposition spectrogram to make room for ui controls
specPos = [0.12 0.3 0.6 0.65];
signalPos = [0.12 0.2 0.6 0.09];
colorbarPos = [0.105 0.3 0.01 0.65];
psdPos = [0.055 0.3 0.04 0.65];
set(handles.haSpectrogram, 'position', specPos);
set(handles.haSignal, 'position', signalPos);
set(handles.haColorbar, 'position', colorbarPos);
set(handles.haPowerSpectralDensity, 'position', psdPos);

%set the time to real, not relative
timelabels = get(handles.haSignal, 'xticklabel');
timelabels = realtime(1) + str2num(timelabels);
set(handles.haSignal, 'xticklabel', num2str(timelabels));

%context menu
cmenu = uicontextmenu;
uimenu(cmenu, 'label', 'save plot', 'callback', @savePlot);
uimenu(cmenu, 'label', 'change colormap', 'callback', @changeCmap);
uimenu(cmenu, 'label', 'change Db range', 'callback', @changeDb);
set(handles.haSpectrogram, 'uicontextmenu', cmenu);
set(get(handles.haSpectrogram,'children'), 'uicontextmenu', cmenu)

% ==========================================================================
% GUI controls
% ==========================================================================
set(fig, 'visible', 'on'); 
bcolor = [0.8 0.8 0.8];

% rangeButton = uicontrol('style', 'pushbutton', 'units', 'normalized', 'position', [0.18 0.1 0.1 0.04], 'string', 'Time Range', 'tooltipstring', 'define time sub range', 'callback', @rangeFun);
sectionButton = uicontrol('style', 'pushbutton', 'units', 'normalized', 'position', [0.18 0.1 0.1 0.04], 'string', 'Section', 'tooltipstring', 'divide spectrogram into sections', 'callback', @defineSections);
% cropButton = uicontrol('style', 'pushbutton', 'units', 'normalized', 'position', [0.05 0.04 0.1 0.04], 'string', 'crop', 'tooltipstring', 'crop to current selection','enable', 'off', 'callback', @newSpecFun);
cutButton = uicontrol('style', 'pushbutton', 'units', 'normalized', 'position', [0.05 0.1 0.1 0.04], 'string', 'cut', 'tooltipstring', 'cut out middle sections', 'callback', @cutoutSections);

panHandle = uipanel(fig, 'units', 'normalized', 'position', [0.75 0.22 0.22 0.75], 'BackgroundColor', bcolor);

checkPanHandle = uipanel(fig, 'units', 'normalized', 'position', [0.5 0.01 0.45 0.15], 'BackgroundColor', bcolor, 'title', 'Plot:');
checkHandles = CreateCheckboxPanel(checkPanHandle, inputs);

uicontrol(fig, 'style', 'text','units', 'normalized',  'position', [0.75 0.17 0.1 0.02],'string', 'State Type:', 'BackgroundColor', bcolor, 'horizontalalignment', 'left');
stateTypeH = uicontrol(fig, 'style', 'popupmenu','units', 'normalized', 'position', [0.81 0.175 0.1 0.02],'string', ['Filtered|Smoothed|Predicted'], 'value', inputs.stateIndex);

booyaButton = uicontrol(fig, 'style', 'pushbutton', 'units', 'normalized', 'position', [0.32 0.02 0.15 0.1], 'string', 'Generate Plots', 'callback', @generateFun);      

%==========================================================================
% Load saved project drawing, save handles
%==========================================================================

if exist('inHandles', 'var')
    if ~isempty(inHandles.sectPoints)
        handles.sectLines = [];
        handles.sectPoints = inHandles.sectPoints;
        ymax = get(handles.haSpectrogram, 'ylim');
        for l = 2:length(handles.sectPoints)-1
            sectTime = handles.sectPoints(l)/sampleRate;
            thisLine = line([sectTime sectTime],ymax, 'color', 'white', 'parent', handles.haSpectrogram);
            handles.sectLines = [handles.sectLines thisLine];
        end
        nSections = length(handles.sectPoints) -1;
        inputHandles = CreatePanelInputs(panHandle, inputs, nSections, 'edit');
        inputs.values(2:end) = [];
    else
        handles.sectPoints = [];        %vector for section boundry points
        inputHandles = CreatePanelInputs(panHandle, inputs, 1, 'edit');
    end
else
    handles.sectPoints = [];        %vector for section boundry points
    inputHandles = CreatePanelInputs(panHandle, inputs, 1, 'edit');
end
inputHandles.stateTypeH = stateTypeH;
handles.rangeLines = [];
handles.cutLines = [];
guidata(fig,handles);
% Functions================================================================

    function h = newSpecFun(h, includeSections)
        %Generates a new spectrogram based on the range points picked by the user
%         h = guidata(fig);
%         newRange = h.specRange;
%         if diff(newRange) == 0
%             errordlg('No Data in the selected range, choose lines inside of spectrogram')
%         end
        
%         if newRange(2) < length(cropSignal)
%             cropSignal = cropSignal(newRange(1):newRange(2));
%         end
        delete([h.haSpectrogram h.haSignal h.haPowerSpectralDensity h.haColorbar]);
        sectPoints = h.sectPoints;
        [spectrogram,t,f, h] = NonparametricSpectrogram(cropSignal, ...
                                             sampleRate, ...
                                             'nTimes',plotTimeSamples, ...
                                             'nFrequencies',plotFreqSamples, ...
                                             'windowDuration',windowDuration, ...
                                             'frequencyRange', frequencyRange, ...
                                             'plotType',2,...
                                             'spectrogramType',2, ...
                                             'dbRange', dbRange, ...
                                             'colormap', inputs.colormap);
        xmax = get(h.haSpectrogram, 'xlim');
        numSamples = length(cropSignal);
        if includeSections
             h.sectLines = [];
             h.sectPoints = sectPoints;
             if ~isempty(sectPoints)
                 for ipoint = sectPoints(2:length(sectPoints)-1);
                     thisLine = line([ipoint*(xmax(2)/numSamples) ipoint*(xmax(2)/numSamples)],ymax, 'color', 'white', 'parent', h.haSpectrogram);
                    h.sectLines = [h.sectLines thisLine];
                 end
             end
         else  
            h.sectPoints = [];
            h.rangeLines = [];
            h.specRange = [1 length(cropSignal)];
        end
        nSections = abs(length(h.sectPoints) - 1);
        inputHandles = CreatePanelInputs(panHandle, inputs, nSections, 'edit', inputHandles);
        drawnow;
        set(h.haSpectrogram, 'position', specPos);
        set(h.haSignal, 'position', signalPos);
        set(h.haColorbar, 'position', colorbarPos);
        set(h.haPowerSpectralDensity, 'position', psdPos);
        set(get(h.haSpectrogram,'children'), 'uicontextmenu', cmenu)
        %set the time to real, not relative
        timelabels = get(h.haSignal, 'xticklabel');
        timelabelstep = (cropTime(2)-cropTime(1))/(length(timelabels)-1);
        timelabels = cropTime(1):timelabelstep:cropTime(2);
        set(h.haSignal, 'xticklabel', timelabels);
        guidata(fig, h);
    end

    function resetSpecFun(hObject, eventdata)
        %Generates a spectrogram using the original range of data
        cropSignal = signal;
        cropTime = realtime;
        h = guidata(hObject);
        delete([h.haSpectrogram h.haSignal h.haPowerSpectralDensity h.haColorbar]);
        if ~isempty(h.sectPoints)
            inputHandles = CreatePanelInputs(panHandle, inputs, 1, 'edit', inputHandles);
        end
        [spectrogram,t,f, h] = NonparametricSpectrogram(signal, ...
                                             sampleRate, ...
                                             'nTimes',plotTimeSamples, ...
                                             'nFrequencies',plotFreqSamples, ...
                                             'windowDuration',windowDuration, ...
                                             'plotType',2,...
                                             'spectrogramType',2, ...
                                             'dbRange', dbRange, ...
                                             'colormap', inputs.colormap);...'frequencyRange',frequency_range); 
        set(h.haSpectrogram, 'position', specPos);
        set(h.haSignal, 'position', signalPos);
        set(h.haColorbar, 'position', colorbarPos);
        set(h.haPowerSpectralDensity, 'position', psdPos);
        set(get(h.haSpectrogram,'children'), 'uicontextmenu', cmenu);
        %set the time to real, not relative
        timelabels = get(h.haSignal, 'xticklabel');
        timelabelstep = (realtime(2)-realtime(1))/(length(timelabels)-1);
        timelabels = realtime(1):timelabelstep:realtime(2);
        set(h.haSignal, 'xticklabel', timelabels);
        h.sectPoints = [];
        h.rangeLines = [];
        h.specRange = [1 signalLength];
        guidata(hObject, h);
    end
       
    function defineSections(hObject, eventdata)
       %Allows the user to define different sections of the spectrogram to
       %input different paramters to
%        set(hObject, 'enable', 'inactive');
       %delete lines from previous sectioning
       h = guidata(hObject);        
       if isfield(h, 'sectLines') && any(ishandle(h.sectLines))
           delete(h.sectLines);
       end           
       h.sectLines = [];        %handles for lines on plot
       h.sectPoints = [];
       xmax = get(h.haSpectrogram, 'xlim');
       ymax = get(h.haSpectrogram, 'ylim');
       numSamples = length(cropSignal);
       
       cursorLine = line([1 1],ymax, 'color', 'white', 'parent', h.haSpectrogram);
       set(fig, 'WindowButtonMotionFcn', @drawLines);
       set(fig, 'WindowKeyPressFcn', @cancelDraw);
       set(fig, 'WindowButtonDownFcn', @commitLines);
       
       setptr(fig, 'crosshair');
        
        function drawLines(src, evtData)
            %moves the line whenever the mouse moves in the axes window
            cp = get(h.haSpectrogram, 'currentpoint');
            set(cursorLine, 'xdata', [cp(1,1) cp(1,1)]);
        end
        
        function commitLines(src, evtData)
            %draws a static line when the user clicks
            cp = get(h.haSpectrogram, 'currentpoint');
            thisLine = line([cp(1) cp(1)],ymax, 'color', 'white', 'parent', h.haSpectrogram);
            h.sectLines = [h.sectLines thisLine];
            h.sectPoints = [h.sectPoints round(cp(1)/(xmax(2)/numSamples))];
            guidata(hObject, h);
        end
        
        function cancelDraw(src, evtData)
            %when user cancels cropping with esc key
            if strcmp(evtData.Key, 'escape')
                set(fig, 'WindowButtonMotionFcn', '');
                set(fig, 'WindowKeyPressFcn', '');
                set(fig, 'WindowButtonDownFcn', '');
                set(hObject, 'enable', 'on');
                delete(cursorLine);
                setptr(fig, 'arrow');
                for l = h.sectLines
                    if ishandle(l)
                        delete(l);
                    end
                end
                h.sectLines = [];
                h.sectPoints = [];
                inputHandles = CreatePanelInputs(panHandle, inputs, 1, 'edit', inputHandles);
            elseif strcmp(evtData.Key, 'return')
                delete(cursorLine);
                setptr(fig, 'arrow');
                if(all(0 < h.sectPoints) && all(h.sectPoints < numSamples))
                    h.sectPoints = [1 h.sectPoints numSamples];
                    h.sectPoints = sort(h.sectPoints);
                    nSections = length(h.sectPoints) -1;
                else
                    msgbox('Selection out of range, must select a point in the figure')
                    h.sectPoints = [];
                    nSections = 1;
                    for l = h.sectLines
                        if ishandle(l)
                            delete(l);
                        end
                    end
                end
                set(fig, 'WindowButtonMotionFcn', '');
                set(fig, 'WindowKeyPressFcn', '');
                set(fig, 'WindowButtonDownFcn', '');
                set(hObject, 'enable', 'on');                
                inputHandles = CreatePanelInputs(panHandle, inputs, nSections, 'edit', inputHandles);
            end
            guidata(hObject, h);
        end       
       guidata(hObject, h);
    end

    function cutoutSections(hObject, eventdata)
       %Allows the user to define different sections of the spectrogram to
       %cut out
       
%        set(hObject, 'enable', 'inactive');
       %delete lines from previous sectioning
       h = guidata(hObject);        
       if isfield(h, 'cutLines') && any(ishandle(h.cutLines))
           delete(h.cutLines);
       end        
       h.cutLines = [];        %handles for lines on plot
       h.cutPoints = [];

       xmax = get(h.haSpectrogram, 'xlim');
       ymax = get(h.haSpectrogram, 'ylim');
       numSamples = length(cropSignal);
       
       cursorLine = line([1 1],ymax, 'color', 'red', 'parent', h.haSpectrogram);
       
       set(fig, 'WindowButtonMotionFcn', @drawLines);
       set(fig, 'WindowButtonDownFcn', @commitLines);
       set(fig, 'WindowKeyPressFcn', @cancelDraw);

       setptr(fig, 'crosshair');
%        set(fig, 'DefaultLineParent', h.haSpectrogram);
        
       function drawLines(src, evtData)
           %moves the line whenever the mouse moves in the axes window
           cp = get(h.haSpectrogram, 'currentpoint');
           set(cursorLine, 'xdata', [cp(1,1) cp(1,1)]);
       end
        
       function commitLines(src, evtData)
           %draws a static line when the user clicks
           cp = get(h.haSpectrogram, 'currentpoint');
           if cp(1) < xmax(1)
               cp(1) = xmax(1);
           elseif cp(1) > xmax(2)
               cp(1) = xmax(2);
           end
           thisLine = line([cp(1) cp(1)],ymax, 'color', 'red', 'parent', h.haSpectrogram);
           h.cutLines = [h.cutLines thisLine];
           thisPoint = round(cp(1)/(xmax(2)/numSamples));
           if thisPoint == 0
               thisPoint = 1;
           end
           h.cutPoints = [h.cutPoints thisPoint];
           if mod(length(h.cutLines),2) == 0
               prevx = get(h.cutLines(end-1), 'xdata');
               if prevx > cp(1)
                   h.cutPoints(end-1:end) = [h.cutPoints(end) h.cutPoints(end-1)];
               end
               xline = line([prevx(1) cp(1)],ymax, 'color', 'red', 'parent', h.haSpectrogram);
               h.cutLines = [h.cutLines xline];
               xline = line([cp(1) prevx(1)],ymax, 'color', 'red', 'parent', h.haSpectrogram);
               h.cutLines = [h.cutLines xline];
           end
%            guidata(hObject, h);
       end
       
       function cancelDraw(src, evtData)
            %when user cancels cropping with esc key
            if strcmp(evtData.Key, 'escape')
                set(fig, 'WindowButtonMotionFcn', '');
                set(fig, 'WindowKeyPressFcn', '');
                set(fig, 'WindowButtonDownFcn', '');
                set(hObject, 'enable', 'on');
                delete(cursorLine);
                setptr(fig, 'arrow');
                for l = h.cutLines
                    if ishandle(l)
                        delete(l);
                    end
                end
                h.cutLines = [];
                h.cutPoints = [];
                nSections = 1;
                inputHandles = CreatePanelInputs(panHandle, inputs, nSections, 'edit', inputHandles);
            elseif strcmp(evtData.Key, 'return')
                delete(cursorLine);
                setptr(fig, 'arrow');
                set(fig, 'WindowButtonMotionFcn', '');
                set(fig, 'WindowKeyPressFcn', '');
                set(fig, 'WindowButtonDownFcn', '');
                set(hObject, 'enable', 'on');
                if mod(length(h.cutPoints),2) ~= 0 %remove partnerless line
                    h.cutPoints(end) = [];
                end
                if(all(0 < h.cutPoints) && all(h.cutPoints <= numSamples))
%                     h.cutPoints = [1 h.cutPoints numSamples];
%                     h.cutPoints = sort(h.cutPoints);
                else
                    msgbox('Selection out of range, must select a point in the figure')
                    h.cutPoints = [];
                    set(hObject, 'enable', 'on');
                   return 
                end
%                 nSections = length(h.cutPoints) -1;
%                 inputHandles = CreatePanelInputs(panHandle, inputs, nSections, 'edit', inputHandles);
                
                prevLength = length(cropSignal);
                inds2clear = [];
                for ctpoint = 1:2:length(h.cutPoints)
                    inds2clear = [inds2clear h.cutPoints(ctpoint):h.cutPoints(ctpoint+1)];
                end
                cropSignal(inds2clear) = [];
                if h.cutPoints(1) == 1
                    %add rise time
                    cropSignal = AddRiseFall(cropSignal, sampleRate, 'risetime', 0.001, 'falltime', 0);
                end
                if h.cutPoints(end) == prevLength
                    %add fall time
                    cropSignal = AddRiseFall(cropSignal, sampleRate, 'risetime', 0, 'falltime', 0.001);
                end
                
                newTimeDuration = length(cropSignal)/sampleRate;
                if any(h.cutPoints ==1) %removed chunk from beginning, adjust time labels to show
                    idx = h.cutPoints(h.cutPoints==1);
                    startcut = length(h.cutPoints(idx):h.cutPoints(idx+1));
                    cropTime(1) = cropTime(1)+(startcut/sampleRate);
                end
                cropTime(2) = cropTime(1)+newTimeDuration;
                h = newSpecFun(h,0);
            end
%             guidata(hObject, h);
        end       
       guidata(hObject, h);
    end

    function generateFun(hObj, evtData)
%        set(hObj, 'enable', 'off')
       drawnow;
       h = guidata(hObj);
       userInputs = RetrieveValues(inputHandles, checkHandles);
       userInputs.colormap = inputs.colormap;
%        userInputs.values.dbRange = inputs.values.dbRange;
       GenerateSynthFigures(cropSignal, sampleRate, h.sectPoints, userInputs)
       set(hObj, 'enable', 'on')
    end

    function savePlot(ho, evtData)
        %saves spectrogram by itself: creates a new figure, moves
        %the axes to it, saves, then moves the axes back and
        %closes the figure.
        handles = guidata(fig);
        saveFormats = {'*.fig'; '*.jpg'; '*.tiff'; '*.pdf'; '*.png'}; 
        [saveName savePath filterIndex] = uiputfile(saveFormats);  
        if filterIndex ~= 0
            axs = findobj(fig, 'type', 'axes'); %the axes that make up the spectrogram
            saveFig = figure('visible', 'off'); %create a temporary figure for saving without the uicontrols
            set(axs, 'parent', saveFig); %move to new figure
            %reposition the spectrogram so that it uses all space
            set(handles.haSpectrogram, 'position', specPosBig);
            set(handles.haSignal, 'position', signalPosBig);
            set(handles.haColorbar, 'position', colorbarPosBig);
            set(handles.haPowerSpectralDensity, 'position', psdPosBig);
            colormap(inputs.colormap);
            saveas(saveFig, [savePath saveName]); %save figure 
            %reposition spectrogram to fit back on the figure with buttons
            set(handles.haSpectrogram, 'position', specPos);
            set(handles.haSignal, 'position', signalPos);
            set(handles.haColorbar, 'position', colorbarPos);
            set(handles.haPowerSpectralDensity, 'position', psdPos);
            set(axs, 'parent', fig); %move back
            close(saveFig)
        end
    end

    function outputWav(ho, evtData)
        [outName outPath] = uiputfile('.wav', 'save location');
        if ~isempty(outName)
            wavwrite(cropSignal, sampleRate, [outPath outName])
        end
    end

    function saveProj(hObj, evtData)
        [saveFile, savePath] = uiputfile('*.mat');
        if ~isempty(saveFile)
            h = guidata(hObj);
            userInputs = RetrieveValues(inputHandles, checkHandles);
            userInputs.colormap = inputs.colormap;
            plotdata.cropSignal = cropSignal;
            plotdata.cropTime = cropTime;        
            save([savePath saveFile], 'userInputs', 'h', 'filePath', 'plotdata');
        end
    end

    function loadNewSpec(hObj, evtData)
        [fileName, filePath, filter] = uigetfile(fileTypes);
        if filter ~= 0
            filePath = [filePath fileName];
            close(fig)
            MainSpecFigure(filePath)
        end
    end
     
    function changeCmap(hObject, ed)
        h = guidata(hObject);
        selectedAxes = findall(fig, 'type', 'axes');
        cmapList = {'jet', 'bone', 'hsv', 'hot', 'cool', 'gray'};
%         cmapIdx = find(strcmp(cmapList, inputs.colormap));
%         [cmapIdx, ok] = listdlg('liststring', cmapList, 'selectionmode', 'single', 'listsize', [150 105], 'name', 'Select colomap');
        [cmapIdx, invCmap] = ColormapDlg('liststring', cmapList,'name', 'Select Colormap');  
        if ~isempty(cmapIdx)
            if invCmap
                inputs.colormap = flipud(eval(cmapList{cmapIdx}));
            else
                inputs.colormap = cmapList{cmapIdx};
            end
            colormap(inputs.colormap);
        end
        guidata(hObject,h)
    end

    function changeDb(hObject,~)
        h = guidata(hObject); 
        dbRange = inputdlg('New Db range: ', 'Db Range', 1, {num2str(dbRange)});
        if ~isempty(dbRange)
            set(inputHandles.dbRangeH, 'string', dbRange{1})
            dbRange = str2num(dbRange{1});
            inputs.values.dbRange = dbRange;
            h = newSpecFun(h,1);
        end
    
    end
end