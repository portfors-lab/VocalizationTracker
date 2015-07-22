function ManipulateSynthFigures(plotHandles, nSections, inputs, stateEstimates)
%creates the figure that displays the plots of the state estimates.  places
%plots in tabs.  Allows modifcation of synthesized signal spectrogram and
%creation of synthesized vocalization audio file

%modified Amy Boyle 3/24/11

    bcolor = [0.8 0.8 0.8];
    
    %these two values are hard-coded constant(at least at the time this 
    %comment was written) througout application, does not affect data, 
    %only presentation.
    plotFreqSamples = 2^10;
    plotTimeSamples = 2^10;
    
    plotFigure = figure('units', 'normalized', 'position', [0.1 0.1 0.8 0.8]);

    warning('off', 'MATLAB:uitab:DeprecatedFunction');
    warning('off', 'MATLAB:uitab:MigratingFunction');
    warning('off', 'MATLAB:uitabgroup:DeprecatedFunction');
    warning('off', 'MATLAB:uitabgroup:MigratingFunction');
    
    hTabGroup = uitabgroup('v0', 'units', 'normalized', 'position', [0 0.2 0.7 0.8]);
    %need to have this renderer to use freezeColors
    set(gcf, 'renderer', 'zbuffer')
    
    for plotHandle = plotHandles
        name = get(plotHandle, 'name');
        warn = warning('off');
        thisTab = uitab('v0',hTabGroup, 'title', name); %v0 argument gets rid of warning, there seems to be no current version for R2010a
        warning(warn)
        kids = get(plotHandle, 'children'); %all components from plot figure
        set(kids, 'parent', thisTab) %place them on a tab
        %set and freeze colormap
        cmap = get(plotHandle, 'colormap');
%         if inputs.bools.invertColor
%             colormap(flipud(eval(inputs.colormap)));
%         else
%             colormap(inputs.colormap)
%         end
%         colormap(cmap);
%         for kid = kids'
%             if strcmp(get(kid, 'tag'), 'Colorbar')
%                  cbfreeze(kid);
%             else
%                freezeColors(kid);              
%             end
%         end     
        AddTabContextMenu(thisTab);
        close(plotHandle);
    end
    
    colormap(inputs.colormap)
    panHandle = uipanel(plotFigure, 'units', 'normalized', 'position', [0.75 0.22 0.22 0.75], 'BackgroundColor', bcolor);
    CreatePanelInputs(panHandle, inputs, nSections, 'text');
    components = findobj(panHandle, 'type', 'uicontrol');
    set(components, 'fontsize', 9, 'backgroundcolor', bcolor); %shrinks all fonts on figure to fit(hopefully; depends on screensize)
    
    boxLength = 0.1;
    boxHeight = 0.03;
    txtLength = 0.15;
    txtHeight = 0.05;
    
    uicontrol(plotFigure, 'style', 'text', 'units', 'normalized', 'position', [0.05 0.13 txtLength boxHeight], 'string', 'Frequency Multiplier', 'backgroundcolor', bcolor, 'horizontalAlignment', 'left');
    freqMultiplierH = uicontrol(plotFigure, 'style', 'edit', 'units', 'normalized', 'position', [0.05 0.1 boxLength boxHeight], 'string', '1');
    
    uicontrol(plotFigure, 'style', 'text', 'units', 'normalized', 'position', [0.05 0.04 txtLength boxHeight], 'string', 'Length Multiplier', 'backgroundcolor', bcolor, 'horizontalAlignment', 'left');
    lenMultiplierH = uicontrol(plotFigure, 'style', 'edit', 'units', 'normalized', 'position', [0.05 0.01 boxLength boxHeight], 'string', '1');
    
    uicontrol(plotFigure, 'style', 'text', 'units', 'normalized', 'position', [0.2 0.04 txtLength boxHeight], 'string', 'Harmonic Selection', 'backgroundcolor', bcolor, 'horizontalAlignment', 'left');
    harmSelectionH = uicontrol(plotFigure, 'style', 'edit', 'units', 'normalized', 'position', [0.2 0.01 boxLength boxHeight]);
    
    uicontrol(plotFigure, 'style', 'text', 'units', 'normalized', 'position', [0.2 0.13 boxLength boxHeight], 'string', 'State Type', 'backgroundcolor', bcolor);
    stateTypeH = uicontrol(plotFigure, 'style', 'popupmenu', 'units', 'normalized', 'position', [0.2 0.1 boxLength boxHeight], 'string', ['Filtered|Smoothed|Predicted'], 'value', inputs.stateIndex);
    
    freqModH = uicontrol(plotFigure, 'style', 'checkbox', 'units', 'normalized', 'position', [0.35  0.13 txtLength+0.05 txtHeight], 'string', 'Remove Frequency Modulations','BackgroundColor', bcolor, 'value', 0);
    ampModH = uicontrol(plotFigure, 'style', 'checkbox', 'units', 'normalized', 'position', [0.35  0.09 txtLength+0.05 txtHeight], 'string', 'Remove Amplitude Modulations','BackgroundColor', bcolor, 'value', 0);
    randPhaseH = uicontrol(plotFigure, 'style', 'checkbox', 'units', 'normalized', 'position', [0.35  0.05 txtLength+0.05 txtHeight], 'string', 'Randomize Phase','BackgroundColor', bcolor, 'value', 0);
    revSigH = uicontrol(plotFigure, 'style', 'checkbox', 'units', 'normalized', 'position', [0.35  0.01 txtLength+0.05 txtHeight], 'string', 'Reverse Signal','BackgroundColor', bcolor, 'value', 0);
    
    makeButton = uicontrol(plotFigure, 'style', 'pushbutton', 'units', 'normalized', 'position', [0.8 0.1 boxLength+0.01 boxHeight], 'string', 'make Audio file', 'enable', 'off', 'callback', @createSynthAudio);
    uicontrol(plotFigure, 'style', 'pushbutton', 'units', 'normalized', 'position', [0.8 0.05 boxLength+0.01 boxHeight], 'string', 'apply modifacations', 'callback', @viewMorphedSpectrogram);
    
    shiftHarmPanel = uipanel(plotFigure, 'units', 'normalized', 'position', [0.55 0.01 0.2 0.17], 'title', 'Move harmonic', 'backgroundcolor', bcolor);
    uicontrol(shiftHarmPanel, 'style', 'text', 'units', 'normalized', 'position', [0.01 0.8 0.2 0.1], 'string', 'start', 'backgroundcolor', bcolor);
    harmStartH = uicontrol(shiftHarmPanel, 'style', 'edit', 'units', 'normalized', 'position', [0.25 0.75 0.6 0.15]);
    uicontrol(shiftHarmPanel, 'style', 'text', 'units', 'normalized', 'position', [0.01 0.55 0.2 0.1], 'string', 'stop', 'backgroundcolor', bcolor);
    harmStopH = uicontrol(shiftHarmPanel, 'style', 'edit', 'units', 'normalized', 'position', [0.25 0.5 0.6 0.15]);
    uicontrol(shiftHarmPanel, 'style', 'text', 'units', 'normalized', 'position', [0.01 0.3 0.2 0.1], 'string', 'shift', 'backgroundcolor', bcolor);
    harmShiftH = uicontrol(shiftHarmPanel, 'style', 'edit', 'units', 'normalized', 'position', [0.25 0.25 0.6 0.15]);
    set([harmStartH harmStopH harmShiftH], 'tooltipstring', 'UNDER CONSTRUCTION');
    
    morphedSignal = [];
    
    function viewMorphedSpectrogram(hObj, evtData)
        modifiers.freqMultiplier = str2double(get(freqMultiplierH, 'string'));
        modifiers.lenMuliplier = str2double(get(lenMultiplierH, 'string'));
        modifiers.harmSelection = str2num(get(harmSelectionH, 'string'));
        modifiers.stateIndex = get(stateTypeH, 'value');
        modifiers.stateList = get(stateTypeH, 'string');
        modifiers.stateType = strtrim(modifiers.stateList(modifiers.stateIndex,:)); %trims spaces
        modifiers.freqMod = get(freqModH, 'value');
        modifiers.ampMod = get(ampModH, 'value');
        modifiers.randPhase = get(randPhaseH, 'value');
        modifiers.reverseSignal = get(revSigH, 'value');
        harmInterval = [str2double(get(harmStartH, 'string')) str2double(get(harmStopH,'string')) str2double(get(harmShiftH,'string'))];
        harmInterval = round(harmInterval*stateEstimates.ModelParameters.sampleRate);
        modifiers.moveHarm = harmInterval;
        
        [morphedStates morphedSignal] = MorphHarmonicStates(stateEstimates, ...
                                                'frequencyMultiplier', modifiers.freqMultiplier, ...
                                                'lengthMultiplier', modifiers.lenMuliplier, ...
                                                'harmonicSelection', modifiers.harmSelection, ...
                                                'stateType', modifiers.stateType, ...
                                                'removefrequencyModulation', modifiers.freqMod, ...
                                                'removeamplitudeModulation', modifiers.ampMod, ...
                                                'randomizePhase', modifiers.randPhase, ...
                                                'moveHarmonic', modifiers.moveHarm);
%                                                 'invertsignal', reverseSignal);
        if modifiers.reverseSignal
            morphedSignal = fliplr(morphedSignal);
        end
        [s, t ,f, specHandles] = NonparametricSpectrogram(morphedSignal, ...
                                             stateEstimates.ModelParameters.sampleRate, ...
                                             'nTimes',plotTimeSamples, ...
                                             'nFrequencies',plotFreqSamples, ...
                                             'windowDuration',inputs.windowDuration, ...
                                             'plotType',1,...
                                             'spectrogramType',2, ...
                                             'dbRange', inputs.values(1).dbRange,...
                                             'frequencyRange', inputs.values(1).frequencyRange, ...
                                             'visible', 'off');
        title(specHandles.haSpectrogram,['Morphed Signal (' modifiers.stateType ')'])
        %set up tab and context menus for this plot
        warn = warning('off');
        tab = uitab('v0', hTabGroup, 'title', 'morphed signal');
        warning(warn)
        kids = get(specHandles.figure, 'children');
        set(kids, 'parent', tab)
        cmap = get(specHandles.figure, 'colormap');
        figure(plotFigure);
        colormap(cmap);
        for kid = kids'
            if strcmp(get(kid, 'tag'), 'Colorbar')
                 cbfreeze(kid);
            else
               freezeColors(kid);              
            end
        end     
        AddTabContextMenu(tab);
        close(specHandles.figure)
        setappdata(tab, 'modifiers', modifiers);
        set(makeButton, 'enable', 'on');
    end    
    
    function createSynthAudio(hObj, evtData)
        filters = {'*.call'; '*.call1';'*.wav'; '*.kanwal'};
        [fileName filePath filterIndex] = uiputfile(filters);
%         fileName = 'testSynthAudio5.wav';
%         filePath = '/home/amy/Desktop/';
        tabs = (get(hTabGroup, 'children'));
        selectedTab = tabs(get(hTabGroup, 'SelectedIndex'));
        modifiers = getappdata(selectedTab, 'modifiers');
        if isempty(modifiers)
            msgbox('first select the tab of the modified signal you wish to synthesize')
        else
            [morphedStates morphedSignal] = MorphHarmonicStates(stateEstimates, ...
                                            'frequencyMultiplier', modifiers.freqMultiplier, ...
                                            'lengthMultiplier', modifiers.lenMuliplier, ...
                                            'harmonicSelection', modifiers.harmSelection, ...
                                            'stateType', modifiers.stateType, ...
                                            'removefrequencyModulation', modifiers.freqMod, ...
                                            'removeamplitudeModulation', modifiers.ampMod, ...
                                            'randomizePhase', modifiers.randPhase);
            if modifiers.reverseSignal
                morphedSignal = fliplr(morphedSignal);
            end
            if fileName ~= 0
                %add the stateType option to morphHarmonicStates...            
                WriteAudioData(morphedSignal, stateEstimates.ModelParameters.sampleRate, [filePath fileName]);
                success = 1;
                if success
                    msgbox('Audio file created successfully');
                end
            end
        end
    end

end