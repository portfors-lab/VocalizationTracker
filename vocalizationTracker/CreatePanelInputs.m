function pHandles = CreatePanelInputs(panelHandle, inputs, ntabs, txtStyle, pHandles)
% creates/displays the user input controls for the main panel

%tooltip strings for the controls
freqminStr = sprintf('lowest frequency to display');
freqmaxStr = sprintf('higest frequency to display');
dbRangeStr = sprintf('allowable intensity range');
winDurStr = sprintf('Time resolution & frequency resolution tradeoff.\nNarrow window has fine time, but coarse frequency,\n and visa versa.');
nharmStr = sprintf('Number of harmonics present in the spectrogram. \n specify all even if you wish to later exclude harmonics');
mvarStr = sprintf('amount of variance present, \nglobal control on spectrogram variance');
freqcoStr = sprintf('values less than 1 bias the frequency toward \n the specified mean frequency');
ampcoStr = sprintf('values less than 1 bias the amplitude towared a constant value');
freqmeanStr  = sprintf('observed mean frequency of fundamental harmonic');
initphaseStr = sprintf('');
initfreqStr = sprintf('starting frequency of the fundamental harmonic');
initampStr = sprintf('');
varphaseStr = sprintf('');
varfreqStr = sprintf('Increase for a spectrogram that changes rapidly, \n decrease for one that changes slowly over time');
varampStr = sprintf('');

    column1 = 0.0686;
    column2 = 0.5379;
    column3 = 305;
    bcolor = [0.8 0.8 0.8];
    boxLength = 0.36;
    boxHeight = 0.045;
    txtHeight = 0.04;
%     freqMultiplier = 1000;

freqMultiplier = 1;

    mouseIcon = imread('grey_pointer.jpeg');
    
    if ~exist('pHandles', 'var') %if pHandles exists, it just means we are updating tabs below
        prevTabs = 0;
        
        uicontrol(panelHandle, 'style', 'text', 'units', 'normalized', 'position', [column1 0.94 boxLength*2 txtHeight], 'string', 'Frequency Range','BackgroundColor', bcolor, 'horizontalalignment', 'left');
        uicontrol(panelHandle, 'style', 'text', 'units', 'normalized', 'position', [column1 0.91 boxLength txtHeight], 'string', 'min:','BackgroundColor', bcolor, 'horizontalalignment', 'left');
        pHandles.freqminH = uicontrol(panelHandle, 'units', 'normalized', 'style', txtStyle, 'position', [column1 0.8750 boxLength boxHeight], 'string', inputs.values(1).frequencyRange(1)/freqMultiplier, 'tooltipstring', freqminStr);
        uicontrol(panelHandle, 'style', 'text',  'units', 'normalized','position', [column2 0.91 boxLength txtHeight], 'string', 'max:','BackgroundColor', bcolor, 'horizontalalignment', 'left');
        pHandles.freqmaxH = uicontrol(panelHandle, 'style', txtStyle, 'units', 'normalized','position', [column2 0.8750 boxLength boxHeight], 'string', inputs.values(1).frequencyRange(2)/freqMultiplier,'tooltipstring', freqmaxStr);    
        uicontrol(panelHandle, 'style', 'text', 'units', 'normalized', 'position', [column1 0.82 boxLength txtHeight], 'string', 'Decibel Range:','BackgroundColor', bcolor, 'horizontalalignment', 'left');
        pHandles.dbRangeH = uicontrol(panelHandle, 'style', 'text', 'units', 'normalized', 'position', [column1 0.7817 boxLength boxHeight], 'string', inputs.values(1).dbRange, 'tooltipstring', dbRangeStr);  
        uicontrol(panelHandle, 'style', 'text', 'units', 'normalized', 'position', [column2 0.82 boxLength+0.06 txtHeight], 'string', 'Window Duration:','BackgroundColor', bcolor, 'horizontalalignment', 'left');
        pHandles.winDurH = uicontrol(panelHandle, 'style', txtStyle, 'units', 'normalized', 'position', [column2 0.7817 boxLength boxHeight], 'string', inputs.windowDuration(1),'tooltipstring', winDurStr);  

        %click selection of frequency from spectrogram
        if strcmp(txtStyle, 'edit')
            buttonSize = [0.055 0.025];
            freqminpos = get(pHandles.freqminH, 'position');
            freqmaxpos = get(pHandles.freqmaxH, 'position');
            pb1pos = [freqminpos(1)+freqminpos(3)-buttonSize(1) freqminpos(2) buttonSize];
            pb2pos = [freqmaxpos(1)+freqmaxpos(3)-buttonSize(1) freqmaxpos(2) buttonSize];
            uicontrol(panelHandle, 'style', 'pushbutton', 'units', 'normalized', 'position', pb1pos, 'cdata', mouseIcon, 'tooltipstring', 'select with mouse', 'callback', {@pickFrequency, pHandles.freqminH});
            uicontrol(panelHandle, 'style', 'pushbutton', 'units', 'normalized', 'position', pb2pos, 'cdata', mouseIcon, 'tooltipstring', 'select with mouse', 'callback', {@pickFrequency, pHandles.freqmaxH});
        end
        uicontrol(panelHandle, 'style', 'text', 'units', 'normalized', 'position', [column1 0.72 boxLength txtHeight], 'string', 'nHarmonics:','BackgroundColor', bcolor, 'horizontalalignment', 'left');
        pHandles.nharmH = uicontrol(panelHandle, 'style', txtStyle, 'units', 'normalized', 'position', [column1 0.6791 boxLength boxHeight], 'string', inputs.values(1).nHarmonics,'tooltipstring', nharmStr);
        
    %   These Parameters are listed in the PlotHarmonicStates function, but not implemented, therefore I have commmented them out 
    %   pHandles.plotPartH = uicontrol(panelHandle, 'style', 'checkbox', 'position', [column3 130 boxLength txtHeight], 'string', 'Particles','BackgroundColor', bcolor, 'value', inputs.bools.plotPart);
    %   pHandles.freqCompH = uicontrol(panelHandle, 'style', 'checkbox', 'position', [column3 90 boxLength txtHeight], 'string', 'Frequency Comparison','BackgroundColor', bcolor, 'value', inputs.bools.freqComp);
    else
        prevTabs = length(pHandles.mvarH);
    end
    
    warning('off', 'MATLAB:uitab:DeprecatedFunction');
    warning('off', 'MATLAB:uitab:MigratingFunction');
    warning('off', 'MATLAB:uitabgroup:DeprecatedFunction');
    warning('off', 'MATLAB:uitabgroup:MigratingFunction');
    
    if ntabs == prevTabs
%         disp('do nothing')
            return
    elseif prevTabs > 0 && ntabs > prevTabs
        %disp('appending tabs')
        inputs.values(1:prevTabs) = inputs.values;
        for a = (prevTabs +1): ntabs
            inputs.values(a) = inputs.values(1);
        end
        panelKids = get(panelHandle, 'children');
        hTabGroup = panelKids(strcmp('uitabgroup', get(panelKids, 'type')));
    elseif ntabs < prevTabs
        %disp('destroying tabs')
        panelKids = get(panelHandle, 'children');
        hTabGroup = panelKids(strcmp('uitabgroup', get(panelKids, 'type')));
        tabs = get(hTabGroup, 'children');
        tabs(end-(prevTabs-ntabs-1):end); %this is causeing error for some reason, trying to delete more than supposed to
        inds2delete = ntabs+1:prevTabs;
        delete(tabs(inds2delete));
        pHandles.mvarH(inds2delete) = [];
        return
    elseif prevTabs==0
       hTabGroup = uitabgroup('parent', panelHandle, 'units', 'normalized', 'Position', [0 0 1 0.55], 'BackgroundColor', bcolor); 
    else
        error('invalid tab creation input')
    end

%     tabPanel = uitabpanel('parent', panelHandle, 'units', 'normalized', 'Position', [0 0 1 0.55], 'title', titles, 'panelbackgroundcolor', bcolor,'framebackgroundcolor', bcolor,'CreateFcn', @createTabs);
    
    smallBoxLength = 0.25;
    boxHeight = 0.08;
    txtHeight = 0.06;
    column2a = 0.37;
    column3a = 0.69;
    for a=prevTabs+1:ntabs
        %inputs
        tab(a) = uitab(hTabGroup, 'title', ['Section' num2str(a)]);
        uicontrol(tab(a), 'style', 'text', 'units', 'normalized', 'position', [column1 0.84 boxLength boxHeight+0.05], 'string', 'Measurement Variance:','BackgroundColor', bcolor, 'horizontalalignment', 'left');
        pHandles.mvarH(a) = uicontrol(tab(a), 'style', txtStyle, 'units', 'normalized', 'position', [column1 0.76 boxLength boxHeight], 'string', inputs.values(a).measurementVariance,'tooltipstring', mvarStr);

        uicontrol(tab(a), 'style', 'text', 'units', 'normalized', 'position', [column1 0.59 boxLength boxHeight+0.05], 'string', 'Frequency Coefficient:','BackgroundColor', bcolor, 'horizontalalignment', 'left');
        pHandles.freqcoH(a) = uicontrol(tab(a), 'style', txtStyle, 'units', 'normalized', 'position', [column1 0.51 boxLength boxHeight], 'string', inputs.values(a).frequencyCoefficient,'tooltipstring', freqcoStr);
        uicontrol(tab(a), 'style', 'text', 'units', 'normalized', 'position', [column2 0.84 boxLength boxHeight+0.05], 'string', 'Amplitude Coefficient:','BackgroundColor', bcolor, 'horizontalalignment', 'left');
        pHandles.ampcoH(a) = uicontrol(tab(a), 'style', txtStyle, 'units', 'normalized', 'position', [column2 0.76 boxLength boxHeight], 'string', inputs.values(a).amplitudeCoefficient,'tooltipstring', ampcoStr);
        uicontrol(tab(a), 'style', 'text', 'units', 'normalized', 'position', [column2 0.59 boxLength boxHeight+0.05], 'string', 'Frequency Mean:','BackgroundColor', bcolor, 'horizontalalignment', 'left');
        pHandles.freqmeanH(a) = uicontrol(tab(a), 'style', txtStyle, 'units', 'normalized', 'position', [column2 0.51 boxLength boxHeight], 'string', inputs.values(a).frequencyMean/freqMultiplier,'tooltipstring', freqmeanStr);

        uicontrol(tab(a), 'style', 'text', 'units', 'normalized', 'position', [column1 0.42 boxLength txtHeight], 'string', 'Initial State','BackgroundColor', bcolor, 'horizontalalignment', 'left');
        uicontrol(tab(a), 'style', 'text', 'units', 'normalized', 'position', [column1 0.37 smallBoxLength txtHeight], 'string', 'phase:','BackgroundColor', bcolor, 'horizontalalignment', 'left');
        pHandles.initphaseH(a) = uicontrol(tab(a), 'style', txtStyle, 'units', 'normalized', 'position', [column1 0.29 smallBoxLength boxHeight], 'string', inputs.values(a).StateInitial.phase,'tooltipstring', initphaseStr);
        uicontrol(tab(a), 'style', 'text', 'units', 'normalized', 'position', [column2a 0.37 smallBoxLength txtHeight], 'string', 'frequency:','BackgroundColor', bcolor, 'horizontalalignment', 'left');
        pHandles.initfreqH(a) = uicontrol(tab(a), 'units', 'normalized', 'style', txtStyle, 'position', [column2a 0.29 smallBoxLength boxHeight], 'string', inputs.values(a).StateInitial.frequency/freqMultiplier,'tooltipstring', initfreqStr);
        uicontrol(tab(a), 'style', 'text', 'units', 'normalized', 'position', [column3a 0.37 smallBoxLength+0.02 txtHeight], 'string', 'amplitudes:','BackgroundColor', bcolor, 'horizontalalignment', 'left');
        pHandles.initampH(a) = uicontrol(tab(a), 'style', txtStyle, 'units', 'normalized', 'position', [column3a 0.29 smallBoxLength boxHeight], 'string', inputs.values(a).StateInitial.amplitudes(length(inputs.values(a).StateInitial.amplitudes)),'tooltipstring', initampStr);

        uicontrol(tab(a), 'style', 'text', 'units', 'normalized', 'position', [column1 0.20 boxLength txtHeight], 'string', 'State Variance','BackgroundColor', bcolor, 'horizontalalignment', 'left');
        uicontrol(tab(a), 'style', 'text', 'units', 'normalized', 'position', [column1 0.15 smallBoxLength txtHeight], 'string', 'phase:','BackgroundColor', bcolor, 'horizontalalignment', 'left');
        pHandles.varphaseH(a) = uicontrol(tab(a), 'style', txtStyle, 'units', 'normalized', 'position', [column1 0.05 smallBoxLength boxHeight], 'string', inputs.values(a).StateVariance.phase,'tooltipstring', varphaseStr);
        uicontrol(tab(a), 'style', 'text', 'units', 'normalized', 'position', [column2a 0.15 smallBoxLength txtHeight], 'string', 'frequency:','BackgroundColor', bcolor, 'horizontalalignment', 'left');
        pHandles.varfreqH(a) = uicontrol(tab(a), 'style', txtStyle, 'units', 'normalized', 'position', [column2a 0.05 smallBoxLength boxHeight], 'string', inputs.values(a).StateVariance.frequency/freqMultiplier,'tooltipstring', varfreqStr);
        uicontrol(tab(a), 'style', 'text', 'units', 'normalized', 'position', [column3a 0.15 smallBoxLength+0.02 txtHeight], 'string', 'amplitudes:','BackgroundColor', bcolor, 'horizontalalignment', 'left');
        pHandles.varampH(a) = uicontrol(tab(a), 'style', txtStyle, 'units', 'normalized', 'position', [column3a 0.05 smallBoxLength boxHeight], 'string', inputs.values(a).StateVariance.amplitudes(1),'tooltipstring', varampStr);

        if strcmp(txtStyle, 'edit')
            buttonSize = [0.055 0.045];
            freqmpos = get(pHandles.freqmeanH(a), 'position');
            freqipos = get(pHandles.initfreqH(a), 'position');
            pb1pos = [freqmpos(1)+freqmpos(3)-buttonSize(1) freqmpos(2) buttonSize];
            pb2pos = [freqipos(1)+freqipos(3)-buttonSize(1) freqipos(2) buttonSize];
            pb1 = uicontrol(tab(a), 'style', 'pushbutton', 'units', 'normalized', 'position', pb1pos, 'cdata', mouseIcon, 'tooltipstring', 'select with mouse', 'callback', {@pickFrequency, pHandles.freqmeanH(a)});
            pb2 = uicontrol(tab(a), 'style', 'pushbutton', 'units', 'normalized', 'position', pb2pos, 'cdata', mouseIcon, 'tooltipstring', 'select with mouse', 'callback', {@pickFrequency, pHandles.initfreqH(a)});
        end
    end

    function pickFrequency(hObj, eventData, inputBoxHandle)
    % select a value from the orignal spectrogram axis? and use this value
    % to set the corresponding edit box
    caxes = gca;
    setptr(gcf, 'crosshair');
    waitforbuttonpress;
    setptr(gcf, 'arrow');
    point = get(caxes, 'currentPoint');
    freq = round(point(1,2));
    ylims = get(caxes, 'ylim');
    if freq < ylims(1)
        freq = ylims(1);
    elseif freq > ylims(2)
        freq = round(ylims(2));
    end
    if (inputs.sampleRate/2) > 40e3
        freq = freq*1000;
    end
    set(inputBoxHandle, 'string', num2str(freq))    
    end
end