function GenerateSynthFigures(signal, sampleRate, sectPoints, inputs)

    multipleHarms = false;

    waitbar(0,'creating plots...');
    ntabs = length(inputs.values);
    
    if multipleHarms
        for a = 1:ntabs
            nHarmArray(a) = inputs.values(a).nHarmonics;
        end
        maxHarm = max(nHarmArray);
    else
        maxHarm = inputs.values(1).nHarmonics; %all the same anyways
    end
    
    %section of the signal for the current parameters
    if ntabs > 1
        for a = 1:ntabs
            if a == ntabs
                b=0;    %final section, include last point
            else
                b=1;    %minus 1 from section so no dulpicate points
            end                     
            trackModelParameters = inputs.values(a);
            sectSignal = signal(sectPoints(a):sectPoints(a+1)-b);
            segmented_EKFStateEstimates{a} = TrackHarmonicsEKF(sectSignal,sampleRate,trackModelParameters, maxHarm,'smoother',true);
        end
        EKFStateEstimates = MergeEKFStateEstimates(segmented_EKFStateEstimates);
    else
        trackModelParameters = inputs.values;
        EKFStateEstimates = TrackHarmonicsEKF(signal,sampleRate,trackModelParameters, trackModelParameters.nHarmonics,'smoother',true);
    end
    
    duration = length(signal)/sampleRate;

    if trackModelParameters.frequencyMean*trackModelParameters.nHarmonics>sampleRate/2, warning('Sampling theorem violated.'); end;
  
    windowDuration = 70/trackModelParameters.frequencyMean;

    bools = inputs.bools;
    fn= fieldnames(bools);
    waitbar(1/(length(fn)+2));
    
%     EKFStateEstimates.Filtered.frequency = abs(EKFStateEstimates.Filtered.frequency)
%         EKFStateEstimates.Predicted.frequency = abs(EKFStateEstimates.Predicted.frequency)
%     EKFStateEstimates.Smoothed.frequency = abs(EKFStateEstimates.Smoothed.frequency)

    try
        if multipleHarms
        EKFStateEstimates.nHarmSections = nHarmArray;    
        PlotHarmonicStates(EKFStateEstimates,...
                           'spectrogramType',2, ...
                           'dbRange', trackModelParameters.dbRange, ...
                           'stateType', inputs.stateType,...
                           'trackedSignal',signal,...
                           'colormap',inputs.colormap ,...
                           'windowDuration',inputs.windowDuration,...
                           'frequencyRange',trackModelParameters.frequencyRange,...
                           'plotSyntheticSpectrogram',bools.synSpec,...
                           'plotSyntheticSpectrogramHarmonics',bools.synSpecHarm,...
                           'plotSignalCompare',bools.signalComp,...
                           'plotAmplitudesAndPhases',bools.ampPhase,...
                           'plotSignalSpectrogram',bools.sigSpec,...
                           'plotSignalSpectrogramHarmonics',bools.sigSpecHarm,...
                           'plotResiduals',bools.plotRes,...
                           'plotMultiSections', sectPoints);
        else
        [plotHandlesStruct plotHandlesVect synthSignal] = PlotHarmonicStates(EKFStateEstimates,...
                           'spectrogramType',2, ...
                           'dbRange', trackModelParameters.dbRange, ...
                           'stateType', inputs.stateType,...
                           'trackedSignal',signal,...
                           'colormap', inputs.colormap,...
                           'windowDuration', inputs.windowDuration,...
                           'frequencyRange',trackModelParameters.frequencyRange,...
                           'plotSyntheticSpectrogram',bools.synSpec,...
                           'plotSyntheticSpectrogramHarmonics',bools.synSpecHarm,...
                           'plotSignalCompare',bools.signalComp,...
                           'plotAmplitudesAndPhases',bools.ampPhase,...
                           'plotSignalSpectrogram',bools.sigSpec,...
                           'plotSignalSpectrogramHarmonics',bools.sigSpecHarm,...
                           'plotResiduals',bools.plotRes);       
        end
%         plotFigure = figure('position', [0 0 700 500]);
%         hTabGroup = uitabgroup('v0');
%         for plotHandle = plotHandlesVect
%             name = get(plotHandle, 'name');
%             thisTab = uitab('v0',hTabGroup, 'title', name);
%             kids =get(plotHandle, 'children');
%             for kid = kids
%                 set(kid, 'parent', thisTab)
%             end
%             close(plotHandle);
%         end
    ManipulateSynthFigures(plotHandlesVect, ntabs, inputs, EKFStateEstimates)
    catch e %version 10
        errordlg({e.message;' '; 'Plotting error, try checking input values'}) %version 10
        %errordlg('Plotting error, try checking input values') %version 7
        rethrow(e)
    end
end