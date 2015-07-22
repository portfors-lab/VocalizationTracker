function mergedVocalization = MergeVocalizations(vocalizationCellArray, sampleRate)

defaultRiseFallTime = 0.0002; % 0.2 milliseconds

nVocalizations = length(vocalizationCellArray);

mergedVocalization = [];

for iVocalization = 1:nVocalizations
    vocalizationComponent = vocalizationCellArray{iVocalization};
    % Remove the mean
    vocalizationComponent = vocalizationComponent - mean(vocalizationComponent);
    if iVocalization == 1
        vocalizationComponent = AddRiseFall(vocalizationComponent, sampleRate, 'fallTime',defaultRiseFallTime);
    elseif iVocalization == nVocalizations
        vocalizationComponent = AddRiseFall(vocalizationComponent, sampleRate, 'riseTime',defaultRiseFallTime);
    else
        vocalizationComponent = AddRiseFall(vocalizationComponent, sampleRate, 'riseFallTime',defaultRiseFallTime);
    end
    mergedVocalization = [mergedVocalization vocalizationComponent'];
end