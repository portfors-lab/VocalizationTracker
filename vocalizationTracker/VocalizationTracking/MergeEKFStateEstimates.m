function EKFStateEstimates = MergeEKFStateEstimates(segmented_EKFStateEstimates)

nSegments = length(segmented_EKFStateEstimates);

EKFStateEstimates = segmented_EKFStateEstimates{1};

risefallTime = 0.0005;

EKFStateEstimates.Smoothed.amplitudesCos = AddRiseFall(EKFStateEstimates.Smoothed.amplitudesCos, EKFStateEstimates.ModelParameters.sampleRate, 'riseTime', risefallTime, 'fallTime', risefallTime);
EKFStateEstimates.Smoothed.amplitudesSin = AddRiseFall(EKFStateEstimates.Smoothed.amplitudesSin, EKFStateEstimates.ModelParameters.sampleRate, 'riseTime', risefallTime, 'fallTime', risefallTime);

EKFStateEstimates.Filtered.amplitudesCos = AddRiseFall(EKFStateEstimates.Filtered.amplitudesCos, EKFStateEstimates.ModelParameters.sampleRate, 'riseTime', risefallTime, 'fallTime', risefallTime);
EKFStateEstimates.Filtered.amplitudesSin = AddRiseFall(EKFStateEstimates.Filtered.amplitudesSin, EKFStateEstimates.ModelParameters.sampleRate, 'riseTime', risefallTime, 'fallTime', risefallTime);

EKFStateEstimates.Predicted.amplitudesCos = AddRiseFall(EKFStateEstimates.Predicted.amplitudesCos, EKFStateEstimates.ModelParameters.sampleRate, 'riseTime', risefallTime, 'fallTime', risefallTime);
EKFStateEstimates.Predicted.amplitudesSin = AddRiseFall(EKFStateEstimates.Predicted.amplitudesSin, EKFStateEstimates.ModelParameters.sampleRate, 'riseTime', risefallTime, 'fallTime', risefallTime);

for iSegment = 2:nSegments
    nextEstimates = segmented_EKFStateEstimates{iSegment};
    
    nextFilteredAmplitudesCos = nextEstimates.Filtered.amplitudesCos;
    nextFilteredAmplitudesCos = AddRiseFall(nextFilteredAmplitudesCos, EKFStateEstimates.ModelParameters.sampleRate, 'riseTime', risefallTime, 'fallTime', risefallTime);
    nextFilteredAmplitudesSin = nextEstimates.Filtered.amplitudesSin;
    nextFilteredAmplitudesSin = AddRiseFall(nextFilteredAmplitudesSin, EKFStateEstimates.ModelParameters.sampleRate, 'riseTime', risefallTime, 'fallTime', risefallTime);
    
    EKFStateEstimates.Filtered.y = [EKFStateEstimates.Filtered.y ; nextEstimates.Filtered.y];
    EKFStateEstimates.Filtered.error = [EKFStateEstimates.Filtered.error ; nextEstimates.Filtered.error];
    EKFStateEstimates.Filtered.phase = [EKFStateEstimates.Filtered.phase ; nextEstimates.Filtered.phase];
    EKFStateEstimates.Filtered.frequency = [EKFStateEstimates.Filtered.frequency ; nextEstimates.Filtered.frequency];
    EKFStateEstimates.Filtered.amplitudesCos = [EKFStateEstimates.Filtered.amplitudesCos ; nextFilteredAmplitudesCos];
    EKFStateEstimates.Filtered.amplitudesSin = [EKFStateEstimates.Filtered.amplitudesSin ; nextFilteredAmplitudesSin];
    
    nextPredictedAmplitudesCos = nextEstimates.Predicted.amplitudesCos;
    nextPredictedAmplitudesCos = AddRiseFall(nextPredictedAmplitudesCos, EKFStateEstimates.ModelParameters.sampleRate, 'riseTime', risefallTime, 'fallTime', risefallTime);
    nextPredictedAmplitudesSin = nextEstimates.Predicted.amplitudesSin;
    nextPredictedAmplitudesSin = AddRiseFall(nextPredictedAmplitudesSin, EKFStateEstimates.ModelParameters.sampleRate, 'riseTime', risefallTime, 'fallTime', risefallTime);
    
    EKFStateEstimates.Predicted.y = [EKFStateEstimates.Predicted.y ; nextEstimates.Predicted.y];
    EKFStateEstimates.Predicted.error = [EKFStateEstimates.Predicted.error ; nextEstimates.Predicted.error];
    EKFStateEstimates.Predicted.phase = [EKFStateEstimates.Predicted.phase ; nextEstimates.Predicted.phase];
    EKFStateEstimates.Predicted.frequency = [EKFStateEstimates.Predicted.frequency ; nextEstimates.Predicted.frequency];
    EKFStateEstimates.Predicted.amplitudesCos = [EKFStateEstimates.Predicted.amplitudesCos ; nextPredictedAmplitudesCos];
    EKFStateEstimates.Predicted.amplitudesSin = [EKFStateEstimates.Predicted.amplitudesSin ; nextPredictedAmplitudesSin];
    
    nextSmoothedAmplitudesCos = nextEstimates.Smoothed.amplitudesCos;
    nextSmoothedAmplitudesCos = AddRiseFall(nextSmoothedAmplitudesCos, EKFStateEstimates.ModelParameters.sampleRate, 'riseTime', risefallTime, 'fallTime', risefallTime);
    nextSmoothedAmplitudesSin = nextEstimates.Smoothed.amplitudesSin;
    nextSmoothedAmplitudesSin = AddRiseFall(nextSmoothedAmplitudesSin, EKFStateEstimates.ModelParameters.sampleRate, 'riseTime', risefallTime, 'fallTime', risefallTime);
    
    EKFStateEstimates.Smoothed.y = [EKFStateEstimates.Smoothed.y ; nextEstimates.Smoothed.y];
    EKFStateEstimates.Smoothed.error = [EKFStateEstimates.Smoothed.error ; nextEstimates.Smoothed.error];
    EKFStateEstimates.Smoothed.phase = [EKFStateEstimates.Smoothed.phase ; nextEstimates.Smoothed.phase];
    EKFStateEstimates.Smoothed.frequency = [EKFStateEstimates.Smoothed.frequency ; nextEstimates.Smoothed.frequency];
    EKFStateEstimates.Smoothed.amplitudesCos = [EKFStateEstimates.Smoothed.amplitudesCos ; nextSmoothedAmplitudesCos];
    EKFStateEstimates.Smoothed.amplitudesSin = [EKFStateEstimates.Smoothed.amplitudesSin ; nextSmoothedAmplitudesSin];
    EKFStateEstimates.Smoothed.jacobianFunction = [EKFStateEstimates.Smoothed.jacobianFunction ; nextEstimates.Smoothed.jacobianFunction];
    EKFStateEstimates.Smoothed.ErrorCovariance = [EKFStateEstimates.Smoothed.ErrorCovariance ; nextEstimates.Smoothed.ErrorCovariance];
end