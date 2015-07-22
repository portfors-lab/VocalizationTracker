function [y yPink States] = GenerateModulatedHarmonicModelStates(sampleRate, duration, ModelParameters)
%ModulatedHarmonicModelStates: Generate the harmonic model states consistent with
%                              the specified modulated harmonic model parameters.
%                              Similar to HarmonicModelStates, but includes a
%                              frequency modulation.
%
%   [y,State] = ModulatedHarmonicModelStates(sampleRate,duration,ModelParameters)
%   
% Input Paramters
%   sampleRate         Sample rate (Hz).
%   duration           The duration of the harmonic signal (s).
%   ModelParameters    Structure containing all of the parameters that 
%                      describe the statistical model.
%
% Output
%   y                 The signal generated from the harmonic model states
%   State             Structure containing all of the parameters
%                     of the harmonic model
%
%   Example: To be written.
%
%   Reference: To be completed.

%====================================================================
% Preprocessing
%====================================================================
nSamples = round(duration*sampleRate);
nHarmonics = ModelParameters.nHarmonics;
sampleInterval = 1/sampleRate;

% phases        = zeros(nSamples,1);
% frequencies   = zeros(nSamples,1);
% amplitudesCos = zeros(nSamples,nHarmonics);
% amplitudesSin = zeros(nSamples,nHarmonics);

%====================================================================
% Main Loop
%====================================================================
y = zeros(nSamples,1);
States.True = struct(...
    'phase'        ,nan(nSamples,1),...
    'frequency'    ,nan(nSamples,1),...
    'amplitudesCos',nan(nSamples,nHarmonics),...
    'amplitudesSin',nan(nSamples,nHarmonics),...
    'measurementNoise',nan(nSamples,1),...
    'modPhase'    ,nan(nSamples,1),...
    'modFrequency',nan(nSamples,1),...
    'modAmplitudeCos',nan(nSamples,1),...
    'modAmplitudeSin',nan(nSamples,1)...
    );

States.label = 'True Modulated Harmonic Model States';

frequencyMean         = ModelParameters.frequencyMean;
frequencySquashStart  = ModelParameters.frequencySquashStart;
frequencySquashAmount = ModelParameters.frequencySquashAmount;
frequencyAlpha        = ModelParameters.frequencyCoefficient;
amplitudeAlpha        = ModelParameters.amplitudeCoefficient;
phase                 = ModelParameters.StateInitial.phase;
frequency             = ModelParameters.StateInitial.frequency;
amplitudesCos         = ModelParameters.StateInitial.amplitudes;
amplitudesSin         = ModelParameters.StateInitial.amplitudes;
% Modulation state variables
modFrequencyMean      = ModelParameters.modFrequencyMean;
modFrequencyAlpha     = ModelParameters.modFrequencyCoefficient;
modAmplitudeMean      = ModelParameters.modAmplitudeMean;
modAmplitudeAlpha     = ModelParameters.modAmplitudeCoefficient;
modPhase              = ModelParameters.StateInitial.modPhase;
modFrequency          = ModelParameters.StateInitial.modFrequency;
modAmplitudeCos       = ModelParameters.StateInitial.modAmplitudeCos;
modAmplitudeSin       = ModelParameters.StateInitial.modAmplitudeSin;

for c1=1:nSamples
    k = (1:nHarmonics).';
% Use these to have the amplitudes vary independently
%     amplitudesCos = amplitudeAlpha*amplitudesCos ...
%         + sqrt(sampleInterval*ModelParameters.StateVariance.amplitudes)*randn(nHarmonics,1);
%     amplitudesSin = amplitudeAlpha*amplitudesSin ...
%         + sqrt(sampleInterval*ModelParameters.StateVariance.amplitudes)*randn(nHarmonics,1);

% Use these to have the amplitudes co-vary
    amplitudesCos = amplitudeAlpha*amplitudesCos ...
        + sqrt(sampleInterval*ModelParameters.StateVariance.amplitudes)'.*(randn(1,1)*ones(nHarmonics,1));
    amplitudesSin = amplitudeAlpha*amplitudesSin ...
        + sqrt(sampleInterval*ModelParameters.StateVariance.amplitudes)'.*(randn(1,1)*ones(nHarmonics,1));
    frequency = frequencyMean + frequencyAlpha*(frequency-frequencyMean) ...
        + sqrt(sampleInterval*ModelParameters.StateVariance.frequency)*randn;
    modPhase = modPhase + sampleInterval*modFrequency ...
        + sqrt(sampleInterval*ModelParameters.StateVariance.modPhase)*randn;
    
    modFrequency = modFrequencyMean + modFrequencyAlpha*(modFrequency-modFrequencyMean) ...
        + sqrt(sampleInterval*ModelParameters.StateVariance.modFrequency)*randn;
    modAmplitudeCos = modAmplitudeMean + modAmplitudeAlpha * (modAmplitudeCos-modAmplitudeMean) ...
        + sqrt(sampleInterval*ModelParameters.StateVariance.modAmplitudes)*randn;
    modAmplitudeSin = modAmplitudeMean + modAmplitudeAlpha * (modAmplitudeSin-modAmplitudeMean) ...
        + sqrt(sampleInterval*ModelParameters.StateVariance.modAmplitudes)*randn;
    
    modulatedFrequency = frequency ...
        + modAmplitudeCos*cos(modPhase) ...
        + modAmplitudeSin*sin(modPhase);
    
    if modulatedFrequency > frequencySquashStart
        squashRange = modulatedFrequency - frequencySquashStart;
        modulatedFrequency = frequencySquashStart + frequencySquashAmount*tanh(squashRange/frequencySquashAmount);
    end
    
    phase = phase + 2*pi*sampleInterval*(modulatedFrequency) ...
        + sqrt(sampleInterval*ModelParameters.StateVariance.phase)*randn;
    phase = mod(phase, 2*pi);
    
    States.True.phase(c1,:) = phase;
    States.True.frequency(c1,:) = modulatedFrequency;
    States.True.amplitudesCos(c1,:) = amplitudesCos.';
    States.True.amplitudesSin(c1,:) = amplitudesSin.';
    States.True.modPhase(c1,:) = modPhase;
    States.True.modFrequency(c1,:) = modFrequency;
    States.True.modAmplitudeCos(c1,:) = modAmplitudeCos;
    States.True.modAmplitudeSin(c1,:) = modAmplitudeSin;
end

y = zeros(1,nSamples);
k = (1:nHarmonics).';
for c1=1:nSamples
    amplitudesCos = States.True.amplitudesCos(c1,:);
    amplitudesSin = States.True.amplitudesSin(c1,:);
    phase         = States.True.phase(c1);
    y(c1)        = amplitudesCos*cos(k*phase) + amplitudesSin*sin(k*phase);
end



pinkNoiseTransient = 1430; %See Pinkify()
%Generate white noise, with extra samples for the pink noise generation
measurementNoise = sqrt(ModelParameters.measurementVariance)*randn(nSamples+pinkNoiseTransient,1);
[pinkNoise transient] = Pinkify(measurementNoise);
assert(transient == pinkNoiseTransient);
measurementNoise = measurementNoise(pinkNoiseTransient+1:end);
pinkNoise = pinkNoise(pinkNoiseTransient+1:end);

if isfield(ModelParameters,'SNR') && ~isnan(ModelParameters.SNR)  

    
    signalRMS = sqrt(mean(y.^2));
    noiseRMS = sqrt(mean(measurementNoise.^2));
    targetRatio = exp(ModelParameters.SNR*log(10)/20);
    multiplier = signalRMS/(noiseRMS*targetRatio);
    measurementNoise = multiplier*measurementNoise;
    
    pinkNoiseRMS = sqrt(mean(pinkNoise.^2));
    multiplier = signalRMS/(pinkNoiseRMS*targetRatio);
    pinkNoise = multiplier*pinkNoise;
    
    % For Debuggins
%     signalPower = mean(y.^2);
%     noisePower = mean(measurementNoise.^2);
%     pinkNoisePower = mean(pinkNoise.^2);
%     targetSNR = ModelParameters.SNR
%     actualSNR = 10*log10(signalPower./noisePower)
%     actualPinkSNR = 10*log10(signalPower./pinkNoisePower)
end


%Add the noise into the generated signal
yPink = y + pinkNoise';
y = y + measurementNoise';

States.True.measurementNoise = measurementNoise;
States.True.pinkMeasurementNoise = pinkNoise;


