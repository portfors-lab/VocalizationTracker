function [y,States] = GenerateHarmonicModelStates(sampleRate,duration,ModelParameters)
%HarmonicModelStates: Generate the harmonic model states consistent with
%                     the specified harmonic model parameters.
%
%   [y,State] = HarmonicModelStates(sampleRate,duration,ModelParameters)
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
nSamples       = round(duration*sampleRate);
nHarmonics     = ModelParameters.nHarmonics;
sampleInterval = 1/sampleRate;
intermittency  = false;

%====================================================================
% Main Loop
%====================================================================
y = zeros(nSamples,1);
States.True = struct(...
    'phase'        ,nan(nSamples,1),...
    'frequency'    ,nan(nSamples,1),...
    'amplitudesCos',nan(nSamples,nHarmonics),...
    'amplitudesSin',nan(nSamples,nHarmonics),...
    'measurementNoise',sqrt(ModelParameters.measurementVariance)*randn(nSamples,1)...
    );
States.label = 'True Harmonic Model States';

frequencyMean  = ModelParameters.frequencyMean;
frequencyAlpha = ModelParameters.frequencyCoefficient;
% amplitudeAlpha = ModelParameters.amplitudeCoefficient;
amplitudeAlpha = 1;
phase          = ModelParameters.StateInitial.phase;
frequency      = ModelParameters.StateInitial.frequency;
amplitudesCos  = ModelParameters.StateInitial.amplitudes;
amplitudesSin  = ModelParameters.StateInitial.amplitudes;

k = (1:nHarmonics).';
for c1=1:nSamples 
    if intermittency && c1>nSamples*1/3 && c1<nSamples*2/3
        amplitudesCos = zeros(nHarmonics,1);
        amplitudesSin = zeros(nHarmonics,1);
        frequency     = 0;
    end
    amplitudesCos = amplitudeAlpha*amplitudesCos                              + sqrt(sampleInterval*ModelParameters.StateVariance.amplitudes).*randn(nHarmonics,1);
    amplitudesSin = amplitudeAlpha*amplitudesSin                              + sqrt(sampleInterval*ModelParameters.StateVariance.amplitudes).*randn(nHarmonics,1);
    frequency     = frequencyMean + frequencyAlpha*(frequency-frequencyMean)  + sqrt(sampleInterval*ModelParameters.StateVariance.frequency)*randn;     
    phase         = phase + 2*pi*sampleInterval*frequency                     + sqrt(sampleInterval*ModelParameters.StateVariance.phase)*randn;  
    y(c1)         = amplitudesCos'*cos(k*phase) + amplitudesSin'*sin(k*phase) + States.True.measurementNoise(c1,1);  
    
    States.True.phase        (c1,:) = phase;
    States.True.frequency    (c1,:) = frequency;
    States.True.amplitudesCos(c1,:) = amplitudesCos.';
    States.True.amplitudesSin(c1,:) = amplitudesSin.';
    if intermittency && c1==nSamples*2/3-1
            frequency      = ModelParameters.StateInitial.frequency;
            amplitudesCos  = ModelParameters.StateInitial.amplitudes;
            amplitudesSin  = ModelParameters.StateInitial.amplitudes;
    end
end

ModelParameters.sampleRate = sampleRate;
States.ModelParameters     = ModelParameters;