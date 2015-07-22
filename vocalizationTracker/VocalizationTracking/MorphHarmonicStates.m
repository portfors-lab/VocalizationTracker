function [MorphedHarmonicStates y] = MorphHarmonicStates(HarmonicStates, varargin)
%MorphHarmonicStates: Modified ("Morphs") the states contained in the supplied
%                     harmonic states
%
%   [MorphedHarmonicStates y] = MorphHarmonicStates(HarmonicStates, varargin)
%   
% Input Paramters
%   HarmonicStates         A structure of states conforming with a harmonic model structure
%
% Optional Parameters
%   'frequencyMultiplier'  The amount to multiply the frequency by.
%                          Default: 1
%   'lengthMultiplier'     The amount to multiply the length of the signal by.
%                          Default: 1
%   'harmonicSelection'    The harmonics to keep in the morphed states.
%                          Format example: [1 3 4], to keep the 1st, 3rd and 4th harmonic
%                          Default: []=All
%   'stateType'            The state type in the harmonic model to use.
%                          Default: 'Filtered'
%   'randomizePhase'       Specifies whether to add a random phase offset to each
%                          harmonic.
%                          Default: false
%
% Output
%   MorphedHarmonicStates  The morphed harmonic state structure
%   y                      The signal generated from the morphed harmonic model states
%
%   Example: To be written.
%
%   Reference: To be completed.

%==============================================================================
% Error Checking
%==============================================================================
if nargin<1
    help MorphHarmonicStates;
    return;
end

%==============================================================================
% Process Function Arguments
%==============================================================================
frequencyMultiplier = 1;
lengthMultiplier = 1;
harmonicSelection = [];
stateType = 'Filtered';
randomizePhase = false;
distortion = 'none';
removeAmplitudeModulation = false;
removeFrequencyModulation = false;
invertSignal = false;
moveHarm = false;

nMandatoryArguments = 1;
if nargin>nMandatoryArguments
    if ~isstruct(varargin{1})
        if rem(length(varargin),2)~=0, error('Optional input arguments must be in name-value pairs.'); end;
        Parameters = struct;
        for c1=1:2:length(varargin)-1
            if ~ischar(varargin{c1}), error(['Error parsing arguments: Expected property name string at argument ' num2str(c1+1)]); end        
            Parameters.(varargin{c1}) = varargin{c1+1};
        end
    else
        Parameters = varargin{1};
    end
    
    parameterNames = fieldnames(Parameters);
    for c1 = 1:length(parameterNames)
        parameterName  = parameterNames{c1};
        parameterValue = Parameters.(parameterName);
        switch lower(parameterName)
            case lower('frequencyMultiplier'), frequencyMultiplier = parameterValue;
            case lower('lengthMultiplier'),    lengthMultiplier    = parameterValue;
            case lower('harmonicSelection'),   harmonicSelection   = parameterValue;    
            case lower('stateType'),           stateType           = parameterValue;   
            case lower('randomizePhase'),      randomizePhase      = parameterValue;
            case lower('distortion'),          distortion          = parameterValue;
            case lower('removeAmplitudeModulation'), removeAmplitudeModulation  = parameterValue;
            case lower('removeFrequencyModulation'), removeFrequencyModulation  = parameterValue;
            case 'invertsignal', invertSignal = parameterValue;
            case 'moveharmonic'
                if ~isempty(parameterValue) && ~any(isnan(parameterValue))
                    moveHarm = true; 
                    moveInterval = parameterValue;
                end
            otherwise,                         error(['Unrecognized property: ''' varargin{c1} '''']);
        end
    end
end
                                                 
%==============================================================================
% Preprocessing
%==============================================================================
if ~isfield(HarmonicStates,stateType)
    error(['There is no state type named' stateType]);
end

morphedStates = getfield(HarmonicStates,stateType);
nHarmonics = size(morphedStates.amplitudesCos,2);

sampleRate = HarmonicStates.ModelParameters.sampleRate;
sampleInterval = 1/sampleRate;
nSamples = length(morphedStates.frequency);

amplitudesCos = morphedStates.amplitudesCos;   
amplitudesSin = morphedStates.amplitudesSin; 
phase = morphedStates.phase; 
frequency = morphedStates.frequency; 

if isfield(morphedStates,'harmonicMultipliers')
    harmonicMultipliers = morphedStates.harmonicMultipliers;
else
    harmonicMultipliers = (1:nHarmonics)';
end

%==============================================================================
% Process Remove Frequency Modulation
%==============================================================================
if (removeFrequencyModulation)
    frequency = ones(length(frequency),1) * mean(frequency);
end

%==============================================================================
% Process Remove Amplitude Modulation
%==============================================================================
if (removeAmplitudeModulation)
    c = amplitudesCos + j*amplitudesSin;
    cAmplitudeEstimate = abs(c);
    meanAmplitudes = mean(cAmplitudeEstimate);
    constantAmplitudes = ones(size(cAmplitudeEstimate,1),size(cAmplitudeEstimate,2)) .* repmat(meanAmplitudes,nSamples,1);
%     constantAmplitudes(:,1) = zeros(size(cAmplitudeEstimate,1),1);
%     cPhaseEstimate     = angle(c);
    %Fixing phase relationship results in less "noise" in synthesized signal
    relativePhaseOffsets = pi/2*[0:nHarmonics-1];
    cPhaseEstimate     = ones(size(cAmplitudeEstimate,1),size(cAmplitudeEstimate,2)) .* repmat(relativePhaseOffsets,nSamples,1);
    amplitudesCos = constantAmplitudes.*cos(cPhaseEstimate);
    amplitudesSin = constantAmplitudes.*sin(cPhaseEstimate);
end

%==============================================================================
% Process Frequency Shifting
%==============================================================================
if length(frequencyMultiplier) > 1
    if length(frequencyMuliplier ~= harmonicSelection)
        error('Frequency muliplier vector does not equal number of harmonics selected for processing');
    end
    harmonicMultipliers(harmonicSelection);
    harmonicMultipliers(harmonicSelection) = harmonicMultipliers(harmonicSelection)*frequencyMultiplier;
end
if isempty(harmonicSelection)
    harmonicMultipliers = harmonicMultipliers*frequencyMultiplier;
else
    harmonicMultipliers(harmonicSelection) = harmonicMultipliers(harmonicSelection)*frequencyMultiplier;
end

%==============================================================================
% Process Length Shifting
%==============================================================================
if lengthMultiplier ~= 1
    nominator = round(lengthMultiplier*10);
    denominator = 10;
    amplitudesCos = resample(amplitudesCos,nominator,denominator);
    amplitudesSin = resample(amplitudesSin,nominator,denominator);
    frequency     = resample(frequency,nominator,denominator);
    %phase         = resample(phase,nominator,denominator);
    nSamples      = length(frequency);
end 
    
%==============================================================================
% Process Harmonic Selection
% If harmonicSelection is chosen along with a frequency multiplier or vector
% of multipliers, then apply the frequency shift to the selected harmonics
% and do not delete the others.
%==============================================================================
if ~isempty(harmonicSelection) && frequencyMultiplier == 1
    harmonicSelection = unique(harmonicSelection);
    if max(harmonicSelection) > nHarmonics || min(harmonicSelection) < 1
        error('Selected Harmonics Are Out Of Range');
    end
    unselectedHarmonics = setdiff(1:nHarmonics,harmonicSelection);
    amplitudesCos(:,unselectedHarmonics) = 0;
    amplitudesSin(:,unselectedHarmonics) = 0;
end

%==============================================================================
% Process Phase Randomization
%==============================================================================
if randomizePhase
    c = amplitudesCos + j*amplitudesSin;
    cAmplitudeEstimate = abs(c);
    cPhaseEstimate     = angle(c);
    randomPhases = cPhaseEstimate + repmat(mod(2*pi*rand(1,nHarmonics),2*pi),nSamples,1);
    amplitudesCos = cAmplitudeEstimate.*cos(randomPhases);
    amplitudesSin = cAmplitudeEstimate.*sin(randomPhases);
end

if invertSignal
    amplitudesCos = flipud(amplitudesCos);
    amplitudesSin = flipud(amplitudesSin);
    frequency = flipud(frequency);
    phase = flipud(phase);
end

if moveHarm
    %Under construction
    intStart = moveInterval(1);
    intStop = moveInterval(2);
    shift = moveInterval(3);
    newStop = intStart+shift+(intStop-intStart);
    if newStop < intStart
        clearStart = intStart;
        clearDur = intStop-intStart+1;
    else
        clearStart = newStop+1;
        clearDur = intStop-clearStart+1;
    end
    amplitudesCos(intStart+shift:newStop,2) = amplitudesCos(intStart:intStop,2);
    amplitudesCos(clearStart:intStop,2) = zeros(clearDur,1);
    amplitudesSin(intStart+shift:newStop,2) = amplitudesSin(intStart:intStop,2);
    amplitudesSin(clearStart:intStop,2) = zeros(clearDur,1);
%     harmonicMultipliers(2) = 1.5;
end

switch lower(distortion)
    %FIXME. Introduces a lot of phase related distortion
    case {'difference'}
        if nHarmonics ~= 2
            error('Difference distortion only applicable to signal models with 2 harmonics');
        end
        c = amplitudesCos + j*amplitudesSin;
        cAmplitudeEstimate = abs(c);
        cPhaseEstimate     = angle(c);
%         newAmplitudes      = diff(cAmplitudeEstimate,1,2);
        newAmplitudes      = min(cAmplitudeEstimate,[],2);
%         newPhases          = diff(cPhaseEstimate,1,2);
        newPhases          = cPhaseEstimate(:,1);
        newCos             = newAmplitudes.*cos(newPhases);
        newSin             = newAmplitudes.*sin(newPhases);
        
        amplitudesSin(:,3) = newSin;
        amplitudesCos(:,3) = newCos;
        harmonicMultipliers(3) = 1.5;
        nHarmonics = 3;
end

%==============================================================================
% Use the morphed model parameters to generate new signal.
%==============================================================================
tmpPhase = zeros(1,nHarmonics)';
y = zeros(1,nSamples);
phase2 = [];
% for c1=1:nSamples
%     tmpAmplitudesCos = amplitudesCos(c1,:);
%     tmpAmplitudesSin = amplitudesSin(c1,:);
%     tmpPhase      = tmpPhase + 2*pi*sampleInterval*frequency(c1)*harmonicMultipliers;
%     y(c1)         = tmpAmplitudesCos*cos(tmpPhase) + tmpAmplitudesSin*sin(tmpPhase);   
%     phase(c1) = tmpPhase(1);
% end
%adapt above method to use the mod of the phase
for c1=1:nSamples
    tmpAmplitudesCos = amplitudesCos(c1,:);
    tmpAmplitudesSin = amplitudesSin(c1,:);
    tmpPhase      = mod(tmpPhase + 2*pi*sampleInterval*frequency(c1)*harmonicMultipliers,2*pi);
%     tmpPhase      = mod(tmpPhase,2*pi);
    y(c1)         = tmpAmplitudesCos*cos(tmpPhase) + tmpAmplitudesSin*sin(tmpPhase);   
    phase(c1) = tmpPhase(1);
end

%--------------------------------------------------------------------------
%another way, from plotHarmonic states...... uses phase caluclated in
%TrackHarmonicsEKF, which takes the mod(prevPhase, 2pi) - of the previous
%phase and 2*pi, to get the instantaneous phase.
% nHarmonics = size(amplitudesCos,2);
% nSamples = length(frequency);
% % Generate synthesized signal from harmonic states
% y = zeros(1,nSamples);
% k = (1:nHarmonics).';
% for c1=1:nSamples
%     tmpAmplitudesCos = amplitudesCos(c1,:);
%     tmpAmplitudesSin = amplitudesSin(c1,:);
%     tmpPhase         = phase(c1);  
%     y(c1)        = tmpAmplitudesCos*cos(k*tmpPhase) + tmpAmplitudesSin*sin(k*tmpPhase);
% end

%--------------------------------------------------------------------------

morphedStates.amplitudesCos = amplitudesCos;   
morphedStates.amplitudesSin = amplitudesSin; 
morphedStates.phase         = phase; 
morphedStates.frequency     = frequency; 
morphedStates.y             = y;
morphedStates.harmonicMultipliers = harmonicMultipliers;

HarmonicStates = setfield(HarmonicStates,stateType, morphedStates);

MorphedHarmonicStates = HarmonicStates;

