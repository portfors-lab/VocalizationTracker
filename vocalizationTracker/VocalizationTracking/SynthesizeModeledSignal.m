function y = SynthesizeModeledSignal(HarmonicStates,varargin)
%SynthesizeModeledSignal: Synthesizes a signal using the provided harmonic states
%
%   [y sampleRate] = GenerateModeledSignal(HarmonicStates,varargin)
%   
% Input Paramters
%   HarmonicStates   A structure of states conforming with a harmonic model
%   structure
%
% Optional Parameters
%   'stateType'      The state type in the harmonic model to use.
%                    Default: 'Filtered'
%
% Output
%   y                The synthesized signal
%   sampleRate       The sample rate of the synthesized signal
%
%   Example: To be written.
%
%   Reference: To be completed.
%
%   See also NonparametricSpectrogram.

%   Tags: Filter, Estimator

%==============================================================================
% Error Checking
%==============================================================================
if nargin<1
    help States;
    return;
end

%==============================================================================
% Process Function Arguments
%==============================================================================
stateType = 'Filtered';

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
            case lower('stateType'),           stateType           = parameterValue;   
            otherwise,                         error(['Unrecognized property: ''' varargin{c1} '''']);
        end
    end
end
                                                 
%==============================================================================
% Preprocessing
%==============================================================================
if ~isfield(HarmonicStates,'True')
    error(['There is no state type named ' 'stateType']);
end

generateStates = getfield(HarmonicStates,'True');
nHarmonics = size(generateStates.amplitudesCos,2);

nSamples = length(generateStates.frequency);
y = zeros(1,nSamples);
k = (1:nHarmonics).';
for c1=1:nSamples
    amplitudesCos = generateStates.amplitudesCos(c1,:);
    amplitudesSin = generateStates.amplitudesSin(c1,:);
    phase         = generateStates.phase(c1);    
    y(c1)         = amplitudesCos*cos(k*phase) + amplitudesSin*sin(k*phase);    
end