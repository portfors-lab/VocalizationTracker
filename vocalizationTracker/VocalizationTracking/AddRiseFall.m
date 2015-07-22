function signal = AddRiseFall(y, sampleRate,varargin)

%==============================================================================
% Error Checking
%==============================================================================
if nargin<2
    help AddRiseFall;
    return;
end

if (size(y,1) < size(y,2))
    y = y';
end

%==============================================================================
% Process Function Arguments
%==============================================================================
riseTime = 0.002; %s
fallTime = 0.002; %s

nMandatoryArguments = 2;
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
            case lower('riseTime'),     riseTime = parameterValue;    
            case lower('fallTime'),     fallTime = parameterValue;   
            case lower('riseFallTime'),     riseTime = parameterValue; fallTime = parameterValue;  
            otherwise,               error(['Unrecognized property: ''' varargin{c1} '''']);
        end
    end
end

riseSamples = round(sampleRate * riseTime);
if riseSamples > length(y)
    error('Rise time specified is longer than the duration of the signal.');
end
fallSamples = round(sampleRate * fallTime);
if fallSamples > length(y)
    error('Fall time specified is longer than the duration of the signal.');
end


rise_envelope = (1-cos(linspace(0,pi,riseSamples)))/2;
y(1:riseSamples,:) = y(1:riseSamples,:) .* repmat(rise_envelope',1,size(y,2));
fall_envelope = fliplr((1-cos(linspace(0,pi,fallSamples)))/2);
y(end-fallSamples+1:end,:) = y(end-fallSamples+1:end,:) .* repmat(fall_envelope',1,size(y,2));

signal = y;