function [inputs]= RetrieveValues(inputHandles, checkHandles)

multipleHarms = false;
% freqScaler = 1000; %frequency is in kHz (bad idea to change this)
freqScaler =1;
try
            dbRange = eval(get(inputHandles.dbRangeH, 'string'));               
            freqMin = eval(get(inputHandles.freqminH, 'string'))*freqScaler;
            freqMax = eval(get(inputHandles.freqmaxH, 'string'))*freqScaler;
            windowDuration = eval(get(inputHandles.winDurH, 'string'));

            
            inputs = nan;        
            ntabs = length(inputHandles.mvarH);

            %this is necessary in order to have sections with differenct
            %nHarmonics, maxHarm must be used inorder to keep matrix
            %dimensions the same so that they may be merged.
            if multipleHarms
                for a = 1:ntabs
                    nHarmArray(a) = eval(get(inputHandles.nharmH(a), 'string'));
                end
                maxHarm = max(nHarmArray);
            else
                nHarmonics = eval(get(inputHandles.nharmH, 'string'));
                maxHarm = nHarmonics;
            end
            
            for a=1:ntabs
                if multipleHarms
                    nHarmonics = eval(get(inputHandles.nharmH(a), 'string'));
                end
                frequencyMean = eval(get(inputHandles.freqmeanH(a), 'string'))*freqScaler; %Hz
                frequencyCoefficient = eval(get(inputHandles.freqcoH(a), 'string'));
                amplitudeCoefficient = eval(get(inputHandles.ampcoH(a), 'string'));
                measurementVariance = eval(get(inputHandles.mvarH(a), 'string'));
                initialPhase = eval(get(inputHandles.initphaseH(a), 'string'));
                initialFrequency = eval(get(inputHandles.initfreqH(a), 'string'))*freqScaler;
                initialAmplitude = eval(get(inputHandles.initampH(a), 'string'));
                variancePhase = eval(get(inputHandles.varphaseH(a), 'string'));
                frequencyVariance = eval(get(inputHandles.varfreqH(a), 'string'))*freqScaler;
                amplitudeVariance = eval(get(inputHandles.varampH(a), 'string'));

                amplitudeVariance = amplitudeVariance*ones(maxHarm,1);

                TrackModelParameters = struct(...
                    'nHarmonics'          ,nHarmonics        ,...
                    'measurementVariance' ,measurementVariance               ,...
                    'frequencyCoefficient',frequencyCoefficient             ,...          % Frequency first-order autoregressive coefficient (pole)
                    'amplitudeCoefficient',amplitudeCoefficient             ,...
                    'frequencyMean'       ,frequencyMean     ,...          % A priori mean frequency (Hz)
                    'frequencyRange'      ,[freqMin freqMax]        ,...          % A priori range of possible frequencies (Hz)
                    'dbRange'             ,dbRange ,...
                    'StateInitial',struct(...
                        'phase'           ,initialPhase                 ,...
                        'frequency'       ,initialFrequency     ,...          % (Hz)
                        'amplitudes'      ,linspace(initialAmplitude/2,initialAmplitude,maxHarm)'       ...0.1*randn(nHarmonics,1)...
                        ),...
                    'StateVariance',struct(...                             % Continuous-time equivalent
                        'phase'           ,variancePhase                 ,...
                        'frequency'       ,frequencyVariance              ,...          % (Hz/second)^2/second
                        'amplitudes'      ,amplitudeVariance...       % (amplitude scale)^2/second 
                        ),...
                    'StateVarianceInitial',struct(...                      % Continuous-time equivalent
                        'phase'           ,variancePhase      ,...
                        'frequency'       ,frequencyVariance,.../(1-frequencyCoefficient^2)                ,...          % (Hz/second)^2/second
                        'amplitudes'      ,amplitudeVariance.../(1-amplitudeCoefficient^2)...         % (amplitude scale)^2/second 
                        )...        
                    );
                warn = warning('off');
                inputs.values(a) = TrackModelParameters; 
                warning(warn)
           end
        catch e %--version 10
            errordlg({e.message;'';'Invalid input, use numbers or valid expressions only'}) %--version 10
           % errordlg('Invalid Input, use numbers or valid mathmatical expressions only') %--version 7
            rethrow(e)
            %lasterr %--version 7
    end

        stateIndex = get(inputHandles.stateTypeH, 'value');
        stateList = get(inputHandles.stateTypeH, 'string');
        stateType = strtrim(stateList(stateIndex,:)); %trims spaces
        
        bools.synSpec = get(checkHandles.synSpecH, 'value');
        bools.synSpecHarm = get(checkHandles.synSpecHarmH, 'value');   
        bools.signalComp = get(checkHandles.signalCompH, 'value');
        bools.ampPhase = get(checkHandles.ampPhaseH, 'value');
        bools.sigSpec = get(checkHandles.sigSpecH, 'value');
        bools.sigSpecHarm = get(checkHandles.sigSpecHarmH, 'value');
        bools.plotRes = get(checkHandles.plotResH, 'value');
        %bools.plotPart = get(checkHandles.plotPartH, 'value');
        %bools.freqComp = get(checkHandles.freqCompH, 'value');
%         bools.invertColor = get(checkHandles.invertColorH, 'value');
        
        inputs.bools = bools;
        inputs.stateIndex = stateIndex;
        inputs.stateType = stateType;
        inputs.windowDuration = windowDuration;
        inputs.stateIndex = stateIndex;
end