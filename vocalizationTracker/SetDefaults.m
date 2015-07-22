function defaults = SetDefaults(ntabs, sampleRate)
%Defines the default values for the user defined parameters

%values =zeros(1,ntabs);
for a=1:ntabs

                values(a) = struct(...
                    'nHarmonics'          ,1        ,...
                    'measurementVariance' ,0.001               ,...
                    'frequencyCoefficient',0.999999             ,...          % Frequency first-order autoregressive coefficient (pole)
                    'amplitudeCoefficient',0.9999            ,...
                    'frequencyMean'       ,50000     ,...          % A priori mean frequency (Hz)
                    'frequencyRange'      ,[0 sampleRate/2]        ,...          % A priori range of possible frequencies (Hz)
                    'dbRange'             ,50 ,...
                    'StateInitial',struct(...
                        'phase'           ,0                 ,...
                        'frequency'       ,38000     ,...          % (Hz)
                        'amplitudes'      ,0       ...0.1*randn(nHarmonics,1)...
                        ),...
                    'StateVariance',struct(...                             % Continuous-time equivalent
                        'phase'           ,1                 ,...
                        'frequency'       ,10^8              ,...          % (Hz/second)^2/second
                        'amplitudes'      ,1 ...       % (amplitude scale)^2/second 
                        ),...
                    'StateVarianceInitial',struct(...                      % Continuous-time equivalent
                        'phase'           ,1      ,...
                        'frequency'       ,10^8,.../(1-frequencyCoefficient^2)                ,...          % (Hz/second)^2/second
                        'amplitudes'      ,3 .../(1-amplitudeCoefficient^2)...         % (amplitude scale)^2/second 
                        )...        
                    );
end


bools.synSpec = true;
bools.synSpecHarm = false;   
bools.signalComp = false;
bools.ampPhase = false;
bools.sigSpec = true;
bools.sigSpecHarm = true;
bools.plotRes = true;
bools.plotPart = false;
bools.freqComp = false;
bools.invertColor = false;

defaults.windowDuration = 1024/sampleRate;
defaults.stateIndex = 1;
defaults.colormap = 'jet';
defaults.values = values;
defaults.bools = bools;

end