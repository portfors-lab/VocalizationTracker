function [handlesStructure handlesVector] = PlotHarmonicStates(States,varargin)
%PlotHarmonicStates: Plot the states contained in the Harmonic Model structure 
%
%   [] =              PlotHarmonicStates(ModelParameters,
%                                        'ParameterName',value,...);
% 
% Input Parameters
%   States            Structure containing all of the parameters that 
%                     describe the statistical model.
%
% Optional Parameters
%   'stateType'       The state type in the harmonic model to use.
%                     Default: 'Filtered'
%   'frequencyRange'  The frequency range of the spectrographic plots
%   'trackedSignal'   The original signal that the harmonic model
%   'spectrogramType' Specifies the format of the plotted spectrograms.
%                     0 = square root of the power spectral density (default)
%                     1 = power spectral density
%                     2 = abs(power spectral density) on a decibel scale
%   'dbRange'         The decibel range of the spectrogram when plotted
%                     on a decibel scale.
%                     Default: 50 dB
%   'windowDuration'  The duration, in seconds, of the window used in the
%                     spectrogram calculation.
%                     Default: 10/frequencyMean of the model.
%   'colorMap'        The colorMap to use in the spectrogram calculations,
%                     passed in as a string.
%                     Default: 'jet'
%   'invertColorMap'  Inverts the RGB intensities of the colormap.
%                     Default: false
%   'trueStates'      Pass is model states for the true, synthetic signal.
%                     By default, this true frequency is plotted against
%                     the frequency provided in States and the NMSE is
%                     displayed.
%
%  The following options enable/disable plotting options. All are 'true'
%  by default.
%
%   'plotSyntheticSpectrogram'            Plotting of the synthetic signal's
%                                         spectrogram
%   'plotSyntheticSpectrogramHarmonics'   Plotting of the synthetic signal's
%                                         harmonics, overlaid on the spectrogram
%   'plotSignalCompare'                   Plots y and the synthetic signal
%                                         on the same plot
%   'plotAmplitudesAndPhases'             Plots the amplitudes and phasees of
%                                         the model states.
%   'plotSignalSpectrogram'               Plots the spectrogram of y
%   'plotSignalSpectrogramHarmonics'      Plots the synthetic signal's harmonics
%                                         overlaid on the spectrogram
%   'plotResiduals'                       Plots the spectrogram of the difference
%                                         between y and the synthetic signal
%   'plotParticles'                       Plots the particle estimates for
%                                         estimates provided by a particle filter
%   'plotFrequencyComparison'             Plots the true frequency 

%   Example: To be written.
%
%   Reference: To be completed.
%
%   See also NonparametricSpectrogram.

%   Tags: Filter, Estimator

% Last Modified Jan 2011 Amy Boyle

handlesVector = [];
%==============================================================================
% Error Checking
%==============================================================================
if nargin<1
    help PlotHarmonicStates;
    return;
end

%==============================================================================
% Process Function Arguments
%==============================================================================
stateType = 'Filtered';
sampleRate = States.ModelParameters.sampleRate;
frequencyRange = [0 sampleRate/2];
y = [];
trueStates = [];
colorMap = 'jet';
invertColorMap = false;
windowDuration = 10/States.ModelParameters.frequencyMean;
dbRange = 50;
spectrogramType = 0;

plotSyntheticSpectrogram = true;
plotSyntheticSpectrogramHarmonics = true;
plotSignalCompare = true;
plotAmplitudesAndPhases = true;
plotSignalSpectrogram = true;
plotSignalSpectrogramHarmonics = true;
plotResiduals = true;
plotParticles = true;
plotTrueSignal = true;
plotParticles = true;

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
            case lower('stateType'),        stateType       = parameterValue;     
            case lower('frequencyRange'),   frequencyRange  = parameterValue;
            case lower('trackedSignal'),    y               = parameterValue;
            case lower('colorMap'),         colorMap        = parameterValue;
            case lower('trueStates'),       trueStates      = parameterValue;
            case lower('windowDuration'),   windowDuration  = parameterValue;
            case lower('spectrogramType'),  spectrogramType = parameterValue;
            case lower('dbRange'),          dbRange         = parameterValue;
            case lower('invertColorMap'),   invertColorMap  = parameterValue;
            case lower('plotSyntheticSpectrogram'),            plotSyntheticSpectrogram  = parameterValue;
            case lower('plotSyntheticSpectrogramHarmonics'),   plotSyntheticSpectrogramHarmonics  = parameterValue;
            case lower('plotSignalCompare'),                   plotSignalCompare  = parameterValue;
            case lower('plotAmplitudesAndPhases'),             plotAmplitudesAndPhases  = parameterValue;
            case lower('plotSignalSpectrogram'),               plotSignalSpectrogram  = parameterValue;
            case lower('plotSignalSpectrogramHarmonics'),      plotSignalSpectrogramHarmonics  = parameterValue;
            case lower('plotResiduals'),                       plotResiduals  = parameterValue;
            case lower('plotParticles'),                       plotParticles  = parameterValue;
            case lower('plotTrueSignal'),                      plotTrueSignal  = parameterValue;
            case lower('plotParticles'),                       plotParticles  = parameterValue;
            case lower('plotFrequencyComparison'),             plotFrequencyComparison = parameterValue;
            case lower('plotMultiSections'),                   plotMultiSections = parameterValue;
            otherwise,                      error(['Unrecognized property: ' parameterName]);
        end
    end
end

%==============================================================================
% Preprocessing
%==============================================================================
if ~isfield(States,stateType)
    error(['There is no state type named ' stateType]);
end

if invertColorMap
    colormapString = ['colormap(1-' colorMap ')'];
else
    colormapString = ['colormap(' colorMap ')'];
end

plotStates = getfield(States,stateType);
nHarmonics = size(plotStates.amplitudesCos,2);
nSamples = length(plotStates.frequency);
%Generate synthesized signal from harmonic states
yp = zeros(1,nSamples);
k = (1:nHarmonics).';
for c1=1:nSamples
    amplitudesCos = plotStates.amplitudesCos(c1,:);
    amplitudesSin = plotStates.amplitudesSin(c1,:);
    phase         = plotStates.phase(c1);  
    yp(c1)        = amplitudesCos*cos(k*phase) + amplitudesSin*sin(k*phase);
end

%Generate "true" signal from true states
if ~isempty(trueStates)
    y = zeros(1,nSamples);
    k = (1:nHarmonics).' ;
    for c1=1:nSamples
        amplitudesCos = trueStates.True.amplitudesCos(c1,:);
        amplitudesSin = trueStates.True.amplitudesSin(c1,:);
        phase         = trueStates.True.phase(c1);             
        % Note: These states are assumed to be generated by GenerateHarmonicModelStates(),
        % which stores the measurement noise in a seperate field.
        y(c1)        = amplitudesCos*cos(k*phase) + amplitudesSin*sin(k*phase) + trueStates.True.measurementNoise(c1);  
  
    end 
end

y = y(:);
yp = yp(:);
n          = 1:nSamples;
t          = (n-0.5)/sampleRate;
plotFrequencies     = plotStates.frequency;
if frequencyRange(2) > 40e3
    %The spectrogram plots in kHz if the frequency range is greater than 40 kHz
    %which necessitates this scaling
    plotFrequencies = plotFrequencies/1000;
end
label               = States.label;

%Parameters for spectrograph plotting
plotFreqSamples = 2^10;
plotTimeSamples = 2^10;

%==============================================================================
% Plot the true and re-synthesized signals
%==============================================================================

if plotSignalCompare && ~isempty(y)
    signalNMSE = mean((y-yp).^2)/var(y);

    statesLabel           = States.label;
    
    fh = figure('numbertitle', 'off', 'name', 'Signal Comparison');
    handlesStructure.signalCompare = fh;
    handlesVector = [handlesVector fh];
    FigureSet;
    CleanUpToolbar(fh)
    colors = DistinctColors(5);
    h = plot(t,y,t,yp);
    for c1=1:length(h)
        set(h(c1),'LineWidth',1.5);
        set(h(c1),'Color',colors(c1,:));
    end
    xlim([0 nSamples/sampleRate]);
    % ylim([-max(y)*2 max(y)*2]);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(['Signal Comparison For ' statesLabel ': NMSE = ' num2str(signalNMSE)]);
    box off;
    AxisSet;
    legend(h,'Original Signal','Synthesized Signal');
end

%==============================================================================
% Amplitudes
%==============================================================================
if plotAmplitudesAndPhases
    fh = figure('numbertitle', 'off', 'name', 'Amplitudes and Phases');
    handlesStructure.amplitudesPhases = fh;
    handlesVector = [handlesVector fh];
    FigureSet;
    CleanUpToolbar(fh);
    colors = DistinctColors(nHarmonics);
    legendHandles = nan(nHarmonics,1);
    legendLabels  = cell(nHarmonics,1);
    for c1=1:nHarmonics
        subplot(2,1,1);       
            amplitudeCosEstimate = plotStates.amplitudesCos(:,c1);
            amplitudeSinEstimate = plotStates.amplitudesSin(:,c1);
            c = amplitudeCosEstimate + j*amplitudeSinEstimate;
            cAmplitudeEstimate = abs(c);
            cPhaseEstimate     = angle(c);

            h = plot(t,cAmplitudeEstimate);
            set(h,'Color',colors(c1,:));
            set(h,'LineWidth',1.5);
            legendHandles(c1) = h;
            legendLabels{c1} = sprintf('%d',c1);        
            hold on;
            ylabel('Amplitudes');
            xlabel('Time (s)');
            title([label ' Amplitude and Phase Estimates']);
            xlim([0 t(end)]);        
        subplot(2,1,2);
            h = plot(t,cPhaseEstimate);
            set(h,'Color',colors(c1,:));
            set(h,'LineWidth',1.5);
            hold on;
            ylabel('Phases (degrees)');
            xlabel('Time (s)');
            xlim([0 t(end)]);
            AxisLines;
    end
    hold off;
    subplot(2,1,1);
    hold off;
    legend(legendHandles,legendLabels,'Location','NorthWest');
    AxisSet;
end

%==============================================================================
% Nonparametric Spectrogram of the tracked signal with the frequency estimates
% of the harmonics overlayed.
%==============================================================================
if plotSignalSpectrogramHarmonics || plotSyntheticSpectrogramHarmonics
    if ~exist('plotMultiSections', 'var')
        plotMultiSections = false;
    else
      nHarmSections = States.nHarmSections;  
    end
end

if (plotSignalSpectrogram || plotSignalSpectrogramHarmonics)  && ~isempty(y)
    NonparametricSpectrogram(y, ...
                             sampleRate, ...
                             'nTimes',plotTimeSamples, ...
                             'nFrequencies',plotFreqSamples, ...
                             'windowDuration',windowDuration, ...
                             'plotType',1,...
                             'spectrogramType',spectrogramType, ...
                             'dbRange', dbRange, ...
                             'frequencyRange', frequencyRange);
    fh= gcf;
    handlesStructure.signalSpectrogram = fh;
    handlesVector = [handlesVector fh];
    %edited part
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotSignalSpectrogramHarmonics
        hold on;
        if strcmp(colorMap,'bone') || strcmp(colorMap,'gray')
            if plotMultiSections
                for a = 1:length(plotMultiSections)-1;
                    tmpPlotFrequencies = plotFrequencies(plotMultiSections(a):plotMultiSections(a+1));
                    tmpT = t(plotMultiSections(a):plotMultiSections(a+1));
                    for c1=1:nHarmSections(a) 
                        h = plot(tmpT,c1*tmpPlotFrequencies,'r');
                        set(h,'LineWidth',1.5);                       
                    end
                end
            else
                for c1=1:nHarmonics
                    h = plot(t,c1*plotFrequencies,'r');
                    set(h,'LineWidth',1.5);
                end
            end
        else
            if plotMultiSections
                for a = 1:length(plotMultiSections)-1;
                    tmpPlotFrequencies = plotFrequencies(plotMultiSections(a):plotMultiSections(a+1));
                    tmpT = t(plotMultiSections(a):plotMultiSections(a+1));
                    for c1=1:nHarmonics(a)                    
                        h = plot(tmpT,c1*tmpPlotFrequencies,'k');
                        set(h,'LineWidth',1.5);
                        h = plot(tmpT,c1*tmpPlotFrequencies,'w');
                        set(h,'LineWidth',1.0);                 
                    end
                end
            else
                for c1=1:nHarmonics
                    h = plot(t,c1*plotFrequencies,'k');
                    set(h,'LineWidth',1.5);
                    h = plot(t,c1*plotFrequencies,'w');
                    set(h,'LineWidth',1.0);
                end
            end
        end
        hold off;
        title([label ' Harmonic Frequency Estimates']);
        set(fh, 'name', 'Harmonic Frequency Estimates')
    else
        title('Original Signal');
        set(fh, 'name', 'Orginial Signal')
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CleanUpToolbar(fh)
    AxisSet;
    FigureSet;
    eval(colormapString);
    colorAxis = caxis;
end

%==============================================================================
% Nonparametric Spectrogram of the synthesized signal
%==============================================================================

if plotSyntheticSpectrogram || plotSyntheticSpectrogramHarmonics
    NonparametricSpectrogram(yp, ...
                             sampleRate, ...
                             'nTimes',plotTimeSamples, ...
                             'nFrequencies',plotFreqSamples, ...
                             'windowDuration',windowDuration, ...
                             'plotType',1,...
                             'spectrogramType',spectrogramType, ...
                             'dbRange', dbRange, ...
                             'frequencyRange', frequencyRange);
    fh = gcf;
    handlesStructure.synthSpectrogram = fh;
    handlesVector = [handlesVector fh];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                         
    if plotSyntheticSpectrogramHarmonics
        hold on;
        if strcmp(colorMap,'bone') || strcmp(colorMap,'gray')
            if plotMultiSections
                for a = 1:length(plotMultiSections)-1;
                    tmpPlotFrequencies = plotFrequencies(plotMultiSections(a):plotMultiSections(a+1));
                    tmpT = t(plotMultiSections(a):plotMultiSections(a+1));
                    for c1=1:nHarmSections(a) 
                        h = plot(tmpT,c1*tmpPlotFrequencies,'r');
                        set(h,'LineWidth',1.5);                       
                    end
                end
            else
                for c1=1:nHarmonics
                    h = plot(t,c1*plotFrequencies,'r');
                    set(h,'LineWidth',1.5);
                end
            end
        else
            if plotMultiSections
                for a = 1:length(plotMultiSections)-1;
                    tmpPlotFrequencies = plotFrequencies(plotMultiSections(a):plotMultiSections(a+1));
                    tmpT = t(plotMultiSections(a):plotMultiSections(a+1));
                    for c1=1:nHarmonics(a)                    
                        h = plot(tmpT,c1*tmpPlotFrequencies,'k');
                        set(h,'LineWidth',1.5);
                        h = plot(tmpT,c1*tmpPlotFrequencies,'w');
                        set(h,'LineWidth',1.0);                 
                    end
                end
            else
                for c1=1:nHarmonics
                    h = plot(t,c1*plotFrequencies,'k');
                    set(h,'LineWidth',1.5);
                    h = plot(t,c1*plotFrequencies,'w');
                    set(h,'LineWidth',1.0);
                end
            end
        end
        hold off;
        title([label ' Harmonic Frequency Estimates']);
        set(fh, 'name', 'Harmonic Frequency Estimates')
    else
        title('Original Signal');
        set(fh, 'name', 'Orginial Signal')
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    fh = gcf;
    set(fh, 'name', 'Synthesized Signal Spectrogram')
    CleanUpToolbar(fh)
%     eval(colormapString);
    AxisSet;
    FigureSet;
end

%==============================================================================
% Plot the Residuals
%==============================================================================
if plotResiduals && ~isempty(y)
    S = NonparametricSpectrogram(y, ...
                                 sampleRate, ...
                                 'nTimes',plotTimeSamples, ...
                                 'nFrequencies',plotFreqSamples, ...
                                 'windowDuration',windowDuration);
    dbRef = max(max(abs(S)));

    NonparametricSpectrogram(y-yp, ...
                             sampleRate, ...
                             'nTimes',plotTimeSamples, ...
                             'nFrequencies',plotFreqSamples, ...
                             'windowDuration',windowDuration, ...
                             'plotType',1,...
                             'spectrogramType',spectrogramType, ...
                             'dbRange', dbRange, ...
                             'dbReference', dbRef, ...
                             'frequencyRange', frequencyRange);
    title([label ' Residuals']);
    fh = gcf;
    handlesStructure.residuals = fh;
    handlesVector = [handlesVector fh];
    set(fh, 'name', 'Residuals')
    CleanUpToolbar(fh)
    colorAxis = caxis;
    AxisSet;
    FigureSet;
    eval(colormapString);
end

%------------------------------------------------------------------------------
% If present, Plot the Particle Trajectories
%------------------------------------------------------------------------------
if plotParticles && isfield(States,'Particles')
    fh = figure('numbertitle', 'off', 'name', 'Particles');
    handlesStructure.particles = fh;
    handlesVector = [handlesVector fh];
    FigureSet;
    CleanUpToolbar(fh)
    nParticles = size(States.Particles.frequency,2);
    k          = 1:nSamples;
    t          = (k-0.5)/sampleRate;
    yMaximum   = max(max(States.Particles.frequency));
    iResampled = find(States.Particles.bResample);
    h = plot([1;1]*t(iResampled),[0;yMaximum]*ones(1,length(iResampled)),'r');
    hold on;
        h = plot(t,States.Particles.frequency,'b');
        set(h,'LineWidth',0.1);
        
%         h = plot(t,nParticlesEffective/nParticles,'m');
%         set(h,'LineWidth',2);

    hold off;
    xlim([0 nSamples/sampleRate]);
    ylim([0 yMaximum]);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    sTitle = sprintf('No. Particles:%d  No. Samples:%d',nParticles,nSamples);
    title(sTitle);
    AxisSet;
    box off;
end

%==============================================================================
% Plot the true and estimated fundamental frequencies
%==============================================================================

if plotTrueSignal && ~isempty(trueStates)
    trueFrequencies          = trueStates.True.frequency(:);
    estimatedFrequencies     = plotStates.frequency(:);

    frequencyNMSE = mean((trueFrequencies-estimatedFrequencies).^2)...
                    /var(trueFrequencies);

    fh = figure('numbertitle', 'off', 'name', 'Frequency Comparison');
    handlesStructure.frequencyCompare = fh;
    handlesVector = [handlesVector fh];
    FigureSet;
    CleanUpToolbar(fh)
    h = plot(t,trueFrequencies,t,estimatedFrequencies);
    for c1=1:length(h)
        set(h(c1),'LineWidth',1.5);
        set(h(c1),'Color',colors(c1,:));
    end
    xlim([0 nSamples/sampleRate]);
    ylim([0 max(trueFrequencies)*2]);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(['Fundamental Frequency Comparison For ' label ': NMSE = ' num2str(frequencyNMSE)]);
    box off;
    AxisSet;
    legend(h,'True',[label ' ' stateType ' Estimate']);
end
