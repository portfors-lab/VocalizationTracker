function axesHandles = VisualizeTraceData(experiment_data, ...
                                            prefs, ...
                                            test_num, ...
                                            trace_num, ...
                                            plot_flag, ...
                                            output_path, ...
                                            saveFormat, ...
                                            resolution, ...
                                            spike_times, ...
                                            spike_idxs)
%
%function VisualizeTraceData(experiment_data,
%                            prefs,
%                            test_num,
%                            trace_num,
%                            plot_flag,
%                            output_path)
%
%   INPUT ARGUMENTS
%   experiment_data     Bat2Matlab data structure
%   prefs               Bat2Matlab preferences
%   test_num            The number of the test to visualize
%   trace_num           The trace of test test_num to visualize
%   plot_flag           A binary vector indicating the plotting features to
%                       display
%                       Default: [1 0 1 1 0 1 0 1 0 1]
%       flag index 1    Plot the raw MER(micro-electrode recording) signals
%       flag index 2    Plot the smoothed MER power signals (used for spike
%                       detection)
%       flag index 3    Plot the detected spikes
%       flag index 4    Plot the spike histogram
%       flag index 5    Plot the smoothed histogram
%       flag index 6    Plot the stimulus
%       flag index 7    Plot the estimated spike rates
%       flag index 8    Plot the spectrogram of the stimulus
%       flag index 9    Plot the peaks of the spectrogram of the stimulus
%       flag index 10   Plot the dB SPL output of the stimulus
%   output_path         The relative or absolute path, including the base
%                       file name, to print resulting figures to.
%                       If empty, the figure will not be saved.
%                       Default: []

% last modified Amy Boyle 30/5/2011, corrected spike rates plot function
% call to spike times parameter

axesHandles = [];

% pf = [1 1 1 1 1 1 1 1 1 1];
pf = [1 0 1 1 0 1 0 1 0 1];
if exist('plot_flag','var') && ~isempty(plot_flag),
    pf = plot_flag;
end;

if exist('force_analysis','var')
    if isempty(force_analysis)
        force_analysis = false;
    end
else
    force_analysis = false;
end

%Plot flags
plot_signal             = 1;
plot_smoothed_signal    = 2;
plot_detected_spikes    = 3;
plot_histogram          = 4;
plot_smoothed_histogram = 5;
plot_stimulus           = 6;
plot_spike_rates        = 7;
plot_spectrogram        = 8;
plot_peaks              = 9;
plot_dB_SPL             = 10;

plots = zeros(1,4);
if pf(plot_smoothed_signal) || pf(plot_signal) || pf(plot_detected_spikes)
    plots(1) = 1;
end
if pf(plot_histogram) || pf(plot_smoothed_histogram) || pf(plot_stimulus) || pf(plot_spike_rates)
    plots(2) = 1;
end
if pf(plot_spectrogram) || pf(plot_peaks)
    plots(3) = 1;
end
if pf(plot_dB_SPL)
    plots(4) = 1;
end
current_plot = 1;
num_plots = sum(plots);

trace = experiment_data.test(test_num).trace(trace_num);

if plots(1) || plots(2) 
    raw_data = ExtractRawData(experiment_data, ...
                              prefs, ...
                              test_num, ...
                              trace_num);
end

if plots(3) || plots(4)
    [individual_stimulus_signals ...
     summed_stimulus_signals ...
     stimulus_sampling_frequencies] = GenerateStimulus(experiment_data, ...
                                                       prefs, ...
                                                       test_num, ...
                                                       trace_num);
 
    [stimulus_spectrograms ...
     stimulus_spectrogram_time_idxs ...
     stimulus_spectrogram_freq_idxs ...
     stimulus_dB_SPLs] =  GenerateSpectrograms(experiment_data, ...
                                               prefs, ...
                                               test_num, ...
                                               trace_num, ...
                                               summed_stimulus_signals);
end



if exist('output_path','var') && ~isempty(output_path)
%     if ~strcmp(prefs.colormap_name, 'jet')
%         output_path = [output_path ' ' prefs.colormap_name];
%     end
    if prefs.force_plot_visible
        f = figure;
    else
        f = figure('Visible','off');
    end
else
    f = figure;
end
name = [prefs.cell_id 'Test_' int2str(test_num) '_Trace_' int2str(trace_num)];
set(f, 'name', name);
if plots(1) || plots(2)
    if exist('spike_times','var') && ~isempty(spike_times) && exist('spike_idxs','var') && ~isempty(spike_idxs) && ~pf(plot_smoothed_signal)
        WriteStatus(['Using cached spike detections for trace visualization for cell ' prefs.cell_id4_plot]);
    else
        [spike_times ...
         spike_idxs ...
         smoothed_data ...
         peak_thresholds] = CalculateSpikeTimes(experiment_data, ...
                                                prefs, ...
                                                test_num, ...
                                                trace_num, ...
                                                raw_data);
    end
end

title_string = ['Results from cell ' prefs.cell_id4_plot ', Test ' int2str(test_num) ', Trace ' int2str(trace_num)];

if plots(3)
    figure(f)
    if num_plots > 1
        subplot(num_plots,1,current_plot);
    end
%     if current_plot == 1
%         title(title_string);
%     end
    if pf(plot_spectrogram)
        spectrogram = stimulus_spectrograms{test_num,trace_num};
        if prefs.spectrogram_absolute_scaling
            dB_ref = prefs.dB_SPL_ref; %Use this to make the scaling absolute for comparisons
        else
            dB_ref = max(max(spectrogram)); %Use this to maximize color contrast in plot
        end
        %Convert to dB scale for plotting, using the maximum value of the
        %spectrogram as the reference.
        warning off
        spectrogram = 10*log10(spectrogram/dB_ref);
        warning on
        spec_dbMax = max(max(spectrogram));
        if spec_dbMax == -Inf
            spectrogram = spectrogram*0;
        else
            if prefs.spectrogram_absolute_scaling
                spectrogram = max(spectrogram,0); %Use this to make the scaling absolute for comparisons
            else
                %Confine range of spectrogram to dbRange decibels.
%                 spectrogram = max(spectrogram,-prefs.dbRange); %Use this to maximize color contrast in plot
                spectrogram = max(spectrogram,-40);
%                 spectrogram = max(spectrogram,-50);
%                 spectrogram = max(spectrogram,-60);
            end
        end
        figure(f);
        hold on
        imagesc(stimulus_spectrogram_time_idxs{test_num,trace_num}*1000,...
                stimulus_spectrogram_freq_idxs{test_num,trace_num}./1000,...
                spectrogram);
        hold off
        set(gca,'YDir','normal'); 
        % Add vocalization file name as title
        if isfield(experiment_data.test(test_num).trace(trace_num).stimulus,'vocal_call_file') && ...
           ~isempty(experiment_data.test(test_num).trace(trace_num).stimulus.vocal_call_file)
            title(strrep(experiment_data.test(test_num).trace(trace_num).stimulus.vocal_call_file,'_',' '));
        end
        if prefs.spectrogram_absolute_scaling
            caxis([0 70]);
        end
        h = colorbar('East');
        set(h,'YAxisLocation','right');
        color_bar_position = get(h,'Position');
        color_bar_position(1) = 0.855;
        set(h,'Position',color_bar_position);
        set(h,'YTickLabelMode','manual');
        set(h,'YTickMode','manual');
        colorbar_tick_labels = get(h,'YTickLabel');
        new_tick_label = [int2str(str2num(colorbar_tick_labels(end,:))) ' dB'];
        colorbar_tick_labels(end,1:length(new_tick_label)) = new_tick_label;
        set(h,'YTickLabel',colorbar_tick_labels);

    end
    if pf(plot_peaks)
        [stimulus_spectrogram_harmonic_vals ...
         stimulus_spectrogram_harmonic_idx] = FindHarmonics(experiment_data, ...
                                                            prefs, ...
                                                            test_num, ...
                                                            trace_num, ...
                                                            stimulus_spectrograms, ...
                                                            stimulus_spectrogram_freq_idxs);
                                                        

        peak_frequencies = stimulus_spectrogram_harmonic_idx{test_num,trace_num};
        freq_idx = stimulus_spectrogram_freq_idxs{test_num,trace_num}./1000;
        time_idx = stimulus_spectrogram_time_idxs{test_num,trace_num}*1000;
        num_peaks = size(peak_frequencies,1);
        figure(f) %make sure that the desired figure is the current figure
        hold on
        for peak = 1:num_peaks
            harmonic_peaks = peak_frequencies(peak,:);
            non_zero_peak_idx = find(harmonic_peaks);
            non_zero_peak_vals = harmonic_peaks(non_zero_peak_idx);
            non_zero_freq_vals = freq_idx(non_zero_peak_vals);
            h = plot(time_idx(non_zero_peak_idx),non_zero_freq_vals,'r.');
            set(h,'LineWidth',2.0);
        end
        hold off    
    end
    ylabel('Frequency (kHz)');
    xlabel('Time (ms)');
    xlim([0 max(stimulus_spectrogram_time_idxs{test_num,trace_num})*1000]);
    ylim([min(stimulus_spectrogram_freq_idxs{test_num,trace_num})./1000 max(stimulus_spectrogram_freq_idxs{test_num,trace_num})./1000]);
    box on
    if (pf(plot_spectrogram) && ~pf(plot_peaks))
        colormap(prefs.colormap);
%         set(h,'YColor',[1 1 1]);
%         set(h,'XColor',[1 1 1]);
    else
        cmap = colormap(bone);
        cmap = flipud(cmap);
        colormap(cmap);
    end
    current_plot = current_plot + 1;
    axesHandles = [axesHandles gca];
end

if plots(1)
    figure(f)
    if num_plots > 1
        subplot(num_plots,1,current_plot);
    end
    data = raw_data{test_num,trace_num};
    fs = trace.samplerate_ad;
    signal_index = 1:length(data);
    time_index = (signal_index-1)/fs*1000;
    hold on
    if pf(plot_smoothed_signal)
        smoothed_data = smoothed_data{test_num,trace_num};
        h = plot(time_index,smoothed_data');
        for line_num = 1:size(h,1)
            set(h(line_num),'Color',0.7*[1 1 1]);
        end
    end
    figure(f)
    if pf(plot_signal)
        plot(time_index,data');
    end
    if pf(plot_detected_spikes)
        num_realizations = size(data,1);
        trace_spike_idx = spike_idxs{test_num,trace_num};
        trace_spike_times = spike_times{test_num,trace_num};
        if pf(plot_smoothed_signal)
            trace_peak_threshold = peak_thresholds{test_num,trace_num};
        end
        for realization_num = 1:num_realizations
            realization_spike_idx = trace_spike_idx{realization_num};
            realization_spike_times = trace_spike_times{realization_num};
            figure(f)
            if pf(plot_signal)
                realization_data = data(realization_num,:);
                plot(realization_spike_times,realization_data(realization_spike_idx),'*k');
            elseif pf(plot_smoothed_signal)
                realization_data = smoothed_data(realization_num,:);
                plot(spike_times{test_num,trace_num}{realization_num},realization_data(realization_spike_idx),'*k');
                plot(time_index,ones(1,length(time_index))*trace_peak_threshold,'k');
            else
                plot(realization_spike_times,zeros(1,length(realization_spike_times)),'*k');
            end
        end
    end
    hold off
    xlabel('Time (ms)');
    ylabel('MER Recording');
    box on
    current_plot = current_plot + 1;
    axesHandles = [axesHandles gca];
end
    
if plots(2)
    figure(f)
    if num_plots > 1
        subplot(num_plots,1,current_plot);
    end
    legend_strs = {};
    if pf(plot_histogram) || pf(plot_smoothed_histogram) || pf(plot_spike_rates)
        if ~exist('spike_times','var')
            [spike_times ...
             spike_idxs ...
             smoothed_data ...
             peak_thresholds] = CalculateSpikeTimes(experiment_data, ...
                                                    prefs, ...
                                                    test_num, ...
                                                    trace_num, ...
                                                    raw_data);
        end
        [histograms ...
         sweep_histograms ...
         histogram_bin_widths ...
         histogram_bin_centers] = GenerateHistograms(experiment_data, ...
                                                     prefs, ...
                                                     test_num, ...
                                                     trace_num, ...
                                                     spike_times);

        histogram_bin_centers = histogram_bin_centers{test_num,trace_num};
        histogram = histograms{test_num,trace_num};
        %Normalize the histogram to average spikes per bin
        histogram = histogram / size(sweep_histograms{test_num,trace_num},1);
    end
    if pf(plot_histogram)
        figure(f)
        hold on
        bar(histogram_bin_centers, histogram,'hist');
        hold off
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor',0.5*[1 1 1],'EdgeColor',0.5*[1 1 1]);
        legend_strs{1,length(legend_strs)+1} = ['PSTH (mean spikes / ' num2str(prefs.histogram_bin_width,2) ' ms bin)'];
    end

    if pf(plot_smoothed_histogram)
        histogram_area = sum(histogram)*histogram_bin_widths{test_num,trace_num};
        y = SmoothSeries(histogram_bin_centers,histogram,histogram_bin_centers,2);
        figure(f)
        hold on
        plot(histogram_bin_centers,y,'b');
        hold off
        legend_strs{1,length(legend_strs)+1} = 'Smoothed PSTH';
    end
    
    if pf(plot_spike_rates)
        if ~exist('spike_idxs','var')
            [spike_times ...
             spike_idxs ...
             smoothed_data ...
             peak_thresholds] = CalculateSpikeTimes(experiment_data, ...
                                                    prefs, ...
                                                    test_num, ...
                                                    trace_num, ...
                                                    raw_data);
        end
        [spike_rates, ...
         spike_rate_averages, ...
         spike_rate_sampling_frequencies] = GenerateSpikeRates(experiment_data, ...
                                                               prefs, ...
                                                               test_num, ...
                                                               trace_num, ...
                                                               spike_times);

        signal_duration = trace.record_duration; %In milliseconds
        spike_rate = spike_rate_averages{test_num,trace_num};  
        spike_rate = spike_rate / 1000; %In spikes per millisecond
        figure(f)
        hold on;
        plot_idx = linspace(0,signal_duration,length(spike_rate));
        plot(plot_idx,spike_rate,'r');
        hold off;
        legend_strs{1,length(legend_strs)+1} = 'Spike Rate (spikes / ms)';
    end
    
    ylabel('Response');
    
    if pf(plot_stimulus) && ~isempty(trace.stimulus)
        if ~exist('stimulus_signals','var')
            [individual_stimulus_signals ...
             summed_stimulus_signals ...
             stimulus_sampling_frequencies] = GenerateStimulus(experiment_data, ...
                                                               prefs, ...
                                                               test_num, ...
                                                               trace_num);
        end
        figure(f)
        hold on;
        if (pf(plot_histogram) || pf(plot_smoothed_histogram) || pf(plot_spike_rates)) 
            max_hist = max(max(histogram),1);
            offset = max_hist + (0.1 * max_hist);
            for stim_num = 1:size(trace.stimulus,2)
                stimulus_begin = trace.stimulus(stim_num).delay;
                if ~isempty(individual_stimulus_signals{test_num,trace_num,stim_num})
                    stim4plot = individual_stimulus_signals{test_num,trace_num,stim_num};
                    stimulus_end = stimulus_begin + ...
                                   length(stim4plot)/stimulus_sampling_frequencies{test_num,trace_num}*1000;
                    stimulus_end = min(stimulus_end,trace.record_duration);
                    plot_idx = linspace(stimulus_begin,stimulus_end,size(stim4plot,2));
                    %Normalize stimulus to 1
                    stim4plot = stim4plot / max(stim4plot);
                    %Shrink range of stimulus for scale of histogram
                    stim4plot = stim4plot * 0.04 * max_hist;
                    %Offset stimulus in plot
                    stim4plot = stim4plot + offset;
                    plot(plot_idx,stim4plot,'k');
                else
                    stimulus_end = stimulus_begin + trace.stimulus(stim_num).duration -1;
                    stim4plot = zeros(1,trace.stimulus(stim_num).duration);
                    %Offset stimulus in plot
                    stim4plot = stim4plot + offset;
                    h = plot(stimulus_begin:stimulus_end,stim4plot,'k');
                    set(h,'lineWidth',5);
                end
                offset = offset + (0.1 * max_hist);
            end
                ylim([0 offset]);
        else
            offset = 0;
            for stim_num = 1:size(trace.stimulus,2)
                stimulus_begin = trace.stimulus(stim_num).delay;
                stimulus_end = stimulus_begin + trace.stimulus(stim_num).duration - 1;
                if ~isempty(individual_stimulus_signals{test_num,trace_num,stim_num})
                    stim4plot = individual_stimulus_signals{test_num,trace_num,stim_num};
                    plot_idx = linspace(stimulus_begin,stimulus_end,size(stim4plot,2));
                    %Normalize stimulus to 0.9
                    stim4plot = (stim4plot / max(stim4plot))*0.9;
                    %Offset stimulus in plot
                    stim4plot = stim4plot + offset;
                    h = plot(plot_idx,stim4plot,'k');
                else
                    stimulus_end = stimulus_begin + trace.stimulus(stim_num).duration -1;
                    stim4plot = zeros(1,trace.stimulus(stim_num).duration);
                    %Offset stimulus in plot
                    stim4plot = stim4plot + offset;
                    h = plot(stimulus_begin:stimulus_end,stim4plot,'k');
                    set(h,'lineWidth',5);
                end
                offset = offset + 2;
            end
            xlim([0 trace.record_duration]);
            ylim([-1 (offset-1)]);
            set(gca,'YTick',[]);
            ylabel('Normalized Stimulus');
        end
        hold off;
        if isfield(experiment_data.test(test_num).trace(trace_num).stimulus,'vocal_call_file') && ...
                ~isempty(experiment_data.test(test_num).trace(trace_num).stimulus.vocal_call_file)
            legend_strs{1,length(legend_strs)+1} = strrep(experiment_data.test(test_num).trace(trace_num).stimulus.vocal_call_file,'_',' ');
        else
            legend_strs{1,length(legend_strs)+1} = 'Stimulus';
        end
    end
    figure(f)
    xlabel('Time (ms)');
    legend(legend_strs,'Location','East');
    box on
    current_plot = current_plot + 1;
    axesHandles = [axesHandles gca];
end

if plots(4)
    figure(f)
    if num_plots > 1
        subplot(num_plots,1,current_plot);
    end
    hold on
    plot(stimulus_spectrogram_time_idxs{test_num,trace_num}*1000,stimulus_dB_SPLs{test_num,trace_num});
    hold off
    ylabel('Stimulus Loudness (dB SPL)');
    xlabel('Time (ms)');
    xlim([0 max(stimulus_spectrogram_time_idxs{test_num,trace_num}*1000)]);
    ylim([prefs.dbMin prefs.dbMax]);
    box on
    current_plot = current_plot + 1;
    axesHandles = [axesHandles gca];
end

% suptitle(title_string);

if exist('output_path','var') && ~isempty(output_path)

    path = [output_path name];
%     PrintFigure(output_path,[],10,8,[],f);
%     PrintFigure(output_path,'jpeg',5,6,[],f);
%     PrintFigure(output_path,'epsc',6,5,[],f);
    if exist('saveFormat', 'var')
        image_format = saveFormat;
    else
        image_format = 'pdf'; %default
    end
    if ~exist('resolution')
        resolution = 300;
    end
    PrintFigure(path,image_format,5,7,resolution,f);
    close(f);
%     saveas(f,output_path,'fig');
end

