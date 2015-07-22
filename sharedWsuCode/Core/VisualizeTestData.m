function [unique_frequencies, ...
          unique_attenuations, ...
          test_data, ...
          tune_classes] = VisualizeTestData(experiment_data, ...
                                            prefs, ...
                                            test_nums, ...
                                            plot_flag, ...
                                            output_path, ...
                                            saveFormat, ...
                                            resolution, ...
                                            colorRange, ...
                                            spike_times, ...
                                            spike_idxs, ...
                                            primary_tuning_range)
%
%function VisualizeTestData(experiment_data, ...
%                           prefs, ...
%                           test_nums, ...
%                           plot_flag)
%
%   INPUT ARGUMENTS
%   experiment_data     Bat2Matlab data structure
%   prefs               Bat2Matlab preferences
%   test_num            The number of the tests to visualize
%   plot_flag           A binary vector indicating the plotting features to
%                       display
%                       Default: [1 1 0 1 0]
%       flag index 1    Plot spike counts for each trace where the
%                       independent axes are frequency and attenuation of
%                       the stimulus
%       flag index 2    Plot PSTHs. Each plot is for a constant
%                       attenuation and a range of frequencies.
%       flag index 3    Plot the results as a surface plot
%       flag index 4    Plot the results as a contour plot
%       flag index 5    Plot the results as a image map
%       output_path     The relative or absolute path, including the base
%                       file name, to print resulting figures to.
%                       If multiple figures are generated, each file will
%                       be appended with '_F' where N is the figure number.
%                       If empty, the figure will not be saved.
%                       Default: []
%   
%Visualize the results of a single or two tone Batlab test. The results can
%be viewed as single and two tone frequency tuning curves. They can also be
%viewed as PSTHs for multiple frequency tests on a single plot. Both plot
%types can be output as surface plots, contour plots, or image maps.


%last modifited Amy Boyle 5/20/11, changed command display colors

warnColor = 'orange';

pf = [1 1 0 1 0]; %Default
pf = [1 0 0 1 0];
if exist('plot_flag','var') && ~isempty(plot_flag)
    pf = plot_flag;
end

%Plot flags
plot_spike_counts   = 1;
plot_histograms     = 2;
plot_surface        = 3;
plot_contour        = 4;
plot_imagemap       = 5;

plot_loudness       = true; %false plots attenuations
plot_percentage_increase_for_2_tone_tests = true; %False plots spike count difference relative to excitory signal
plot_statistical_differences = false;
alpha = 0.01;

% image_format = 'jpeg'; %Default
% image_format = 'epsc2'; %Default
if exist('output_path', 'var')
    if ~exist('resolution', 'var')
        resolution = 300;
    end
    if exist('saveFormat', 'var')
        image_format = saveFormat;
    else
        image_format = 'pdf'; %default
    end
end
calculate_tuning_type = true;

% figure_exists = get(0,'CurrentFigure');
% cm = colormap(prefs.colormap);
% if isempty(figure_exists)
%     close;
% end

if ~exist('test_nums','var')
    test_nums = [];
end

if isempty(test_nums)
    tests_2_process = 1:size(experiment_data.test,2);
else
    tests_2_process = test_nums;
end

if exist('spike_times','var') && ~isempty(spike_times) && exist('spike_idxs','var') && ~isempty(spike_idxs)
    calculate_spike_times = false;
else
    calculate_spike_times = true;
end

%Supported test types and identifiers
%'tone'     = 1;
%'twotone'  = 2;

colorbar_label = ['Mean spikes per presentation'];
if ~exist('output_path','var')
    output_path = [];
end
base_output_path = output_path;

all_unique_frequencies = {};
all_unique_attenuations = {};
all_test_data = {};
all_tune_classes = {};

test_idx = 0;

for test_num = tests_2_process
    test_idx = test_idx + 1;
    
    if calculate_spike_times        
        raw_data = ExtractRawData(experiment_data,prefs,test_num,[]);
    
        [spike_times ...
         spike_idxs] = CalculateSpikeTimes(experiment_data, ...
                                           prefs, ...
                                           test_num, ...
                                           [], ...
                                           raw_data);
        clear raw_data;                          
    end
    
    
    [histograms ...
     sweep_histograms ...
     histogram_bin_widths ...
     histogram_bin_centers] = GenerateHistograms(experiment_data, ...
                                                 prefs, ...
                                                 test_num, ...
                                                 [], ...
                                                 spike_times);
    
    num_histogram_bins = size(histograms{test_num,1},2);
                                             
    test = experiment_data.test(test_num);
    traces = test.trace;
    num_traces = length(traces);
    if strcmp(test.testtype,'tone')
        testtype = 1;
        test_label = 'Single Tone';
        
    elseif strcmp(test.testtype,'twotone')
        testtype = 2;
        if num_traces < 3
%             display(['Test ' int2str(test_num) ' does not have standard two-tone test structure. Continuing.']);
            WriteStatus(['Test ' int2str(test_num) ' does not have standard two-tone test structure. Continuing.'], warnColor);
            continue;
        end
        %Find frequency and attenuation of excitory signal in twotone tests
        last_trace = traces(num_traces);
        second_trace = traces(2);
        third_trace = traces(3);
        if length(second_trace.stimulus) ~= 1
%             display(['Second trace of two tone test ' int2str(test_num) ' does not have single tone.']);
            WriteStatus(['Second trace of two tone test ' int2str(test_num) ' does not have single tone.'], warnColor);
            continue;
        end
        if length(third_trace.stimulus) ~= 1
%             display(['Third trace of two tone test ' int2str(test_num) ' does not have single tone.']);
            WriteStatus(['Third trace of two tone test ' int2str(test_num) ' does not have single tone.'], warnColor);
            continue;
        end
        if length(last_trace.stimulus) ~= 2
%             display(['Last trace of two tone test ' int2str(test_num) ' does not have two tones.']);
            WriteStatus(['Last trace of two tone test ' int2str(test_num) ' does not have two tones.'], warnColor);
            continue;
        end
        
        %Locate the excitory signal test
        if last_trace.stimulus(1).frequency == second_trace.stimulus(1).frequency && last_trace.stimulus(1).attenuation == second_trace.stimulus(1).attenuation
            excitory_frequency = second_trace.stimulus(1).frequency;
            excitory_attenuation = second_trace.stimulus(1).attenuation;
            two_tone_changing_idx = 2;
        elseif last_trace.stimulus(2).frequency == third_trace.stimulus(1).frequency && last_trace.stimulus(2).attenuation == third_trace.stimulus(1).attenuation
            excitory_frequency = third_trace.stimulus(1).frequency;
            excitory_attenuation = third_trace.stimulus(1).attenuation;
            two_tone_changing_idx = 1;
        else
%             display(['No excitory signal found for two tone test ' int2str(test_num)]);
            WriteStatus(['No excitory signal found for two tone test ' int2str(test_num) ''], warnColor);
            continue;
        end
        test_label = 'Two Tone';
%*PDR***************begin inset (020408)
    elseif strcmp(test.testtype,'threetone')
        testtype = 3;
        if num_traces < 3
%             display(['Test ' int2str(test_num) ' does not have standart tree-tone test structure. Continuing.']);
            WriteStatus(['Test ' int2str(test_num) ' does not have standart tree-tone test structure. Continuing.'], warnColor);
            continue;
        end
        %Find frequency and attenuation of excitory signal in twotone tests
        last_trace = traces(num_traces);
        second_trace = traces(2);
        third_trace = traces(3);
        fourth_trace = traces(4);
        if length(second_trace.stimulus) ~= 1
%             display(['Second trace of three tone test ' int2str(test_num) ' does not have single tone.']);
            WriteStatus(['Second trace of three tone test ' int2str(test_num) ' does not have single tone.'], warnColor);
            continue;
        end
        if length(third_trace.stimulus) ~= 1
%             display(['Third trace of three tone test ' int2str(test_num) ' does not have single tone.']);
            WriteStatus(['Third trace of three tone test ' int2str(test_num) ' does not have single tone.'], warnColor);
            continue;
        end
        if length(fourth_trace.stimulus) ~= 1
%             display(['Fourth trace of three tone test ' int2str(test_num) ' does not have single tone.']);
            WriteStatus(['Fourth trace of three tone test ' int2str(test_num) ' does not have single tone.'], warnColor);
            continue;
        end
        if length(last_trace.stimulus) ~= 3
%             display(['Last trace of three tone test ' int2str(test_num) ' does not have three tones.']);
            WriteStatus(['Last trace of three tone test ' int2str(test_num) ' does not have three tones.'], warnColor);
            continue;
        end
        
        %Locate the changing signal test
        if last_trace.stimulus(1).frequency ~= second_trace.stimulus(1).frequency
            two_tone_changing_idx = 1;
        elseif last_trace.stimulus(2).frequency ~= third_trace.stimulus(1).frequency
            two_tone_changing_idx = 2;
        elseif last_trace.stimulus(3).frequency ~= fourth_trace.stimulus(1).frequency
            two_tone_changing_idx = 3;
        else
%             display(['No changing signal found for three tone test ' int2str(test_num)]);
            WriteStatus(['No changing signal found for three tone test ' int2str(test_num) ], warnColor);
            continue;
        end
        %Locate the excitory signal test
        if two_tone_changing_idx ~= 1
            excitory_frequency = second_trace.stimulus(1).frequency;
            excitory_attenuation = second_trace.stimulus(1).attenuation;
        elseif two_tone_changing_idx == 1 && two_tone_changing_idx ~= 2
            excitory_frequency = third_trace.stimulus(1).frequency;
            excitory_attenuation = third_trace.stimulus(1).attenuation;
        elseif two_tone_changing_idx ~= 1 && two_tone_changing_idx ~= 2
            excitory_frequency2 = third_trace.stimulus(1).frequency;
            excitory_attenuation2 = third_trace.stimulus(1).attenuation;
        elseif two_tone_changing_idx ~= 3
            excitory_frequency2 = fourth_trace.stimulus(1).frequency;
            excitory_attenuation2 = fourth_trace.stimulus(1).attenuation;
        else
%             display(['Not enough signals found for three tone test ' int2str(test_num)]);
             WriteStatus(['Not enough signals found for three tone test ' int2str(test_num)], warnColor);
            continue;
        end
        test_label = 'Three Tone';
        suppression_duration = last_trace.stimulus(two_tone_changing_idx).delay+last_trace.stimulus(two_tone_changing_idx).duration;
%*PDR***************end inset (020408)
    else
%         display(['Test ' int2str(test_num) ' has an unsupported test type for visualization. Only single and two tone tests are supported.']);
        WriteStatus(['Test ' int2str(test_num) ' has an unsupported test type for visualization. Only single and two tone tests are supported.'], warnColor);
        continue;
    end
    trace_num_spikes = zeros(num_traces,1);
    trace_frequencies = zeros(num_traces,1);
    trace_attenuations = zeros(num_traces,1);
    trace_histograms = zeros(num_traces,num_histogram_bins);
    trace_num_sweep_spikes = cell(num_traces,1);
    for trace_num = 1:num_traces
        trace = traces(trace_num);
        if ~isempty(trace.stimulus) || strcmp(trace.is_control,'True')
%*PDR***************begin inset (022008)
            if testtype == 3
                spike_time = spike_times{test_num,trace_num};
                num_sweeps = length(spike_time);
                total_num_spikes = 0;
                total_num_sweep_spikes = zeros(1,num_sweeps);
                for sweep_num = 1:num_sweeps
                    late_spks = find(spike_time{sweep_num} > suppression_duration);
                    total_num_sweep_spikes(sweep_num) = length(late_spks);
                    total_num_spikes = total_num_spikes + length(late_spks);
                end
            else
                spike_time = spike_times{test_num,trace_num};
                num_sweeps = length(spike_time);
                total_num_spikes = 0;
                total_num_sweep_spikes = zeros(1,num_sweeps);
                for sweep_num = 1:num_sweeps
                    total_num_sweep_spikes(sweep_num) = length(spike_time{sweep_num});
                    total_num_spikes = total_num_spikes + length(spike_time{sweep_num});
                end
            end
%*PDR***************end inset (022008)
            average_num_spikes = total_num_spikes/num_sweeps;
            if strcmp(trace.is_control,'True')
                if testtype == 1
                    num_spontaneous_spikes = total_num_sweep_spikes;
                    spontaneous_firing_rate = average_num_spikes;
                    spontaneous_histogram = histograms{test_num,trace_num}/num_sweeps;
                end
            elseif (testtype == 2 || testtype == 3) && length(trace.stimulus) == 1
                if trace.stimulus(1).frequency == excitory_frequency && trace.stimulus(1).attenuation == excitory_attenuation
                    excitory_num_sweep_spikes = total_num_sweep_spikes;
                    excitory_num_spikes = average_num_spikes;
                    excitory_histogram = histograms{test_num,trace_num}/num_sweeps;
                end
            else
                if testtype == 1
                    stim_idx = 1;
                else
                    stim_idx = two_tone_changing_idx;
                end
                trace_num_sweep_spikes{trace_num} = total_num_sweep_spikes;
                trace_num_spikes(trace_num) = average_num_spikes;
                trace_frequencies(trace_num) = trace.stimulus(stim_idx).frequency;
                trace_attenuations(trace_num) = trace.stimulus(stim_idx).attenuation;
                trace_histograms(trace_num,:) = histograms{test_num,trace_num}/num_sweeps;
            end
        end
    end
    if testtype == 1
%         %Subtract spontaneous spike count from spike totals
%         trace_num_spikes = trace_num_spikes - spontaneous_firing_rate;
%         trace_histograms = trace_histograms - repmat(spontaneous_histogram, size(trace_histograms,1), 1);

        significant_differences = ones(num_traces)*0.5;
        for trace_num = 1:num_traces
            if ~isempty(trace_num_sweep_spikes{trace_num}) && exist('num_spontaneous_spikes')
                warning off
                significant_difference = ttest2(num_spontaneous_spikes,trace_num_sweep_spikes{trace_num},alpha);
                warning on
                if ~isnan(significant_difference) && significant_difference
                    if mean(num_spontaneous_spikes) > mean(trace_num_sweep_spikes{trace_num})
                        significant_differences(trace_num) = 0;
                    else
                        significant_differences(trace_num) = 1;
                    end
                end
            end
        end
    else
        if plot_percentage_increase_for_2_tone_tests
            %Subtract excitory signal's spike count from spike totals
            if excitory_num_spikes == 0
                trace_num_spikes = trace_num_spikes;
            else
                trace_num_spikes = trace_num_spikes/excitory_num_spikes;
            end
        else
            %Subtract excitory signal's spike count from spike totals
            trace_num_spikes = trace_num_spikes - excitory_num_spikes;
        end
        trace_histograms = trace_histograms - repmat(excitory_histogram, size(trace_histograms,1), 1);
        
        significant_differences = ones(num_traces)*0.5;
        for trace_num = 1:num_traces
            if ~isempty(trace_num_sweep_spikes{trace_num})
                warning off
                significant_difference = ttest2(excitory_num_sweep_spikes,trace_num_sweep_spikes{trace_num},alpha);
                warning on
                if ~isnan(significant_difference) && significant_difference
                    if mean(excitory_num_sweep_spikes) > mean(trace_num_sweep_spikes{trace_num})
                        significant_differences(trace_num) = 0;
                    else
                        significant_differences(trace_num) = 1;
                    end
                end
            end
        end
    end
    %Convert to kHz
    trace_frequencies = trace_frequencies/1000;
    unique_frequencies = unique(trace_frequencies);
    %Remove zero frequency entries
    unique_frequencies = unique_frequencies(find(unique_frequencies));
    num_frequencies = length(unique_frequencies);
    unique_attenuations = unique(trace_attenuations);
    %Remove zero attenuations entries
    unique_attenuations = unique_attenuations(find(unique_attenuations));
    num_attenuations = length(unique_attenuations);
    
    test_data = zeros(num_frequencies,num_attenuations);
    test_histogram_data = zeros(num_frequencies,num_attenuations,num_histogram_bins);
    test_stat_data = zeros(num_frequencies,num_attenuations);
    for trace_num = 1:num_traces
        trace = traces(trace_num);
        if ~isempty(trace.stimulus)
            frequency = trace_frequencies(trace_num);
            attenuation = trace_attenuations(trace_num);
            freq_idx = find(frequency == unique_frequencies);
            attn_idx = find(attenuation == unique_attenuations);
            if ~isempty(freq_idx) && ~isempty(attn_idx)
                test_data(freq_idx,attn_idx) = trace_num_spikes(trace_num);
                test_histogram_data(freq_idx,attn_idx,:) = trace_histograms(trace_num,:);
                test_stat_data(freq_idx,attn_idx) = significant_differences(trace_num);
            end
        end
    end
    
    tuneClass = 0;
    if calculate_tuning_type && length(unique_frequencies) > 1
        %Calculate the tuning class
        freqResponse = sum(test_data, 2);
        maxResponse = max(freqResponse);
        tuneIndicies = find(freqResponse > 0.3*maxResponse);
        maxTune = max(tuneIndicies);
        minTune = min(tuneIndicies);
        thresholdIndex = find(unique_frequencies>35);
        if isempty(thresholdIndex) || unique_frequencies(maxTune) < unique_frequencies(thresholdIndex(1))
            if unique_frequencies(maxTune)-unique_frequencies(minTune) < 8
                tuneClass = 1;
%                 if ~isempty(output_path)
%                     output_path = [base_output_path prefs.cell_id '_narrow-low'];
%                 end
                WriteStatus('Tuning type: narrow-low');
            else
                tuneClass = 2;
%                 output_path = [base_output_path prefs.cell_id '_broad-low'];
                WriteStatus('Tuning type: broad-low');
            end
        else
            if unique_frequencies(minTune) > unique_frequencies(thresholdIndex(1))
                tuneClass = 3;
%                 if ~isempty(output_path)
%                     output_path = [base_output_path prefs.cell_id '_high'];
%                 end
                WriteStatus('Tuning type: high')
            else
                tuneClass = 4;
%                 if ~isempty(output_path)
%                     output_path = [base_output_path prefs.cell_id '_multi'];
%                 end
                WriteStatus('Tuning type: multi')
            end
        end
    else
%         output_path = base_output_path;
    end
    
    % Store return variables
    all_unique_frequencies{test_idx} = unique_frequencies;
    all_unique_attenuations{test_idx} = unique_attenuations;
    all_test_data{test_idx} = test_data;
    all_tune_classes{test_idx} = tuneClass;
    
    
    if plot_statistical_differences
        if (testtype == 2 || testtype == 3) || (testtype == 1 && exist('num_spontaneous_spikes'))
            primary_freq_idx = find(unique_frequencies >= primary_tuning_range(1) & unique_frequencies <= primary_tuning_range(2));
            test_stat_data(primary_freq_idx,:) = 0.5;
            exication_found = false;
            inhibition_found = false;
            stat_file_name = '';%['Test ' test_num];
            if max(max(test_stat_data)) == 1
                stat_file_name = [stat_file_name ' Excitation Found'];
            end
            if min(min(test_stat_data)) == 0
                stat_file_name = [stat_file_name ' Inhibition Found'];
            end
            stat_output_path = [output_path stat_file_name];
            save(stat_output_path,'test_stat_data');
        end
    end
    
    if plot_loudness
%         unique_attenuations = flipud(100 - unique_attenuations);
        plot_xlabel = 'Intensity (dB SPL)';
    else %plot attenuations
        plot_xlabel = 'Attenuation (dB SPL)';
    end
    
    if pf(plot_spike_counts)
        figure_title = [test_label ' Tuning Curve for Test ' int2str(test_num)];
        if ( testtype == 2 || testtype == 3 )
            figure_title = [figure_title '. Excitory Tone: ' num2str(excitory_frequency/1000) ' kHz / ' num2str(excitory_attenuation) ' dB Attenuation'];
        end
        if length(unique_frequencies) == 1            
            if length(unique_attenuations) > 1
%                 disp(['Warning: only single unique frequency, cannot plot test ' num2str(test_num) ' as tuning curve']);
                WriteStatus(['Only single unique frequency, cannot plot test ' num2str(test_num) ' as tuning curve '], warnColor);
                if exist('output_path','var') && ~isempty(output_path)
                    if ~strcmp(prefs.colormap_name, 'jet')
%                         output_path = [output_path ' ' prefs.colormap_name];
                    end
                    if prefs.force_plot_visible
                        f = figure;
                    else
                        f = figure('Visible','off');
                    end
                else
                    f = figure;
                end
                plot(unique_attenuations,squeeze(test_data));
                xlim([min(unique_attenuations) max(unique_attenuations)]);
                xlabel('Attenuation (dB SPL)');
                ylabel('Mean spikes per sweep');
                set(gca,'XTick',unique_attenuations);
                %title(figure_title);
                figureName = [prefs.cell_id num2str(test_num) '_spike_count_vs_loudness'];
                set(f, 'name',  figureName);
                if exist('output_path','var') && ~isempty(output_path)
                    PrintFigure([output_path figureName],image_format,5,4,resolution,f);
                    close(f);
                end
            else
                cprintf(warnColor, ['Only single unique frequency, cannot plot test ' num2str(test_num) ' as tuning curve']);
            end
        elseif length(unique_attenuations) == 1
%             disp(['Warning: only single unique attenuation, cannot plot test ' num2str(test_num) ' as tuning curve']);
            WriteStatus(['Only single unique attenuation, cannot plot test ' num2str(test_num) ' as tuning curve'], warnColor);
            if exist('output_path','var') && ~isempty(output_path)
                if ~strcmp(prefs.colormap_name, 'jet')
%                     output_path = [output_path ' ' prefs.colormap_name];
                end
                if prefs.force_plot_visible
                    f = figure;
                else
                    f = figure('Visible','off');
                end
            else
                f = figure;
            end
            plot(unique_frequencies,test_data);
            xlim([min(unique_frequencies) max(unique_frequencies)]);
            xlabel('Frequency (kHz)');
            ylabel('Mean spikes per sweep');
            %title(figure_title);
            figureName = [prefs.cell_id num2str(test_num) '_spike_count_vs_freq'];
            set(f, 'name',  figureName);
            if exist('output_path','var') && ~isempty(output_path)
                PrintFigure([output_path figureName],image_format,5,4,resolution,f);
                close(f);
            end
        else
            [X Y] = meshgrid(unique_frequencies,unique_attenuations);
            if pf(plot_surface)
                if exist('output_path','var') && ~isempty(output_path)
                    if ~strcmp(prefs.colormap_name, 'jet')
%                         output_path = [output_path ' ' prefs.colormap_name];
                    end
                    if prefs.force_plot_visible
                        f = figure;
                    else
                        f = figure('Visible','off');
                    end
                else
                    f = figure;
                end
                surf(X,Y,test_data');
                view(2);
                shading interp
                if prefs.invert_color
                    colormap(flipud(prefs.colormap));
                else
                    colormap(prefs.colormap);
                end
                xlim([min(unique_frequencies) max(unique_frequencies)]);
                ylim([min(unique_attenuations) max(unique_attenuations)]);
                xlabel('Frequency (kHz)');
                ylabel('Intensity (dB SPL)');
                if length(unique_frequencies) <= 10
                    set(gca,'XTick',unique_frequencies);
                end
                if length(unique_attenuations) <= 10
                    set(gca,'YTick',unique_attenuations);
                    if plot_loudness
                        set(gca,'YTickLabel',num2str(100 - unique_attenuations));
                    end
                end
                %title(figure_title);
                figureName = [prefs.cell_id num2str(test_num) '_tuning_curve_surface'];
                set(f, 'name',  figureName);
                set(gca,'YDir','reverse');
                if exist('colorRange', 'var') && ~isempty(colorRange)
                    SetColorBar(colorbar_label, colorRange);
                elseif strcmp(test.testtype,'tone')
                    SetColorBar(colorbar_label,2);
                elseif plot_percentage_increase_for_2_tone_tests
                    SetColorBar(colorbar_label,0);
                else
                    SetColorBar(colorbar_label,1);
                end
                if (testtype == 2 || testtype == 3)
                    hold on
                    h = plot(excitory_frequency/1000,excitory_attenuation,'.');
                    set(h,'MarkerSize',20);
                    hold off
                end
                if exist('output_path','var') && ~isempty(output_path)
                    PrintFigure([output_path figureName],image_format,5,4,resolution,f);
                    close(f);
                end
            end
            if pf(plot_contour)
                if exist('output_path','var') && ~isempty(output_path)
                    if ~strcmp(prefs.colormap_name, 'jet')
%                         output_path = [output_path ' ' prefs.colormap_name];
                    end
                    if prefs.force_plot_visible
                        f = figure;
                    else
                        f = figure('Visible','off');
                    end
                else
                    f = figure;
                end
                contourf(X,Y,test_data');
                view(2);
                if prefs.invert_color
                    colormap(flipud(prefs.colormap));
                else
                    colormap(prefs.colormap);
                end
                xlim([min(unique_frequencies) max(unique_frequencies)]);
                ylim([min(unique_attenuations) max(unique_attenuations)]);
                xlabel('Frequency (kHz)');
                ylabel('Intensity (dB SPL)');
                set(gca,'YDir','reverse');
                if length(unique_frequencies) <= 10
                    set(gca,'XTick',unique_frequencies);
                end
                if length(unique_attenuations) <= 10
                    set(gca,'YTick',unique_attenuations);
                    if plot_loudness
                        set(gca,'YTickLabel',num2str(100 - unique_attenuations));
                    end
                end
                %title(figure_title);
                figureName = [prefs.cell_id num2str(test_num) '_tuning_curve_contour'];
                set(f, 'name',  figureName);
                if exist('colorRange', 'var') && ~isempty(colorRange)
                    SetColorBar(colorbar_label, colorRange);
                elseif strcmp(test.testtype,'tone')
                    SetColorBar(colorbar_label,2);
                elseif plot_percentage_increase_for_2_tone_tests
                    SetColorBar(colorbar_label,0);
                else
                    SetColorBar(colorbar_label,1);
                end
                if (testtype == 2 || testtype == 3)
                    hold on
                    h = plot(excitory_frequency/1000,excitory_attenuation,'.');
                    set(h,'MarkerSize',20);
                    hold off
                end
                if exist('output_path','var') && ~isempty(output_path)
                    PrintFigure([output_path figureName],image_format,5,4,resolution,f);
                    close(f);
                end
            end
            if pf(plot_imagemap)
                if exist('output_path','var') && ~isempty(output_path)
                    if ~strcmp(prefs.colormap_name, 'jet')
%                         output_path = [output_path ' ' prefs.colormap_name];
                    end
                    if prefs.force_plot_visible
                        f = figure;
                    else
                        f = figure('Visible','off');
                    end
                else
                    f = figure;
                end
%                 surf(X,Y,test_data','EdgeAlpha',0);
                imagesc(unique_frequencies,unique_attenuations,test_data');
                view(2);
                if prefs.invert_color
                    colormap(flipud(prefs.colormap));
                else
                    colormap(prefs.colormap);
                end
                xlim([min(unique_frequencies) max(unique_frequencies)]);
                ylim([min(unique_attenuations) max(unique_attenuations)]);
                xlabel('Frequency (kHz)');
                ylabel('Intensity (dB SPL)');
                if length(unique_frequencies) <= 10
                    set(gca,'XTick',unique_frequencies);
                end
                if length(unique_attenuations) <= 10
                    set(gca,'YTick',unique_attenuations);
                    if plot_loudness
                        set(gca,'YTickLabel',num2str(100 - unique_attenuations));
                    end
                end
                %title(figure_title);
                figureName = [prefs.cell_id num2str(test_num) '_tuning_curve_imagemap'];
                set(f, 'name',  figureName);
                set(gca,'YDir','reverse');
                if exist('colorRange', 'var') && ~isempty(colorRange)
                    SetColorBar(colorbar_label, colorRange);                
                elseif strcmp(test.testtype,'tone')
                    SetColorBar(colorbar_label,2);
                elseif plot_percentage_increase_for_2_tone_tests
                    SetColorBar(colorbar_label,0);
                else
                    SetColorBar(colorbar_label,1);
                end
                if testtype == 2
                    hold on
                    h = plot(excitory_frequency/1000,excitory_attenuation,'.');
                    set(h,'MarkerSize',20);
                    hold off
                end
                if exist('output_path','var') && ~isempty(output_path)
                    PrintFigure([output_path figureName],image_format,5,4,resolution,f);
                    close(f);
                end
            end
        end
    end
    
    if pf(plot_histograms)
        cm = flipud(bone);
        for attenuation = unique_attenuations'
            num_freqs = length(unique_frequencies);
            if exist('output_path','var') && ~isempty(output_path)
                hist_output_path = [output_path int2str(attenuation) 'dB'];
            end
            attn_idx = find(attenuation == unique_attenuations);
            figure_title = [test_label ' Multi-frequency Histogram for Test ' int2str(test_num) ' at ' int2str(attenuation) ' dB Attenuation'];
            if num_freqs == 1
                WriteStatus(['Only one frequency found while generating ' figure_title]);
%                 f = figure;
%                 bar(histogram_bin_centers{test_num,trace_num}, squeeze(test_histogram_data(1,1,:)),'hist');
%                 h = findobj(gca,'Type','patch');
%                 set(h,'FaceColor',0.5*[1 1 1],'EdgeColor','w');
%                 xlabel('Time (ms)');
%                 ylabel(['Mean spikes per ' num2str(prefs.histogram_bin_width,2) ' ms histogram bin']);
%                 title(figure_title);
            else
                hist_data = zeros(num_freqs,num_histogram_bins);
                count = 0;
                for frequency = unique_frequencies'
                    count = count + 1;
                    freq_idx = find(frequency == unique_frequencies);
                    attn_idx = find(attenuation == unique_attenuations);
                    hist_data(count,:) = test_histogram_data(freq_idx,attn_idx,:);
                end

                [X Y] = meshgrid(histogram_bin_centers{test_num,trace_num},unique_frequencies);
                if pf(plot_surface)
                    if exist('output_path','var') && ~isempty(output_path)
                        if ~strcmp(prefs.colormap_name, 'jet')
%                             output_path = [output_path ' ' prefs.colormap_name];
                        end
                        if prefs.force_plot_visible
                            f = figure;
                        else
                            f = figure('Visible','off');
                        end
                    else
                        f = figure;
                    end
                    surf(X,Y,hist_data);
                    view(2);
                    shading interp
                    colormap(cm);
                    ylim([min(unique_frequencies) max(unique_frequencies)]);
                    xlabel('Time (ms)');
                    ylabel('Frequency (kHz)');
                    if length(unique_frequencies) <= 10
                        set(gca,'YTick',unique_frequencies);
                    end
                    %title(figure_title);
                    figureName = [prefs.cell_id num2str(test_num) '_psth_surface_' int2str(attenuation) 'dB'];
                    set(f, 'name',  figureName);
                    if exist('colorRange', 'var') && ~isempty(colorRange)
                        SetColorBar(colorbar_label, colorRange);                    
                    elseif strcmp(test.testtype,'tone')
                        SetColorBar(colorbar_label,4);
                    else
                        SetColorBar(colorbar_label,5);
                    end
                    if exist('output_path','var') && ~isempty(output_path)
                        PrintFigure([output_path figureName],image_format,5,4,resolution,f);
                        close(f);
                    end
                end
                if pf(plot_contour)
                    if exist('output_path','var') && ~isempty(output_path)
                        if ~strcmp(prefs.colormap_name, 'jet')
%                             output_path = [output_path ' ' prefs.colormap_name];
                        end
                        if prefs.force_plot_visible
                            f = figure;
                        else
                            f = figure('Visible','off');
                        end
                    else
                        f = figure;
                    end
                    contourf(X,Y,hist_data);
                    view(2);
                    colormap(cm);
                    ylim([min(unique_frequencies) max(unique_frequencies)]);
                    xlabel('Time (ms)');
                    ylabel('Frequency (kHz)');
                    if length(unique_frequencies) <= 10
                        set(gca,'YTick',unique_frequencies);
                    end
                    %title(figure_title);
                    figureName = [prefs.cell_id num2str(test_num) '_psth_contour' int2str(attenuation) 'dB'];
                    set(f, 'name',  figureName);                    
                    if exist('colorRange', 'var') && ~isempty(colorRange)
                        SetColorBar(colorbar_label, colorRange);
                    elseif strcmp(test.testtype,'tone') 
                        SetColorBar(colorbar_label,4);
                    else
                        SetColorBar(colorbar_label,5);
                    end
                    if exist('output_path','var') && ~isempty(output_path)
                        PrintFigure([output_path figureName],image_format,5,4,resolution,f);
                        close(f);
                    end
                end
                if pf(plot_imagemap)
                    if exist('output_path','var') && ~isempty(output_path)
                        if ~strcmp(prefs.colormap_name, 'jet')
%                             output_path = [output_path ' ' prefs.colormap_name];
                        end
                        if prefs.force_plot_visible
                            f = figure;
                        else
                            f = figure('Visible','off');
                        end
                    else
                        f = figure;
                    end
                    surf(X,Y,hist_data,'EdgeAlpha',0);
                    view(2);
                    colormap(cm);
                    xlim([min(histogram_bin_centers{test_num,trace_num}) max(histogram_bin_centers{test_num,trace_num})]);
                    ylim([min(unique_frequencies) max(unique_frequencies)]);
                    xlabel('Time (ms)');
                    ylabel('Frequency (kHz)');
                    if length(unique_frequencies) <= 10
                        set(gca,'YTick',unique_frequencies);
                    end
                    %title(figure_title);
                    figureName = [prefs.cell_id num2str(test_num) '_psth_imagemap' int2str(attenuation) 'dB'];
                    set(f, 'name',  figureName);
                    if exist('colorRange', 'var') && ~isempty(colorRange)
                        SetColorBar(colorbar_label, colorRange);                    
                    elseif strcmp(test.testtype,'tone')
                        SetColorBar(colorbar_label,4);
                    else
                        SetColorBar(colorbar_label,5);
                    end
                    if exist('output_path','var') && ~isempty(output_path)
                        PrintFigure([output_path figureName],image_format,5,4,resolution,f);
                        close(f);
                    end
                end
            end
        end
    end
end

unique_frequencies = all_unique_frequencies;
unique_attenuations = all_unique_attenuations;
test_data = all_test_data;
tune_classes = all_tune_classes;

        
