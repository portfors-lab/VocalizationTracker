function [input_data ...
          target_data] = GenerateModelData(experiment_data, ...
                                           prefs, ...
                                           input_idx, ...
                                           AR_idx, ...
                                           test_nums, ...
                                           trace_nums, ...
                                           model_mode, ...
                                           spectrogram_intervals, ...
                                           spike_times, ...
                                           spike_idxs)
%
%function [input_data
%          target_data] = GenerateModelData(experiment_data,
%                                           input_idx,
%                                           test_nums,
%                                           trace_nums,
%                                           model_mode,
%                                           spectrogram_intervals)
%
%   INPUT ARGUMENTS
%   experiment_data     Bat2Matlab data structure
%   prefs               Bat2Matlab preferences
%   input_idx           The input indices for
%   trace_nums          The numbers of the traces to calculate spike times
%                       for.
%                       Default: all traces
%   bin_width           The width of the histogram bins (in milliseconds)
%                       Default = 1
%   spike_times         Pre-calculated spike rates, stored in a cell array 
%                       indexed by {test_num,trace_num}. If empty, the
%                       spike rates are extracted on the fly.
%                       Default = []
%
%   OUTPUT ARGUMENTS
%   input_data          The calculated input data matrix, in the
%                       Bat2Matlab data cell format
%   output_data         The calculated output data vector, in the
%                       Bat2Matlab data cell format

if ~exist('sampled_frequencies','var')
    sampled_frequencies = [];
end
if ~exist('two_tone_excitory_info','var')
    two_tone_excitory_info = [];
end
if ~exist('two_tone_excitory_level','var')
    two_tone_excitory_level = [];
end
if ~exist('excitory_frequency_band_idx','var')
    excitory_frequency_band_idx = [];
end
if exist('spike_times','var') && ~isempty(spike_times) && exist('spike_idxs','var') && ~isempty(spike_idxs)
    calculate_spike_times = false;
else
    calculate_spike_times = true;
end

force_binary_spike_train = false;
include_bias = false;
num_data_rows_per_cell = prefs.model_num_data_rows_per_cell;

%Set collection parameters based on model mode
%0 : Collect data for model training
%1 : Collect data for model testing
if model_mode == 0
    combine_sweep_data = false;
elseif model_mode == 1
    combine_sweep_data = true;
else
    error('Unknown modeling mode specified in GenerateModelData()');
end


%Specify model input feature type
%0 : Spectrogram
%1 : Speech
model_input_type = 0;

%Specify model output feature type
%0 : Spike rate
%1 : Binary spike train
model_output_type = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% PREALLOCATE MODEL DATA VARIABLES %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Count total number of input data rows for preallocating variables
rms_pressure_power = cell(max(test_nums),100);
total_input_rows = 0;
%Use the specified model time resolution when calculating the spectrograms
prefs.spectrogram_time_samples_per_millisecond = prefs.model_time_samples_per_millisecond;
prefs.histogram_bin_width = 1/prefs.model_time_samples_per_millisecond;
for test_num = test_nums
    test = experiment_data.test(test_num);

    traces = test.trace;
    if isempty(trace_nums)
        traces_2_process = 1:size(traces,2);
    else
        traces_2_process = trace_nums;
    end
    for trace_num = traces_2_process
        trace = traces(trace_num);
        cache_filepath = [prefs.cache_dir filesep 'Test ' int2str(test_num) ' Trace ' int2str(trace_num) '.mat'];
        if exist(cache_filepath,'file')
            display(['Reading spectrogram for Test ' int2str(test_num) ', Trace ' int2str(trace_num) ' out of the cache']);
            load(cache_filepath);
        else
            display(['Calculating spectrogram for test ' int2str(test_num) ', trace ' int2str(trace_num)]);
            spectrograms = GenerateSpectrograms(experiment_data,prefs,test_num,trace_num);
            spectrogram = spectrograms{test_num,trace_num};
            save(cache_filepath,'spectrogram');
        end
        if model_input_type == 0
            rms_pressure_power{test_num,trace_num} = IntegrateSpectralDensity(spectrogram,prefs.spectrogram_range,spectrogram_intervals,prefs.model_spectral_integration);
        else
            error('Handle this scenario');
        end
        if ~combine_sweep_data && isfield(trace,'num_samples') && ~isempty(trace.num_samples)
            total_input_rows = total_input_rows + size(spectrogram,2) * trace.num_samples;
        else
            total_input_rows = total_input_rows + size(spectrogram,2);
        end
    end
end
if model_input_type == 0
    input_block_rows = size(spectrogram_intervals,1);
    input_block_columns = length(input_idx);
    total_num_features = (size(spectrogram_intervals,1))*input_block_columns+length(AR_idx);
    if include_bias
        total_num_features = total_num_features + 1;
    end
else
    error('Unsupported input type for preallocation.');
end
num_data_cells = ceil(total_input_rows / num_data_rows_per_cell);
input_data = cell(num_data_cells,1);
target_data = cell(num_data_cells,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Generate spike times %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if calculate_spike_times
    display('Generating Spike Times');
    [spike_times ...
     spike_idxs] = CalculateSpikeTimes(experiment_data, ...
                                       prefs, ...
                                       test_nums, ...
                                       trace_nums);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Generate histograms %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Generating Histograms');
bin_width = 1/prefs.model_time_samples_per_millisecond;
[histograms ...
 sweep_histograms ...
 histogram_bin_widths ...
 histogram_bin_centers] = GenerateHistograms(experiment_data, ...
                                             prefs, ...
                                             test_nums, ...
                                             trace_nums, ...
                                             spike_times);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Generate spike rates, if necessary %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampling_frequency = prefs.model_time_samples_per_millisecond * 1000;
if model_output_type == 0
    display('Generating Spike Rates');
    [spike_rates, ...
     spike_rate_averages, ...
     spike_rate_sampling_frequencies] = GenerateSpikeRates(experiment_data, ...
                                                           prefs, ...
                                                           test_nums, ...
                                                           trace_nums, ...
                                                           spike_idxs);
end
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Collect/Generate Model Data %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_row_idx = 0;
for test_num = test_nums
    test = experiment_data.test(test_num);
    traces = test.trace;
    if isempty(trace_nums)
        traces_2_process = 1:size(traces,2);
    else
        traces_2_process = trace_nums;
    end
    
    for trace_num = traces_2_process
        trace = traces(trace_num);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%% LOAD INPUT DATA %%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if model_input_type == 0
            warning off %Suppress log(0) warnings
            dB_SPL = 10*log10(rms_pressure_power{test_num,trace_num}/prefs.dB_SPL_ref);
            warning on
            dB_SPL = max(dB_SPL,prefs.dbMin);
            dB_SPL(find(dB_SPL < prefs.model_data_dbMin)) = 0;
            input_features = dB_SPL;
            [num_features num_time_samples] = size(input_features);
        elseif model_input_type == 1
            peak_vals = stimulus_spectrogram_harmonic_vals{test_num,trace_num};
            peak_indices = stimulus_spectrogram_harmonic_idx{test_num,trace_num};
            input_features = [peak_vals ; peak_indices];
            [num_features num_time_samples] = size(peak_vals);
        else
            error('Unkown model input feature type');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%% LOAD TARGET DATA %%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if model_output_type == 0 %Spike Rates
            if combine_sweep_data
                target_features = spike_rate_averages{test_num,trace_num};
            else
                target_features = spike_rates{test_num,trace_num};
            end
        elseif model_output_type == 1 %Binary Spike Trains
            trace_sweep_histograms = sweep_histograms{test_num,trace_num};
            if force_binary_spike_train && max(max(trace_sweep_histograms)) > 1
                warning('Multiple spikes within time bin for binary spike train calculation. Forcing binary spike train.');
                trace_sweep_histograms(find(trace_sweep_histograms > 1)) = 1;
            end
            num_sweeps = size(trace_sweep_histograms,1);
            if combine_sweep_data
                %Turn into histogram
                target_features = mean(trace_sweep_histograms);
            else
                target_features = trace_sweep_histograms;
            end
        else
            error('Unkown model input feature type');
        end

        for time_sample = 1:num_time_samples
            %Both input_idx and AR_idx correspond to indices that end
            %up in the model's input matrices. input_idx corresponds to
            %zeros and AR_idx correspond to poles for a pole-zero
            %model.
            candidate_input_idx = time_sample + input_idx;
            candidate_AR_idx = time_sample + AR_idx;
            candidate_input = zeros(input_block_rows,input_block_columns);
            in_bound_idx = find(candidate_input_idx > 0 & candidate_input_idx <= num_time_samples);
            AR_bound_idx = find(candidate_AR_idx > 0 & candidate_AR_idx <= num_time_samples);
            candidate_input(:,in_bound_idx) = input_features(:,candidate_input_idx(in_bound_idx));
            candidate_input = candidate_input(:)';
            if include_bias
                candidate_input = [candidate_input 1];
            end
            
            for sweep_num = 1:size(target_features,1)
                candidate_target = target_features(sweep_num,time_sample);
                data_row_idx = data_row_idx + 1;
                data_cell_num = ceil(data_row_idx / num_data_rows_per_cell);
                data_cell_row = rem(data_row_idx,num_data_rows_per_cell);
                if data_cell_row == 0
                    data_cell_row = num_data_rows_per_cell;
                end
                if data_cell_row == 1
                    display(['Generating data cell ' int2str(data_cell_num) ' of ' int2str(num_data_cells)]);
                    if data_cell_num > 1
                        input_data{data_cell_num-1} = sparse(input_data{data_cell_num-1});
                        target_data{data_cell_num-1} = sparse(target_data{data_cell_num-1});
                    end
                    if data_cell_num == num_data_cells
                        if rem(total_input_rows,(num_data_rows_per_cell)) == 0
                            %Total rows is exact multiple of rows per data cell
                            input_data{data_cell_num} = zeros(num_data_rows_per_cell,total_num_features);
                            target_data{data_cell_num} = zeros(num_data_rows_per_cell,1);
                        else
                            input_data{data_cell_num} = zeros(rem(total_input_rows,(num_data_rows_per_cell)),total_num_features);
                            target_data{data_cell_num} = zeros(rem(total_input_rows,(num_data_rows_per_cell)),1);
                        end
                    else
                        input_data{data_cell_num} = zeros(num_data_rows_per_cell,total_num_features);
                        target_data{data_cell_num} = zeros(num_data_rows_per_cell,1);
                    end
                end
                if include_bias
                    candidate_input(end) = 1;
                end
                candidate_AR_input = zeros(1,length(AR_idx));
                candidate_AR_input(AR_bound_idx) = target_features(sweep_num,candidate_AR_idx(AR_bound_idx));
                if size(input_data{data_cell_num},1) < data_cell_row
                    error('Error in preallocation code: Not enough data rows');
                elseif size(input_data{data_cell_num},2) ~= length(candidate_input)+ length(candidate_AR_input)
                    error('Error in preallocation code: Not enough feature columns');
                end
                input_data{data_cell_num}(data_cell_row,:) = [candidate_input candidate_AR_input];
                target_data{data_cell_num}(data_cell_row,:) = candidate_target;
            end
        end
    end
end

input_data{data_cell_num} = sparse(input_data{data_cell_num});
target_data{data_cell_num} = sparse(target_data{data_cell_num});

if data_cell_row ~= size(input_data{num_data_cells},1) ||...
   data_cell_num ~= num_data_cells
   error('Index counting error');
end

