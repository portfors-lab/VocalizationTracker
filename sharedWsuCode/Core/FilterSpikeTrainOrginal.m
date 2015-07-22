function [out_filtered_spikes, ...
          out_filtered_spikes_averages] = FilterSpikeTrain(experiment_data, ...
                                                           prefs, ...
                                                           test_nums, ...
                                                           trace_nums, ...
                                                           spike_times, ...
                                                           shuffle_spikes)
%
%function [out_spike_rates,
%          out_spike_rate_averages,
%          out_spike_rate_sampling_frequencies] = GenerateSpikeRates(experiment_data,
%                                                                    test_nums,
%                                                                    trace_nums,
%                                                                    sampling_frequency,
%                                                                    spike_times)
%
%   INPUT ARGUMENTS
%   experiment_data         Bat2Matlab data structure
%   prefs                   Bat2Matlab preferences
%   test_nums               The numbers of the tests to calculate spike rates
%                           for.
%   trace_nums              The numbers of the traces to calculate spike rates
%                           for.
%                           Default: all traces
%   sampling_frequency      Sampling frequency for rate calculation
%                           Default: Sampling rate of MER recording
%   spike_times             Pre-calculated spike rates, stored in a cell array 
%                           indexed by {test_num,trace_num}. If empty, the
%                           spike rates are extracted on the fly.
%                           Default = []
%
%   OUTPUT ARGUMENTS      
%   out_spike_rates                     The calculated spike times returned in a cell array 
%                                       indexed by {test_num,trace_num}
%   out_spike_rate_averages             The calculated spike rate averagest returned in a cell 
%                                       array indexed by {test_num,trace_num}
%   out_spike_rate_sampling_frequencies The calculated spike rate sampling frequencies 
%                                       returned in a cell array indexed by {test_num,trace_num}

if ~exist('trace_nums','var')
    trace_nums = [];
end

if ~exist('sampling_frequency','var')
    sampling_frequency = [];
end

if ~exist('shuffle_spikes','var')
    shuffle_spikes = false;
end

if exist('spike_times','var')
    if isempty(spike_times)
        spike_times = CalculateSpikeTimes(experiment_data, ...
                                          prefs, ...
                                          test_nums, ...
                                          trace_nums);
    end
else
    spike_times = CalculateSpikeTimes(experiment_data, ...
                                      prefs, ...
                                      test_nums, ...
                                      trace_nums);
end

% Generate the kernel for smoothing.
kernel = GenerateNormalKernel(prefs);

out_filtered_spikes = cell(max(test_nums),max(trace_nums));
out_filtered_spikes_averages = cell(max(test_nums),max(trace_nums));

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

        if isnan(prefs.spike_rate_cutoff_frequency)
            sampling_frequency = trace.samplerate_ad;
        else
            sampling_frequency = prefs.spike_rate_sampling_frequency;
        end
        
        %Generate the reconstituded spike trains from the spike time data
        spike_trains = GenerateSpikeTrains(experiment_data, prefs, test_num, trace_num, spike_times{test_num,trace_num}, shuffle_spikes);

        %Filter these using the gaussian kernel
        filtered_spike_trains = FilterSpikeTrains(kernel, spike_trains, spike_times{test_num,trace_num});
        
        out_filtered_spikes{test_num,trace_num} = filtered_spike_trains;
        out_filtered_spikes_averages{test_num,trace_num} =  mean(filtered_spike_trains);
    end
end

function spike_trains = GenerateSpikeTrains(experiment_data, prefs, test_num, trace_num, these_spike_times, shuffle_spikes)
%Reconstitute the spike trains as signals from the spike times
num_sweeps = length(these_spike_times);
duration = experiment_data.test(test_num).trace(trace_num).record_duration; %in Ms
spike_trains = zeros(num_sweeps,duration*prefs.filtered_gaussian_sample_frequency/1000);
for sweep_num = 1:num_sweeps
    if shuffle_spikes
        num_spikes = length(these_spike_times{sweep_num});
        sweep_spikes = rand(1,num_spikes)*duration;
    else
        sweep_spikes = these_spike_times{sweep_num};
    end
    num_spikes = length(sweep_spikes);
    for spike_num = 1:num_spikes
        spike_time = sweep_spikes(spike_num);
        spike_idx = round(spike_time * prefs.filtered_gaussian_sample_frequency/1000)+1;
        spike_trains(sweep_num, spike_idx) = 1;
    end
end

function filtered_spike_trains = FilterSpikeTrains(kernel, spike_trains, these_spike_times)
[num_sweeps num_samples] = size(spike_trains);
kernel_num_samples = length(kernel);
filtered_spike_trains = zeros(num_sweeps, num_samples + kernel_num_samples - 1);
for sweep_num = 1:num_sweeps
    if  ~isempty(these_spike_times{sweep_num})
        filtered_spike_trains(sweep_num,:) = conv(kernel,spike_trains(sweep_num,:));
    end
end

function kernel = GenerateNormalKernel(prefs)
%Convolve the spike trains with the normal kernel
kernel_stdev = prefs.filtered_gaussian_stdev/1000; %In seconds
kernel_length = 10 * kernel_stdev; %In seconds
kernel_num_samples = ceil(kernel_length * prefs.filtered_gaussian_sample_frequency);
if mod(kernel_num_samples,2) == 0
    kernel_num_samples = kernel_num_samples + 1;
end
kernel_time_idx = linspace(-kernel_length/2,kernel_length/2,kernel_num_samples); %In seconds
kernel = exp(-kernel_time_idx.^2/(2*kernel_stdev^2));

