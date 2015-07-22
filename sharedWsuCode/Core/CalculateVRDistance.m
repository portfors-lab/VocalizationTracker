function [mean_distance ...
          stdev_distance ...
          distances ...
          filtered_spike_trains1 ...
          filtered_spike_trains2 ] = CalculateVRDistance(experiment_data, ...
                                                         prefs, ...
                                                         test_num1, ...
                                                         trace_num1, ...
                                                         test_num2, ...
                                                         trace_num2, ...
                                                         spike_times, ...
                                                         shuffle_spikes)

find_minimum_distance = false;
normalize_distance = true;

if  length(test_num1) ~= 1 ...
    ||  length(trace_num1) ~= 1 ...
    ||  length(test_num2) ~= 1 ...
    ||  length(trace_num2) ~= 1
    error('Test and trace numbers must be specified and must be of length 1.');
end

if ~exist('shuffle_spikes','var')
    shuffle_spikes = false;
end
                                                     
if exist('spike_times','var')
    if isempty(spike_times)
        [spike_times spike_idx] = CalculateSpikeTimes(experiment_data,prefs,test_num1,trace_num1);
        [spike_times2 spike_idx2] = CalculateSpikeTimes(experiment_data,prefs,test_num2,trace_num2);
        spike_times{test_num2,trace_num2} = spike_times2{test_num2,trace_num2};
        spike_idx{test_num2,trace_num2} = spike_idx2{test_num2,trace_num2};
    end
else
    [spike_times spike_idx] = CalculateSpikeTimes(experiment_data,prefs,test_num1,trace_num1);
    [spike_times2 spike_idx2] = CalculateSpikeTimes(experiment_data,prefs,test_num2,trace_num2);
    spike_times{test_num2,trace_num2} = spike_times2{test_num2,trace_num2};
    spike_idx{test_num2,trace_num2} = spike_idx2{test_num2,trace_num2};
end

spike_times1 = spike_times{test_num1,trace_num1};
spike_times2 = spike_times{test_num2,trace_num2};

%Check to see if at least one cell response hits threshold (50% of presentations result in spike)
threshold = 0.5;
fraction_of_spikes_present_in_each_sweep_1 = sum((sum(spike_times1,2)>0)) / size(spike_times1,1);
fraction_of_spikes_present_in_each_sweep_2 = sum((sum(spike_times2,2)>0)) / size(spike_times2,1);
if fraction_of_spikes_present_in_each_sweep_1 < threshold ...
   && fraction_of_spikes_present_in_each_sweep_2 < threshold
   mean_distance = nan;
   stdev_distance = nan;
   distances = nan;
   return;
end

%Generate the reconstituded spike trains from the spike time data
spike_trains1 = GenerateSpikeTrains(experiment_data, prefs, test_num1, trace_num1, spike_times1, shuffle_spikes);
spike_trains2 = GenerateSpikeTrains(experiment_data, prefs, test_num2, trace_num2, spike_times2, shuffle_spikes);

% Generate the kernel for smoothing.
kernel = GenerateExponentialKernel(prefs);

%Filter these using the heavyside exponential decay kernel
filtered_spike_trains1 = FilterSpikeTrains(kernel, spike_trains1, spike_times{test_num1,trace_num1});
filtered_spike_trains2 = FilterSpikeTrains(kernel, spike_trains2, spike_times{test_num2,trace_num2});

%Calculate the summed squared distance between all pairwise combinations of the two filtered spike trains
distance = 0;
num_sweeps1 = length(spike_times1);
num_sweeps2 = length(spike_times2);
for sweep_num1 = 1:num_sweeps1
    for sweep_num2 = 1:num_sweeps2
        filtered1 = filtered_spike_trains1(sweep_num1,:);
        filtered2 = filtered_spike_trains2(sweep_num2,:);
        num_samples1 = length(filtered1);
        num_samples2 = length(filtered2);
        
        if find_minimum_distance
            comparison_length = num_samples1 + num_samples2;
            zero_pad = comparison_length - num_samples1;
            filtered1 = [filtered1 zeros(1,zero_pad)];
            sweep_distance = inf;
            for compare_start = 1:(comparison_length - num_samples2 + 1)
                this_distance = sum((filtered1(compare_start:(compare_start + num_samples1 - 1)) - filtered2).^2);
                if this_distance < sweep_distance
                    sweep_distance = this_distance;
                end
            end
        else
            if num_samples1 > num_samples2
                zero_pad = num_samples1 - num_samples2;
                filtered2 = [filtered2 zeros(1,zero_pad)];
            elseif num_samples2 > num_samples1
                zero_pad = num_samples2 - num_samples1;
                filtered1 = [filtered1 zeros(1,zero_pad)];
            end
            sweep_distance = sum((filtered1 - filtered2).^2);
        end
        if normalize_distance && (sum(filtered1.^2) + sum(filtered2.^2)) ~= 0
            %Normalize the comparison by the maximum possible distance between the two filtered spike trains
            sweep_distance = sweep_distance / (sum(filtered1.^2) + sum(filtered2.^2));
        end
        distance = distance + sweep_distance;
    end
end
%Find the mean distance over all of the sweep combinations
num_combinations = num_sweeps1 * num_sweeps2;
distance = distance / num_combinations;


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

function kernel = GenerateExponentialKernel(prefs)
%Convolve the spike trains with the normal kernel
kernel_decay = prefs.exponential_decay/1000; %In seconds
kernel_length = 10 * kernel_decay; %In seconds
kernel_num_samples = ceil(kernel_length * prefs.filtered_exponential_sample_frequency);
if mod(kernel_num_samples,2) == 0
    %Make it symmetric
    kernel_num_samples = kernel_num_samples + 1;
end
kernel_time_idx = linspace(0,kernel_length,kernel_num_samples); %In seconds
kernel = exp(-kernel_time_idx/kernel_decay);
    
    
    
