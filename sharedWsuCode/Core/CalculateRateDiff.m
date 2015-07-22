function [mean_diff ...
          stdev_diffs ...
          diffs ...
          mean_spike_rates1 ...
          mean_spike_rates2] = CalculateRateDiff(experiment_data, ...
                                                 prefs, ...
                                                 test_num1, ...
                                                 trace_num1, ...
                                                 test_num2, ...
                                                 trace_num2, ...
                                                 spike_times)
                                             
%Modified Amy Boyle 10/3/2011; changed spike_times1 & 2 variables to cell
%from matrix

use_mean_spike_rates = false;
use_threshold = false;
threshold = 0.5;

if  length(test_num1) ~= 1 ...
    ||  length(trace_num1) ~= 1 ...
    ||  length(test_num2) ~= 1 ...
    ||  length(trace_num2) ~= 1
    error('Test and trace numbers must be specified and must be of length 1.');
end

if exist('spike_times','var')
    if isempty(spike_times)
        spike_times = CalculateSpikeTimes(experiment_data,prefs,test_num1,trace_num1);
        spike_times2 = CalculateSpikeTimes(experiment_data,prefs,test_num2,trace_num2);
        spike_times{test_num2,trace_num2} = spike_times2{test_num2,trace_num2};
    end
else
    spike_times = CalculateSpikeTimes(experiment_data,prefs,test_num1,trace_num1);
    spike_times2 = CalculateSpikeTimes(experiment_data,prefs,test_num2,trace_num2);
    spike_times{test_num2,trace_num2} = spike_times2{test_num2,trace_num2};
end
spike_times1 = spike_times{test_num1,trace_num1};
spike_times2 = spike_times{test_num2,trace_num2};

if ~iscell(spike_times1)
    %convert matrix to cell and remove zeros
   spike_times1 = mat2cell(spike_times1, ones(1,size(spike_times1,1)), size(spike_times1, 2));
   spike_times2 = mat2cell(spike_times2, ones(1,size(spike_times2,1)), size(spike_times2, 2));

   spike_times1 = cellfun(@(x) x(x~=0),spike_times1, 'uniformoutput', false);
   spike_times2 = cellfun(@(x) x(x~=0),spike_times2, 'uniformoutput', false);
end

%Generate the mean spike rate for each sweep. Arbitrarily, calculate the rate
%in spikes per second. The units will drop out in the normalization.
duration1 = experiment_data.test(test_num1).trace(trace_num1).record_duration; %in Ms;
spike_rates1 = zeros(1,length(spike_times1));
for sweep_idx = 1:length(spike_times1)
    spike_rates1(sweep_idx) = length(spike_times1{sweep_idx}) * 1000/duration1;
end
duration2 = experiment_data.test(test_num2).trace(trace_num2).record_duration; %in Ms;
spike_rates2 = zeros(1,length(spike_times2));
for sweep_idx = 1:length(spike_times2)
    spike_rates2(sweep_idx) = length(spike_times2{sweep_idx}) * 1000/duration2;
end

mean_spike_rates1 = mean(spike_rates1);
mean_spike_rates2 = mean(spike_rates2);

if use_threshold
    %Check to see if at least one cell response hits threshold (50% of presentations result in spike)
    threshold = 0.5;
    num_with_spikes = 0;
    for sweep_idx = 1:length(spike_times1)
        if ~isempty(spike_times1{sweep_idx})
            num_with_spikes = num_with_spikes + 1;
        end
    end
    fraction_of_spikes_present_in_each_sweep_1 = num_with_spikes / length(spike_times1);
    num_with_spikes = 0;
    for sweep_idx = 1:length(spike_times2)
        if ~isempty(spike_times2{sweep_idx})
            num_with_spikes = num_with_spikes + 1;
        end
    end
    fraction_of_spikes_present_in_each_sweep_2 = num_with_spikes / length(spike_times2);
    if fraction_of_spikes_present_in_each_sweep_1 < threshold ...
            && fraction_of_spikes_present_in_each_sweep_2 < threshold
        mean_diff = nan;
        stdev_diffs = nan;
        diffs = nan;
        mean_spike_rates1 = 0;
        mean_spike_rates2 = 0;
        return;
    end
end

if use_mean_spike_rates
    diffs = (mean_spike_rates1-mean_spike_rates2)/ (mean_spike_rates1+mean_spike_rates2);
    stdev_diffs = 0;
    mean_diff = diffs;
    return;
end

%Calculate the correlation between the two rate filtered spike trains.
diffs = [];
num_sweeps1 = length(spike_rates1);
num_sweeps2 = length(spike_rates2);
for sweep_num1 = 1:num_sweeps1
    for sweep_num2 = sweep_num1:num_sweeps2
        if (test_num1 == test_num2 && trace_num1 == trace_num2 && sweep_num1 == sweep_num2)
            % Correlation between same signals, so skip computation and set to 1
            sweep_diff = NaN;
        else
            sweep_diff = spike_rates1(sweep_num1)-spike_rates2(sweep_num2);
            sweep_normalization = spike_rates1(sweep_num1)+spike_rates2(sweep_num2);
            %This happens when niether sweep has spikes.
            if sweep_normalization == 0
                %         sweep_diff = NaN;
                sweep_diff = 0;
            else
                sweep_diff = sweep_diff / sweep_normalization;
            end
            diffs = [diffs sweep_diff];
        end
    end
end
% Ignore NaNs
notNanIdx = find(~isnan(diffs));
diffs = diffs(notNanIdx);
mean_diff = mean(diffs);
stdev_diffs = std(diffs);


    
    
    
    
