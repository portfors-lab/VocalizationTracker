function [mean_correlation ...
          stdev_correlation ...
          correlations ...
          spike_rates1 ...
          spike_rates2 ] = CalculateSpikeCorr(experiment_data, ...
                                              prefs, ...
                                              test_num1, ...
                                              trace_num1, ...
                                              test_num2, ...
                                              trace_num2, ...
                                              spike_rates, ...
                                              sweep_index1)
%to find Rx

find_maximum_correlation = false;
use_mean_spike_rates = false;

if  length(test_num1) ~= 1 ...
    ||  length(trace_num1) ~= 1 ...
    ||  length(test_num2) ~= 1 ...
    ||  length(trace_num2) ~= 1
    error('Test and trace numbers must be specified and must be of length 1.');
end
                                                     
if exist('spike_rates','var')
    if isempty(spike_rates)
        spike_rates = GenerateSpikeRates(experiment_data,prefs,test_num1,trace_num1);
        spike_rates2 = GenerateSpikeRates(experiment_data,prefs,test_num2,trace_num2);
        spike_rates{test_num2,trace_num2} = spike_rates2{test_num2,trace_num2};
    end
else
    spike_rates = GenerateSpikeRates(experiment_data,prefs,test_num1,trace_num1);
    spike_rates2 = GenerateSpikeRates(experiment_data,prefs,test_num2,trace_num2);
    spike_rates{test_num2,trace_num2} = spike_rates2{test_num2,trace_num2};
end
spike_rates1 = spike_rates{test_num1,trace_num1};
spike_rates2 = spike_rates{test_num2,trace_num2};

if ~exist('sweep_index1','var')
    sweep_index1 = [];
end

if use_mean_spike_rates
    spike_rates1 = mean(spike_rates1);
    spike_rates2 = mean(spike_rates2);
end

%Calculate the correlation between the two rate filtered spike trains.
correlations = [];
num_sweeps1 = size(spike_rates1,1);
num_sweeps2 = size(spike_rates2,1);

% If a comparison is desired between a single sweep and another set of spikes,
% sweep_index1 can be set
if ~isempty(sweep_index1)
    sweep_num1_range = sweep_index1;
else
    sweep_num1_range = 1:num_sweeps1; %modified by Amy Boyle from 1:num_sweeps2- caused error
end

for sweep_num1 = sweep_num1_range
    if ~isempty(sweep_index1)
        sweep_num2_range = 1:num_sweeps2;
    else
        sweep_num2_range = sweep_num1:num_sweeps2;
    end
    for sweep_num2 = sweep_num2_range
        sweep_correlation = NaN; %Initialize
        if (test_num1 == test_num2 && trace_num1 == trace_num2 && sweep_num1 == sweep_num2)
            % Correlation between same signals, so skip computation and set to 1
            sweep_correlation = NaN;
        else
            % Using abs(*) here since the spike rate resulting from RateFilter may have negative ringing.
            filtered1 = abs(spike_rates1(sweep_num1,:));
            filtered2 = abs(spike_rates2(sweep_num2,:));
            num_samples1 = length(filtered1);
            num_samples2 = length(filtered2);

            if find_maximum_correlation
                error('Check logic. Many changes have happened to code below without modifications here.');
                %If this flag is set, then we are going to compare all time offsets of the
                %signal to determine the maximum possible correlation
                comparison_filtered_signal_base = [zeros(1,num_samples1 - 1) zeros(1,num_samples2) zeros(1,num_samples1 - 1)];
                comparison_filtered_signal = comparison_filtered_signal_base;
                comparison_filtered_signal(num_samples1:num_samples1+num_samples2-1) = filtered2;

                max_sweep_correlation = 0;
                sweep_correlation = 0;
                for compare_idx = 1:num_samples1+num_samples2-1
                    sliding_filtered_signal = comparison_filtered_signal_base;
                    sliding_filtered_signal(compare_idx:compare_idx+num_samples2-1) = filtered1;
                    sweep_correlation = sum(sliding_filtered_signal.*comparison_filtered_signal);
                    sweep_normalization = sum(abs(sliding_filtered_signal).*abs(comparison_filtered_signal));
                    if sweep_normalization == 0
                        sweep_normalized_correlation = 0;
                    else
                        sweep_normalized_correlation = sweep_correlation / sweep_normalization;
                    end
                    if sweep_normalized_correlation > max_sweep_correlation
                        max_sweep_correlation = sweep_normalized_correlation;
                    end
                end
                sweep_correlation = max_sweep_correlation;
            else
                %Otherwise, just find the correlation of the signals as they are
                if num_samples1 > num_samples2
                    zero_pad = num_samples1 - num_samples2;
                    filtered2 = [filtered2 zeros(1,zero_pad)];
                elseif num_samples2 > num_samples1
                    zero_pad = num_samples2 - num_samples1;
                    filtered1 = [filtered1 zeros(1,zero_pad)];
                end
                sweep_correlation = sum(filtered1.*filtered2);
                sweep_normalization = sqrt(sum(filtered1.^2))*sqrt(sum(filtered2.^2));
                %This happens when at least one of the recordings has no spikes
                if sweep_normalization == 0
                    if (sum(filtered1) + sum(filtered2)) > 0
                        %At least one of the recording has a spike, so call it uncorrelated
                        sweep_correlation = 0;
                    else
                        %Neither recording has spikes. Again call it
                        %uncorrelated (grey area here)
                        sweep_correlation = 0;
                    end
                else
                    sweep_correlation = sweep_correlation / sweep_normalization;
                end
            end
        end
        correlations = [correlations sweep_correlation];
    end
end
% Ignore NaNs
notNanIdx = find(~isnan(correlations));
correlations = correlations(notNanIdx);
if length(correlations) == 0
    mean_correlation = 0;
    stdev_correlation = 0;
else
    mean_correlation = mean(correlations);
    %otherMethod = sum(correlations)*(2/(num_sweeps1*(num_sweeps1-1)));
    stdev_correlation = std(correlations);
end

end


    
    
    
    
