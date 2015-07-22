function [model ...
          spontaneous_rate] = CreateLMSModel(input_data, ...
                                             target_data, ...
                                             sampled_frequencies, ...
                                             num_poles)
%
%function [model
%          spontaneous_rate] = CreateLMSModel(input_data,
%                                             target_data,
%                                             sampled_frequencies,
%                                             num_poles)
%
%   INPUT ARGUMENTS
%   input_data              The (N X M) input data matrix with N input rows
%                           of M parameters each
%   target_data             The (N X 1) target data matrix with N input rows
%   sampled_frequencies     The center frequencies of each frequency band
%                           of the model.
%   num_poles               The number of poles in the model.
%
%   OUTPUT ARGUMENTS
%   model                   The trained model.
%   spontaneous_rate        The DC offset of the model output.  
%
%CreateLMSModel creates the best linear model of the data for each
%frequency band independently. 

total_num_features = size(input_data{1},2);
num_frequency_bands = length(sampled_frequencies);

if exist('num_poles','var')
    if isempty(num_poles)
        num_poles = 0;
    end
else
    num_poles = 0;
end

if rem(total_num_features - num_poles, num_frequency_bands) ~= 0
    error('Total number of model parameters is not a multiple of the number of frequency intervals');
end

freq_model_order = (total_num_features - num_poles) / num_frequency_bands;

num_data_cells = length(input_data);
spontaneous_firing_sum = 0;
total_data_rows = 0;
total_zero_rows = 0;
%Calculate the spontaneous firing rate, here defined as the spikes that
%cannot be accounted for by stimulus inputs (zero rows of the stimulus
%component input data, not the autoregressive component of the input
%data.
for data_cell_num = 1:num_data_cells
    zero_rows = find(sum(input_data{data_cell_num},2) == 0);
    spontaneous_firing_sum = spontaneous_firing_sum + sum(target_data{data_cell_num}(zero_rows));
    total_zero_rows = total_zero_rows + length(zero_rows);
    total_data_rows = total_data_rows + length(target_data{data_cell_num});
end

spontaneous_rate = spontaneous_firing_sum / total_zero_rows;
model = zeros(total_num_features - num_poles,1);
AR_model = zeros(num_poles,1);
for frequency_band_num = 1:num_frequency_bands
    display(['Creating model ' int2str(frequency_band_num) ' of ' int2str(num_frequency_bands) '. Frequency Band : ' int2str(sampled_frequencies(frequency_band_num)) ' Hz']);
    input_freq_data = [];
    AR_data = [];
    target_freq_data = [];
    freq_idx = frequency_band_num + (0:num_frequency_bands:num_frequency_bands*(freq_model_order-1));
    for data_cell_num = 1:num_data_cells
        freq_data = input_data{data_cell_num}(:,freq_idx);
        non_zero_rows = find(sum(freq_data,2));
        if num_poles > 0
            AR_data = input_data{data_cell_num}(non_zero_rows,(total_num_features-num_poles+1:total_num_features));
            input_freq_data = [input_freq_data ; input_data{data_cell_num}(non_zero_rows,freq_idx) AR_data];
        else
            input_freq_data = [input_freq_data ; input_data{data_cell_num}(non_zero_rows,freq_idx)];
        end
        target_freq_data = [target_freq_data ; target_data{data_cell_num}(non_zero_rows,1)];
    end
    target_freq_data = target_freq_data - spontaneous_rate;
    freq_model = FitData(input_freq_data,target_freq_data);
    if num_poles > 0
        %Sum all of the model parameters associated with the AR component
        %of the model for averaging below.
        AR_model = AR_model + freq_model(end-num_poles+1:end);
        freq_model = freq_model(1:end-num_poles);
    end
    model(freq_idx) = freq_model;
end

if num_poles > 0
    AR_model = AR_model / num_frequency_bands;
    model = [model ; AR_model];
end
    
    
    