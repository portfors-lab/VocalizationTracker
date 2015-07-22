function [constrained_model ...
          LMS_model] = CreateModel(experiment_data, ...
                                   prefs, ...
                                   train_data, ...
                                   input_idx, ...
                                   model_label, ...
                                   spike_times, ...
                                   spike_idxs, ...
                                   output_path, ...
                                   input_model, ...
                                   num_epochs)
%
%function [constrained_model
%          LMS_model] = CreateModel(experiment_data,
%                                   prefs,
%                                   train_data,
%                                   input_idx,
%                                   model_label,
%                                   plot_flag)
%
%   INPUT ARGUMENTS
%   experiment_data         Bat2Matlab data structure
%   prefs                   Bat2Matlab preferences
%   train_data              The numbers of the tests used to train the
%                           model.
%   input_indx              The delay indices to use in the model.
%                           Default: -10:0
%   model_label             The label for the model.
%                           Default: 'Constrained'
%   num_epochs              The number of training epochs for the
%                           gradient descent algorithm
%                           Default: 3000
%
%   OUTPUT ARGUMENTS
%   constrained_model       The model fit on the train data and trained
%                           using the constrained optimization procedure.
%   LMS_model               The model fit on the train data and trained
%                           using the independent frequency band LMS 
%                           optimization procedure.

plot_results = true;

trace_nums = [];
AR_idx = [];

if exist('model_label','var')
    if isempty(model_label)
        model_label = 'Constrained';
    end
else
    model_label = 'Constrained';
end

if exist('data_label','var')
    if isempty(data_label)
        model_label = '';
    end
else
    data_label = '';
end

if exist('plot_flag','var')
    if isempty(plot_flag)
        plot_flag = 0;
    end
else
    plot_flag = 0;
end

if exist('input_idx','var')
    if isempty(plot_flag)
        input_idx = -10:0;
    end
else
    input_idx = -10:0;
end

if exist('num_epochs','var')
    if isempty(num_epochs)
        num_epochs = 3000;
    end
else
    num_epochs = 3000;
end

if ~exist('spike_times','var')
    spike_times = [];
end

if ~exist('spike_idxs','var')
    spike_idxs = [];
end

if ~exist('input_model','var')
    input_model = [];
end

[spectrogram_intervals ...
 sampled_frequencies ...
 two_tone_excitory_info] = CalculateSpectrogramIntervals(experiment_data, ...
                                                         prefs, ...
                                                         train_data, ...
                                                         trace_nums, ...
                                                         prefs.max_interval_width);
                                                     
% spectrogram_intervals = prefs.default_intervals;
% sampled_frequencies = prefs.default_sampled_frequencies;

if exist('input_model','var') && ~isempty(input_model)
    if size(unique(two_tone_excitory_info(:,1)),1) == 0
        error('Fixing row of model that is not a two-tone model');
    elseif size(unique(two_tone_excitory_info(:,1)),1) > 1
        error('Multiple excitation frequencies present in two-tone tests');
    end
    freq_idx = find(two_tone_excitory_info(1) == input_model.sampled_frequencies);
    if isempty(freq_idx)
        error('Two-tone excitory tone not found in input model');
    end
    num_frequency_bands = length(input_model.sampled_frequencies);
    num_features_per_frequency = length(input_model.input_idx);
    fix_idx = freq_idx + (0:num_frequency_bands:num_frequency_bands*(num_features_per_frequency-1));
    spectrogram_intervals = input_model.spectrogram_intervals;
    sampled_frequencies = input_model.sampled_frequencies;
else
    fix_idx = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate test data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_mode = 1;
[input_data target_data] = GenerateModelData(experiment_data, ...
                                             prefs, ...
                                             input_idx, ...
                                             AR_idx, ...
                                             train_data, ...
                                             trace_nums, ...
                                             model_mode, ...
                                             spectrogram_intervals, ...
                                             spike_times, ...
                                             spike_idxs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate Filterbank LMS model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
if isempty(input_model)
    [LMS_model spontaneous_rate] = CreateLMSModel(input_data,target_data,sampled_frequencies,length(AR_idx));

    tmpModel = real(LMS_model);
    tmpModel(isnan(tmpModel)) = 0;
    LMS_model = [];
    LMS_model.parameters = tmpModel;
    LMS_model.spontaneous_rate = spontaneous_rate;
    LMS_model.input_idx = input_idx;
    LMS_model.AR_idx = AR_idx;
    LMS_model.spectrogram_intervals = spectrogram_intervals;
    LMS_model.sampled_frequencies = sampled_frequencies;
    LMS_model.label = 'LMS Filterbank';
else
    LMS_model = input_model;
    spontaneous_rate = LMS_model.spontaneous_rate;
end
   
%Compare the difference between the model's output and the actual output
output = ProcessModel(LMS_model, input_data);
target_data_4_plot = Cell2Array(target_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Optimize Filterbank LMS model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
if isempty(input_model)
    learning_rate = [   3000     1000    1000 ; ...
                        0.001   0.0005  0.0001];
%     learning_rate = [   3     1    1 ; ...
%                         0.001   0.0005  0.0001];
else
    learning_rate = [   1000    1000 ; ...
                       0.0005  0.0001];
end

[constrained_model constrained_spontaneous_rate errors] = ConstrainedFit(input_data, ...
                                                                         target_data, ...
                                                                         spontaneous_rate, ...
                                                                         LMS_model.parameters, ...
                                                                         num_epochs, ...
                                                                         learning_rate, ...
                                                                         fix_idx);

tmpModel = constrained_model;
constrained_model = [];
constrained_model.parameters = tmpModel;
constrained_model.spontaneous_rate = spontaneous_rate;
constrained_model.input_idx = input_idx;
constrained_model.AR_idx = AR_idx;
constrained_model.spectrogram_intervals = spectrogram_intervals;
constrained_model.sampled_frequencies = sampled_frequencies;
constrained_model.label = model_label;

constrained_output = ProcessModel(constrained_model, input_data);

plot_flag = true;
if plot_flag
    if ~isempty(errors)
        %Plot the training results
        if exist('output_path','var') && ~isempty(output_path)
            f1 = figure('Visible','off');
        else
            f1 = figure;
        end
        semilogy(errors);
        xlabel('Training Epoch');
        ylabel('MSE');
        title(['Results of gradient descent optimization for ' model_label ' for cell ' prefs.cell_id4_plot]);
    end
    
    %Compare the difference between the constrained model's output and the
    %actual output
    if exist('output_path','var') && ~isempty(output_path)
        f2 = figure('Visible','off');
    else
        f2 = figure;
    end
    hold on
    plot(target_data_4_plot,'b');
    plot(output,'r');
    plot(constrained_output,'k');
    hold off
    title([model_label ' Model Performance on ' model_label ' Train Data for cell ' prefs.cell_id4_plot]);
    LMS_mse = mean((target_data_4_plot - max(output,0)).^2)/mean(target_data_4_plot.^2);
    Constrained_mse = mean((target_data_4_plot - max(constrained_output,0)).^2)/mean(target_data_4_plot.^2);

    legend('Physiological Data',...
           ['Filterbank LMS Model Output. MSE = ' num2str(LMS_mse)],...
           ['Constrained Model Output. MSE = ' num2str(Constrained_mse)]);
       
    if exist('output_path','var') && ~isempty(output_path)
        printFigure([output_path prefs.cell_id ' ' model_label ' Training'],[],10,8,[],f1);
        printFigure([output_path prefs.cell_id ' ' model_label ' MSE on Training Data'],[],10,8,[],f2);
        VisualizeModel(prefs, constrained_model,[model_label ' Model'], [output_path prefs.cell_id ' ' model_label]);
    else
        VisualizeModel(prefs, constrained_model,[model_label ' Model']);
    end
end



