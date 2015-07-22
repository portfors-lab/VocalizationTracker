function [model_errors ...
          model_output, axesHandles]= VisualizeTracePredictions(experiment_data, ...
                                                   prefs, ...
                                                   models, ...
                                                   test_num, ...
                                                   trace_num, ...
                                                   output_path, ...
                                                   spike_times, ...
                                                   spike_idxs)
%
%function model_errors = VisualizeTracePredictions(experiment_data, ...
%                                                  prefs, ...
%                                                  models, ...
%                                                  test_num, ...
%                                                  trace_num, ...
%                                                  output_path)
%
%   INPUT ARGUMENTS
%   experiment_data     Bat2Matlab data structure
%   prefs               Bat2Matlab preferences
%   models              A cell vector containing Bat2Matlab model
%                       structures
%   test_num            The number of the Batlab test to test the 
%                       model on
%   trace_num           The number of the Batlab trace to test the 
%                       model on
%   output_path         The relative or absolute path, including the base
%                       file name, to print resulting figures to.
%                       If empty, the figure will not be saved.
%                       Default: []
%
%   OUTPUT ARGUMENTS
%   model_errors        The errors for each model on the test data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate model input for plotting purposes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Put model into cell array if only one passed in
if ~iscell(models)
    model = models;
    clear models;
    models = cell(1,1);
    models{1} = model;
end

if ~exist('spike_times','var')
    spike_times = [];
end

if ~exist('spike_idxs','var')
    spike_idxs = [];
end

spectrogram_plot_type = 1; %NOT IMPLEMENTED. 0 : input data , 1 : full spectrogram 
visualize_data = true;

if exist('output_path','var') && ~isempty(output_path) && ~strcmp(prefs.colormap_name, 'jet')
	output_path = [output_path ' ' prefs.colormap_name];
end

model_mode = 1;
%Plot the data input associated with the last model in the list,
%somewhat arbitrarily.
model = models{end};
input_data_4_plot = GenerateModelData(experiment_data, ...
                                         prefs, ...
                                         0, ...
                                         [], ...
                                         test_num, ...
                                         trace_num, ...
                                         model_mode, ...
                                         model.spectrogram_intervals, ...
                                         spike_times, ...
                                         spike_idxs);
                                     
[input_data_4_prediction ...
 output_data_4_plot] = GenerateModelData(experiment_data, ...
                                               prefs, ...
                                               model.input_idx, ...
                                               model.AR_idx, ...
                                               test_num, ...
                                               trace_num, ...
                                               model_mode, ...
                                               model.spectrogram_intervals, ...
                                               spike_times, ...
                                               spike_idxs);
                                     
output_data_4_plot = Cell2Array(output_data_4_plot);

num_models = size(models,2);
model_output = zeros(num_models,length(output_data_4_plot));
model_errors = zeros(num_models,1);
%Iterate over the models in the model list
for model_num = 1:num_models
    model = models{model_num};                                     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Generate model input for each model
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     input_data = GenerateModelData(experiment_data, ...
%                                    prefs, ...
%                                    model.input_idx, ...
%                                    model.AR_idx, ...
%                                    test_num, ...
%                                    trace_num, ...
%                                    model_mode, ...
%                                    model.spectrogram_intervals, ...
%                                    spike_times, ...
%                                    spike_idxs);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Generate the output for the constrained one tone model contribution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     output = ProcessPoleZeroModel(model, input_data);  
    output = ProcessModel(model, input_data_4_prediction);  
    model_output(model_num,:) = output;
    model_errors(model_num) = mean((output_data_4_plot - max(output,0)).^2)/mean(output_data_4_plot.^2);
end

if visualize_data
    prefs.histogram_bin_width = 1/prefs.model_time_samples_per_millisecond;
    axesHandles = VisualizeTraceData(experiment_data, ...
                                       prefs, ...
                                       test_num, ...
                                       trace_num, ...
                                       [0 0 0 1 0 1 0 1 0 0]);


    time_idx = linspace(1,size(output_data_4_plot,1) / prefs.model_time_samples_per_millisecond,size(output_data_4_plot,1));

    subplot(2,1,2)
    hold all
%     axes(axesHandles);
    set(gca,'YLimMode','auto');
    % plot(time_idx,output_data_4_plot);
    for model_num = 1:num_models
        plot(time_idx,model_output(model_num,:));
    end
    % plot(time_idx,target_data{1},'Color',[0.6 0.6 0.6]);
    hold off
    xlabel('Time (ms)');
    ylabel(['Spike Rate (mean spikes / ' num2str(prefs.histogram_bin_width,2) ' ms bin)']);
    if strcmp(experiment_data.test(test_num).testtype,'twotone')
        legend_str = {['PSTH (mean spikes / ' num2str(prefs.histogram_bin_width,2) ' ms bin)'],'Stimulus 1 Waveform','Stimulus 2 Waveform'};
    else
        legend_str = {['PSTH (mean spikes / ' num2str(prefs.histogram_bin_width,2) ' ms bin)'],'Stimulus Waveform'};
    end
    entry_num = length(legend_str)+1;
    for model_num = 1:num_models
        model = models{model_num};
        legend_str{entry_num} = [model.label ' Model Prediction. NMSE: ' num2str(model_errors(model_num),3)];
        entry_num = entry_num + 1;
    end
    legend(legend_str,'Location','NorthEast');

%*PDR*    legend(legend_str,'Location','SouthEast');

    if spectrogram_plot_type == 0
        subplot(2,1,2)
        [X Y] = meshgrid(time_idx,model.sampled_frequencies);
        data4plot = full(Cell2Array(input_data_4_plot)');
        surf(X,Y,data4plot,'EdgeAlpha',0);
        view(2);
        xlim([min(time_idx) max(time_idx)]);
        ylim([min(model.sampled_frequencies) max(model.sampled_frequencies)]);
        % imagesc(time_idx,model.sampled_frequencies,input_data{1}');
        % set(gca,'YDir','normal');
        h = colorbar('East');    
        colormap(jet);
        set(h,'YColor',[1 1 1]);
        set(h,'XColor',[1 1 1]);
        title('Model Input Data');
        xlabel('Time (ms)');
        ylabel('Frequency (Hz)');
    end
end


if exist('output_path','var') && ~isempty(output_path)
    %Append file extensions because the decimal in the error makes it so
    %that printFigure() doesn't add the extension.
%     output_variance = mean(var(model_output));
%     jpg_output_path = [output_path '_' num2str(output_variance,3) '.jpg'];
%     fig_output_path = [output_path '_' num2str(output_variance,3) '.fig'];
    output_jpg = true;
    output_tiff = false;
    output_epsc = false;
    output_fig = false;
    save_data = true;
    if visualize_data
        if output_jpg
            jpg_output_path = [output_path '.jpg'];
            printFigure(jpg_output_path, 'jpeg', 11, 9, 150);
%             PrintFigure(jpg_output_path, 'jpeg', 6, 5);
        end
        if output_epsc
            epsc_output_path = [output_path '.eps'];
%             printFigure(epsc_output_path, 'epsc', 11, 9);
%             PrintFigure(epsc_output_path, 'epsc', 6, 5);
        end
        if output_tiff
            tiff_output_path = [output_path '.tiff'];
            PrintFigure(tiff_output_path, 'tiff', 11, 9);
        end
        if output_fig
            fig_output_path = [output_path '.fig'];
            saveas(gca,fig_output_path,'fig');
        end
    end
    if save_data
        save(output_path, 'model_output', 'model_errors','output_data_4_plot');
    end
end





