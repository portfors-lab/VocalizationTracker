function VisualizeModel(prefs, ...
                        models, ...
                        plot_title, ...
                        output_path, ...
                        primary_tuning_range)
%
%function VisualizeModel(model)
%
%   INPUT ARGUMENTS
%   model               The model to visualize
%   title               Optional plot title argument
%                       Default = 'Model Visualization'
%   output_path         The relative or absolute path, including the base
%                       file name, to print resulting figures to.
%                       If empty, the figure will not be saved.
%                       Default: []

if ~iscell(models)
    model = models;
    models = cell(1);
    models{1} = model;
end

num_models = length(models);

if exist('output_path','var') && ~isempty(output_path)
    if ~strcmp(prefs.colormap_name, 'jet')
        output_path = [output_path ' ' prefs.colormap_name];
    end
    if prefs.force_plot_visible
        f = figure;
    else
        f = figure('Visible','off');
    end
else
    f = figure;
end

plot_together = false;
% output_format = 'jpeg';
output_format = 'epsc2';

for model_num = 1:num_models
    model = models{model_num};

    total_num_features = size(model.parameters,1);
    num_frequency_bands = length(model.sampled_frequencies);
    num_features_per_frequency = length(model.input_idx);
    num_poles = length(model.AR_idx);

    if rem(total_num_features - num_poles, num_frequency_bands) ~= 0
        error('Total number of model parameters is not a multiple of the number of frequency intervals');
    end

    if exist('plot_title','var')
        if isempty(plot_title)
            plot_title = 'Model Visualization';
        end
    else
        plot_title = 'Model Visualization';
    end

    reshaped_model = zeros(num_frequency_bands,num_features_per_frequency);
    for frequency_band_num = 1:num_frequency_bands
        freq_idx = frequency_band_num + (0:num_frequency_bands:num_frequency_bands*(num_features_per_frequency-1));
        reshaped_model(frequency_band_num,:) = model.parameters(freq_idx);
    end


%     if num_poles > 0
%         subplot(2,1,2);
%         AR_model = model(end-num_poles+1:end);
%         stem(model.AR_idx,AR_model,'filled');
%         xlim([min(model.AR_idx)-0.5 max(model.AR_idx)+0.5]);
%         xlabel('AR Coefficients');
%         ylabel('Coefficient Amplitude');
%         subplot(2,1,1);
%     end
    time_idx = model.input_idx*(1/prefs.model_time_samples_per_millisecond);
    freq_idx = model.sampled_frequencies/1000; %Convert to Hz
    
    %This is all to be able to use imagesc instead of surf for postscript output
    resamples = 512;
    resampled_freq_idx = linspace(min(freq_idx),max(freq_idx),resamples);
    resampled_model = zeros(resamples,size(reshaped_model,2));
    for column_idx = 1:size(reshaped_model,2)
        resampled_model(:,column_idx) = interp1(freq_idx,reshaped_model(:,column_idx)',resampled_freq_idx,'nearest');
    end
    freq_idx = resampled_freq_idx;
    reshaped_model = resampled_model;

    if num_models > 1 && plot_together
        subplot(1,num_models,model_num);
    else
        if model_num > 1
            if exist('output_path','var') && ~isempty(output_path)
                g = figure('Visible','off');
%                 g = figure;
            else
                g = figure;
            end
            %Add the figure to the list for plotting
            f = [f g];
        end
    end
    if min(size(reshaped_model)) == 1
        stem(model.input_idx,reshaped_model);
        xlabel('Time (ms)');
        ylabel('Coefficient Amplitude');
        title(['Calculated Model for ' int2str(freq_idx) ' Hz Frequency Band']);
    else
%         [X Y] = meshgrid(time_idx,freq_idx);
% %         surf(X,Y,reshaped_model.*1000,'EdgeAlpha',0);
%         view(2);
        imagesc(time_idx,freq_idx,reshaped_model.*1000);
        ylim([min(freq_idx) max(freq_idx)]);
        xlim([min(time_idx) max(time_idx)]);
%         caxis([-7e-3 7e-3]);
%         caxis([-7 7]);
        caxis([-3 3]);
        colorbar;
        
        colormap(prefs.colormap)
        set(gca,'YDir','normal');   

        if num_models == 1 || ceil(num_models/2) == model_num || ~plot_together
            xlabel('Time (ms)');
%             title(plot_title);
        end
        if model_num == 1 || ~plot_together
            ylabel('Frequency (Hz)');
        else
            set(gca,'YTick',[]);
        end
%         title(model.label)
        if ~isempty(strfind(model.label,'BF'))
            lower_freq = model.sampled_frequencies(find(model.sampled_frequencies/1000 >= max(min(freq_idx),primary_tuning_range(1))))/1000;
            lower_freq = lower_freq(1);
            upper_freq = model.sampled_frequencies(find(model.sampled_frequencies/1000 <= min(max(freq_idx),primary_tuning_range(2))))/1000;
            upper_freq = upper_freq(end);
            ylim([lower_freq upper_freq]);
        end
    end
end

% if num_models > 1
%     subplot(1,num_models+1,num_models+1);
%     axis off;
% % % %     caxis([-7e-3 7e-3]);
% % % %     c = colorbar('East');
% % % %     set(c,'YTick',[]);
% % % %     set(c,'YTickLabel',[]);
% % % %     position = get(c,'OuterPosition');
% % % %     position(1) = 0.92;
% % % %     set(c,'OuterPosition',position);
% else
%     caxis([-7e-3 7e-3]);
%     colorbar
% end

if exist('output_path','var') && ~isempty(output_path)
%     printFigure(output_path,[],3*num_models,3,[],f);
    if plot_together
        printFigure(output_path,output_format,4*num_models,3,600,f);
    else
        printFigure(output_path,output_format,4,3,600,f);
    end
%     saveas(f,output_path,'fig')
end
