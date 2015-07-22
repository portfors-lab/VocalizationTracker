function diff = TuningCurveDifference(experiment_data, prefs, test1, test2)
% element-for-element difference of two tuning curves

%Amy Boyle Jan 2011
scale_colorbar = true;

    %get tuning curves
    [freqs, attenuations, data, classes] = VisualizeTestData(experiment_data, prefs, [test1 test2], [1 0 0 1 0], [], [], []);
    if scale_colorbar
        close(gcf)
        close(gcf)
        max_spikes = max(max([data{1} data{2}]));
        colorRange = [0 max_spikes];
    else
        colorRange =[];
    end
    
    [freqs, attenuations, data, classes] = VisualizeTestData(experiment_data, prefs, [test1 test2], [1 0 0 1 0], [], [], [],colorRange);
    %find difference
    diff = abs(data{1} - data{2});
%     diff2 = data{1} - data{2};
    %plot difference
    plot_loudness = 1;
    testtype = 'diff';
    colorbar_label = 'Mean spikes per presentation difference';
    
    unique_frequencies = freqs{1};
    unique_attenuations = attenuations{1};
    test_data = diff;
    [X Y] = meshgrid(unique_frequencies,unique_attenuations);
    f = figure;
    contourf(X,Y,test_data');
    view(2);
    colormap(prefs.colormap);
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
    
    figureName = [prefs.cell_id num2str(test1) '_' num2str(test2) '_tuning_curve_difference'];
    set(f, 'name',  figureName);
    if exist('colorRange', 'var') && ~isempty(colorRange)
        SetColorBar(colorbar_label, colorRange);
    elseif strcmp(testtype,'diff')
        SetColorBar(colorbar_label,2);
    else
        SetColorBar(colorbar_label,1);
    end