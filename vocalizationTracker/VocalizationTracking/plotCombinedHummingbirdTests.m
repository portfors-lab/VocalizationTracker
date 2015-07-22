this_EKS_nmse = [];
this_spectrogram_nmse = [];
this_pink_EKS_nmse = [];
this_pink_spectrogram_nmse = [];
workspaces = dir('workspace*');
for workspaceNum = 1:length(workspaces)
    workspace = workspaces(workspaceNum);
    workspace_name = workspace.name
    load(workspace_name);
    this_EKS_nmse = [this_EKS_nmse EKS_nmse];
    this_spectrogram_nmse = [this_spectrogram_nmse spectrogram_nmse];
    this_pink_EKS_nmse = [this_pink_EKS_nmse pink_EKS_nmse];
    this_pink_spectrogram_nmse = [this_pink_spectrogram_nmse pink_spectrogram_nmse];
end

EKS_nmse = this_EKS_nmse;
spectrogram_nmse = this_spectrogram_nmse;
pink_EKS_nmse = this_pink_EKS_nmse;
pink_spectrogram_nmse = this_pink_spectrogram_nmse;
num_simulations = size(EKS_nmse,2);

% [a b] = find(spectrogram_nmse<2);
% num_bad_spectrogram_tracks = length(find(spectrogram_nmse > 2))
% EKS_nmse = EKS_nmse(:,a);
% spectrogram_nmse = spectrogram_nmse(:,a);

percent = .95;
percent = .80;

figure;
hold on
legend_handles = [];
legend_labels = {};
legend_index = 1;
if include_EKS
%     mean_EKS_nmse = mean(EKS_nmse,2);
%     stdev_EKS_nmse = std(EKS_nmse,0,2);
%     EKS_upper_confidence = mean_EKS_nmse + 1.96*stdev_EKS_nmse/sqrt(num_simulations);
%     EKS_lower_confidence = mean_EKS_nmse - 1.96*stdev_EKS_nmse/sqrt(num_simulations);
%     h1 = plot_with_confidence_intervals(10*log10(mean_signal_power./synthetic_noise),mean_EKS_nmse,EKS_upper_confidence,EKS_lower_confidence,'b');
    
    median_EKS_nmse = median(EKS_nmse,2);
    pct_99_EKS = zeros(num_params,1);
    pct_01_EKS = zeros(num_params,1);
    for i = 1:num_params
        pct_99_EKS(i,1) = Percentile(EKS_nmse(i,:),percent);
        pct_01_EKS(i,1) = Percentile(EKS_nmse(i,:),1-percent);
    end
    EKS_upper_confidence = median_EKS_nmse + pct_99_EKS;
    EKS_lower_confidence = median_EKS_nmse - pct_01_EKS;     
    h1 = plot_with_confidence_intervals(10*log10(mean_signal_power./synthetic_noise),median_EKS_nmse,EKS_upper_confidence,EKS_lower_confidence,'r');
        
    legend_handles = [legend_handles h1];
    legend_labels{legend_index} = 'EKS Tracker';
    legend_index = legend_index + 1;
end
if include_spectrogram
%     mean_spectrogram_nmse = mean(spectrogram_nmse,2);
%     stdev_spectrogram_nmse = std(spectrogram_nmse,0,2);
%     spectrogram_upper_confidence = mean_spectrogram_nmse + 1.96*stdev_spectrogram_nmse/sqrt(num_simulations);
%     spectrogram_lower_confidence = mean_spectrogram_nmse - 1.96*stdev_spectrogram_nmse/sqrt(num_simulations);
%     h2 = plot_with_confidence_intervals(10*log10(mean_signal_power./synthetic_noise),mean_spectrogram_nmse,spectrogram_upper_confidence,spectrogram_lower_confidence,'r');
    
    median_spectrogram_nmse = median(spectrogram_nmse,2);
    pct_99_EKS = zeros(num_params,1);
    pct_01_EKS = zeros(num_params,1);
    for i = 1:num_params
        pct_99_spectrogram(i,1) = Percentile(spectrogram_nmse(i,:),percent);
        pct_01_spectrogram(i,1) = Percentile(spectrogram_nmse(i,:),1-percent);
    end
    spectrogram_upper_confidence = median_spectrogram_nmse + pct_99_spectrogram;
    spectrogram_lower_confidence = median_spectrogram_nmse - pct_01_spectrogram;  
    h2 = plot_with_confidence_intervals(10*log10(mean_signal_power./synthetic_noise),median_spectrogram_nmse,spectrogram_upper_confidence,spectrogram_lower_confidence,'--b');
    
    legend_handles = [legend_handles h2];
    legend_labels{legend_index} = 'Spectrogram Tracker';
    legend_index = legend_index + 1;
end
if include_PF
    mean_PF_nmse = mean(PF_nmse,2);
    stdev_PF_nmse = std(PF_nmse,0,2);
    PF_upper_confidence = mean_PF_nmse + 1.96*stdev_PF_nmse/sqrt(num_simulations);
    PF_lower_confidence = mean_PF_nmse - 1.96*stdev_PF_nmse/sqrt(num_simulations);

    h2 = plot_with_confidence_intervals(10*log10(mean_signal_power./synthetic_noise),mean_PF_nmse,PF_upper_confidence,spectrogram_lower_confidence,'g');
    legend_handles = [legend_handles h3];
    legend_labels{legend_index} = 'PF Tracker';
    legend_index = legend_index + 1;
end

plot(10*log10(mean_signal_power./synthetic_noise),EKS_upper_confidence,'k');

xlabel('SNR');
ylabel('NMSE');
% set(gca,'XScale','log');
set(gca,'yScale','log');
legend(legend_handles,legend_labels);
% xlim([min(10*log10(mean_signal_power./synthetic_noise)) max(10*log10(mean_signal_power./synthetic_noise))]);
xlim([0 35]);

title('Fundamental Frequency Tracking Performance');
hold off
AxisSet;
ylim([10^-3.2 10^0.6])

%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------

EKS_nmse = pink_EKS_nmse;
spectrogram_nmse = pink_spectrogram_nmse;

figure;
hold on
legend_handles = [];
legend_labels = {};
legend_index = 1;
if include_EKS
%     mean_EKS_nmse = mean(EKS_nmse,2);
%     stdev_EKS_nmse = std(EKS_nmse,0,2);
%     EKS_upper_confidence = mean_EKS_nmse + 1.96*stdev_EKS_nmse/sqrt(num_simulations);
%     EKS_lower_confidence = mean_EKS_nmse - 1.96*stdev_EKS_nmse/sqrt(num_simulations);
%     h1 = plot_with_confidence_intervals(10*log10(mean_signal_power./synthetic_noise),mean_EKS_nmse,EKS_upper_confidence,EKS_lower_confidence,'b');
    
    median_EKS_nmse = median(EKS_nmse,2);
    pct_99_EKS = zeros(num_params,1);
    pct_01_EKS = zeros(num_params,1);
    for i = 1:num_params
        pct_99_EKS(i,1) = Percentile(EKS_nmse(i,:),percent);
        pct_01_EKS(i,1) = Percentile(EKS_nmse(i,:),1-percent);
    end
    EKS_upper_confidence = median_EKS_nmse + pct_99_EKS;
    EKS_lower_confidence = median_EKS_nmse - pct_01_EKS;     
    h1 = plot_with_confidence_intervals(10*log10(mean_signal_power./synthetic_noise),median_EKS_nmse,EKS_upper_confidence,EKS_lower_confidence,'r');
        
    legend_handles = [legend_handles h1];
    legend_labels{legend_index} = 'EKS Tracker';
    legend_index = legend_index + 1;
end
if include_spectrogram
%     mean_spectrogram_nmse = mean(spectrogram_nmse,2);
%     stdev_spectrogram_nmse = std(spectrogram_nmse,0,2);
%     spectrogram_upper_confidence = mean_spectrogram_nmse + 1.96*stdev_spectrogram_nmse/sqrt(num_simulations);
%     spectrogram_lower_confidence = mean_spectrogram_nmse - 1.96*stdev_spectrogram_nmse/sqrt(num_simulations);
%     h2 = plot_with_confidence_intervals(10*log10(mean_signal_power./synthetic_noise),mean_spectrogram_nmse,spectrogram_upper_confidence,spectrogram_lower_confidence,'r');
    
    median_spectrogram_nmse = median(spectrogram_nmse,2);
    pct_99_EKS = zeros(num_params,1);
    pct_01_EKS = zeros(num_params,1);
    for i = 1:num_params
        pct_99_spectrogram(i,1) = Percentile(spectrogram_nmse(i,:),percent);
        pct_01_spectrogram(i,1) = Percentile(spectrogram_nmse(i,:),1-percent);
    end
    spectrogram_upper_confidence = median_spectrogram_nmse + pct_99_spectrogram;
    spectrogram_lower_confidence = median_spectrogram_nmse - pct_01_spectrogram;  
    h2 = plot_with_confidence_intervals(10*log10(mean_signal_power./synthetic_noise),median_spectrogram_nmse,spectrogram_upper_confidence,spectrogram_lower_confidence,'--b');
    
    legend_handles = [legend_handles h2];
    legend_labels{legend_index} = 'Spectrogram Tracker';
    legend_index = legend_index + 1;
end
if include_PF
    mean_PF_nmse = mean(PF_nmse,2);
    stdev_PF_nmse = std(PF_nmse,0,2);
    PF_upper_confidence = mean_PF_nmse + 1.96*stdev_PF_nmse/sqrt(num_simulations);
    PF_lower_confidence = mean_PF_nmse - 1.96*stdev_PF_nmse/sqrt(num_simulations);

    h2 = plot_with_confidence_intervals(10*log10(mean_signal_power./synthetic_noise),mean_PF_nmse,PF_upper_confidence,spectrogram_lower_confidence,'g');
    legend_handles = [legend_handles h3];
    legend_labels{legend_index} = 'PF Tracker';
    legend_index = legend_index + 1;
end

plot(10*log10(mean_signal_power./synthetic_noise),EKS_upper_confidence,'k');

xlabel('SNR');
ylabel('NMSE');
% set(gca,'XScale','log');
set(gca,'yScale','log');
legend(legend_handles,legend_labels);
% xlim([min(10*log10(mean_signal_power./synthetic_noise)) max(10*log10(mean_signal_power./synthetic_noise))]);
xlim([0 35]);

title('Fundamental Frequency Tracking Performance');
hold off
AxisSet;
ylim([10^-3.2 10^0.6])