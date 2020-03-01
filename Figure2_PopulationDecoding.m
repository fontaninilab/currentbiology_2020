% add the path to the NDT so add_ndt_paths_and_init_rand_generator can be called
clear
toolbox_basedir_name = 'D:\MATLAB\Ephys-analysis\ndt.1.0.4\'
addpath(toolbox_basedir_name);
 
% add the NDT paths using add_ndt_paths_and_init_rand_generator
add_ndt_paths_and_init_rand_generator

%% start population decoding analysis using the NDT
mkdir decoding_results % create a folder to store data

%
specific_binned_labels_names = 'stimulus_ID';  % use object identity labels to decode which object was shown
num_cv_splits = 10;     % use cross-validation splits;
binned_data_file_name = 'all_100ms_bins_100ms_sampled30repeats.mat'; % use the data that was previously binned

% the name of where to save the results
save_file_name = 'decoding_results/population_decoding_results';

% create the basic objects needed for decoding
ds = basic_DS(binned_data_file_name, specific_binned_labels_names,  num_cv_splits); % create the basic datasource object
the_feature_preprocessors{1} = zscore_normalize_FP;  % create a feature preprocess that z-score normalizes each feature
the_classifier = max_correlation_coefficient_CL;  % select a classifier
% the_classifier = libsvm_CL;
the_cross_validator = standard_resample_CV(ds, the_classifier, the_feature_preprocessors);

fprintf('Currently running regular decoding results\n')
% if running the regular results, create the regular cross-validator
the_cross_validator = standard_resample_CV(ds, the_classifier, the_feature_preprocessors);
the_cross_validator.num_resample_runs = 10; % only repeat it 10 times to get the variance
% the name of where to save the results for regular (non-shuffled) decoding results as before
save_file_name = 'decoding_results/population_decoding_results';
% we will also greatly speed up the run-time of the analysis by not creating a full TCT matrix 
% (i.e., we are only training and testing the classifier on the same time bin)
the_cross_validator.test_only_at_training_times = 1;  
 
% run the decoding analysis and save the results
DECODING_RESULTS = the_cross_validator.run_cv_decoding; 
 
% save the results
save(save_file_name, 'DECODING_RESULTS');
%%

%% 
load('population_decoding_results.mat') % load the population_decoding_results.mat generated
performance = squeeze(mean(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.decoding_results,2));

% mean decoding during sampling
de_sampling = mean(performance(:,16:20),2); % The taste detection starts from 16th bins
de_delay1   = mean(performance(:,21:30),2); % Early delay
de_delay2   = mean(performance(:,31:40),2); % Later delay

bar_plot_multi([de_sampling,de_delay1, de_delay2])
ylim([0,1])
ylabel('Decoding accuracy')
names = {'','Sampling', 'Delay1','Delay2',''}
xticks([0,1,2,3,4])
xticklabels(names)
[~,~,stats] = anova1([de_sampling,de_delay1, de_delay2]);
% [~,~,stats] = anova1([de_sampling,de_delay1, de_delay2, de_sampling2,de_delay12, de_delay22]);
figure
[c,~,~,gnames] = multcompare(stats);
%% get the confusion matrix
confusion = DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.confusion_matrix_results.confusion_matrix;
conf_sampling = mean(confusion(:,:,16:20),3);
conf_delay1   = mean(confusion(:,:,21:30),3);
conf_delay2   = mean(confusion(:,:,31:40),3);

% plot the confusion matrix for sampling period
figure; imagesc(conf_sampling)
caxis([5,75])
colormap(jet)
xticks([1,2,3,4])
xticklabels(DECODING_RESULTS.DS_PARAMETERS.label_names_to_use)
yticks([1,2,3,4])
yticklabels(DECODING_RESULTS.DS_PARAMETERS.label_names_to_use)
title('Confusion matrix during sampling')

% plot the confusion matrix for delay1 period
figure; imagesc(conf_delay1)
caxis([5,75])
colormap(jet)
xticks([1,2,3,4])
xticklabels(DECODING_RESULTS.DS_PARAMETERS.label_names_to_use)
yticks([1,2,3,4])
yticklabels(DECODING_RESULTS.DS_PARAMETERS.label_names_to_use)
title('Confusion matrix during delay1')

% plot the confusion matrix for delay2 period
figure; imagesc(conf_delay2./mean(sum(conf_delay2)).*100)
caxis([5,75])
colormap(jet)
xticks([1,2,3,4])
xticklabels(DECODING_RESULTS.DS_PARAMETERS.label_names_to_use)
yticks([1,2,3,4])
yticklabels(DECODING_RESULTS.DS_PARAMETERS.label_names_to_use)
title('Confusion matrix during delay2')

