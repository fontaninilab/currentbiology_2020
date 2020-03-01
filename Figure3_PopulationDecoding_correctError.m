% add the path to the NDT so add_ndt_paths_and_init_rand_generator can be called
clear
toolbox_basedir_name = 'D:\MATLAB\Ephys-analysis\ndt.1.0.4\'
addpath(toolbox_basedir_name);
 
% add the NDT paths using add_ndt_paths_and_init_rand_generator
add_ndt_paths_and_init_rand_generator

%% start population decoding analysis using the NDT
mkdir L_decoding_results % create a folder to store data
mkdir L_decoding_results\shuff_results % create a folder to store data

mkdir R_decoding_results
mkdir R_decoding_results\shuff_results % create a folder to store data

filename = 'LeftCued_CorrectError_100ms_bins_100ms_sampled10repeats.mat';
for i = 0:5
    run_basic_decoding_shuff(i,filename,'L')
end

filename = 'RightCued_CorrectError_100ms_bins_100ms_sampled10repeats.mat';
for i = 0:5
    run_basic_decoding_shuff(i,filename,'R')
end
%% Plot the decoding performance
figure;
result_names{1} = 'R_decoding_results/population_decoding_results';
plot_obj = plot_standard_results_object(result_names);
 
 
% create the names of directories that contain the shuffled data for creating null distributions
% (this is a cell array so that multiple p-values are created when comparing results)
pval_dir_name{1} = 'R_decoding_results/shuff_results/';
plot_obj.p_values = pval_dir_name;
 
% use data from all time bins when creating the null distribution
plot_obj.collapse_all_times_when_estimating_pvals = 1;
plot_obj.p_value_alpha_level = 0.001;
% plot the results as usual
plot_obj.plot_results;
% title('taste Prep R')
%% 
load('population_decoding_results.mat') % load the population_decoding_results.mat generated
performance_L = squeeze(mean(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.decoding_results,2));

load('population_decoding_results.mat') % load the population_decoding_results.mat generated
performance_R = squeeze(mean(DECODING_RESULTS.ZERO_ONE_LOSS_RESULTS.decoding_results,2));

performance_LR = [performance_L; performance_R];
t = -4.45:0.1:1.45;
figure;
boundedline(t,mean(performance,1), 2.807* std(performance_LR)./sqrt(size(performance_LR,1))) % plot 99.5% confidence interval
