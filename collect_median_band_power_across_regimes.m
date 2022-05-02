% ===================================================================
% Get median power in time regimes  across trials for all animals


% Author ::
% Brijesh Modi
% Cortical Microcircuits Lab, 
% European Brain Research Institute, 
% Roma, Italia.
% brijeshmodi12@gmail.com

% For :: 
% Memory Dynamics Lab,
% Donders Centre for Brain, Cognition and Behavior
% Nijmegen, the Netherlands

% Date Created: 29 Oct 2020
% ===================================================================

clc; clear all ; close all; tic;

% % Root Path - 
	root_directory = 'D:\matteo_datasets\Open Field Datasets\';
	% root_directory = 'D:\Matteo Datasets\Open Field Datasets';

	
	animal_folders_list = ["AH2","AH3","AH4","AH5"];

	sub_directory_name1 = 'csd_manifold_awake_sleep';

	sub_directory_name2 = 'quantify_movements_awake_lfp'

	file_to_get = 'global_regime_median_power.mat';

	total_folders = length(animal_folders_list);



	% awake
	awake_bands_median_values_across_trials_across_animals = {[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]};

	% sleep
	sleep_bands_median_values_across_trials_across_animals = {[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]};
	

% Core Loop

	% Core Loop

	for animal_iterator = 1:total_folders

		animal_id = char(animal_folders_list(animal_iterator))
		full_filepath = fullfile(root_directory, animal_id, sub_directory_name1, sub_directory_name2, file_to_get);

		try
			load(full_filepath);
		catch
			fprintf('File not found\n');
			return
		end

		total_trials = size(global_regime_median_power,2);


		for trial_iterator = 1:total_trials

			if trial_iterator == 1 | trial_iterator == 4
				output_variable = sleep_bands_median_values_across_trials_across_animals;
			else
				output_variable = awake_bands_median_values_across_trials_across_animals;
			end


			total_bands = size(global_regime_median_power{1},2);

			% bands for current trials
			ct_bands = global_regime_median_power{trial_iterator};

			
			for band_iterator = 1:18


				temp1 = output_variable{band_iterator};

				temp1 = [temp1 ct_bands(:, band_iterator) ];

				output_variable{band_iterator} = temp1;

			end % End for band_iterator


			if trial_iterator == 1 | trial_iterator == 4
				sleep_bands_median_values_across_trials_across_animals = output_variable;
			else
				awake_bands_median_values_across_trials_across_animals = output_variable;
			end

			% return
 			
		end % End for trial_iterator


	end % End for animal_iterator


return