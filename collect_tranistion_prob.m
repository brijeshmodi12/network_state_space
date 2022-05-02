% ===================================================================
% Get Change in Transition Probabilities across trials

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

% Date Created: 28 March 2021
% ===================================================================

clc; clear all ; close all; tic;
fprintf('Collecting Transition Probabilities\n')
% Root Path - 
	% root_directory = 'C:\Users\Matteo\Desktop\Brijesh\Open Field Dataset\';
	root_directory = 'D:\matteo_datasets\Open Field Datasets\';

	
	animal_folders_list = ["AH2","AH3","AH4","AH5"];

	sub_directory_name1 = 'network_state_space_bin200';

	lfp_folder = 'lfp_terrain';
	awake_folder = 'bin_10';
	sleep_folder = 'bin_3';

	total_folders = length(animal_folders_list);

	awake_data = [];
	sleep_data = [];


	trial_labels = ["Pre-Sleep","Arena1","Arena2", "Post-Sleep"]
% Core Loop

	for animal_iterator = 1:total_folders

		animal_id = char(animal_folders_list(animal_iterator))
		fname = sprintf('%s_resultant_diff.mat', animal_id);

		full_filepath = fullfile(root_directory, animal_id, sub_directory_name1, lfp_folder, awake_folder, fname);

		try
			load(full_filepath);
		catch
			fprintf('File not found\n');
			return
		end
		

		awake_data = [awake_data; [resultant_incoming_change, resultant_outgoing_change]] ;

		full_filepath = fullfile(root_directory, animal_id, sub_directory_name1, lfp_folder, sleep_folder, fname);

		try
			load(full_filepath);
		catch
			fprintf('File not found\n');
			return
		end


		sleep_data = [sleep_data; [resultant_incoming_change, resultant_outgoing_change]] ;


		

		

	end % End for animal_iterator

	






	

	% output_global_filename = fullfile(output_folder,'subspace_occupancy_across_animals.mat');
	% save(output_global_filename,'subspace_occupancy_across_animals');


