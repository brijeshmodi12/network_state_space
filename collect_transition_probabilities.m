% ===================================================================
% Collect transition probabilities for specific sleep states

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

% Date Created: 25 Feb, 2022
% ===================================================================

clc; clear all ; close all; tic;


% root_directory = 'C:\Users\Matteo\Desktop\Brijesh\Open Field Dataset\';
	root_directory = 'D:\matteo_datasets\Open Field Datasets\';

	animal_folders_list = ["AH2","AH3","AH4","AH5"];

	sub_directory_name1 = 'network_state_space_bin200';


	total_folders = length(animal_folders_list);


	rem2rem = [];
	rem2nonrem = [];
	rem2int = [];

	nonrem2rem = [];
	nonrem2nonrem = [];
	nonrem2int = [];

	int2rem = [];
	int2nonrem = [];
	int2int = [];
	
	for animal_iterator = 1:total_folders

		animal_id = char(animal_folders_list(animal_iterator))
	
		animal_root_folder = fullfile(root_directory,animal_id,sub_directory_name1);


		bin_status_folder = fullfile(animal_root_folder,'correlate_lfp_in_bins');
		bin_status_filename = fullfile(bin_status_folder, 'corr_vector_across_trials_lfp.mat');

		% transition matrix folder
		tm_folder = fullfile(animal_root_folder,'lfp_terrain','bin_3_Sleep');
		tm_filename = sprintf('%s_global_transition_data_sleep',animal_id);
		tm_filepath = fullfile(tm_folder, tm_filename);


		try
			load(bin_status_filename);
			load(tm_filepath);
		catch
			fprintf('Bin state classification or transition data file not found\n')
			return
		end


	
		diff_matrix = abs(new_tm_across_trials(:,:,2) - new_tm_across_trials(:,:,1));

		sleep1_indices = 1:9;
	
		sleep2_indices = 28:36;

		rem1_indices = global_isrem_status(sleep1_indices);
		rem2_indices = global_isrem_status(sleep2_indices);

		temp_rem = rem1_indices + rem2_indices;
		rem_ids = find(temp_rem == 2);

		nonrem1_indices = global_isnonrem_status(sleep1_indices);
		nonrem2_indices = global_isnonrem_status(sleep2_indices);

		temp_nonrem = nonrem1_indices + nonrem2_indices;
		nonrem_ids = find(temp_nonrem == 2);

		sleep1_status = rem1_indices + nonrem1_indices;
		sleep2_status = rem2_indices + nonrem2_indices;

		sleep_status = sleep1_status + sleep2_status;
		sleep_status(find(isnan(sleep_status))) = 99;
		int_ids = find(sleep_status == 0);


		% Collect change in transition probabilites from rem
			
			for rem_iterator = rem_ids

				% rem_iterator

				temp1 = []; temp2 = []; temp3 = [];
				temp1 = diff_matrix(rem_iterator, rem_ids);
				rem2rem = [rem2rem; temp1' ];

				temp2 = diff_matrix(rem_iterator, int_ids);
				rem2int = [rem2int ; temp2'];
				
				temp3 = diff_matrix(rem_iterator, nonrem_ids);
				rem2nonrem = [rem2nonrem ; temp3'];		

			end % End for bin_iterator



		% Collect change in transition probabilites from nonrem
					
			for nonrem_iterator = nonrem_ids

				% nonrem_iterator

				temp1 = []; temp2 = []; temp3 = [];

				temp1 = diff_matrix(nonrem_iterator, rem_ids);
				nonrem2rem = [nonrem2rem; temp1' ];

				temp2 = diff_matrix(nonrem_iterator, int_ids);
				nonrem2int = [nonrem2int ; temp2'];
			
				temp3 = diff_matrix(nonrem_iterator, nonrem_ids);
				nonrem2nonrem = [nonrem2nonrem ; temp3'];		

			end % End for bin_iterator


		% Collect change in transition probabilites from int
			
			for int_iterator = int_ids

				% int_iterator

				temp1 = []; temp2 = []; temp3 = [];

				temp1 = diff_matrix(int_iterator, rem_ids);
				int2rem = [int2rem; temp1' ];

				temp2 = diff_matrix(int_iterator, nonrem_ids);
				int2nonrem = [int2nonrem ; temp2'];
			
				temp3 = diff_matrix(int_iterator, int_ids);
				int2int = [int2int ; temp3'];		

			end % End for bin_iterator



		
	end % End for animal_iterator
