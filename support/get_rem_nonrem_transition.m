%% get_rem_nonrem_transition: function description
function [rem2rem, nonrem2nonrem] = get_rem_nonrem_transition(diff_matrix, animal_id)

	rem2rem = [];
	nonrem2nonrem = [];

	root_directory = 'D:\matteo_datasets\Open Field Datasets\';
	
	sub_directory_name1 = 'network_state_space_bin200';
	
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

			temp1 = []; temp2 = []; temp3 = [];
			temp1 = diff_matrix(rem_iterator, rem_ids);
			rem2rem = [rem2rem; temp1' ];


		end % End for bin_iterator



	% Collect change in transition probabilites from nonrem
				
		for nonrem_iterator = nonrem_ids

			temp1 = []; temp2 = []; temp3 = [];
			
			temp3 = diff_matrix(nonrem_iterator, nonrem_ids);
			nonrem2nonrem = [nonrem2nonrem ; temp3'];		

		end % End for bin_iterator



