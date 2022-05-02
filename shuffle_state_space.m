% ===================================================================
% Get state space features and calculate occupancy

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

% Date Created: 11 April 2022
% ===================================================================

clc; clear all ; close all; tic;

% Root Path - 
	% root_directory = 'C:\Users\Matteo\Desktop\Brijesh\Open Field Dataset\';
	root_directory = 'D:\matteo_datasets\Open Field Datasets\';

	
	animal_folders_list = ["AH2","AH3","AH4","AH5"];

	sub_directory_name1 = 'network_state_space_bin200';

	subspace_occupancy_folder = 'subspace_occupation';


	do_joint_shuffling = 00;
	compute_shuffled_occupancy = 01;

	occupancy_csd_filename = 'global_subspace_occupancy_fraction_csd.mat';
	occupancy_lfp_filename = 'global_subspace_occupancy_fraction_lfp.mat';


	og_state_space_filename ='state_space_data.mat';

	total_folders = length(animal_folders_list);

	color_gray = [134, 138, 145 ] ./ 255;

	animal_colors = viridis(4);

	trial_labels = ["Pre-Sleep","Arena1","Arena2", "Post-Sleep"]
% Core Loop

	global_diff_occupancy = [];

	make_eps_images = 1;

	for animal_iterator = 1:total_folders

		animal_id = char(animal_folders_list(animal_iterator))
		full_filepath = fullfile(root_directory, animal_id, sub_directory_name1, og_state_space_filename);

		try
			load(full_filepath);
		catch
			fprintf('File not found\n');
			return
		end

		output_folder = fullfile(root_directory, animal_id, sub_directory_name1, subspace_occupancy_folder);


		% Shuffle State Space
		total_features = size(global_lfp_states, 2);
		total_states = size(global_lfp_states,1);

		total_random_iterations = 100;

		% Joint Shuffling

		if do_joint_shuffling

			% random_matrix_generator 

			global_random_matrices = [];
			for random_iterator = 1:total_random_iterations
				mat = [];
				for i = 1:18
					mat(:,i) = randperm(total_states);
				end
				global_random_matrices(:,:,random_iterator) = mat;
			end




			for random_iterator = 1:total_random_iterations

				shuffled_state_space = NaN * zeros(size(global_lfp_states));
				shuffled_projection = [];
				shuffled_indices_matrix = [];
				

				for feature_iterator = 1:total_features
					
					current_feature =[]; shuffled_indices = [];

					shuffled_indices = global_random_matrices(:,feature_iterator, random_iterator);

					current_feature = global_lfp_states(:,feature_iterator);

					shuffled_state_space(:,feature_iterator) = current_feature(shuffled_indices);

				end % End for feature_iterator

				% UMAP on Shuffled State Space


				[shuffled_projection, umap_params] = run_umap(shuffled_state_space,  'n_components', 2, 'n_neighbors', 25 ,'min_dist', 0.1, 'metric', 'euclidean' );
				close all
				shuffled_projection_filename = sprintf('shuffled_projection_%d.mat',random_iterator);
				shuffled_projection_filepath = fullfile(output_folder, shuffled_projection_filename);
				save(shuffled_projection_filepath, 'shuffled_projection', 'shuffled_state_space','shuffled_indices_matrix');
				

			end % End for random_iterator

		end % End if do_joint_shuffli01;



		if compute_shuffled_occupancy

			shuffled_awake_occupancy = [];
			shuffled_sleep_occupancy = [];
			diff_occupancy = [];

			for random_iterator = 1:total_random_iterations

				shuffled_projection_filename = sprintf('shuffled_projection_%d.mat',random_iterator);
				shuffled_projection_filepath = fullfile(output_folder, shuffled_projection_filename);

				shuffled_projection = [];

				load(shuffled_projection_filepath)

				global_x = shuffled_projection(:,1);
				global_y = shuffled_projection(:,2);

				% Divide entire umap output into bins
				global_reference_bins_x  = linspace( min(global_x), max(global_x), 100 );
				global_reference_bins_y  = linspace( min(global_y), max(global_y), 100 );

				all_occupancy = hist3([global_x, global_y],'Edges',{global_reference_bins_x, global_reference_bins_y});

				total_occupied_bins = sum(all_occupancy(:) > 0);

				all_pv = 1:global_states_count_lfp(5);

				total_trials = length(trial_labels);
				trial_colors = magma(total_trials);

				occupancy_values = [];
				
				for trial_iterator = 1 : total_trials

					range_start = global_states_count_lfp(trial_iterator) + 1 ;
					range_end = global_states_count_lfp(trial_iterator + 1) ;

					pv_range = range_start : range_end; 

					x = shuffled_projection(pv_range,1);
					y = shuffled_projection(pv_range,2);		
					non_pv_range = setdiff(all_pv, pv_range);	

					% figure(random_iterator)
					% set(gcf, 'Position', get(0, 'Screensize'));
					% subplot(1,4,trial_iterator)
					% scatter(shuffled_projection(non_pv_range,1),shuffled_projection(non_pv_range,2),5 ,color_gray, 'filled')
					% hold on
					% scatter(x,y,5, trial_colors(trial_iterator,:) ,'filled');		
					% pbaspect([1 1 1])
					% title(sprintf('%s', char(trial_labels(trial_iterator))))
					% xlim(lfp_xlim);
					% ylim(lfp_ylim);


					trial_occupancy = hist3([x, y],'Edges',{global_reference_bins_x, global_reference_bins_y});
					trial_occupied_bins = sum(trial_occupancy(:) > 0);

					occupancy_values(trial_iterator) = trial_occupied_bins / total_occupied_bins;

					if ismember(trial_iterator, [2,3])
						shuffled_awake_occupancy = [shuffled_awake_occupancy; occupancy_values(trial_iterator)];
					else 
						shuffled_sleep_occupancy = [shuffled_awake_occupancy; occupancy_values(trial_iterator)];
					end



					
				end % End for trial_iterator

				diff_occupancy = [ diff_occupancy abs(mean(occupancy_values([1,4])) - mean(occupancy_values([2 3]))  )];
					

					

			end % End for random_iterator

			observed_data_filepath = fullfile(output_folder,'global_subspace_occupancy_fraction_lfp.mat')
			load(observed_data_filepath)

			temp = global_subspace_occupancy_fraction_lfp;

			observed_diff_occupancy = abs(mean(temp([1,4])) - mean(temp([2 3])))


			figure(1)
			histogram(diff_occupancy, 50, 'FaceColor', color_gray, 'EdgeColor', color_gray);	
			xlim([min(diff_occupancy)- 0.1 , max(diff_occupancy)+0.1])
			title(sprintf('95 percentile = %0.3f; observed data = %0.3f', prctile(diff_occupancy, 95), observed_diff_occupancy ))	
			filename = sprintf('distribution of diff occupancy.png');
			fullname = fullfile(output_folder,filename);
			saveas(gcf,fullname);

			if make_eps_images	
				epsname = sprintf('%s_distribution_of_diff_occupancy.eps', animal_id);
				fullname_eps = fullfile(output_folder, epsname)	
				set(gcf,'renderer','Painters')
				saveas(gcf,fullname_eps,'epsc')
			end

			return


		end % End if compute_shuffled_occupancy


		

	


	end % End for animal_iterator




	prctile(global_diff_occupancy,[5,95])
	filename = sprintf('fraction_area.png');
	fullname = fullfile(output_folder,filename);
	saveas(gcf,fullname);

	


	if make_eps_images	
		epsname = sprintf('%s_global_state_space.eps', animal_id);
		fullname_eps = fullfile(global_output_folder,epsname)	
		set(gcf,'renderer','Painters')
		saveas(gcf,fullname_eps,'epsc')
	end


	return
	






	

	% output_global_filename = fullfile(output_folder,'subspace_occupancy_across_animals.mat');
	% save(output_global_filename,'subspace_occupancy_across_animals');


