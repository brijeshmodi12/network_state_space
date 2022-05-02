% ===================================================================
% Get Oscillators correlation time frames csd and do UMAP

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

% Date Created: 21 Jan 2021
% ===================================================================

clc; clear all ; close all; tic;

% Root Path - 
	% root_directory = 'C:\Users\Matteo\Desktop\Brijesh\Open Field Dataset\';
	root_directory = 'D:\matteo_datasets\Open Field Datasets\';

	
	animal_folders_list = ["AH2","AH3","AH4","AH5"];

	sub_directory_name1 = 'network_state_space_bin200';

	correlate_csd_folder = 'correlate_csd_time';

	correlate_csd_filename = 'corr_vector_across_trials_csd_over_time.mat';

	total_folders = length(animal_folders_list);

	csd_trials_across_animals = [];
	csd_trials_across_animals = [];

	color_gray = [134, 138, 145 ] ./ 255;
	color_pink = [189, 11, 73] ./ 255;

	animal_colors = plasma(total_folders);

	color_matrix = [];

	lower_tri_matrix = tril(magic(18),-1);
	lower_tri_matrix = lower_tri_matrix > 0;

	output_folder = fullfile(root_directory,'correlate_bands_over_time_across_trials_across_animals');
	mkdir(output_folder);

	global_rem_values = [];
	global_nonrem_values = [];
	sleep_status_across_animals = [];
	csd_trials_across_animals = [];

	trial_time_count = [0];
	animal_point_count = [0];

	awake_sleep_color = viridis(2);

	% Core Loop
	for animal_iterator = 1:total_folders
		
		% csd_trials_across_animals = [];
		
		% global_rem_values = [];
		% global_nonrem_values = [];

		rem_colors = [];
		nonrem_colors = [];
		x = []; y = []; umap_output_csd = [];
		

		animal_id = char(animal_folders_list(animal_iterator));



		% -------------------------------------------------------------------
		% Collect csd correlation matrix
		% -------------------------------------------------------------------	
			full_filepath = fullfile(root_directory, animal_id, sub_directory_name1, correlate_csd_folder, correlate_csd_filename);

			try
				load(full_filepath);
			catch
				fprintf('File not found\n');
				return
			end

			total_trials = size(sorted_mat_across_trials,2);

			trial_colors = viridis(total_trials);
			trial_color_matrix = [];

			animal_points_counter = 0;

			for trial_iterator = 1:total_trials		

				current_trial_matrix = [];

				current_trial_matrix = sorted_mat_across_trials{1,trial_iterator} ;

				total_time_points = size(current_trial_matrix,3);

				% Iterate over each time frame
				for time_iterator = 1 : total_time_points

					current_time_matrix = current_trial_matrix(:,:,time_iterator);

					current_time_matrix = current_time_matrix(lower_tri_matrix);

					current_time_matrix = reshape(current_time_matrix, 1, [] ) ;

					csd_trials_across_animals = [csd_trials_across_animals;  current_time_matrix];

				end % End for time_iterator
				
				global_rem_values = [global_rem_values  global_median_rem{1,trial_iterator} ]	;
				global_nonrem_values = [global_nonrem_values  global_median_nonrem{1,trial_iterator} ]	;
				sleep_status_across_animals = [sleep_status_across_animals global_sleep_status{1,trial_iterator} ];
				

			end % End for trial_iterator

		end	% End for animal_iterator	

			% Assign Colors
			nonrem_range = [min(global_nonrem_values) prctile(global_nonrem_values,95)]	;
			[nonrem_cmap, nonrem_colors] = assign_colors(global_nonrem_values, nonrem_range, 10, 'pink');

			% Assign Colors
			rem_range = [min(global_rem_values) prctile(global_rem_values,95)]	;
			[rem_cmap, rem_colors] = assign_colors(global_rem_values, rem_range, 10, 'magma');

			sleep_status_range = [min(sleep_status_across_animals) max(sleep_status_across_animals)];
			[ss_cmap, sleep_status_colors] = assign_colors(sleep_status_across_animals, sleep_status_range, 2, 'viridis');

		% -------------------------------------------------------------------
		% Do Umap csd
		% -------------------------------------------------------------------	

			n_components = 2;

			[umap_output_csd, umap_params] = run_umap(csd_trials_across_animals, 'n_components', n_components, 'n_neighbors', 5 ,'min_dist', 0.1, 'metric', 'cosine' );
			x = umap_output_csd(:,1);
			y = umap_output_csd(:,2);
			if n_components > 2
				z = umap_output_csd(:,3);
			end

			figure(123)
			set(gcf, 'Position', get(0, 'Screensize'));
			scatter(x(:), y(:), 50, sleep_status_colors , 'filled');
			colormap(ss_cmap)
			h = colorbar ;
			ylabel(h,'Sleep Status')
			caxis(sleep_status_range);

			title(sprintf('All Sleep Status Overlay on CSD CorrMat over time',animal_id))
			pbaspect([1 1 1])
			legend('Sleep')



			fname = sprintf('All_sleep_status_overlay_csd_corr_mat_across_trials_raw.png', animal_id);
			csd_filepath = fullfile(output_folder, fname);
			saveas(gcf, csd_filepath);


			figure(animal_iterator+4)
			set(gcf, 'Position', get(0, 'Screensize'));
			scatter(x(:), y(:), 50,  'k' , 'filled');
			title(sprintf('All csd CorrMat over time',animal_id))
			pbaspect([1 1 1])

			fname = sprintf('All_csd_corr_mat_across_trials_raw.png', animal_id);
			csd_filepath = fullfile(output_folder, fname);
			saveas(gcf, csd_filepath);



			figure(animal_iterator+2)
			set(gcf, 'Position', get(0, 'Screensize'));
			subplot(1,2,1)
			scatter(x(:), y(:), 50,  rem_colors(:,:) , 'filled');
			title(sprintf('REM overlay'))
			
			cbh1 = colorbar;
			colormap(cbh1,'magma')
			ylabel(cbh1,'Median Theta/Delta Ratio')
			caxis([min(global_rem_values) max(global_rem_values)])
			% xlim([min(x) max(x)])
			% ylim([min(y) max(y)])
			pbaspect([1 1 1])

			figure(animal_iterator+2)
			subplot(1,2,2)
			scatter(x(:), y(:), 50,  nonrem_colors(:,:) , 'filled');
			title(sprintf('NonREM overlay'))
			cbh2 = colorbar;
			colormap(cbh2,'pink')
			ylabel(cbh2,'Median Delta x Spindle Product')
			caxis([min(global_nonrem_values) max(global_nonrem_values)])
			% xlim([min(x) max(x)])
			% ylim([min(y) max(y)])
			pbaspect([1 1 1])

			suptitle(sprintf('All REM and NonREM overlay on CSD CorrMat over time',animal_id))

			fname = sprintf('All_csd_corr_mat_across_trials.png', animal_id);
			csd_filepath = fullfile(output_folder, fname);
			saveas(gcf, csd_filepath);
			
			close all
			
			return
		% -------------------------------------------------------------------
		% Plot REM Trajectories
		% -------------------------------------------------------------------	
			% for trial_iterator = 1:total_trials

			% 	ct_rem_peaks = global_rem_peak_locations{1,trial_iterator};

			% 	if isempty(ct_rem_peaks)
			% 		continue;
			% 	end

			% 	total_peaks = length(ct_rem_peaks);

			% 	range_start = global_timeframe_counts(trial_iterator) + 1 ;
			% 	range_end = global_timeframe_counts(trial_iterator + 1) ;

			% 	trial_timeframe_range = [];
			% 	trial_timeframe_range = range_start : range_end; 

			% 	trial_x = umap_output_csd(trial_timeframe_range,1);
			% 	trial_y = umap_output_csd(trial_timeframe_range,2);

			% 	rem_peak_fig_id = 112321;
			% 	figure(rem_peak_fig_id)

			% 	trajectory_colors = viridis(total_peaks);
			% 	for peak_iterator = 1:total_peaks

			% 		current_peak = ct_rem_peaks(peak_iterator);

			% 		trajectory_indices = current_peak-4 : current_peak;
			% 		subplot(2,2,peak_iterator)
			% 		scatter(x,y,30,color_gray,'filled');
			% 		hold on
			% 		plot(x(trajectory_indices), y(trajectory_indices), 'Color',trajectory_colors(peak_iterator,:), 'LineWidth',3);
			% 		hold on
			% 		scatter(x(current_peak), y(current_peak), 30, color_pink, 'filled');

			% 		pbaspect([1 1 1])


			



			% 	end % End for peak_iterator


			% end % End for trial_iterator

			






		
	% end 

	return
