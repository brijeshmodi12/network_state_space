% ===================================================================
% Get Oscillators correlation time frames CSD and do UMAP

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

	correlate_csd_folder = 'correlate_csd_in_bins';

	correlate_csd_filename = 'corr_vector_across_trials_csd.mat';

	total_folders = length(animal_folders_list);


	color_gray = [134, 138, 145 ] ./ 255;
	color_pink = [189, 11, 73] ./ 255;

	animal_colors = plasma(total_folders);

	color_matrix = [];

	lower_tri_matrix = tril(magic(18),-1);
	lower_tri_matrix = lower_tri_matrix > 0;

	output_folder = fullfile(root_directory,'correlate_bands_in_bins_across_animals');
	mkdir(output_folder);

	rem_across_animals = [];
	nonrem_across_animals = [];
	sleep_status_across_animals = [];
	bin_corrmat_across_animals = [];
	
	awake_sleep_color = viridis(2);

	% Core Loop
	for animal_iterator = 1:total_folders
		
		
		animal_id = char(animal_folders_list(animal_iterator));

		% -------------------------------------------------------------------
		% Collect CSD correlation matrix
		% -------------------------------------------------------------------	
			full_filepath = fullfile(root_directory, animal_id, sub_directory_name1, correlate_csd_folder, correlate_csd_filename);

			try
				load(full_filepath);
			catch
				fprintf('File not found\n');
				return
			end

			

			total_bins = numel(global_bin_median_rem);

			for bin_iterator = 1:total_bins

				temp_sorted_matrix = [];
				temp_sorted_matrix = global_bin_sorted_matrix(:,:,bin_iterator);

				temp_sorted_matrix = temp_sorted_matrix(lower_tri_matrix);

				temp_sorted_matrix = reshape(temp_sorted_matrix, 1, [] ) ;

				bin_corrmat_across_animals = [bin_corrmat_across_animals;  temp_sorted_matrix];

			end  % End for bin_iterator

			rem_across_animals = [rem_across_animals global_bin_median_rem];

			nonrem_across_animals = [nonrem_across_animals global_bin_median_nonrem];

			sleep_status_across_animals = [sleep_status_across_animals global_bin_sleep_status]; 

		end	% End for animal_iterator	
		

		% Assign Colors
		nonrem_range = [min(nonrem_across_animals) max(nonrem_across_animals)];
		[nonrem_cmap, nonrem_colors] = assign_colors(nonrem_across_animals, nonrem_range, 10, 'pink');

		% Assign Colors
		rem_range = [min(rem_across_animals) max(rem_across_animals)]	;
		[rem_cmap, rem_colors] = assign_colors(rem_across_animals, rem_range, 10, 'magma');

		sleep_status_range = [min(sleep_status_across_animals) max(sleep_status_across_animals)];
		[ss_cmap, sleep_status_colors] = assign_colors(sleep_status_across_animals, sleep_status_range, 2, 'viridis');


		% -------------------------------------------------------------------
		% Do Umap CSD
		% -------------------------------------------------------------------	

			n_components = 2;

			[umap_output_csd, umap_params] = run_umap(bin_corrmat_across_animals, 'n_components', n_components, 'n_neighbors', 5 ,'min_dist', 0.1, 'metric', 'cosine' );
			x = umap_output_csd(:,1);
			y = umap_output_csd(:,2);
			if n_components > 2
				z = umap_output_csd(:,3);
			end

			figure(123)
			set(gcf, 'Position', get(0, 'Screensize'));
			scatter(x(:), y(:), 100, sleep_status_colors , 'filled');
			colormap(ss_cmap)
			h = colorbar ;
			ylabel(h,'Sleep Status')
			caxis(sleep_status_range);
			% ylim([0 15])
			% xlim([-10 10])

			title(sprintf('All Sleep Status Overlay on CSD CorrMat in Bins',animal_id))
			pbaspect([1 1 1])
			legend('Sleep')



			fname = sprintf('All_sleep_status_overlay_csd_corrmat_bins.png', animal_id);
			csd_filepath = fullfile(output_folder, fname);
			saveas(gcf, csd_filepath);


			figure(animal_iterator+4)
			set(gcf, 'Position', get(0, 'Screensize'));
			scatter(x(:), y(:), 100,  'k' , 'filled');
			% ylim([0 15])
			% xlim([-10 10])
			title(sprintf('All CSD CorrMat in Bins',animal_id))
			pbaspect([1 1 1])


			fname = sprintf('All_csd_corrmat_bins.png', animal_id);
			csd_filepath = fullfile(output_folder, fname);
			saveas(gcf, csd_filepath);



			figure(animal_iterator+2)
			set(gcf, 'Position', get(0, 'Screensize'));
			subplot(1,2,1)
			scatter(x(:), y(:), 100,  rem_colors(:,:) , 'filled');
			title(sprintf('REM overlay'))
			
			cbh1 = colorbar;
			colormap(cbh1,'magma')
			ylabel(cbh1,'Median Theta/Delta Ratio')
			caxis(rem_range)
			% ylim([0 15])
			% xlim([-10 10])
			% xlim([min(x) max(x)])
			% ylim([min(y) max(y)])
			pbaspect([1 1 1])

			figure(animal_iterator+2)
			subplot(1,2,2)
			scatter(x(:), y(:), 100,  nonrem_colors(:,:) , 'filled');
			title(sprintf('NonREM overlay'))
			cbh2 = colorbar;
			colormap(cbh2,'pink')
			ylabel(cbh2,'Median Delta x Spindle Product')
			caxis(nonrem_range)
			% ylim([0 15])
			% xlim([-10 10])
			pbaspect([1 1 1])

			suptitle(sprintf('All REM and NonREM overlay on CSD CorrMat in Bins',animal_id))

			fname = sprintf('All_csd_corr_mat_bins.png', animal_id);
			csd_filepath = fullfile(output_folder, fname);
			saveas(gcf, csd_filepath);
				return
				
			close all
			
			return
		


		
	% end 

	return
