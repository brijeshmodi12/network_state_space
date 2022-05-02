% ===================================================================
% Get fraction area across trials for all animals

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

% Date Created: 21 Jan, 2021
% ===================================================================

clc; clear all ; close all; tic;

% Root Path - 
	% root_directory = 'C:\Users\Matteo\Desktop\Brijesh\Open Field Dataset\';
	root_directory = 'D:\matteo_datasets\Open Field Datasets\';

	
	animal_folders_list = ["AH2","AH3","AH4","AH5"];

	sub_directory_name1 = 'network_state_space_bin200';

	area_coverage_lfp_folder = 'area_coverage_lfp';
	area_coverage_csd_folder = 'area_coverage_csd';

	area_coverage_csd_filename = 'global_coverage_variables_csd.mat';
	area_coverage_lfp_filename = 'global_coverage_variables_lfp.mat';

	total_folders = length(animal_folders_list);

	lfp_trials_across_animals = [];
	csd_trials_across_animals = [];

	animal_colors = plasma(total_folders);

	color_matrix = [];



	output_folder = fullfile(root_directory,'coverage_speed_across_trials')
	mkdir(output_folder);

% Core Loop
	for animal_iterator = 1:total_folders

		animal_id = char(animal_folders_list(animal_iterator))
		full_filepath = fullfile(root_directory, animal_id, sub_directory_name1, area_coverage_lfp_folder, area_coverage_lfp_filename);

		try
			load(full_filepath);
		catch
			fprintf('File not found\n');
			return
		end

		current_animal_average_speed = [];
	
		current_animal_average_speed = cellfun(@median, global_speed_of_coverage);

		current_animal_average_speed = current_animal_average_speed / time_interval_for_area_sec;

		lfp_trials_across_animals = [lfp_trials_across_animals; current_animal_average_speed];


		continue;

		full_filepath = fullfile(root_directory, animal_id, sub_directory_name1, area_coverage_csd_folder, area_coverage_csd_filename);

		try
			load(full_filepath);
		catch
			fprintf('File not found\n');
			return
		end
	

		current_animal_average_speed = [];
	
		current_animal_average_speed = cellfun(@median, global_speed_of_coverage);
		
		current_animal_average_speed = current_animal_average_speed / time_interval_for_area_sec;

		csd_trials_across_animals = [csd_trials_across_animals; current_animal_average_speed];


	

	end % End for animal_iterator

	return


		for animal_iterator = 1:total_folders

			animal_points = (animal_iterator - 1 ) * total_trials + 1 : total_trials * animal_iterator;
			figure(5)
			subplot(2,2,animal_iterator)
			scatter(x(animal_points), y(animal_points),  100,  viridis(total_trials) , 'filled');
			hold on
			plot(x(animal_points), y(animal_points), 'lineWidth', 2, 'color',animal_colors(animal_iterator,:));
			hold on
			title(sprintf('Animal %d',animal_iterator+1))
			cbh = colorbar
			colormap('viridis')
			ylabel(cbh,'Trials')
			caxis([1 4])
		end 

		suptitle(sprintf('LFP'))

		fname = sprintf('lfp_corr_mat_across_trials_across_animals.png', animal_id);
		lfp_filepath = fullfile(output_folder, fname);
		saveas(gcf, lfp_filepath);


	% -------------------------------------------------------------------
	% CSD
	% -------------------------------------------------------------------

		n_components = 2;


		[umap_output_csd, umap_params] = run_umap(csd_trials_across_animals, 'n_components', n_components, 'n_neighbors', 5 ,'min_dist', 0.1, 'metric', 'euclidean' );
		x = umap_output_csd(:,1);
		y = umap_output_csd(:,2);
		if n_components > 2
			z = umap_output_csd(:,3);
		end

		figure(9)
		subplot(1,2,2)
		if n_components > 2
			scatter3(x,y,z, 50,  color_matrix ,'filled');
			str = '3D';
		else
			scatter(x,y, 50,  color_matrix ,'filled');	
			str = '2D';
			hold
		end

		pbaspect([1 1 1]) 
		title(sprintf('CSD Corr Across Trials'));


		% for animal_iterator = 1:total_folders

		% 	animal_points = (animal_iterator - 1 ) * total_trials + 1 : total_trials * animal_iterator;
		% 	figure(5)
		% 	scatter3(x(animal_points), y(animal_points), z(animal_points), 100,  viridis(total_trials) , 'filled');
		% 	hold on
		% 	plot3(x(animal_points), y(animal_points), z(animal_points), 'color', animal_colors(animal_iterator,:) , 'lineWidth', 2);
		% 	hold on
			
		% end 

		for animal_iterator = 1:total_folders

			animal_points = (animal_iterator - 1 ) * total_trials + 1 : total_trials * animal_iterator;
			figure(51)
			subplot(2,2,animal_iterator)
			scatter(x(animal_points), y(animal_points),  100,  viridis(total_trials) , 'filled');
			hold on
			plot(x(animal_points), y(animal_points), 'lineWidth', 2, 'color',animal_colors(animal_iterator,:));
			hold on
			title(sprintf('Animal %d',animal_iterator+1))
			cbh = colorbar
			colormap('viridis')
			ylabel(cbh,'Trials')
			caxis([1 4])

			
		end 

		suptitle(sprintf('CSD'))

		fname = sprintf('CSD_corr_mat_across_trials_across_animals.png', animal_id);
		csd_filepath = fullfile(output_folder, fname);
		saveas(gcf, csd_filepath);



