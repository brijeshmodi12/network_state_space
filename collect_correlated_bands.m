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

	correlate_lfp_folder = 'correlate_lfp';
	correlate_csd_folder = 'correlate_csd';

	correlate_csd_filename = 'corr_vector_across_trials_csd.mat';
	correlate_lfp_filename = 'corr_vector_across_trials_lfp.mat';

	total_folders = length(animal_folders_list);

	lfp_trials_across_animals = [];
	csd_trials_across_animals = [];

	animal_colors = plasma(total_folders);

	color_matrix = [];

	lower_tri_matrix = tril(magic(18),-1);
	lower_tri_matrix = lower_tri_matrix > 0;

	output_folder = fullfile(root_directory,'correlate_bands_across_trials_across_animals')
	mkdir(output_folder);

	temp_color = viridis(4);

	awake_color = [5 190 120] / 255;
	sleep_color = [96 96 96] / 255;

	trial_colors =  [sleep_color; awake_color; awake_color; sleep_color];

% Core Loop
	for animal_iterator = 1:total_folders

		animal_id = char(animal_folders_list(animal_iterator))
		full_filepath = fullfile(root_directory, animal_id, sub_directory_name1, correlate_lfp_folder, correlate_lfp_filename);

		try
			load(full_filepath);
		catch
			fprintf('File not found\n');
			return
		end

		total_trials = size(sorted_mat_across_trials,2);


		for trial_iterator = 1 : total_trials

			current_trial_matrix = [];

			current_trial_matrix = sorted_mat_across_trials{1,trial_iterator} ;

			current_trial_matrix = current_trial_matrix(lower_tri_matrix);

			current_trial_matrix = reshape(current_trial_matrix, 1, [] ) ;

			lfp_trials_across_animals = [lfp_trials_across_animals;  current_trial_matrix];

		end % End for trial_iterator
		


		full_filepath = fullfile(root_directory, animal_id, sub_directory_name1, correlate_csd_folder, correlate_csd_filename);

		try
			load(full_filepath);
		catch
			fprintf('File not found\n');
			return
		end
	

		for trial_iterator = 1 : total_trials

			current_trial_matrix = [];

			current_trial_matrix = sorted_mat_across_trials{1,trial_iterator}; 

			current_trial_matrix = current_trial_matrix(lower_tri_matrix);

			current_trial_matrix = reshape(current_trial_matrix, 1, [] ) ;

			csd_trials_across_animals = [csd_trials_across_animals;  current_trial_matrix];

		end % End for trial_iterator


		% color_matrix = [color_matrix; repmat(animal_colors(animal_iterator,:),total_trials,1) ];

		color_matrix = [color_matrix; trial_colors];
	

	end % End for animal_iterator


	% -------------------------------------------------------------------
	% LFP 
	% -------------------------------------------------------------------

		n_components = 2;


		[umap_output_lfp, umap_params] = run_umap(lfp_trials_across_animals, 'n_components', n_components, 'n_neighbors', 5 ,'min_dist', 0.1, 'metric', 'euclidean' );
		x = umap_output_lfp(:,1);
		y = umap_output_lfp(:,2);
		if n_components > 2
			z = umap_output_lfp(:,3);
		end

		figure(4)
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
		title(sprintf('LFP Corr Across Trials'));
		legend({'Sleep','Awake'})
		xlim([min(x)-1 max(x)+1])
		ylim([min(y)-1 max(y)+1])

	

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

		figure(4)
		subplot(1,2,1)
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
		xlim([min(x)-1 max(x)+1])
		ylim([min(y)-1 max(y)+1])


		fname = sprintf('alltrials.png', animal_id);
		csd_filepath = fullfile(output_folder, fname);
		saveas(gcf, csd_filepath);
		make_eps_images = 1;
		if make_eps_images	
			epsname =sprintf('alltrials.eps', animal_id);
			fullname_eps = fullfile(output_folder,epsname)	
			set(gcf,'renderer','Painters')
			saveas(gcf,fullname_eps,'epsc')
		end



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



