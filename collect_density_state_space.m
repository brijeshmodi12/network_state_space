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

% Date Created: 27 Oct 2020
% ===================================================================

clc; clear all ; close all; tic;

% Root Path - 
	% root_directory = 'C:\Users\Matteo\Desktop\Brijesh\Open Field Dataset\';
	root_directory = 'D:\matteo_datasets\Open Field Datasets\';

	
	animal_folders_list = ["AH2","AH3","AH4","AH5"];

	sub_directory_name1 = 'network_state_space_bin200';

	subspace_occupancy_folder = 'subspace_occupation';

	density_csd_filename = 'global_median_density_csd.mat';
	density_lfp_filename = 'global_median_density_lfp.mat';

	total_folders = length(animal_folders_list);

	density_across_animals_csd = [];
	density_across_animals_lfp = [];

	animal_colors = viridis(4);

	trial_labels = ["Pre-Sleep","Arena1","Arena2", "Post-Sleep"]
% Core Loop

	for animal_iterator = 1:total_folders

		animal_id = char(animal_folders_list(animal_iterator))
		full_filepath = fullfile(root_directory, animal_id, sub_directory_name1, subspace_occupancy_folder, density_lfp_filename);

		try
			load(full_filepath);
		catch
			fprintf('File not found\n');
			return
		end
		

		density_across_animals_lfp = [density_across_animals_lfp ; global_median_density_lfp];

		fig_lfp = 1 ;
		figure(fig_lfp)
		% set(gcf, 'Position', get(0, 'Screensize'));
		subplot(1,2,1)
		plot(global_median_density_lfp, 'Color', animal_colors(animal_iterator,:),'lineWidth', 2,'HandleVisibility','off')
		hold on
		scatter(1:4, global_median_density_lfp, 100, animal_colors(animal_iterator,:), 'filled')
		title('LFP')
		hold on
		xlim([0 total_folders+1])
		ylim([0 1])
		xticks(1:4)
		yticks(0:0.3:1)
		xticklabels(trial_labels)
		xtickangle(45)
		ylabel('Fraction of Subspace')
		xlabel('Trials')
		ax = gca;
		ax.XAxis.FontSize = 14 ;
		ax.YAxis.FontSize = 14 ;
		legend(animal_folders_list)
		pbaspect([1 1 1])
		set(gca,'box','off')





		full_filepath = fullfile(root_directory, animal_id, sub_directory_name1, subspace_occupancy_folder, density_csd_filename);

		try
			load(full_filepath);
		catch
			fprintf('File not found\n');
			return
		end
	

		density_across_animals_csd = [density_across_animals_csd ; global_median_density_csd];

		fig_lfp = 1 ;
		figure(fig_lfp)
		subplot(1,2,2)
		plot(global_median_density_csd, 'Color', animal_colors(animal_iterator,:),'lineWidth', 2,'HandleVisibility','off')
		hold on
		scatter(1:4, global_median_density_csd, 100, animal_colors(animal_iterator,:), 'filled')
		title('CSD')
		hold on
		
		xlim([0 total_folders+1])
		ylim([0 1])
		xticks(1:4)
		yticks(0:0.3:1)
		xticklabels(trial_labels)
		xtickangle(45)
		ylabel('Fraction of Subspace')
		xlabel('Trials')
		ax = gca;
		ax.XAxis.FontSize = 14 ;
		ax.YAxis.FontSize = 14 ;
		legend(animal_folders_list)
		pbaspect([1 1 1])
		set(gca,'box','off')


		

	end % End for animal_iterator

	output_folder = fullfile(root_directory, 'subspace_occupancy_across_animals')
	mkdir(output_folder)
	filename = sprintf('fraction_area.png');
	fullname = fullfile(output_folder,filename);
	saveas(gcf,fullname);

	






	

	% output_global_filename = fullfile(output_folder,'subspace_occupancy_across_animals.mat');
	% save(output_global_filename,'subspace_occupancy_across_animals');


