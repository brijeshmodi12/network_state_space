% ===================================================================
% Get coverage speed v/s rem

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

	lower_tri_matrix = tril(magic(18),-1);
	lower_tri_matrix = lower_tri_matrix > 0;

	output_folder = fullfile(root_directory,'coverage_speed_across_trials')
	mkdir(output_folder);


	global_speed_matrix = [];
	global_rem_ratio_matrix = [];
	global_nonrem_ratio_matrix = [];

	global_power_around_rem = {[], [],[],[],[],[],[]}; 

	npoints_pre = 50;

	global_counter = 1;
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
	

		for trial_iterator = [1 4]

			fprintf('Animal %s Trial : %d ---------------\n', animal_id ,trial_iterator)
			
			ct_coverage_speed = []; ct_rem_ratio = []; ct_rem_indices = []; ct_nonrem_ratio = [];

			ct_coverage_speed = zscore(global_speed_of_coverage{1,trial_iterator});
			ct_rem_ratio = zscore(global_median_rem{1,trial_iterator});
			ct_nonrem_ratio = zscore(global_median_nonrem{1,trial_iterator});

			ct_rem_indices = global_rem_peak_locations{1,trial_iterator};

			if isempty(ct_rem_indices) 
				continue;
			end

			total_rem_bouts = numel(ct_rem_indices);

			diff_between_adjacent_peaks = [100 diff(ct_rem_indices)];

			for rem_iterator = 1:total_rem_bouts

				temp_ind = rem_iterator + 1;
				if diff_between_adjacent_peaks(temp_ind-1) < 10 
					fprintf('Skipped %d\n', rem_iterator)
					continue;
				end

				

				points_around_rem_peak = ct_rem_indices(rem_iterator) - npoints_pre : ct_rem_indices(rem_iterator) + npoints_pre ;

				if max(points_around_rem_peak) > length(ct_coverage_speed) | min(points_around_rem_peak) < 1
					continue;
				end

				fprintf('Accepted %d\n',rem_iterator)

				global_speed_matrix(global_counter,:) = ct_coverage_speed(points_around_rem_peak);

				global_rem_ratio_matrix(global_counter,:) = ct_rem_ratio(points_around_rem_peak);

				global_nonrem_ratio_matrix(global_counter,:) = ct_nonrem_ratio(points_around_rem_peak);


				ct_power_bands = global_median_power{1,trial_iterator};

				band_counter = 1;

				for band_iterator = 1:6
					ct_band = [];

					ct_band = zscore(ct_power_bands(band_iterator,:));

					temp = global_power_around_rem{1,band_counter};

					temp = [temp; ct_band(points_around_rem_peak)];

					global_power_around_rem{1,band_counter} = temp;

					band_counter = band_counter + 1;
				end %
				

				global_counter = global_counter + 1;

			end % End for total_rem_bouts


		end % End for trial_iterator


	end % End for animal_iterator

	% plot(global_speed_matrix')

	figure
	subplot(3,1,1)
	plot(mean(global_speed_matrix))
	subplot(3,1,2)
	plot(mean(global_rem_ratio_matrix))
	subplot(3,1,3)
	plot(mean(global_nonrem_ratio_matrix))

	x_axis = -npoints_pre:npoints_pre;
	x_axis = x_axis * 10;

	output_var_path = fullfile(output_folder, 'speed_vs_rem.mat');
	save(output_var_path, 'global_rem_ratio_matrix', 'global_speed_matrix','global_nonrem_ratio_matrix','x_axis')	
	return



