% ===================================================================
	% State space for current source density
	% from Pyramidal, Radiatum and SLM layer of CA1

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

	% Date Created: 23 Nov 2020

	% AH2

	% CSD Layers:
		% Pyr = 2 
		% Rad = 5 or 6 
		% SLM = 9 

	% for example: Raw_CSD(2,:) to get the pyramidal layer CSD. Etc.

	% !Attention: the signal are already DownSampled! SR = 1000

	% In C_CSD you find the signal already filtered for some bands. But you can directly filter the Raw_CSD if you want. 

	
% ===================================================================

clc; clear all ; close all; tic;

fprintf('------------------------------------------------\n')
fprintf('CSD State Space\n')
fprintf('------------------------------------------------\n')

% -------------------------------------------------------------------
% Plot Variables
% -------------------------------------------------------------------
	% Set the following variables to '1' for plotting
	% Plotting variables
		power_overlay = 1;
		% make_video = 0;	
		% cell_firing_overlay = 0;
		% plot_distance_matrix = 0;
		% rem_detection = 1;
		% nonrem_detection = 1;
		% quantify_angles = 1;
		% quantify_movements = 0;
		% quantify_area = 1;
		% subspace_occupation = 1;
		% time_overlay = 0;
		% quantify_movements_2 = 0;


% -------------------------------------------------------------------
% General Constants
% -------------------------------------------------------------------
	% Paths and Filenames
		root_directory = 'C:\Users\Matteo\Desktop\Brijesh\Open Field Dataset\';
		% root_directory = 'D:\matteo_datasets\Open Field Datasets\';

		animal_folders_list = ["AH2","AH3","AH4","AH5"];

		csd_root_dirname = 'CSD';

		csd_filename = 'C_Raw_CSD.mat';
		lfp_filename = 'C_Raw_LFP.mat';

		total_folders = length(animal_folders_list);

		color_gray = [134, 138, 145 ] ./ 255;

		% Combine all the bands in an single array so you can apply loops 
		bands_strings = ["Delta", "Spindle", "Theta", "Slow gamma", "Medium gamma", "Fast gamma"];

	

		% Indices to access pyramidal, radiatum and slm layers
		csd_layer_indices = [2 6 9];
		lfp_layer_indices = csd_layer_indices + 1;

		layer_strings = ["Pyramidal", "Radiatum", "SLM"];


		% Custom Frequency Bands
		% Frequency range for delta (in Hz)
		delta_band = [ 1 5 ];

		% Frequency range for spindle (in Hz)
		spindle_band = [ 10 20 ];			

		% Frequency range for theta (in Hz)
		theta_band = [ 6 10 ];

		% Frequency range for slow gamma (in Hz)
		sgamma_band = [ 20 45 ];

		% Frequency range for medium gamma (in Hz)
		mgamma_band = [ 60 90 ];

		% Frequency range for fast gamma (in Hz) 
		fgamma_band = [ 100 200 ];
	

		% Combine all the bands in an single array so you can apply loops 
		bands_strings = ["Delta", "Spindle", "Theta", "Slow gamma", "Medium gamma", "Fast gamma"];
		bands_array = [delta_band; spindle_band; theta_band ; sgamma_band; mgamma_band ; fgamma_band];


	% Log/Linearly spaced Frequency Bands
		% start_band = 4;	
		% end_band = 300;
		% total_bands = 15;
		% continous_bands = linspace(start_band,end_band, total_bands+1);
		% continous_bands = logspace(log10(start_band),log10(end_band),total_bands+1);
		% bands_array(1:total_bands,1) = continous_bands(1:end-1);
		% bands_array(1:total_bands,2) = continous_bands(2:end);

		% return

		total_bands = size(bands_array,1);

		total_layers = length(csd_layer_indices);

		downsampling_frequency = 1000;

		% Bin size for state space construction -  seconds
		bin_size_sec = 0.2;

		% Total samples in each bin
		bin_size_lfp = bin_size_sec * downsampling_frequency;

		% For smoothing
		% Gaussian kernel
		smoothing_std = 3;

		smoothing_width = smoothing_std * 6;



% -------------------------------------------------------------------
% Core Loop
% -------------------------------------------------------------------
	% Process each animal
	for animal_iterator = 1:total_folders

		animal_id = char(animal_folders_list(animal_iterator));
		fprintf('------ Animal %s ------\n', animal_id)
		
		
		% -------------------------------------------------------------------
		% Animal specific Constants
		% -------------------------------------------------------------------		
			switch animal_id
				case 'AH2'
					trial_folders = ["preSleep", "1Rectangle", "2Circle", "postSleep"]; % AH2 
					output_folder_string = 'csd_manifold_awake_sleep';

					% trial_folders = [ "1Rectangle", "2Circle"]; % AH2 
					% output_folder_string = 'all_trials_power_manifold_awake';
					% only_awake = 1;

					% trial_folders = ["preSleep", "postSleep"]; % AH2 
					% output_folder_string = 	'all_trials_power_manifold_sleep';
					% only_sleep = 1;

					% lfp_filename = '135_CH42.continuous'; % AH2 

					lfp_cut_off_vector = [12 12 12 12 12 12; 
										  10 10 10 10 10 10;
										  1 10 10 10 10 10;
										  12 12 12 12 12 12;	
										  	];	


				case 'AH3'
					trial_folders = ["preSleep", "1Square", "2Circle", "postSleep"]; % AH3 
					output_folder_string = 	'csd_manifold_awake_sleep';
					
					% trial_folders = [ "1Square", "2Circle" ]; % AH3
					% output_folder_string = 	'all_trials_power_manifold_awake';
					% only_awake = 1;

					% trial_folders = ["preSleep", "postSleep"]; % AH3  
					% output_folder_string = 	'all_trials_power_manifold_sleep';
					% only_sleep = 1;

					lfp_filename = '159_CH46.continuous'; % AH3 

					lfp_cut_off_vector = [12 12 12 12 12 12; 
										  5 6 6 5 10 10;
										  0.5 5 5 5 10 10;
										  12 12 12 12 12 12;	
										  	];	
					

				case 'AH4'	
					trial_folders = ["preSleep", "1Circle", "2Square", "postSleep"]; % AH4 
					output_folder_string = 'csd_manifold_awake_sleep';

					% trial_folders = [ "1Circle", "2Square" ]; % AH4 
					% output_folder_string = 'all_trials_power_manifold_awake';
					% only_awake = 1;

					% trial_folders = ["preSleep", "postSleep"]; % AH4 
					% output_folder_string = 	'all_trials_power_manifold_sleep';
					% only_sleep = 1;

					lfp_filename = '126_CH46.continuous'; % AH4

					lfp_cut_off_vector = [12 12 12 12 12 12; 
										  2 10 10 2 10 2;
										  8 8 8 5 10 1.5;
										  12 12 12 12 12 12;	
										  	];	


				case 'AH5'	
					trial_folders = ["preSleep", "1Circle", "2Square", "postSleep"]; % AH5 
					output_folder_string = 'csd_manifold_awake_sleep';

					% trial_folders = ["1Circle", "2Square"]; % AH5 
					% output_folder_string = 'all_trials_power_manifold_awake';
					% only_awake = 1;

					% trial_folders = ["preSleep", "postSleep"]; % AH5
					% output_folder_string = 	'all_trials_power_manifold_sleep';
					% only_sleep = 1;
					
					lfp_filename = '101_CH42.continuous'; % AH5

					lfp_cut_off_vector = [12 12 12 12 12 12; 
										  1.5 10 10 2 10 10;
										  5 8 8 1 6 8;
										  12 12 12 12 12 12;	
										  	];
					

				otherwise
					fprintf('Incorrect animal id\n')
					return
			end % End switch animal_id

			output_folder = fullfile(root_directory,animal_id,output_folder_string);

			global_output_folder = output_folder;
			mkdir(global_output_folder);


			total_trials = length(trial_folders);


		% -------------------------------------------------------------------
		% Storage Variables
		% -------------------------------------------------------------------		
			global_csd_states = [];

			global_lfp_states = [];

			% Stores no. of states in each trials
			% Used for reading/reloading saved states 
			global_states_count = [0];


		% -------------------------------------------------------------------
		% Process each trial
		% -------------------------------------------------------------------	
			for trial_iterator = 1:total_trials

				current_trial = char(trial_folders(trial_iterator));

				fprintf('Processing: %s %s...\n', animal_id, current_trial);		

				% -------------------------------------------------------------------
				% Load Trial specific CSD/LFP files
				% -------------------------------------------------------------------	

					try
						csd_filepath = fullfile(root_directory, animal_id, csd_root_dirname, current_trial, csd_filename);
						lfp_filepath = fullfile(root_directory, animal_id, csd_root_dirname, current_trial, lfp_filename);
						
						load(csd_filepath)
						load(lfp_filepath)

					catch
						fprintf('Files not found\n')
						return
					end

				% -------------------------------------------------------------------
				% Common noise elimination
				% -------------------------------------------------------------------

					common_noise_lfp = mean(Raw_CSD);

					lfp_layers = Raw_CSD(csd_layer_indices,:);

					for layer_iterator = 1:total_layers

						current_layer = char(layer_strings(layer_iterator));

						temp_signal = [];
						temp_signal = lfp_layers(layer_iterator,1:end) - common_noise_lfp(1:end);

						figure(2)
						subplot(4,1,layer_iterator)
						plot(lfp_layers(layer_iterator,1:end))
						title(sprintf('%s %s %s - Raw CSD' ,animal_id, current_trial, current_layer))
						ylim_vals = [min(lfp_layers(layer_iterator,1:end)) max(lfp_layers(layer_iterator,1:end))];
						% ylim(ylim_vals)

						figure(1)
						subplot(4,1,layer_iterator)
						plot(temp_signal)
						title(sprintf('%s %s %s - CSD' ,animal_id, current_trial, current_layer))
						ylim(ylim_vals)
						

					end % End for layer_iterator

						figure(1)
						subplot(4,1,4)
						plot(common_noise_lfp(1:end),'r')
						title(sprintf('%s %s - Common Noise' ,animal_id, current_trial))



					return	

			

			

			end % End for trial_iterator	

		
		return

	end % End for animal_iterator

toc;
