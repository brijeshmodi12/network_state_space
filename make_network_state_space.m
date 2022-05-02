% ===================================================================
	% State space for current source density and local field potentials
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

% ===================================================================

clc; clear all ; close all; tic;

fprintf('------------------------------------------------\n')
fprintf('CSD / LFP State Space\n')
fprintf('------------------------------------------------\n')

% -------------------------------------------------------------------
% Plot Variables
% -------------------------------------------------------------------
	% Set the following variables to '1' for plotting/computing
	% Plotting variables
		plot_raw_bands = 0;
		plot_raw_lfp_csd = 0;
		plot_raw_acceleration = 0;
		plot_swr_with_layers = 0;

		plot_global_state_space = 0;

		plot_subspace_occupation = 0;
		plot_subspace_density_csd = 0;
		plot_subspace_density_lfp = 0;
		plot_subspace_occupation_time_control = 1;

		

		plot_acceleration_overlay = 0;
		
		feature_overlay_bins_lfp = 0;
		feature_overlay_bins_csd = 0;

		variance_overlay_bins_lfp = 0;
		rem_nonrem_overlay_bins_lfp = 0;

		make_video_csd = 0;
		make_video_lfp = 0;

		correlate_csd = 0;
		correlate_lfp = 0;

		correlate_csd_in_bins = 0;
		correlate_lfp_in_bins = 0;

		area_coverage_lfp = 0;
		area_coverage_csd = 0;


		% Outdated
			quantify_movements_lfp = 0;
			quantify_movements_csd = 0;

			rem_overlay_csd = 0;
			nonrem_overlay_csd = 0;

			rem_overlay_lfp = 0;
			nonrem_overlay_lfp = 0;

			correlate_lfp_time = 0;
			correlate_csd_time = 0;

			lfp_overlay_across_trials = 0;
			csd_overlay_across_trials = 0;

			feature_overlay_csd = 0;
			feature_overlay_lfp = 0;


		% Construct State Space with following settings
			do_umap = 0;

			do_pca = 0;

			do_isomap = 0;

			do_filtering = 0;

			do_rescale_global = 0;

			do_relative_bands = 0;

			do_construction = 1;

			remove_extreme_bins = 1;

			make_pdf_images = 0;

			make_eps_images = 00;

			do_hilbert = 0;


		% Load pre-computed state space data
		load_precomputed_data = 01;
	


% -------------------------------------------------------------------
% General Constants
% -------------------------------------------------------------------	
	% Paths and Filenames
		% root_directory = 'C:\Users\Matteo\Desktop\Brijesh\Open Field Dataset\';
		root_directory = 'D:\matteo_datasets\Open Field Datasets\';
		% root_directory = 'D:\Matteo Datasets\Open Field Datasets\';

		animal_folders_list = ["AH2","AH3","AH4","AH5"];

		csd_root_dirname = 'CSD';

		csd_filename = 'C_Raw_CSD.mat';
		lfp_filename = 'C_Raw_LFP.mat';
		spikedata_filename = 'spike_times.mat';

		total_folders = length(animal_folders_list);

		color_gray = [134, 138, 145 ] ./ 255;
		color_pink = [189, 11, 73] ./ 255;
		color_hist =  [134, 138, 145 ] ./ 255;	
		rem_color = [227, 14, 85] / 255;
		nonrem_color = [26, 130, 150] / 255;
			

		% Combine all the bands in an single array so you can apply loops 
		bands_strings = ["Delta", "Theta", "Spindle", "Slow gamma", "Medium gamma", "Fast gamma"];

		layer_strings = ["Pyramidal", "Radiatum", "SLM"];
		layer_letters =  ["P", "R", "S"];

		% Custom Frequency Bands
		% Frequency range for delta (in Hz)
		delta_band = [ 1 5 ];

		% Frequency range for theta (in Hz)
		theta_band = [ 6 10 ];

		% Frequency range for spindle (in Hz)
		spindle_band = [ 10 20 ];			

		% Frequency range for slow gamma (in Hz)
		sgamma_band = [ 20 45 ];

		% Frequency range for medium gamma (in Hz)
		mgamma_band = [ 60 90 ];

		% Frequency range for fast gamma (in Hz) 
		fgamma_band = [ 100 200 ];
	

		% Combine all the bands in an single array so you can apply loops 
		bands_strings = ["Delta", "Theta", "Spindle", "Slow gamma","Medium gamma","Fast gamma"];
		bands_array = [delta_band; theta_band; spindle_band; sgamma_band; mgamma_band ; fgamma_band];

		bands_strings2 = ["Delta", "Theta", "Spindle", "SlowG", "MediumG", "FastG"];

		delta_index = 1;
		theta_index = 2;
		spindle_index = 3;



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

		total_layers = 3;

		total_features = total_bands * total_layers;

		sf_ephys = 30000;

		downsampling_frequency = 1000;

		downsampling_factor = sf_ephys / downsampling_frequency;

		% Bin size for state space construction -  seconds
		bin_size_sec = 0.2;

		% Total samples in each bin
		bin_size_lfp = bin_size_sec * downsampling_frequency;

		% For smoothing
		% Gaussian kernel
		smoothing_std = 3;

		smoothing_width = smoothing_std;


		% For videos
		plot_speed_vector = [50,50,50,50];
		
		total_freq_bands = length(bands_strings);

		% Time interval between two measurements of area coverage
		% in seconds;
		time_interval_for_area_sec = 10;

		time_interval_for_correlation_evolution_sec = 10;

		time_interval_for_area = fix(time_interval_for_area_sec / bin_size_sec);
		time_interval_for_correlation_evolution = fix(time_interval_for_correlation_evolution_sec / bin_size_sec);
		

		% Make strings for layers
		for layer_iterator = 1:total_layers

			for band_iterator = 1:total_freq_bands

				ind =  (layer_iterator - 1)*total_freq_bands + band_iterator;

				x_axis_string = strcat(layer_letters(layer_iterator),'-',bands_strings2(band_iterator));

				x_ticks_array(ind) = x_axis_string;

				x_axis_string2 = strcat(layer_letters(layer_iterator),num2str(band_iterator));

				x_ticks_array2(ind) = x_axis_string2;

			end % End for band_iterator

		end % End for layer_iterator

		% variance threshold for accelerometer (to classify sleep / awake)
		acc_var_threshold = 0.005;

		sleep_acc_threshold = 0.01;

		% Sorts all delta bands together, all theta bands together and so on.
		sorted_band_indices = reshape(1:18, total_freq_bands, []);
		sorted_band_indices = reshape(sorted_band_indices', [], 1 ) ;


		band_layer_struct.band_strings = bands_strings;
		band_layer_struct.band_strings_short = bands_strings2;
		band_layer_struct.frequency = bands_array;
		band_layer_struct.layer_strings = layer_strings;
		band_layer_struct.layer_letters = layer_letters;
		band_layer_struct.sorted_band_indices = sorted_band_indices;
		band_layer_struct.x_ticks_array = x_ticks_array;
		band_layer_struct.x_ticks_array_short = x_ticks_array2;


		% For extreme bins ()
		extreme_bins_cutoff = 6;




% -------------------------------------------------------------------
% Core Loop
% -------------------------------------------------------------------
	% Process each animal
	for animal_iterator = 4:total_folders

		animal_id = char(animal_folders_list(animal_iterator));
		fprintf('------ Animal %s ------\n', animal_id)
		
		
		% -------------------------------------------------------------------
		% Animal specific Constants
		% -------------------------------------------------------------------		
			switch animal_id
				case 'AH2'	
					trial_folders = ["preSleep", "1Rectangle", "2Circle", "postSleep"]; % AH2 
					output_folder_string = 'network_state_space';


					lfp_cutoff_vector = [6 8 8 ;
										 5 4 5;
										 4 4 3;
										 7 7 7;	
										  	];
					% trials x layers					  	
				  	csd_cutoff_vector = [ 5 5 5;
				  						  3.5 3 5;
				  						  1 2 4;	
				  						  3 5 5;
				  	 ];

				  	% Indices to access pyramidal, radiatum and slm layers
					csd_layer_indices = [1 6 9];
					lfp_layer_indices = csd_layer_indices + 1;
					swr_id = 8;

					% accelerometer filename
					x_acc_filename = '135_AUX1.continuous';
					y_acc_filename = '135_AUX2.continuous';
					z_acc_filename = '135_AUX3.continuous';
				
					csd_xlim = [-15 15];
					csd_ylim = [-8 6];

					lfp_xlim = [-15 15];
					lfp_ylim = [-8 8];
					% x_acc_filename = '135_AUX7.continuous';
					% y_acc_filename = '135_AUX8.continuous';
					% z_acc_filename = '135_AUX9.continuous';


				case 'AH3'	
					trial_folders = ["preSleep", "1Square", "2Circle", "postSleep"]; % AH3 
					output_folder_string = 	'network_state_space';
				

					lfp_cutoff_vector = [7 7 7;
										 5 5 5;
										 4 4 4;
										 7 7 7;	
										  	];
					% trials x layers					  	
				  	csd_cutoff_vector = [ 10 5 5;
				  						  4.5 4 5;
				  						  2.5 4 4;	
				  						  7 5 5;
				  	 ];

				  	% Indices to access pyramidal, radiatum and slm layers
					csd_layer_indices = [1 5 9];
					lfp_layer_indices = csd_layer_indices + 1;
					swr_id = 12;

					csd_xlim = [-15 15];
					csd_ylim = [-8 12];

					lfp_xlim = [-15 10];
					lfp_ylim = [-8 6];

					x_acc_filename = '159_AUX1.continuous';
					y_acc_filename = '159_AUX2.continuous';
					z_acc_filename = '159_AUX3.continuous';
			

				case 'AH4'	
					trial_folders = ["preSleep", "1Circle", "2Square", "postSleep"]; % AH4 
					output_folder_string = 'network_state_space';


					lfp_cutoff_vector = [8 8 8;
										 5 4 5;
										 4.5 4.5 4;
										 8 8 8;	
										  	];
					% trials x layers					  	
				  	csd_cutoff_vector = [ 8 8 8;
				  						  4 4 5;
				  						  4.5 4 4;	
				  						  9 9 9 ;
				  	 ];


					% Indices to access pyramidal, radiatum and slm layers
					csd_layer_indices = [1 5 9];
					lfp_layer_indices = csd_layer_indices + 1;	
					swr_id = 32;

					csd_xlim = [-10 15];
					csd_ylim = [-15 15];

					lfp_xlim = [-15 15];
					lfp_ylim = [-10 10];

					x_acc_filename = '126_AUX1.continuous';
					y_acc_filename = '126_AUX2.continuous';
					z_acc_filename = '126_AUX3.continuous';									  	


				case 'AH5'	
					trial_folders = ["preSleep", "1Circle", "2Square", "postSleep"]; % AH5 
					output_folder_string = 'network_state_space';

				

					lfp_cutoff_vector = [8 8 8;
										 5 4 5;
										 4.5 4.5 4.5;
										 8 8 8;	
										  	];
					% trials x layers					  	
				  	csd_cutoff_vector = [ 8 9 8;
				  						  4 1.5 4;
				  						  4.5 4 4;	
				  						  7 9 7 ;
				  	 ];

					% Indices to access pyramidal, radiatum and slm layers
					csd_layer_indices =  [2 5 10];
					lfp_layer_indices = csd_layer_indices + 1;	
					swr_id = 2;

					csd_xlim = [-10 10];
					csd_ylim = [-8 12];

					lfp_xlim = [-15 10];
					lfp_ylim = [-8 8];

					x_acc_filename = '101_AUX1.continuous';
					y_acc_filename = '101_AUX2.continuous';
					z_acc_filename = '101_AUX3.continuous';
								

				otherwise
					fprintf('Incorrect animal id\n')
					return
			end % End switch animal_id

			total_trials = length(trial_folders);

			% For correlation of lfp bands in bins (for cleaning out noisy bins)
			% Manually specified by visual inspection of state space
			awake_bin_id_lfp = {[2 3 5 6 8 9], [7 8 9], [7 8 9], [1 2 3]};	

			sleep_bin_id_lfp = {[1:9], [ 2 3 5 6 7 8 9], [1:9], [1:9] };


			% awake_bin_id_csd = {[5 8 9], [2 3 5], [2 4 5], [1 2 3 4 7 8 9 ]};	

			% sleep_bin_id_csd = {[1 2 4 5 7 8], [1:8], [2 3 4 5 6 8 9], [1 2 3 4 7 8 9] };


		% -------------------------------------------------------------------
		% Make Output Folder based on plotting parameters
		% -------------------------------------------------------------------		
			output_folder_string = strcat(output_folder_string, '_bin' ,num2str(bin_size_lfp));
	
			if do_pca
				output_folder_string = strcat(output_folder_string,'_pca');
			end
		
			if do_isomap
				output_folder_string = strcat(output_folder_string,'_isomap');
			end

			if do_rescale_global
				output_folder_string = strcat(output_folder_string,'_global_rescale');
			end

			if do_relative_bands
				output_folder_string = strcat(output_folder_string,'_relative_bands');
			end

			if do_hilbert
				output_folder_string = strcat(output_folder_string,'_hilbert');
			end % End if do_hilbert

			output_folder = fullfile(root_directory,animal_id,output_folder_string);

			global_output_folder = output_folder;
			mkdir(global_output_folder);


		% -------------------------------------------------------------------
		% Storage Variables
		% -------------------------------------------------------------------		
			global_csd_states = [];
			global_lfp_states = [];

			% Stores no. of states in each trials
			% Used for reading/reloading saved states 
			global_states_count_csd = [0];
			global_states_count_lfp = [0];


			global_median_power_in_bands_across_trials_csd = [];
			global_median_power_in_bands_across_trials_lfp = [];

			global_regime_median_power_csd = [];
			global_regime_median_power_lfp = [];

			global_artifacts_bins_csd = {};
			global_artifacts_bins_lfp = {};

			global_acceleration_data_lfp = [];
			global_acceleration_data_csd = [];

			global_binned_spikes_csd = [];
			global_binned_spikes_lfp = [];

			global_lfp_for_csd = [];

			global_raw_bin_indices_lfp = [];
			global_raw_bin_indices_csd = [];
		

		% -------------------------------------------------------------------
		% Construction of State Space
		% ------------------------------------------------------------------		
			if size(global_csd_states,1) == load_precomputed_data	

				if do_construction
					
					% -------------------------------------------------------------------
					% Process each trial and Construct State Space
					% -------------------------------------------------------------------	
					for trial_iterator = 1:total_trials

						current_trial = char(trial_folders(trial_iterator));

						fprintf('Processing: %s %s...\n', animal_id, current_trial);		

						% -------------------------------------------------------------------
						% Load Trial specific CSD/LFP files
						% -------------------------------------------------------------------	
							fprintf('--Loading CSD and LFP\n');			

							try
								csd_filepath = fullfile(root_directory, animal_id, csd_root_dirname, current_trial, csd_filename);
								lfp_filepath = fullfile(root_directory, animal_id, csd_root_dirname, current_trial, lfp_filename);
								
								load(csd_filepath)
								load(lfp_filepath)

							catch
								fprintf('Files not found\n')
								return
							end

							trial_sample_count = size(Raw_CSD,2);

							fprintf('--Loading Spikes\n');		

							try
								spikedata_filepath = fullfile(root_directory, animal_id, current_trial, spikedata_filename);
								load(spikedata_filepath)
							catch
								fprintf('Spike data not found\n')
								return
							end


						% -------------------------------------------------------------------
						% Bin Spiketrain for each cell
						% -------------------------------------------------------------------
							total_cells = length(spike_times);

							% total bins for spike train binning
							total_bins = ceil(trial_sample_count / bin_size_lfp);

							% Edges for binning spikes
							bin_edges = linspace(1 , trial_sample_count, total_bins);

							bin_start = bin_edges(1:end-1);
							
							bin_end = bin_edges(2:end);

							bin_edges = [bin_start' bin_end'];

							% in seconds
							trial_duration = trial_sample_count / downsampling_frequency;
			
							global_trial_duration(trial_iterator) = trial_duration;
					
						
							% Pre process spike timings	
								cells_list = spike_times;

								% Adjusting occurence of spike time on lfp samples
								for cell_iterator = 1:total_cells
									cells_list{cell_iterator} = floor(cells_list{cell_iterator} * downsampling_frequency);
								end

						
							% Bin spike times for each cell
					
								trial_binned_spikes = [];

								trial_binned_spikes = create_pv_spiketrain_raw(cells_list, bin_edges);


										
						% -------------------------------------------------------------------
						% Load acceleration files and bin it
						% -------------------------------------------------------------------	

							if plot_raw_acceleration
								acceleration_output_folder = fullfile(global_output_folder,'raw_acceleration');
								mkdir(acceleration_output_folder);
							end

							fprintf('--Loading Accelerometer Data\n');			
							current_trial_path = fullfile(root_directory, animal_id, current_trial);

							acc_matrix = []; binned_acceleration = []; 	mean_acceleration = [];

							x_acc_filepath = fullfile(current_trial_path, x_acc_filename);
							[temp_var, ~, ~] = load_open_ephys_data_faster(x_acc_filepath);
							acc_matrix(1,:) = zscore(downsample(temp_var, downsampling_factor));
							clear temp_var
							
							y_acc_filepath = fullfile(current_trial_path, y_acc_filename);
							[temp_var, ~, ~] = load_open_ephys_data_faster(y_acc_filepath);
							acc_matrix(2,:) = zscore(downsample(temp_var, downsampling_factor));
							clear temp_var

							z_acc_filepath = fullfile(current_trial_path, z_acc_filename);
							[temp_var, ~, ~] = load_open_ephys_data_faster(z_acc_filepath);
							acc_matrix(3,:) = zscore(downsample(temp_var, downsampling_factor));
							clear temp_var
							
						
							mean_acceleration = mean(acc_matrix);

							% mean_acceleration = acc_matrix;

							if numel(mean_acceleration) ~= trial_sample_count
								fprintf('Accelerometer Data Inconsistent\n')
								return
							end
								
							remainder = mod(trial_sample_count, bin_size_lfp);

							adjusted_length = trial_sample_count - remainder;

							
							temp = []; sleep_bin_indices = []; awake_bin_indices = [];
							temp = reshape(mean_acceleration(1:adjusted_length), bin_size_lfp, [] );
							binned_acceleration = mean(temp);
							binned_acceleration = smooth_gaussian(binned_acceleration, smoothing_std, smoothing_width);

							binned_acceleration = reshape(binned_acceleration,[],1);

							movvar_binned_acc = movvar(binned_acceleration,100);

							movvar_binned_acc = smooth_gaussian(movvar_binned_acc,20,40);

							sleep_bin_indices = find(movvar_binned_acc < sleep_acc_threshold);

							awake_bin_indices = find(movvar_binned_acc >= sleep_acc_threshold);
														
							% tt = movvar_binned_acc < sleep_acc_threshold;

							% figure
							% subplot(3,1,1)
							% plot(binned_acceleration)
							% subplot(3,1,2)
							% plot(movvar_binned_acc)
							% subplot(3,1,3)
							% plot(tt)
							% ylim([0 2])
							% title('median')
							% return
							
							
							if plot_raw_acceleration	

								acc_colors = viridis(3);

								for acc_iterator = 1:3	

									raw_acc_fig = 1;
									figure(raw_acc_fig)
									subplot(3,1,acc_iterator)
									plot(acc_matrix(acc_iterator,:), 'Color', acc_colors(acc_iterator,:));
									title(sprintf('Accelerometer Axis %d',acc_iterator));

								end % End for acc_iterator

								figure(raw_acc_fig)
								suptitle(sprintf('%s %s', animal_id, current_trial))
								filename = sprintf('raw_acceleration_%s.png',current_trial);
								fullname = fullfile(acceleration_output_folder,filename);
								saveas(gcf,fullname);
								close(figure(raw_acc_fig))

							end % End if plot_raw_acceleration	

							clear acc_matrix;	
								


						% -------------------------------------------------------------------
						% Filter CSD/LFP into bands
						% -------------------------------------------------------------------
							if do_filtering	
								raw_output_folder = fullfile(global_output_folder,'raw_traces');
								mkdir(raw_output_folder);
								% Initialize to zeros
								% filtered_lfp = zeros(total_layers*total_bands, trial_sample_count);
								% filtered_csd = filtered_lfp;

								filtered_lfp = [];
								filtered_csd = [];

								lfp_artifacts = []; artifacts_indices_lfp = []; 
								csd_artifacts = []; artifacts_indices_csd = [];
								
								lfp_pyramidal = [];

								for layer_iterator = 1:total_layers
									tc = viridis(total_layers);

									fprintf('--Processing Layer: %d\n', layer_iterator)

									current_layer = char(layer_strings(layer_iterator));

									csd_ind = csd_layer_indices(layer_iterator);
									current_csd = Raw_CSD(csd_ind,:);

									lfp_ind = lfp_layer_indices(layer_iterator);
									current_lfp = Raw_LFP(lfp_ind,:);

									if layer_iterator == 1
										lfp_pyramidal = current_lfp;
									end

									cut_off_csd = csd_cutoff_vector(trial_iterator,layer_iterator)*mad(current_csd);
									txcsd = 1:length(current_csd);
									tycsd = cut_off_csd * ones(1,length(current_csd));	

									cut_off_lfp = lfp_cutoff_vector(trial_iterator,layer_iterator)*mad(current_lfp);
									txlfp = 1:length(current_lfp);
									tylfp = cut_off_lfp * ones(1,length(current_lfp));	

									lfp_artifacts(:,layer_iterator) = abs(current_lfp) >= cut_off_lfp;
									csd_artifacts(:,layer_iterator) = abs(current_csd) >= cut_off_csd;

									artifacts_indices_lfp = find(lfp_artifacts);
									artifacts_indices_csd = find(csd_artifacts);

									clean_lfp = current_lfp;
									clean_lfp(artifacts_indices_lfp)  = 0;

									clean_csd = current_csd;
									clean_csd(artifacts_indices_csd)  = 0;

									if plot_raw_lfp_csd
										figure(1)
										set(gcf, 'Position', get(0, 'Screensize'));
										subplot(total_layers,1,layer_iterator)
										plot(current_csd,'Color',tc(:,layer_iterator))

										hold on
										plot(txcsd, tycsd,'r' ,'lineWidth',2)
										hold on 
										plot(txcsd, -tycsd,'r' ,'lineWidth',2)
										ylim([-cut_off_csd*2 cut_off_csd*2])
										title(sprintf('%s %s %s Cutoff %d x m.a.d - CSD' ,animal_id, current_trial, current_layer, csd_cutoff_vector(trial_iterator,layer_iterator)))

										% figure(1)
										% % set(gcf, 'Position', get(0, 'Screensize'));
										% subplot(total_layers,1,layer_iterator+1)
										% plot(clean_csd(1:100*downsampling_frequency),'Color',tc(:,layer_iterator))	
										% ylim([-cut_off_csd*2 cut_off_csd*2])
										% title('clean csd')
										% return
										
										figure(2)
										set(gcf, 'Position', get(0, 'Screensize'));
										subplot(total_layers,1,layer_iterator)
										plot(current_lfp,'Color',tc(:,layer_iterator))

										hold on
										plot(txlfp, tylfp,'r' ,'lineWidth',2)
										hold on 
										plot(txlfp, -tylfp,'r' ,'lineWidth',2)
										ylim([-cut_off_lfp*2 cut_off_lfp*2])
										title(sprintf('%s %s %s Cutoff %d x m.a.d - LFP' ,animal_id, current_trial, current_layer, lfp_cutoff_vector(trial_iterator,layer_iterator)))		

									end % End if plot_raw_lfp_csd
									

									templfp = [];
									tempcsd = [];

									mad_lfp = mad(current_lfp);
									mad_csd = mad(current_csd);

									fprintf('  --Filtering Bands\n')

									for band_iterator = 1:total_bands

										band_index = (layer_iterator - 1) * total_bands + band_iterator;

										% Fiter LFP
											temp_filtered_lfp = [];
											temp_filtered_lfp = eegfilt(current_lfp, downsampling_frequency, bands_array(band_iterator,1), bands_array(band_iterator,2));

											if do_relative_bands
												templfp(band_iterator,:) =  zscore_mod(temp_filtered_lfp, mad_lfp);	
											else
												templfp(band_iterator,:) = temp_filtered_lfp;
											end	


									

										% Fiter CSD
											temp_filtered_csd = [];
											temp_filtered_csd = eegfilt(current_csd, downsampling_frequency, bands_array(band_iterator,1), bands_array(band_iterator,2));

											if do_relative_bands
												tempcsd(band_iterator,:) =  zscore_mod(temp_filtered_csd, mad_csd);	
											else
												tempcsd(band_iterator,:) = temp_filtered_csd;
											end	
												

									end % End For band_iterator


									filtered_lfp = [filtered_lfp; templfp];
									filtered_csd = [filtered_csd; tempcsd];	

												
								end % End for layer_iterator

								fprintf('--Saving Filtered Data\n')
								filtered_data_fname = sprintf('filtered_data.mat');
								filtered_data_fpath = fullfile(root_directory, animal_id, current_trial, filtered_data_fname);
								save(filtered_data_fpath, 'filtered_lfp', 'filtered_csd','lfp_artifacts','csd_artifacts')

								if plot_raw_lfp_csd
									filename = sprintf('%s %s - Raw CSD.png', animal_id, current_trial);	
									fullname = fullfile(raw_output_folder,filename);
									saveas(figure(1),fullname);

									filename = sprintf('%s %s - Raw LFP.png', animal_id, current_trial);	
									fullname = fullfile(raw_output_folder,filename);
									saveas(figure(2),fullname);

									close all
									
								end % End if plot_raw_lfp_csd

							else
								fprintf('--Loading Filtered Data...\n')
								filtered_data_fname = sprintf('filtered_data.mat');
								filtered_data_fpath = fullfile(root_directory, animal_id, current_trial, filtered_data_fname);
								load(filtered_data_fpath)

							end % End if do_filtering
				



						% -------------------------------------------------------------------
						% Detect SWRs and bin it presence vector
						% -------------------------------------------------------------------
							if plot_swr_with_layers
								all_layers_folder = fullfile(global_output_folder,'all_layers_plot');
								mkdir(all_layers_folder)
								fprintf('--Detecting SWR...\n')
							
								% Get Pyramidal LFP
								lfp_ind = lfp_layer_indices(1);
								lfp_pyramidal = Raw_LFP(lfp_ind,:);

								lfp_swr_band = filtered_lfp(6,:);

								[swr_event_indices, swr_presence_vector] = detect_swr(lfp_swr_band, lfp_pyramidal,0);
								% return
								swr_start = swr_event_indices(:,1);
								swr_end = swr_event_indices(:,2);

								swr_peak = swr_start + fix((swr_end - swr_start) / 2);

								total_peaks= min([100, length(swr_end)])
								collected_swr_traces_lfp = [];
								collected_swr_traces_csd = [];

								for peak_iterator = 1:total_peaks
									
									swr_plot_ind = swr_peak(peak_iterator) - 250 : swr_peak(peak_iterator) + 250;

									collected_swr_traces_lfp(:,:, peak_iterator) = Raw_LFP(:,swr_plot_ind);
									collected_swr_traces_csd(:,:, peak_iterator) = Raw_CSD(:,swr_plot_ind);

								end 

								average_profile_around_swr_lfp = mean(collected_swr_traces_lfp,3);
								average_profile_around_swr_csd = mean(collected_swr_traces_csd,3);

								xax = -250:250;	
							
								swr_plot_ind = swr_event_indices(swr_id,1) - 250 : swr_event_indices(swr_id,2) + 250;

								x0 = 10;
								y0 = 10;
								width1 = 500;
								height = 1200;

								% LFP Layers
									total_rows_lfp = size(Raw_LFP,1);
									total_cols = 1;

									layer_palette = repmat(color_gray,total_rows_lfp,1);
									layer_palette(lfp_layer_indices,:) = magma(3);
									% index = reshape(1:total_rows_lfp*total_cols, total_cols,total_rows_lfp).';
									temp = Raw_LFP(:, swr_plot_ind);

									yaxis_range = [min(temp(:)) max(temp(:)) ];
									yaxis_range = [-1000 1000 ];

									for row_iterator  = 1:total_rows_lfp

										layer_color = layer_palette(row_iterator,:);
										figure(985)
										% set(gcf, 'Position',  [100, 100, 500, 400])
										
										set(gcf,'position',[x0,y0,width1,height])
										subplot(total_rows_lfp+1,total_cols,row_iterator)
										temp_row = [];
										% temp_row = Raw_LFP(row_iterator, swr_plot_ind);
										temp_row = average_profile_around_swr_lfp(row_iterator,:);
										temp_row = temp_row-mean(temp_row);
										plot(xax, temp_row,'Color',layer_color,'lineWidth',1.2)
										set(gca, 'visible' ,'off')
										% ylim(yaxis_range)

										if row_iterator == total_rows_lfp
											subplot(total_rows_lfp+1,1,row_iterator+1)
											plot([0; 0], [500; -500], '-k',  [0; 100], [-500; -500], '-k', 'LineWidth', 1.5)
											xlim([0 length(temp_row)]);
											ylim(yaxis_range)
											set(gca, 'visible', 'off')
											suptitle(sprintf('%s LFP Layers', animal_id))
										end

										

										if make_eps_images	&  row_iterator == total_rows_lfp
											epsname = sprintf('%s_lfp_layers.eps', animal_id);
											fullname_eps = fullfile(all_layers_folder,epsname)
											set(gcf,'renderer','Painters')
											saveas(gcf,fullname_eps,'epsc')
										end


										filename = sprintf('%s_lfp_layers.png', animal_id );
										fullname = fullfile(all_layers_folder,filename);
										saveas(gcf,fullname);


									end % End for row_iterator
									
								

								% CSD Layers
									total_rows_csd = size(Raw_CSD,1);
									total_cols = 1;

									layer_palette = repmat(color_gray,total_rows_csd,1);
									layer_palette(csd_layer_indices,:) = magma(3);
									% index = reshape(1:total_rows_csd*total_cols, total_cols,total_rows_csd).';
									temp = Raw_CSD(:, swr_plot_ind);

									yaxis_range = [min(temp(:)) max(temp(:)) ];
									yaxis_range = [-1000 1000 ];

									for row_iterator  = 1:total_rows_csd

										layer_color = layer_palette(row_iterator,:);
										figure(985123)
										% set(gcf, 'Position',  [100, 100, 500, 400])
										
										set(gcf,'position',[x0,y0,width1,height])
										subplot(total_rows_csd+1,total_cols,row_iterator)
										temp_row = [];
										% temp_row = Raw_CSD(row_iterator, swr_plot_ind);
										temp_row = average_profile_around_swr_csd(row_iterator,:);

										temp_row = temp_row-mean(temp_row);
										plot(temp_row,'Color',layer_color,'lineWidth',1.2)
										set(gca, 'visible' ,'off')
										% ylim(yaxis_range)

										if row_iterator == total_rows_csd
											subplot(total_rows_csd+1,1,row_iterator+1)
											plot([0; 0], [500; -500], '-k',  [0; 100], [-500; -500], '-k', 'LineWidth', 1.5)
											xlim([0 length(temp_row)]);
											ylim(yaxis_range)
											set(gca, 'visible', 'off')
											suptitle(sprintf('%s CSD Layers', animal_id))
										end

										

										if make_eps_images	&  row_iterator == total_rows_csd
											epsname = sprintf('%s_csd_layers.eps', animal_id);
											fullname_eps = fullfile(all_layers_folder,epsname)
											set(gcf,'renderer','Painters')
											saveas(gcf,fullname_eps,'epsc')
										end



										filename = sprintf('%s_csd_layers.png', animal_id );
										fullname = fullfile(all_layers_folder,filename);
										saveas(gcf,fullname);


									end % End for row_iterator
						




							return
							end % End if plot_swr_with_layers



						% -------------------------------------------------------------------
						% Trial Storage Variables
						% -------------------------------------------------------------------	
							% To store state space of each trial
							trial_csd_state_space = [];
							trial_lfp_state_space = [];
							trial_lfp_for_csd_state_space = [];

							% Raw state space, without any smoothing or modification
							raw_csd_state_space = [];
							raw_lfp_state_space = [];
							
							artifacts_bins_indices_csd = [];
							artifacts_bins_indices_lfp = [];

							raw_indices = [];



						% -------------------------------------------------------------------
						% Detect Artifacts bins
						% -------------------------------------------------------------------
							fprintf('--Detecting Artifacts\n')

							% Sum Binary artifacts vectors 
							artifacts_across_layers_csd = sum(csd_artifacts,2);
							artifacts_across_layers_lfp = sum(lfp_artifacts,2);

							% bin artifacts vectors as you do with LFP/CSD
							remainder = mod(trial_sample_count, bin_size_lfp);
							adjusted_length = trial_sample_count - remainder;
							
							
							temp = []; binned_artifacts_csd = []; binned_artifacts_lfp = [];

							% All the elements in a bin are now adjusted across a column	
							temp = reshape(artifacts_across_layers_csd(1:adjusted_length), bin_size_lfp , []);	
							binned_artifacts_csd = sum(temp);

							temp = [];
							temp = reshape(artifacts_across_layers_lfp(1:adjusted_length), bin_size_lfp , []);	
							binned_artifacts_lfp = sum(temp);


							artifacts_bins_indices_csd = find(binned_artifacts_csd);
							artifacts_bins_indices_lfp = find(binned_artifacts_lfp);

							to_remove_bins_csd = [];
							to_remove_bins_lfp = [];


							if ismember(trial_iterator,[2,3])
								to_remove_bins_lfp = artifacts_bins_indices_lfp;
								to_remove_bins_csd = artifacts_bins_indices_csd;
							else
								to_remove_bins_lfp = union(artifacts_bins_indices_lfp, awake_bin_indices);
								to_remove_bins_csd = union(artifacts_bins_indices_csd, awake_bin_indices);
							end



						% -------------------------------------------------------------------
						% Compute CSD/LFP states from bands
						% -------------------------------------------------------------------
							total_features = size(filtered_csd,1);

							% Where bin power exceeds above 6 std for an oscillator.
							extreme_bins_artifacts_lfp = [];
							extreme_bins_artifacts_csd = [];

							for feature_iterator = 1:total_features
								% ------------------------------
								% CSD
								% ------------------------------
									current_csd_band = [];
									current_csd_band = filtered_csd(feature_iterator,:);

									current_csd_band = current_csd_band .^2;

									remainder = mod(trial_sample_count, bin_size_lfp);
									adjusted_length = trial_sample_count - remainder;
									
									% All the elements in a bin are now adjusted across a column
									temp = []; binned_csd = [];
									temp = reshape(current_csd_band(1:adjusted_length), bin_size_lfp , []);	
									binned_csd = median(temp);
									
									% Remove Artifacts bins from further processing pipeline
									binned_csd(to_remove_bins_csd) = [];
																	
									raw_csd_state_space = [ raw_csd_state_space  binned_csd'];

									smoothed_csd = smooth_gaussian(binned_csd, smoothing_std, smoothing_width);

									trial_csd_state_space(:,feature_iterator) = smoothed_csd';

							    % ------------------------------
								% LFP
								% ------------------------------
									current_lfp_band = [];
									current_lfp_band = filtered_lfp(feature_iterator,:);

									% if do_hilbert
									% current_lfp_power_hilbert =  abs(hilbert(current_lfp_band)).^2;	
									% else
									current_lfp_power = current_lfp_band .^ 2;
									% end


									% figure
									% subplot(2,1,1)
									% plot(current_lfp_power_hilbert(1:10000))
									% title('Hilbert')
									% subplot(2,1,2)
									% plot(current_lfp_power(1:10000))

									% return
									
									% All the elements in a bin are now adjusted across a column
									temp = []; binned_lfp = []; binned_lfp_for_csd = [];
									temp = reshape(current_lfp_power(1:adjusted_length), bin_size_lfp , []);	
									binned_lfp = median(temp);

									binned_lfp_for_csd = binned_lfp;

									% Remove Artifacts bins from further processing pipeline
									binned_lfp(to_remove_bins_lfp) = [];

									
									
									raw_lfp_state_space = [ raw_lfp_state_space  binned_lfp'];

									smoothed_power = smooth_gaussian(binned_lfp, smoothing_std, smoothing_width);

									trial_lfp_state_space(:, feature_iterator) = smoothed_power';

									
									binned_lfp_for_csd(to_remove_bins_csd) = [];
									smoothed_power2 = smooth_gaussian(binned_lfp_for_csd, smoothing_std, smoothing_width);
									trial_lfp_for_csd_state_space(:, feature_iterator) = smoothed_power2';


								% ------------------------------
								% Raw samples indices
								% ------------------------------
									if feature_iterator == 1
										raw_indices_temp = [];

										raw_indices_temp = 1:adjusted_length;

										raw_indices_temp = reshape(raw_indices_temp, bin_size_lfp , []);

										raw_indices_start = raw_indices_temp(1,:);
										raw_indices_end = raw_indices_temp(bin_size_lfp,:);

										raw_indices = [raw_indices_start' raw_indices_end'];
										
									end  % End if feature_iterator


							end % End for feature_iterator

							clear filtered_lfp filtered_csd



						% -------------------------------------------------------------------
						% Remove Artifacts bins from further processing pipeline
						% -------------------------------------------------------------------
												
							% Accelerometer Data
								acc_for_lfp = binned_acceleration;
								acc_for_lfp(to_remove_bins_lfp) = [];

								acc_for_csd = binned_acceleration;
								acc_for_csd(to_remove_bins_csd) = [];

								global_acceleration_data_lfp = [global_acceleration_data_lfp; acc_for_lfp ];
								global_acceleration_data_csd = [global_acceleration_data_csd; acc_for_csd ];

							% Network State Space Datasets
								global_csd_states = [global_csd_states ; trial_csd_state_space];

								global_lfp_states = [global_lfp_states; trial_lfp_state_space];

								global_states_count_csd = [global_states_count_csd size(global_csd_states,1) ];
								global_states_count_lfp = [global_states_count_lfp size(global_lfp_states,1) ];

								global_artifacts_bins_csd{1,trial_iterator} = to_remove_bins_csd;
								global_artifacts_bins_lfp{1,trial_iterator} = to_remove_bins_lfp;


							% Binned Spike Trains
								trial_pv_csd = trial_binned_spikes;
								trial_pv_csd(to_remove_bins_csd,:) = [];

								trial_pv_lfp = trial_binned_spikes;
								trial_pv_lfp(to_remove_bins_lfp,:) = [];	

								global_binned_spikes_csd = [global_binned_spikes_csd; trial_pv_csd ];
								global_binned_spikes_lfp = [global_binned_spikes_lfp; trial_pv_lfp ] ;

							% LFP bins for CSD Indices
								global_lfp_for_csd = [global_lfp_for_csd; trial_lfp_for_csd_state_space];

							% Raw Indices

								temp_ri_csd = raw_indices;
								temp_ri_lfp = raw_indices;

								temp_ri_csd(to_remove_bins_csd,:) = [];
								temp_ri_lfp(to_remove_bins_lfp,:) = [];

								global_raw_bin_indices_lfp = [global_raw_bin_indices_lfp; temp_ri_lfp];
								global_raw_bin_indices_csd = [global_raw_bin_indices_csd; temp_ri_csd];


						fprintf('------------------------------\n')

					end % End for trial_iterator	

					% -------------------------------------------------------------------
					% Remove Extreme bins from all trials
					% -------------------------------------------------------------------
						if remove_extreme_bins
							band_power_folder = fullfile(global_output_folder,'band_power_distribution');
							mkdir(band_power_folder);

							fprintf('Removing Global Extreme Bins...\n')

							total_features = size(global_lfp_states,2);
							
							extreme_global_bins_csd = [];
							extreme_global_bins_lfp = [];
							band_color = magma(6);
							band_color = repmat(band_color, 3,1);

							for feature_iterator = 1:total_features

								temp_extreme_bins_csd = [];
								temp_extreme_bins_lfp = [];
								current_feature_csd = [];
								current_feature_lfp = []; 

								current_feature_csd = zscore(global_csd_states(:,feature_iterator));
								current_feature_lfp = zscore(global_lfp_states(:,feature_iterator));

								temp_extreme_bins_csd = find(current_feature_csd > extreme_bins_cutoff);
								temp_extreme_bins_lfp = find(current_feature_lfp > extreme_bins_cutoff);

								extreme_global_bins_csd = union(extreme_global_bins_csd, temp_extreme_bins_csd);
								extreme_global_bins_lfp = union(extreme_global_bins_lfp, temp_extreme_bins_lfp);

								hist_temp = [];
								figure(234)
								set(gcf, 'Position', get(0, 'Screensize'));
								subplot(3,6,feature_iterator)
								hist_temp = histogram(current_feature_lfp, 100, 'FaceColor', band_color(feature_iterator,:), 'EdgeColor', band_color(feature_iterator,:));					
								set(gca,'box','off') 
								yaxis_range = [min(hist_temp.Values) max(hist_temp.Values)];
								yaxis_range = ceil(yaxis_range/10)*10;% xaxis_range = ceil(caxis_range/10)*10;
								yticks(yaxis_range)
								yticklabels(yaxis_range)

								hist_temp = [];

								figure(235)
								set(gcf, 'Position', get(0, 'Screensize'));
								subplot(3,6,feature_iterator)
								hist_temp = histogram(current_feature_csd, 100, 'FaceColor', band_color(feature_iterator,:), 'EdgeColor', band_color(feature_iterator,:));					
								set(gca,'box','off') 
								yaxis_range = [min(hist_temp.Values) max(hist_temp.Values)];
								yaxis_range = ceil(yaxis_range/10)*10;% xaxis_range = ceil(caxis_range/10)*10;
								yticks(yaxis_range)
								yticklabels(yaxis_range)
								
							end % End for feature_iterator

							lfp_states_count_new = [0];
							csd_states_count_new = [0];

							% Adjust bin counts for each trial after removal of extreme bins
							for trial_iterator = 1:total_trials
								range_start = global_states_count_lfp(trial_iterator) + 1;
								range_end_lfp = global_states_count_lfp(trial_iterator + 1);
								range_end_csd = global_states_count_csd(trial_iterator + 1);

								extreme_count_upto_trial_lfp = 0;
								extreme_count_upto_trial_csd = 0;
								extreme_count_upto_trial_lfp =  numel(find(extreme_global_bins_lfp <= range_end_lfp));
								extreme_count_upto_trial_csd =  numel(find(extreme_global_bins_csd <= range_end_csd));

								lfp_states_count_new = [lfp_states_count_new  (range_end_lfp - extreme_count_upto_trial_lfp)	];
								csd_states_count_new = [csd_states_count_new  (range_end_csd - extreme_count_upto_trial_csd)	];

							end % End for trial_iterator

							global_states_count_lfp = lfp_states_count_new;
							global_states_count_csd = csd_states_count_new;


							% Remove extreme bins from all global variables
								% Acceleration
								global_acceleration_data_lfp(extreme_global_bins_lfp) = [];
								global_acceleration_data_csd(extreme_global_bins_csd)  = [];

								% Network State Space Datasets
								global_csd_states(extreme_global_bins_csd,:) = [];
								global_lfp_states(extreme_global_bins_lfp,:) = [];

								% Binned Spike Trains
								global_binned_spikes_csd(extreme_global_bins_csd,:) = [];
								global_binned_spikes_lfp(extreme_global_bins_lfp,:) = [];

								% LFP bins for CSD Indices
								global_lfp_for_csd(extreme_global_bins_csd,:) = [];

								% Raw Indices
								global_raw_bin_indices_lfp(extreme_global_bins_lfp,:) = [];
								global_raw_bin_indices_csd(extreme_global_bins_csd,:) = []; 




							figure(234)
							suptitle('Distribution of Median Band Power LFP')
							filename = sprintf('distribution_band_power_lfp.png');
							fullname = fullfile(band_power_folder,filename);
							saveas(gcf,fullname);

							figure(235)
							suptitle('Distribution of Median Band Power CSD')
							filename = sprintf('distribution_band_power_csd.png');
							fullname = fullfile(band_power_folder,filename);
							saveas(gcf,fullname);



	
						end % End if remove_extreme_bins


				else 
					fprintf('Construction Skipped...\n')
					load(fullfile(global_output_folder,'state_space_data.mat'))	

				end % End if do_construction



				if do_rescale_global
					popmin = min(global_csd_states);
					popmax = max(global_csd_states);
					% global_csd_states = rescale(global_csd_states, 0, 1, 'InputMin',popmin,'InputMax',popmax);
					global_csd_states = zscore(global_csd_states);
				end


				% -------------------------------------------------------------------
				%  Run Dimensionality Reduction on global state spaces -- CSD
				% -------------------------------------------------------------------	
					projection_csd = [];	
					
					if do_umap
						fprintf('Performing UMAP on CSD\n')
						projection_csd = [];

						[projection_csd, umap_params_csd] = run_umap(global_csd_states,  'n_components', 2, 'n_neighbors', 25 ,'min_dist', 0.1, 'metric', 'euclidean' );
						x = projection_csd(:,1);
						y = projection_csd(:,2);
						
					end % End if do_umap

					if do_pca
						fprintf('Performing PCA CSD\n')
						% projection_csd = pca scores;
						% same name for code consistency;
						[coeff_csd, scores_csd, latent_csd, tsquared_csd, variance_explained_csd, mu_csd] = pca(global_csd_states);
						projection_csd = scores_csd(:,1:2);
						x = projection_csd(:,1);
						y = projection_csd(:,2);

					end % End if do_pca


					if do_isomap

						fprintf('Performing ISOMAP on CSD\n')
						% pair-wise distance between population vectors
						pairwise_distance = pdist(global_lfp_states(:,:),'cosine');
						
						% Just put the results in matrix form and add some small noise to avoid coincident points
						pairwise_distance = squareform(pairwise_distance + rand(size(pairwise_distance))*0.0001 );
						
						isomap_type = 'k'; % Nearest neighbour
						% % % isomap_type = 'epsilon'; % radial distance

						% % must be integer if isomap_type = k
						distance_value = 20;

						options.dims = [2];

						[isomap_output] = IsoMap(pairwise_distance, isomap_type, distance_value, options); 

						projection_lfp = isomap_output.coords{2,1};
						x = projection_lfp(1,:);
						y = projection_lfp(2,:);
						
					end % End if do_isomap


					
					
			
				% -------------------------------------------------------------------
				%  Run Dimensionality Reduction on global state spaces -- LFP
				% -------------------------------------------------------------------	
					projection_lfp = [];

					global_lfp_states_og = global_lfp_states;

					if do_rescale_global
						popmin = min(global_lfp_states);
						popmax = max(global_lfp_states);
						% global_lfp_states = rescale(global_lfp_states, 0, 1, 'InputMin',popmin,'InputMax',popmax);
						global_lfp_states = zscore(global_lfp_states);
					end

					if do_umap
					
						[projection_lfp, umap_params_lfp] = run_umap(global_lfp_states, 'n_components', 2, 'n_neighbors', 25 ,'min_dist', 0.1, 'metric', 'euclidean' );
						x = projection_lfp(:,1);
						y = projection_lfp(:,2);
						
					end % End if do_umap

					if do_pca
						fprintf('Performing PCA on LFP')
						% projection_lfp = pca scores;
						% same name for code consistency;
						[coeff_lfp, scores_lfp, latent_lfp, tsquared_lfp, variance_explained_lfp, mu_lfp] = pca(global_lfp_states);
						projection_lfp = scores_lfp(:,1:2);
						x = projection_lfp(:,1);
						y = projection_lfp(:,2);

					end % End if do_pca

					


				% -------------------------------------------------------------------
				%  Save global variables into memory for future use.
				% -------------------------------------------------------------------
				
					output_global_filename = fullfile(global_output_folder,'state_space_data.mat')

					if do_umap
						save(output_global_filename,'global_lfp_states_og','global_csd_states','global_lfp_states', 'global_states_count_csd','global_states_count_lfp',...
							'projection_csd', 'projection_lfp', 'bin_size_sec','global_artifacts_bins_lfp',...
							'global_artifacts_bins_csd', 'umap_params_lfp', 'umap_params_csd',...
							'global_acceleration_data_csd','global_acceleration_data_lfp',...
							'global_binned_spikes_csd','global_binned_spikes_lfp','global_lfp_for_csd','global_raw_bin_indices_csd','global_raw_bin_indices_lfp',...
							'csd_layer_indices','band_layer_struct')
					end % End if do_umap


					if do_pca
						save(output_global_filename,'global_lfp_states_og','global_csd_states','global_lfp_states', 'global_states_count_csd','global_states_count_lfp','projection_csd', 'projection_lfp',...
							'coeff_lfp','scores_lfp','latent_lfp','tsquared_lfp','variance_explained_lfp','mu_lfp',...
							'coeff_csd', 'scores_csd','latent_csd','tsquared_csd','variance_explained_csd','mu_csd',...
							'bin_size_sec','global_artifacts_bins_lfp','global_artifacts_bins_csd',...
							'global_acceleration_data_csd','global_acceleration_data_lfp',...
							'global_binned_spikes_csd','global_binned_spikes_lfp','global_lfp_for_csd','global_raw_bin_indices_csd','global_raw_bin_indices_lfp',...
							'csd_layer_indices','band_layer_struct')
					end % End if do_pca


			else

				load(fullfile(global_output_folder,'state_space_data.mat'))	

			end % End loading pre-computed data


		% -------------------------------------------------------------------
		% Plot Global State Space
		% -------------------------------------------------------------------	
			if plot_global_state_space
				fprintf('Plotting Global State Space...\n')

				x_lfp = projection_lfp(:,1);
				y_lfp = projection_lfp(:,2);
				x_csd = projection_csd(:,1);
				y_csd = projection_csd(:,2);

				figure(4)
				subplot(1,2,1)
				scatter(x_csd,y_csd, 5,  color_gray ,'filled');	
				pbaspect([1 1 1]) 
				% title(sprintf('CSD State Space'));
				% xlabel('UMAP Axis 1')
				% ylabel('UMAP Axis 2')
				xlim(csd_xlim);
				ylim(csd_ylim);


				subplot(1,2,2)
				scatter(x_lfp,y_lfp, 5,  color_gray ,'filled');	
				pbaspect([1 1 1]) 
				xlim(lfp_xlim);
				ylim(lfp_ylim);

				% title(sprintf('LFP State Space'));
				% xlabel('UMAP Axis 1')
				% ylabel('UMAP Axis 2')



				if make_pdf_images
					pdfname = sprintf('%s_global_state_space.pdf', animal_id);
					fullname_pdf = fullfile(global_output_folder,pdfname)
					figure(4)
					exportgraphics(gcf,fullname_pdf,'ContentType','image','BackgroundColor','white', 'Resolution',300)
				end

				if make_eps_images	
					epsname = sprintf('%s_global_state_space.eps', animal_id);
					fullname_eps = fullfile(global_output_folder,epsname)	
					set(gcf,'renderer','Painters')
					saveas(gcf,fullname_eps,'epsc')
				end

			end % End if plot_global_state_space


		% -------------------------------------------------------------------
		% Plot acceleration overlay
		% -------------------------------------------------------------------
			if plot_acceleration_overlay
				acceleration_overlay_folder = fullfile(global_output_folder,'acceleration_overlay');
				mkdir(acceleration_overlay_folder)
				fprintf('Accelerometer Overlay\n')

				acc_var =	movvar(global_acceleration_data_lfp,10);

				for trial_iterator = 1 : total_trials

					range_start = global_states_count_lfp(trial_iterator) + 1;
					range_end = global_states_count_lfp(trial_iterator + 1);

					pv_range = range_start : range_end; 

					x = projection_lfp(pv_range,1);
					y = projection_lfp(pv_range,2);
					
					current_trial = char(trial_folders(trial_iterator));
		
					cmin = min(acc_var);
					cmax = prctile(acc_var, 95);

					caxis_range = [cmin, cmax];

					trial_acceleration = acc_var(pv_range);

					[acc_cmap, acc_colors] = assign_colors(trial_acceleration, caxis_range, 10000, 'magma');


					% New fig for each band
					fig_id = 1525 + band_iterator;
					figure(fig_id);
					set(gcf, 'Position', get(0, 'Screensize'));
					subplot(1,4,trial_iterator)						
					scatter(projection_lfp(:,1), projection_lfp(:,2), 5 ,color_gray, 'filled')
					hold on
					scatter(x,y,5, acc_colors ,'filled');
					colormap(acc_cmap)
					colorbar
					caxis(caxis_range)
					h = colorbar;
					ylabel(h, 'Acc Variance')
					title(sprintf('%s', current_trial ))
					pbaspect([1 1 1]) 


					if trial_iterator == total_trials
						suptitle(sprintf('%s Accelerometer Overlay', animal_id))
						filename = sprintf('%s_acceleration_overlay.png', animal_id );
						fullname = fullfile(acceleration_overlay_folder,filename);
						saveas(figure(fig_id),fullname);

						if make_eps_images	
							epsname = sprintf('%s_accelerometer_overlay.eps', animal_id);
							fullname_eps = fullfile(acceleration_overlay_folder,epsname);
							set(gcf,'renderer','Painters')
							set(gca, 'visible' ,'off')
							saveas(gcf,fullname_eps,'epsc')
						end % End if make_eps_images

					end % End if trial_iterator		


				end % End for trial_iterator			

				return


			end % End if plot_acceleration_overlay
	

		% -------------------------------------------------------------------
		%  Plot Subspace Occupancy - CSD
		% -------------------------------------------------------------------		
			% if plot_subspace_occupation
			% 	subspace_occupation_folder = fullfile(global_output_folder,'subspace_occupation');
			% 	mkdir(subspace_occupation_folder)
			% 	trial_colors = magma(total_trials);


			% 	global_x = projection_csd(:,1);
			% 	global_y = projection_csd(:,2);

			% 	% Divide entire umap output into bins
			% 	global_reference_bins_x  = linspace( min(global_x), max(global_x), 100 );
			% 	global_reference_bins_y  = linspace( min(global_y), max(global_y), 100 );

			% 	all_occupancy = hist3([global_x, global_y],'Edges',{global_reference_bins_x, global_reference_bins_y});

			% 	total_occupied_bins = sum(all_occupancy(:) > 0);	

			% 	all_pv = 1:global_states_count_csd(5);

			% 	for trial_iterator = 1 : total_trials			

			% 		range_start = global_states_count_csd(trial_iterator) + 1 ;
			% 		range_end = global_states_count_csd(trial_iterator + 1) ;

			% 		pv_range = range_start : range_end; 

			% 		x = projection_csd(pv_range,1);
			% 		y = projection_csd(pv_range,2);
			% 		% x0 = 10;
			% 		% y0 = 10;
			% 		% widthbox = 800;
			% 		% heightbox = 700;

			% 		non_pv_range = setdiff(all_pv, pv_range);

			% 		figure(4)
			% 		set(gcf, 'Position', get(0, 'Screensize'));
			% 		% set(gcf,'position',[x0,y0,widthbox,heightbox])
			% 		subplot(1,4,trial_iterator)			
			% 		scatter(projection_csd(non_pv_range,1),projection_csd(non_pv_range,2),5 ,color_gray, 'filled')
			% 		hold on
			% 		scatter(x,y,5, trial_colors(trial_iterator,:) ,'filled');
			% 		pbaspect([1 1 1])
			% 		title(sprintf('%s', char(trial_folders(trial_iterator))))
			% 		xlim(csd_xlim);
			% 		ylim(csd_ylim);
			% 		% set(gca,'box','off') 
					


			% 		trial_occupancy = hist3([x, y],'Edges',{global_reference_bins_x, global_reference_bins_y});
			% 		trial_occupied_bins = sum(trial_occupancy(:) > 0);

			% 		global_subspace_occupancy_fraction_csd(trial_iterator) = trial_occupied_bins / total_occupied_bins;

					
			% 	end % End for trial_iterator
			% 	suptitle(sprintf('%s CSD State Space Occpuancy', animal_id))
			% 	filename = sprintf('%s_CSD_subspace_occupancy.png', animal_id);
			% 	fullname = fullfile(subspace_occupation_folder,filename);
			% 	saveas(gcf,fullname);
				

			% 	if make_pdf_images
			% 		pdfname = sprintf('%s_CSD_subspace_occupancy.pdf', animal_id);
			% 		fullname_pdf = fullfile(subspace_occupation_folder,pdfname)
			% 		figure(4)
			% 		exportgraphics(gcf,fullname_pdf,'ContentType','image','BackgroundColor','white', 'Resolution',300)
			% 	end
				
			% 	if make_eps_images	
			% 		epsname = sprintf('%s_CSD_subspace_occupancy.eps', animal_id);
			% 		fullname_eps = fullfile(subspace_occupation_folder,epsname)
			% 		set(gcf,'renderer','Painters')
			% 		set(gca, 'visible' ,'off')
			% 		saveas(gcf,fullname_eps,'epsc')
			% 	end

				
			% 	close all

			% 	output_global_filename = fullfile(subspace_occupation_folder,'global_subspace_occupancy_fraction_csd.mat');
			% 	save(output_global_filename,'global_subspace_occupancy_fraction_csd')

			% end % End if plot_subspace_occupation


		% -------------------------------------------------------------------
		%  Plot Subspace Occupancy - LFP
		% -------------------------------------------------------------------		
			if plot_subspace_occupation
				
				trial_colors = magma(total_trials);
				% trial_colors = trial_colors(end-3:end,:);

				global_x = projection_lfp(:,1);
				global_y = projection_lfp(:,2);

				% Divide entire umap output into bins
				global_reference_bins_x  = linspace( min(global_x), max(global_x), 100 );
				global_reference_bins_y  = linspace( min(global_y), max(global_y), 100 );

				all_occupancy = hist3([global_x, global_y],'Edges',{global_reference_bins_x, global_reference_bins_y});

				total_occupied_bins = sum(all_occupancy(:) > 0);

				all_pv = 1:global_states_count_lfp(5);

				
				for trial_iterator = 1 : total_trials

					range_start = global_states_count_lfp(trial_iterator) + 1 ;
					range_end = global_states_count_lfp(trial_iterator + 1) ;

					pv_range = range_start : range_end; 

					x = projection_lfp(pv_range,1);
					y = projection_lfp(pv_range,2);		
					non_pv_range = setdiff(all_pv, pv_range);	

					figure(5)
					set(gcf, 'Position', get(0, 'Screensize'));
					subplot(1,4,trial_iterator)
					scatter(projection_lfp(non_pv_range,1),projection_lfp(non_pv_range,2),5 ,color_gray, 'filled')
					hold on
					scatter(x,y,5, trial_colors(trial_iterator,:) ,'filled');		
					pbaspect([1 1 1])
					title(sprintf('%s', char(trial_folders(trial_iterator))))
					% xlim(lfp_xlim);
					% ylim(lfp_ylim);
					


					trial_occupancy = hist3([x, y],'Edges',{global_reference_bins_x, global_reference_bins_y});
					trial_occupied_bins = sum(trial_occupancy(:) > 0);

					global_subspace_occupancy_fraction_lfp(trial_iterator) = trial_occupied_bins / total_occupied_bins;


					
				end % End for trial_iterator
				suptitle(sprintf('%s LFP State Space Occpuancy', animal_id))
				filename = sprintf('%s_LFP_subspace_occupancy.png', animal_id);
				fullname = fullfile(subspace_occupation_folder,filename);
				saveas(gcf,fullname);

				if make_pdf_images
					pdfname = sprintf('%s_LFP_subspace_occupancy.pdf', animal_id);
					fullname_pdf = fullfile(subspace_occupation_folder,pdfname)
					figure(5)
					exportgraphics(gcf,fullname_pdf,'ContentType','image','BackgroundColor','white', 'Resolution',300)
				end

				if make_eps_images	
					epsname = sprintf('%s_LFP_subspace_occupancy.eps', animal_id);
					fullname_eps = fullfile(subspace_occupation_folder,epsname)
					set(gcf,'renderer','Painters')
					set(gca, 'visible' ,'off')
					saveas(gcf,fullname_eps,'epsc')
				end

				close all

				output_global_filename = fullfile(subspace_occupation_folder,'global_subspace_occupancy_fraction_lfp.mat');
				save(output_global_filename,'global_subspace_occupancy_fraction_lfp')
			
			end % End if plot_subspace_occupation
			

		% -------------------------------------------------------------------
		%  Plot Subspace Occupancy - LFP
		% -------------------------------------------------------------------		
			if plot_subspace_occupation_time_control

				fprintf('Subspace Occupancy Time Control')

				subspace_occupation_folder = fullfile(global_output_folder,'subspace_occupation');
				mkdir(subspace_occupation_folder)
				
				trial_colors = magma(total_trials);

				global_x = projection_lfp(:,1);
				global_y = projection_lfp(:,2);

				global_x_sub = [];
				global_y_sub = [];
				ind = {};

				for trial_iterator = 1 : total_trials

					x = []; y = [];

					range_start = global_states_count_lfp(trial_iterator) + 1 ;
					range_end = global_states_count_lfp(trial_iterator + 1) ;

					pv_range = range_start : range_end; 

					x = projection_lfp(pv_range,1);
					y = projection_lfp(pv_range,2);		

					random_5k = randperm(numel(x), 4000);

					ind{1,trial_iterator} = random_5k;

					global_x_sub = [global_x_sub; reshape(x(random_5k),[],1)];
					global_y_sub = [global_y_sub; reshape(y(random_5k),[],1)];

				end % End for trial_iterator


				% Divide entire umap output into bins
				global_reference_bins_x  = linspace( min(global_x_sub), max(global_x_sub), 50 );
				global_reference_bins_y  = linspace( min(global_y_sub), max(global_y_sub), 50 );

				all_occupancy = hist3([global_x_sub, global_y_sub],'Edges',{global_reference_bins_x, global_reference_bins_y});

				total_occupied_bins = sum(all_occupancy(:) > 0);

				all_pv = 1:global_states_count_lfp(5);

				
				for trial_iterator = 1 : total_trials

					range_start = global_states_count_lfp(trial_iterator) + 1;
					range_end = global_states_count_lfp(trial_iterator + 1);

					pv_range = range_start : range_end; 

					x = projection_lfp(pv_range,1);
					y = projection_lfp(pv_range,2);	


					x = x(ind{1,trial_iterator});
					y = y(ind{1,trial_iterator});	

					figure(5)
					set(gcf, 'Position', get(0, 'Screensize'));
					subplot(1,4,trial_iterator)
					scatter(global_x_sub, global_y_sub,5 ,color_gray, 'filled')
					hold on
					scatter(x,y,5, trial_colors(trial_iterator,:) ,'filled');		
					pbaspect([1 1 1])
					title(sprintf('%s', char(trial_folders(trial_iterator))))
					% xlim(lfp_xlim);
					% ylim(lfp_ylim);
					


					trial_occupancy = hist3([x, y],'Edges',{global_reference_bins_x, global_reference_bins_y});
					trial_occupied_bins = sum(trial_occupancy(:) > 0);

					global_subspace_occupancy_fraction_lfp(trial_iterator) = trial_occupied_bins / total_occupied_bins;


					
				end % End for trial_iterator


				suptitle(sprintf('%s LFP State Space Occpuancy Time Control', animal_id))
				filename = sprintf('%s_LFP_subspace_occupancy_time_control.png', animal_id);
				fullname = fullfile(subspace_occupation_folder,filename);
				saveas(gcf,fullname);
			

				if make_pdf_images
					pdfname = sprintf('%s_LFP_subspace_occupancy.pdf', animal_id);
					fullname_pdf = fullfile(subspace_occupation_folder,pdfname)
					figure(5)
					exportgraphics(gcf,fullname_pdf,'ContentType','image','BackgroundColor','white', 'Resolution',300)
				end

				if make_eps_images	
					epsname = sprintf('%s_LFP_subspace_occupancy_time_control.eps', animal_id);
					fullname_eps = fullfile(subspace_occupation_folder,epsname)
					set(gcf,'renderer','Painters')
					set(gca, 'visible' ,'off')
					saveas(gcf,fullname_eps,'epsc')
				end

				close all

				output_global_filename = fullfile(subspace_occupation_folder,'occupancy_time_control.mat');
				save(output_global_filename,'global_subspace_occupancy_fraction_lfp')
			
			end % End if plot_subspace_occupation


		% -------------------------------------------------------------------
		%  Plot Subspace Density - CSD
		% -------------------------------------------------------------------		
			if plot_subspace_density_csd
				fprintf('Subspace Density CSD...\n')
				subspace_density_folder = fullfile(global_output_folder,'subspace_occupation');
				mkdir(subspace_density_folder)



				no_xbins = 20;
				no_ybins = 20;
				total_bins = no_xbins * no_ybins;

				global_x = projection_csd(:,1);
				global_y = projection_csd(:,2);

				% Divide entire umap output into bins
				global_edges_x  = linspace( min(global_x), max(global_x), no_xbins );
				global_edges_y  = linspace( min(global_y), max(global_y), no_ybins );

				global_trial_density = []; global_median_density_csd = [];
				
				all_pv = 1:global_states_count_csd(5);
				
				
				% Compute density for all trials
				for trial_iterator = 1 : total_trials			

					range_start = global_states_count_csd(trial_iterator) + 1 ;
					range_end = global_states_count_csd(trial_iterator + 1) ;

					pv_range = range_start : range_end; 

					x = projection_csd(pv_range,1);
					y = projection_csd(pv_range,2);
					trial_duration = numel(x) * bin_size_sec;
 					non_pv_range = setdiff(all_pv, pv_range);
					trial_density = [];

					% Get xindices 
					[xcounts, x_edges, xindices] = histcounts(x, 'BinEdges', global_edges_x);
					[ycounts, y_edges, yindices] = histcounts(y, 'BinEdges', global_edges_y);


					for xbin_iterator = 1:no_xbins
						tempx = []; 
						tempx =	find(xindices == xbin_iterator);

						for ybin_iterator = 1:no_ybins
							tempy = [];
							
							tempy =	find(yindices == ybin_iterator);

							bin_states_indices = intersect(tempx,tempy);
							bsi = bin_states_indices;	

							if numel(bsi) == 0;
								continue;
							end			 

							trial_density(bsi) = numel(bsi) / trial_duration;

						end % End for ybin_iterator

					end % End for xbin_iterator	
					
					global_trial_density = [global_trial_density; reshape(trial_density,[],1)];

				end % End for trial_iterator



				% Plot density on same scale;
				% Compute density for all trials
				for trial_iterator = 1 : total_trials			

					range_start = global_states_count_csd(trial_iterator) + 1 ;
					range_end = global_states_count_csd(trial_iterator + 1) ;

					pv_range = range_start : range_end; 
					non_pv_range = setdiff(all_pv, pv_range);	

					x = []; y = []; trial_density = []; density_colors = [];

					x = projection_csd(pv_range,1);
					y = projection_csd(pv_range,2);
					trial_density = global_trial_density(pv_range);
					global_median_density_csd(trial_iterator) = median(trial_density);
					

					caxis_range = [min(global_trial_density), prctile(global_trial_density,95)];

					[density_cmap, density_colors] = assign_colors(trial_density, caxis_range, 1000, '*purd'); 

					figure(123)
					set(gcf, 'Position', get(0, 'Screensize'));
					subplot(1,4,trial_iterator)
					scatter(projection_csd(non_pv_range,1),projection_csd(non_pv_range,2),5 ,color_gray, 'filled')
					hold on
					scatter(x,y,5, density_colors ,'filled');	
					pbaspect([1 1 1])
					title(sprintf('%s', char(trial_folders(trial_iterator))))
					cbh = colorbar;
					ylabel(cbh,'Norm. Density (No. of States / sec / bin)')
					caxis(caxis_range)
					colormap(density_cmap)
					xlim(csd_xlim);
					ylim(csd_ylim);

 

 				end % End for trial_iterator

				suptitle(sprintf('%s CSD State Space Density', animal_id))
				filename = sprintf('%s_CSD_subspace_density.png', animal_id);
				fullname = fullfile(subspace_density_folder,filename);
				saveas(gcf,fullname);


				if make_pdf_images
					pdfname = sprintf('%s_CSD_subspace_density.pdf', animal_id);
					fullname_pdf = fullfile(subspace_density_folder,pdfname)
					figure(123)
					exportgraphics(gcf,fullname_pdf,'ContentType','image','BackgroundColor','white', 'Resolution',300)
				end

				if make_eps_images	
					epsname = sprintf('%s_CSD_subspace_density.eps', animal_id);
					fullname_eps = fullfile(subspace_density_folder,epsname)
					set(gcf,'renderer','Painters')
					saveas(gcf,fullname_eps,'epsc')
				end

				close all
		

				output_global_filename = fullfile(subspace_density_folder,'global_median_density_csd.mat');
				save(output_global_filename,'global_median_density_csd')

			end % End if plot_subspace_density_csd


		% -------------------------------------------------------------------
		%  Plot Subspace Density - LFP
		% -------------------------------------------------------------------		
			if plot_subspace_density_lfp
				fprintf('Subspace Density LFP...\n')
				subspace_density_folder = fullfile(global_output_folder,'subspace_occupation');
				mkdir(subspace_density_folder)

				lfp_overlay_output_folder = fullfile(global_output_folder,'overlay_lfp_power_bins');
				load(fullfile(lfp_overlay_output_folder,'is_rem_nrem_lfp.mat'))


				global_bin_density_lfp = [];
				global_bin_bsi = {};

				no_xbins = 20;
				no_ybins = 20;
				total_bins = no_xbins * no_ybins;

				global_x = projection_lfp(:,1);
				global_y = projection_lfp(:,2);

				% Divide entire umap output into bins
				global_edges_x  = linspace( min(global_x), max(global_x), no_xbins );
				global_edges_y  = linspace( min(global_y), max(global_y), no_ybins );

				global_trial_density = [];
				global_median_density_lfp = [];

				global_rem_density = [];
				global_nonrem_density = [];

				global_rem_ratio = rescale(global_lfp_states(:,theta_index) ./ global_lfp_states(:,delta_index));
				global_nonrem_ratio = rescale(global_lfp_states(:,delta_index) .* global_lfp_states(:,spindle_index));

				all_pv = 1:global_states_count_lfp(5);


				% Compute density for all trials
				for trial_iterator = 1 : total_trials			

					range_start = global_states_count_lfp(trial_iterator) + 1 ;
					range_end = global_states_count_lfp(trial_iterator + 1) ;

					pv_range = range_start : range_end; 

					x = projection_lfp(pv_range,1);
					y = projection_lfp(pv_range,2);
					trial_duration = numel(x) * bin_size_sec;

					trial_isrem = is_rem(pv_range);
					trial_isnonrem = is_nonrem(pv_range);

					trial_rem_density = [];
					trial_nonrem_density = [];
 
					trial_density = [];

					% Get xindices 
					[xcounts, x_edges, xindices] = histcounts(x, 'BinEdges', global_edges_x);
					[ycounts, y_edges, yindices] = histcounts(y, 'BinEdges', global_edges_y);

					bin_id = 1;
					trial_bin_density = [];
					

					for xbin_iterator = 1:no_xbins
						tempx = []; 
						tempx =	find(xindices == xbin_iterator);

						for ybin_iterator = 1:no_ybins
							tempy = [];
							bin_id = (xbin_iterator - 1) * no_xbins + ybin_iterator;

							
							tempy =	find(yindices == ybin_iterator);

							bin_states_indices = intersect(tempx,tempy);
							bsi = bin_states_indices;	

							if numel(bsi) == 0;
								trial_bin_density(bin_id) = NaN;
								global_bin_bsi{bin_id,trial_iterator} = [];
								continue;
							end		

							bin_density = [];

							bin_density = numel(bsi) / trial_duration;

							trial_density(bsi) = bin_density;

							trial_bin_density(bin_id) = bin_density;

							global_bin_bsi{bin_id,trial_iterator} = bsi;

							% Collect density for rem and nonrem bins
							if ismember(trial_iterator,[1 4])
								half_size = length(bsi) / 2;

								if sum(trial_isrem(bsi)) > half_size
									trial_rem_density = [trial_rem_density bin_density];
								end

								if sum(trial_isnonrem(bsi)) > half_size
									trial_nonrem_density = [trial_nonrem_density bin_density];
								end


							end % End if ismember



						end % End for ybin_iterator

					end % End for xbin_iterator	

					global_bin_density_lfp(:,trial_iterator) = trial_bin_density;
					
					global_trial_density = [global_trial_density; reshape(trial_density,[],1)];

					global_rem_density(trial_iterator) = median(trial_rem_density);

					global_nonrem_density(trial_iterator) = median(trial_nonrem_density);

					
				end % End for trial_iterator

			

				% Compute difference in density for pre and post sleep
				difference_density = global_bin_density_lfp(:,4) - global_bin_density_lfp(:,1);
				
				nan_ind = find(isnan(difference_density));

				non_nan_ind = find(~isnan(difference_density));

				difference_density_colors = zeros(length(difference_density), 3);

				max_diff = max(abs(difference_density));

				differnce_density_range = [-max_diff max_diff];

				[ddcmap, difference_density_colors] = assign_colors(difference_density, differnce_density_range, [] , 'spectral');
				
				total_bins = length(difference_density);

				state_space_colors = [];

				% Plot difference on post sleep state space
					for bin_iterator = 1:total_bins

						bsi = global_bin_bsi{bin_iterator,4};	
						
						if isnan(difference_density(bin_iterator))

							if ~isempty(global_bin_bsi{bin_iterator,4})
								state_space_colors(bsi,:) = repmat(color_gray, length(bsi), 1);
							end % End if ~isempty

						else 
									
							state_space_colors(bsi,:) = repmat(difference_density_colors(bin_iterator,:), length(bsi), 1);	

						end % End if isnan

					end % End for bin_iterator	









				% Plot density on same scale;
				% Compute density for all trials
				for trial_iterator = 1 : total_trials			

					range_start = global_states_count_lfp(trial_iterator) + 1 ;
					range_end = global_states_count_lfp(trial_iterator + 1) ;

					pv_range = range_start : range_end; 
					non_pv_range = setdiff(all_pv, pv_range);

					x = []; y = []; trial_density = []; density_colors = [];

					x = projection_lfp(pv_range,1);
					y = projection_lfp(pv_range,2);
					trial_density = global_trial_density(pv_range);

					global_median_density_lfp(trial_iterator) = median(trial_density);

					caxis_range = [min(global_trial_density), prctile(global_trial_density,95)];

					[density_cmap, density_colors] = assign_colors(trial_density, caxis_range, 1000, 'summer'); 

					figure(123)
					set(gcf, 'Position', get(0, 'Screensize'));
					subplot(1,4,trial_iterator)
					scatter(projection_lfp(non_pv_range,1),projection_lfp(non_pv_range,2),5 ,color_gray, 'filled')
					hold on
					scatter(x,y,5, density_colors ,'filled');	
					pbaspect([1 1 1])
					title(sprintf('%s', char(trial_folders(trial_iterator))))
					cbh = colorbar;
					ylabel(cbh,'Norm. Density (No. of States / sec / bin)')
					caxis(caxis_range)
					colormap(density_cmap)
					xlim(lfp_xlim);
					ylim(lfp_ylim);



					if trial_iterator == 4
						figure(987)
						% set(gcf, 'Position', get(gca,'PropertyName');(0, 'Screensize'));
						scatter(projection_lfp(non_pv_range,1),projection_lfp(non_pv_range,2),5 ,color_gray, 'filled')
						hold on
						scatter(x,y,5, state_space_colors ,'filled');	
						pbaspect([1 1 1])
						title(sprintf('postSleep - presleep'))
						cbh = colorbar;
						ylabel(cbh,'Change in Density')
						caxis(differnce_density_range)
						colormap(ddcmap)
						xlim(lfp_xlim);
						ylim(lfp_ylim);
					
					end

	 

 				end % End for trial_iterator

 				figure(123)
				suptitle(sprintf('%s LFP State Space Density', animal_id))
				filename = sprintf('%s_LFP_subspace_density.png', animal_id);
				fullname = fullfile(subspace_density_folder,filename);
				saveas(gcf,fullname);

				figure(987)
				suptitle(sprintf('%s LFP Difference Density', animal_id))
				filename = sprintf('%s_LFP_density_diff.png', animal_id);
				fullname = fullfile(subspace_density_folder,filename);
				saveas(gcf,fullname);

				if make_pdf_images
					pdfname = sprintf('%s_LFP_subspace_density.pdf', animal_id);
					fullname_pdf = fullfile(subspace_density_folder,pdfname);
					figure(123)
					exportgraphics(gcf,fullname_pdf,'ContentType','image','BackgroundColor','white', 'Resolution',300)
				end


				if make_eps_images	
					epsname = sprintf('%s_LFP_subspace_density.eps', animal_id);
					fullname_eps = fullfile(subspace_density_folder,epsname);
					figure(123)
					set(gcf,'renderer','Painters')
					saveas(gcf,fullname_eps,'epsc')

					epsname = sprintf('%s_LFP_density_diff.eps', animal_id);
					fullname_eps = fullfile(subspace_density_folder,epsname);
					figure(987)
					set(gcf,'renderer','Painters')
					saveas(gcf,fullname_eps,'epsc')
				end

				return
				close all
			
				output_global_filename = fullfile(subspace_density_folder,'global_median_density_lfp.mat');
				save(output_global_filename,'global_median_density_lfp', 'global_rem_density','global_nonrem_density')
			end % End if plot_subspace_density_lfp


		% -------------------------------------------------------------------
		% LFP Overlay (Median in Bins)
		% -------------------------------------------------------------------		
			if feature_overlay_bins_lfp

				fprintf('Power Overlay in Bins...\n')

				lfp_overlay_output_folder = fullfile(global_output_folder,'overlay_lfp_power_bins');
				mkdir(lfp_overlay_output_folder)

				no_xbins = 20;
				no_ybins = 20;
				total_bins = no_xbins * no_ybins;

				global_x = projection_lfp(:,1);
				global_y = projection_lfp(:,2);

				% Divide entire umap output into bins
				global_edges_x  = linspace( min(global_x), max(global_x), no_xbins+1 );
				global_edges_y  = linspace( min(global_y), max(global_y), no_ybins+1 );

				global_bin_power_all_bands = [];

				total_features = size(global_lfp_states,2);
				
				% Compute density for all trials
				for trial_iterator = 1 : total_trials			

					range_start = global_states_count_lfp(trial_iterator) + 1 ;
					range_end = global_states_count_lfp(trial_iterator + 1) ;

					pv_range = range_start : range_end; 

					x = projection_lfp(pv_range,1);
					y = projection_lfp(pv_range,2);
					trial_duration = numel(x) * bin_size_sec;

					trial_lfp_states = global_lfp_states(pv_range,:);
 
					trial_median_bin_power = [];
					trail_bin_variance = [];

					% Get xindices 
					[xcounts, x_edges, xindices] = histcounts(x, 'BinEdges', global_edges_x);
					[ycounts, y_edges, yindices] = histcounts(y, 'BinEdges', global_edges_y);


					for xbin_iterator = 1:no_xbins
						tempx = []; 
						tempx =	find(xindices == xbin_iterator);

						for ybin_iterator = 1:no_ybins
							tempy = [];
							
							tempy =	find(yindices == ybin_iterator);

							bin_states_indices = intersect(tempx,tempy);
							bsi = bin_states_indices;	

							if numel(bsi) < 1;
								continue;
							end			 
							median_power_bin_all_bands = [];
							median_power_bin_all_bands = median(trial_lfp_states(bsi,:));

							trial_median_bin_power(bsi,1:total_features) = repmat(median_power_bin_all_bands, length(bsi), 1);
							
						end % End for ybin_iterator

					end % End for xbin_iterator	
					
					global_bin_power_all_bands = [global_bin_power_all_bands; trial_median_bin_power];

				end % End for trial_iterator

				band_counter = 1;
				global_median_power_in_trials = [];
				band_color = magma(6);
				all_pv = 1:global_states_count_lfp(5);
				

				for feature_iterator = 1:total_features
					fprintf('Processing Feature %d...\n', feature_iterator)

					caxis_range = [];
					caxis_range = [min(global_bin_power_all_bands(:,feature_iterator)), prctile(global_bin_power_all_bands(:,feature_iterator),95	)];

					% Plot density on same scale;
					% Compute density for all trials
					for trial_iterator = 1 : total_trials			

						range_start = global_states_count_lfp(trial_iterator) + 1 ;
						range_end = global_states_count_lfp(trial_iterator + 1) ;

						pv_range = range_start : range_end; 

						x = []; y = []; trial_band_power = []; trial_band_colors = [];

						x = projection_lfp(pv_range,1);
						y = projection_lfp(pv_range,2);
						trial_band_power = global_bin_power_all_bands(pv_range, feature_iterator);


						global_median_power_in_trials(feature_iterator, trial_iterator) = median(trial_band_power);

						[power_cmap, trial_band_colors] = assign_colors(trial_band_power, caxis_range, 1000, 'viridis'); 

						non_pv_range = setdiff(all_pv, pv_range);

						plot_overlay_maps = 1;
						if plot_overlay_maps 
							figure(123)
							set(gcf, 'Position', get(0, 'Screensize'));
							subplot(1,4,trial_iterator)
							scatter(projection_lfp(non_pv_range,1),projection_lfp(non_pv_range,2),5 ,color_gray, 'filled')
							hold on
							scatter(x,y,5, trial_band_colors ,'filled');	
							pbaspect([1 1 1])
							title(sprintf('%s', char(trial_folders(trial_iterator))))
							% if trial_iterator == 4
							% 	cbh = colorbar;
							% 	ylabel(cbh,'Median Power in Bins')
							% 	caxis(caxis_range)

							% end
							colormap(power_cmap)
							set(gca,'XColor','none','YColor','none')
							xlim([min(projection_lfp(:,1)) max(projection_lfp(:,1)) ])
							ylim([min(projection_lfp(:,2)) max(projection_lfp(:,2)) ])
						end % End if plot_overlay_maps

						plot_distribution = 0;
						if trial_iterator == 4 & plot_distribution
							figure(234)
							set(gcf, 'Position', get(0, 'Screensize'));
							subplot(3,6,feature_iterator)
							hist_temp = histogram(zscore(global_bin_power_all_bands(:,feature_iterator)), 100, 'FaceColor', band_color(band_counter,:), 'EdgeColor', band_color(band_counter,:));
							title(sprintf('%s %s', animal_id, x_ticks_array(feature_iterator)))
							caxis_range(2) = prctile(global_bin_power_all_bands(:,feature_iterator), 99);
							
							set(gca,'box','off') 
							xticks(caxis_range)
							xticklabels(fix(caxis_range))
							yaxis_range = [min(hist_temp.Values) max(hist_temp.Values)];

							yaxis_range = ceil(yaxis_range/10)*10;
							xaxis_range = ceil(caxis_range/10)*10;

							yticks(yaxis_range)
							yticklabels(yaxis_range)
							xticks(xaxis_range)
							xticklabels(xaxis_range)

							xlim(xaxis_range)
						end % End if trial_iterator

						if feature_iterator == 1 & trial_iterator == 1
							figure(185)
							set(gcf, 'Position', get(0, 'Screensize'));
							subplot(1,4,trial_iterator)
							scatter(projection_lfp(:,1),projection_lfp(:,2),5 ,'k', 'filled')
							xticks(global_edges_x)
							yticks(global_edges_y)
							grid on
							ax = gca;
							% ax.GridColor = [0 0 0];
							ax.GridLineStyle = '-';
							ax.GridColor = 'k';
							ax.GridAlpha = 1; % maximum line opacity
							set (gca, 'xticklabel' , {[]});
							set (gca, 'yticklabel' , {[]});
							pbaspect([1 1 1])
							xlim([min(projection_lfp(:,1)) max(projection_lfp(:,1)) ])
							ylim([min(projection_lfp(:,2)) max(projection_lfp(:,2)) ])
														
						end % End if feature_iterator
 

	 				end % End for trial_iterator



	 				if plot_overlay_maps 
		 				figure(123)
						suptitle(sprintf('%s %s ( %d-%d Hz ) Median Power on LFP', animal_id, x_ticks_array(feature_iterator), bands_array(band_counter,1), bands_array(band_counter,2)))
						filename = sprintf('%s_%s_LFP_power_overlay.png', animal_id, x_ticks_array(feature_iterator));;
						fullname = fullfile(lfp_overlay_output_folder,filename);
						saveas(gcf,fullname);

						if make_eps_images
							epsname = sprintf('%s_%s_LFP_power_overlay.eps', animal_id, x_ticks_array(feature_iterator));;
							fullname_eps = fullfile(lfp_overlay_output_folder,epsname)	
							set(gcf,'renderer','Painters')
							saveas(gcf,fullname_eps,'epsc')
						end % End if make_eps_images
				

						close(figure(123))

					end % End if plot_overlay_maps
				

					if band_counter == 6
						band_counter = 1;
					else
						band_counter = band_counter + 1;
					end


				end % End for feature_iterator

				if plot_distribution

					% figure(234)
					% suptitle(sprintf('%s Distribution of Band Power LFP', animal_id ))
					% filename = sprintf('%s_LFP_power_distribution.png', animal_id);
					% fullname = fullfile(lfp_overlay_output_folder,filename);

					% saveas(gcf,fullname);
					% % exportgraphics(gcf,fullname,'Resolution',300)

					% if make_pdf_images
					% 	pdfname = sprintf('%s_LFP_power_distribution.pdf', animal_id);
					% 	fullname_pdf = fullfile(lfp_overlay_output_folder,pdfname)
					% 	figure(234)
					% 	exportgraphics(gcf,fullname_pdf,'ContentType','image','BackgroundColor','none', 'Resolution',300)
					% end
				end % End if plot_distribution
					
				output_global_filename = fullfile(lfp_overlay_output_folder,'global_median_power_lfp.mat');
				save(output_global_filename,'global_median_power_in_trials','global_bin_power_all_bands')



			end % End if feature_overlay_bins_lfp


		% -------------------------------------------------------------------
		% LFP Variance Overlay (in Bins)
		% -------------------------------------------------------------------
			if variance_overlay_bins_lfp

				fprintf('Variance Overlay in Bins...\n')

				lfp_variance_output_folder = fullfile(global_output_folder,'overlay_lfp_variance_bins');
				mkdir(lfp_variance_output_folder)

				no_xbins = 20;
				no_ybins = 20;
				total_bins = no_xbins * no_ybins;

				global_x = projection_lfp(:,1);
				global_y = projection_lfp(:,2);

				% Divide entire umap output into bins
				global_edges_x  = linspace( min(global_x), max(global_x), no_xbins+1 );
				global_edges_y  = linspace( min(global_y), max(global_y), no_ybins+1 );

				global_bin_variance = [];

				total_features = size(global_lfp_states,2);
				
				% Compute density for all trials
				for trial_iterator = 1 : total_trials			

					range_start = global_states_count_lfp(trial_iterator) + 1 ;
					range_end = global_states_count_lfp(trial_iterator + 1) ;

					pv_range = range_start : range_end; 

					x = projection_lfp(pv_range,1);
					y = projection_lfp(pv_range,2);
					trial_duration = numel(x) * bin_size_sec;

					trial_lfp_states = global_lfp_states(pv_range,:);
 
					trail_bin_variance = [];

					% Get xindices 
					[xcounts, x_edges, xindices] = histcounts(x, 'BinEdges', global_edges_x);
					[ycounts, y_edges, yindices] = histcounts(y, 'BinEdges', global_edges_y);


					for xbin_iterator = 1:no_xbins
						tempx = []; 
						tempx =	find(xindices == xbin_iterator);

						for ybin_iterator = 1:no_ybins
							tempy = [];
							
							tempy =	find(yindices == ybin_iterator);

							bin_states_indices = intersect(tempx,tempy);
							bsi = bin_states_indices;	

							if numel(bsi) < 1;
								continue;
							end			 
							bin_variance_all_bins = [];
							bin_variance_all_bins = var(trial_lfp_states(bsi,:));

							trail_bin_variance(bsi,1:total_features) = repmat(bin_variance_all_bins, length(bsi), 1);
							
						end % End for ybin_iterator

					end % End for xbin_iterator	
					
					global_bin_variance = [global_bin_variance; trail_bin_variance];

				end % End for trial_iterator

				band_counter = 1;
				global_median_variance = [];
				band_color = magma(6);
				all_pv = 1:global_states_count_lfp(5);
				

				for feature_iterator = 1:total_features
					fprintf('Processing Feature %d...\n', feature_iterator)

					caxis_range = [];
					caxis_range = [min(global_bin_variance(:,feature_iterator)), prctile(global_bin_variance(:,feature_iterator),95 ) ];

					% Plot density on same scale;
					% Compute density for all trials
					for trial_iterator = 1 : total_trials			

						range_start = global_states_count_lfp(trial_iterator) + 1 ;
						range_end = global_states_count_lfp(trial_iterator + 1) ;

						pv_range = range_start : range_end; 

						x = []; y = []; trial_band_power = []; trial_band_colors = [];

						x = projection_lfp(pv_range,1);
						y = projection_lfp(pv_range,2);
						trial_band_power = global_bin_variance(pv_range, feature_iterator);


						global_median_variance(feature_iterator, trial_iterator) = median(trial_band_power);

						[variance_cmap, trial_band_colors] = assign_colors(trial_band_power, caxis_range, 1000, 'magma'); 

						non_pv_range = setdiff(all_pv, pv_range);

						plot_overlay_maps = 1;
						if plot_overlay_maps 
							figure(123)
							set(gcf, 'Position', get(0, 'Screensize'));
							subplot(1,4,trial_iterator)
							scatter(projection_lfp(non_pv_range,1),projection_lfp(non_pv_range,2),5 ,color_gray, 'filled')
							hold on
							scatter(x,y,5, trial_band_colors ,'filled');	
							pbaspect([1 1 1])
							title(sprintf('%s', char(trial_folders(trial_iterator))))
							cbh = colorbar;
							ylabel(cbh,'Variance')
							caxis(caxis_range)

						
							colormap(variance_cmap)
							set(gca,'XColor','none','YColor','none')
							xlim([min(projection_lfp(:,1)) max(projection_lfp(:,1)) ])
							ylim([min(projection_lfp(:,2)) max(projection_lfp(:,2)) ])
						end % End if plot_overlay_maps

						plot_distribution = 0;
						if trial_iterator == 4 & plot_distribution
							figure(234)
							set(gcf, 'Position', get(0, 'Screensize'));
							subplot(3,6,feature_iterator)
							hist_temp = histogram(zscore(global_bin_variance(:,feature_iterator)), 100, 'FaceColor', band_color(band_counter,:), 'EdgeColor', band_color(band_counter,:));
							title(sprintf('%s %s', animal_id, x_ticks_array(feature_iterator)))
							caxis_range(2) = prctile(global_bin_variance(:,feature_iterator), 99);
							
							set(gca,'box','off') 
							xticks(caxis_range)
							xticklabels(fix(caxis_range))
							yaxis_range = [min(hist_temp.Values) max(hist_temp.Values)];

							yaxis_range = ceil(yaxis_range/10)*10;
							xaxis_range = ceil(caxis_range/10)*10;

							yticks(yaxis_range)
							yticklabels(yaxis_range)
							xticks(xaxis_range)
							xticklabels(xaxis_range)

							xlim(xaxis_range)
						end % End if trial_iterator

						if feature_iterator == 1 & trial_iterator == 1
							figure(185)
							set(gcf, 'Position', get(0, 'Screensize'));
							subplot(1,4,trial_iterator)
							scatter(projection_lfp(:,1),projection_lfp(:,2),5 ,'k', 'filled')
							xticks(global_edges_x)
							yticks(global_edges_y)
							grid on
							ax = gca;
							% ax.GridColor = [0 0 0];
							ax.GridLineStyle = '-';
							ax.GridColor = 'k';
							ax.GridAlpha = 1; % maximum line opacity
							set (gca, 'xticklabel' , {[]});
							set (gca, 'yticklabel' , {[]});
							pbaspect([1 1 1])
							xlim([min(projection_lfp(:,1)) max(projection_lfp(:,1)) ])
							ylim([min(projection_lfp(:,2)) max(projection_lfp(:,2)) ])
														
						end % End if feature_iterator
 

	 				end % End for trial_iterator



	 				if plot_overlay_maps 
		 				figure(123)
						suptitle(sprintf('%s %s ( %d-%d Hz ) Variance  on LFP', animal_id, x_ticks_array(feature_iterator), bands_array(band_counter,1), bands_array(band_counter,2)))
						filename = sprintf('%s_%s_LFP_power_variance_overlay.png', animal_id, x_ticks_array(feature_iterator));;
						fullname = fullfile(lfp_variance_output_folder,filename);
						saveas(gcf,fullname);

						if make_eps_images
							epsname = sprintf('%s_%s_LFP_power_variance_overlay.eps', animal_id, x_ticks_array(feature_iterator));;
							fullname_eps = fullfile(lfp_variance_output_folder,epsname)	
							set(gcf,'renderer','Painters')
							saveas(gcf,fullname_eps,'epsc')
						end % End if make_eps_images
				

						close(figure(123))

					end % End if plot_overlay_maps
				

					if band_counter == 6
						band_counter = 1;
					else
						band_counter = band_counter + 1;
					end


				end % End for feature_iterator

				if plot_distribution

					% figure(234)
					% suptitle(sprintf('%s Distribution of Band Power LFP', animal_id ))
					% filename = sprintf('%s_LFP_power_distribution.png', animal_id);
					% fullname = fullfile(lfp_variance_output_folder,filename);

					% saveas(gcf,fullname);
					% % exportgraphics(gcf,fullname,'Resolution',300)

					% if make_pdf_images
					% 	pdfname = sprintf('%s_LFP_power_distribution.pdf', animal_id);
					% 	fullname_pdf = fullfile(lfp_variance_output_folder,pdfname)
					% 	figure(234)
					% 	exportgraphics(gcf,fullname_pdf,'ContentType','image','BackgroundColor','none', 'Resolution',300)
					% end
				end % End if plot_distribution
					
				output_global_filename = fullfile(lfp_variance_output_folder,'global_bin_variance_lfp.mat');
				save(output_global_filename,'global_median_variance','global_bin_variance')



			end % End if variance_overlay_bins_lfp


		% -------------------------------------------------------------------
		% REM / Non-REM Overlay (Median in Bins)
		% -------------------------------------------------------------------
			if rem_nonrem_overlay_bins_lfp

				fprintf('REM/Non-REM Overlay in Bins...\n')

				lfp_overlay_output_folder = fullfile(global_output_folder,'overlay_lfp_power_bins');

				load(fullfile(lfp_overlay_output_folder,'global_median_power_lfp.mat'))


				
				all_theta = global_bin_power_all_bands(:,theta_index);
				all_delta = global_bin_power_all_bands(:,delta_index);
				all_spindle = global_bin_power_all_bands(:,spindle_index);

				global_rem_ratio = zscore(all_theta ./ all_delta);
				global_nonrem_ratio = zscore(all_spindle .* all_delta);


				% For Non-REM
					nonrem_threshold_std = prctile(global_nonrem_ratio,70);

					is_nonrem = global_nonrem_ratio > nonrem_threshold_std;

					nonrem_axis = [min(global_nonrem_ratio) prctile(global_nonrem_ratio,95)];
					nonrem_colors = [];
					[nonrem_cmap, nonrem_colors] = assign_colors(global_nonrem_ratio, nonrem_axis, 10, 'pink' );

					nonrem_state_indices = find(is_nonrem == 1);
					nonrem_remaining_states = find(is_nonrem == 0);

					is_nonrem_colors(nonrem_state_indices,:) = repmat(nonrem_color, length(nonrem_state_indices), 1);
					is_nonrem_colors(nonrem_remaining_states,:) = repmat(color_gray, length(nonrem_remaining_states), 1);

				% For REM
					rem_threshold_std = prctile(global_rem_ratio,80);

					is_rem = global_rem_ratio > rem_threshold_std;

					rem_axis = [min(global_rem_ratio) prctile(global_rem_ratio,95)];
					rem_colors = [];
					[rem_cmap, rem_colors] = assign_colors(global_rem_ratio, rem_axis, 10, 'magma' );

					rem_state_indices = find(is_rem == 1);
					rem_remaining_states = find(is_rem == 0);

					is_rem_colors(rem_state_indices,:) = repmat(rem_color, length(rem_state_indices), 1);
					is_rem_colors(rem_remaining_states,:) = repmat(color_gray, length(rem_remaining_states), 1);

					combined_rem_nonrem = sum([reshape(is_rem,[],1) reshape(is_nonrem,[],1)], 2);
			
					remaining_states_combined = find(combined_rem_nonrem == 0);

					combined_colors(rem_state_indices,:) = repmat(rem_color, length(rem_state_indices), 1);
					combined_colors(nonrem_state_indices,:) = repmat(nonrem_color, length(nonrem_state_indices), 1); 
					combined_colors(remaining_states_combined,:) = repmat(color_gray, length(remaining_states_combined), 1);

					all_pv = 1:global_states_count_lfp(5);


					ratio_id = 1;
					isrem_nonrem_id = 2;
										 

					for trial_iterator = [1,4]

						current_trial = char(trial_folders(trial_iterator));
						range_start = global_states_count_lfp(trial_iterator) + 1;
						range_end = global_states_count_lfp(trial_iterator + 1);
						pv_range = range_start : range_end; 

						
						non_pv_range = setdiff(all_pv, pv_range);

						x = projection_lfp(pv_range,1);
						y = projection_lfp(pv_range,2);
						

						trial_rem_colors = [] ; trial_nonrem_colors = [] ;
						trial_rem_colors = rem_colors(pv_range,:);
						trial_nonrem_colors = nonrem_colors(pv_range,:);

				

						figure(ratio_id)
						set(gcf, 'Position', get(0, 'Screensize'));
						subplot(2,4,trial_iterator)			
						scatter(projection_lfp(non_pv_range,1), projection_lfp(non_pv_range,2),5 ,color_gray, 'filled')
						hold on
						scatter(x,y,5, trial_rem_colors ,'filled');
						colormap(gca, rem_cmap)
						caxis(gca, rem_axis)
						h = colorbar;
						ylabel(h, 'Theta/Delta')
						title(sprintf('%s',current_trial))
						pbaspect([1 1 1]) 

						figure(ratio_id)
						set(gcf, 'Position', get(0, 'Screensize'));
						subplot(2,4,trial_iterator+4)			
						scatter(projection_lfp(non_pv_range,1), projection_lfp(non_pv_range,2),5 ,color_gray, 'filled')
						hold on
						scatter(x,y,5, trial_nonrem_colors ,'filled');
						colormap(gca,nonrem_cmap)
						caxis(gca,nonrem_axis)
						h = colorbar;
						ylabel(h, 'Delta x spindle')
						title(sprintf('%s',current_trial))
						pbaspect([1 1 1]) 
					


						figure(isrem_nonrem_id)
						set(gcf, 'Position', get(0, 'Screensize'));
						subplot(3,4,trial_iterator)			
						scatter(projection_lfp(non_pv_range,1), projection_lfp(non_pv_range,2),5 ,color_gray, 'filled')
						hold on
						scatter(x,y,5, is_rem_colors(pv_range,:) ,'filled');
						title(sprintf('%s',current_trial))
						pbaspect([1 1 1]) 

						subplot(3,4,trial_iterator+4)			
						scatter(projection_lfp(non_pv_range,1), projection_lfp(non_pv_range,2),5 ,color_gray, 'filled')
						hold on
						scatter(x,y,5, is_nonrem_colors(pv_range,:) ,'filled');
						title(sprintf('%s',current_trial))
						pbaspect([1 1 1]) 

						subplot(3,4,trial_iterator+8)			
						scatter(projection_lfp(non_pv_range,1), projection_lfp(non_pv_range,2),5 ,color_gray, 'filled')
						hold on
						scatter(x,y,5, combined_colors(pv_range,:) ,'filled');
						title(sprintf('%s',current_trial))
						pbaspect([1 1 1]) 


						
						if trial_iterator == total_trials
							figure(ratio_id)	
							suptitle(sprintf('%s Distribution of bands ratio on LFP', animal_id))
							filename = sprintf('%s_bands_ratio_on_LFP.png', animal_id );
							fullname = fullfile(lfp_overlay_output_folder,filename);
							saveas(gcf,fullname);

							figure(isrem_nonrem_id)
							suptitle(sprintf('%s Is REM / Non-REM LFP', animal_id))
							filename = sprintf('%s_isrem_on_LFP.png', animal_id );
							fullname = fullfile(lfp_overlay_output_folder,filename);
							saveas(gcf,fullname);

							if make_eps_images

								figure(ratio_id)
								epsname = sprintf('%s_bands_ratio_on_LFP.eps', animal_id );
								fullname_eps = fullfile(lfp_overlay_output_folder,epsname)
								set(gcf,'renderer','Painters')
								saveas(gcf,fullname_eps,'epsc')


								figure(isrem_nonrem_id)
								epsname = sprintf('%s_is_rem_overlay_lfp.eps', animal_id);
								fullname_eps = fullfile(lfp_overlay_output_folder,epsname)
								set(gcf,'renderer','Painters')
								saveas(gcf,fullname_eps,'epsc')

							end % End if make_eps_images


						end % End if trial_iterator


					end % End for trial_iterator	


					output_global_filename = fullfile(lfp_overlay_output_folder,'is_rem_nrem_lfp.mat');
					save(output_global_filename,'is_rem','is_nonrem')

			end % End if rem_nonrem_overlay_bins_lfp


		% -------------------------------------------------------------------
		% CSD Overlay (Median in Bins)
		% -------------------------------------------------------------------		
			if feature_overlay_bins_csd

				fprintf('Power Overlay in Bins...\n')

				csd_overlay_output_folder = fullfile(global_output_folder,'overlay_csd_power_bins');
				mkdir(csd_overlay_output_folder)

				no_xbins = 50;
				no_ybins = 50;
				total_bins = no_xbins * no_ybins;

				global_x = projection_csd(:,1);
				global_y = projection_csd(:,2);

				% Divide entire umap output into bins
				global_edges_x  = linspace( min(global_x), max(global_x), no_xbins+1 );
				global_edges_y  = linspace( min(global_y), max(global_y), no_ybins+1 );

				global_bin_power_all_bands = [];

				total_features = size(global_csd_states,2);
				
				% Compute density for all trials
				for trial_iterator = 1 : total_trials			

					range_start = global_states_count_csd(trial_iterator) + 1 ;
					range_end = global_states_count_csd(trial_iterator + 1) ;

					pv_range = range_start : range_end; 

					x = projection_csd(pv_range,1);
					y = projection_csd(pv_range,2);
					trial_duration = numel(x) * bin_size_sec;

					trial_csd_states = global_csd_states(pv_range,:);
 
					trial_median_bin_power = [];

					% Get xindices 
					[xcounts, x_edges, xindices] = histcounts(x, 'BinEdges', global_edges_x);
					[ycounts, y_edges, yindices] = histcounts(y, 'BinEdges', global_edges_y);


					for xbin_iterator = 1:no_xbins
						tempx = []; 
						tempx =	find(xindices == xbin_iterator);

						for ybin_iterator = 1:no_ybins
							tempy = [];
							
							tempy =	find(yindices == ybin_iterator);

							bin_states_indices = intersect(tempx,tempy);
							bsi = bin_states_indices;	

							if numel(bsi) < 1;
								continue;
							end			 
							median_power_bin_all_bands = [];
							median_power_bin_all_bands = median(trial_csd_states(bsi,:));

							trial_median_bin_power(bsi,1:total_features) = repmat(median_power_bin_all_bands, length(bsi), 1);
						end % End for ybin_iterator

					end % End for xbin_iterator	
					
					global_bin_power_all_bands = [global_bin_power_all_bands; trial_median_bin_power];

				end % End for trial_iterator

				band_counter = 1;
				global_median_power_in_trials = [];
				band_color = magma(6);

				for feature_iterator = 1:total_features
					fprintf('Processing Feature %d...\n', feature_iterator)

					caxis_range = [];
					caxis_range = [min(global_bin_power_all_bands(:,feature_iterator)), prctile(global_bin_power_all_bands(:,feature_iterator),95	)];

					% Plot density on same scale;
					% Compute density for all trials
					for trial_iterator = 1 : total_trials			

						range_start = global_states_count_csd(trial_iterator) + 1 ;
						range_end = global_states_count_csd(trial_iterator + 1) ;

						pv_range = range_start : range_end; 

						x = []; y = []; trial_band_power = []; trial_band_colors = [];

						x = projection_csd(pv_range,1);
						y = projection_csd(pv_range,2);
						trial_band_power = global_bin_power_all_bands(pv_range, feature_iterator);


						global_median_power_in_trials(feature_iterator, trial_iterator) = median(trial_band_power);

			
						[power_cmap, trial_band_colors] = assign_colors(trial_band_power, caxis_range, 1000, 'viridis'); 

						plot_overlay_maps = 0;
						if plot_overlay_maps 
							figure(123)
							set(gcf, 'Position', get(0, 'Screensize'));
							subplot(1,4,trial_iterator)
							scatter(projection_csd(:,1),projection_csd(:,2),5 ,color_gray, 'filled')
							hold on
							scatter(x,y,5, trial_band_colors ,'filled');	
							pbaspect([1 1 1])
							title(sprintf('%s', char(trial_folders(trial_iterator))))
							% if trial_iterator == 4
							% 	cbh = colorbar;
							% 	ylabel(cbh,'Median Power in Bins')
							% 	caxis(caxis_range)

							% end
							colormap(power_cmap)
							set(gca,'XColor','none','YColor','none')
							xlim([min(projection_csd(:,1)) max(projection_csd(:,1)) ])
							ylim([min(projection_csd(:,2)) max(projection_csd(:,2)) ])
						end % End if plot_overlay_maps

						plot_distribution = 0;
						if trial_iterator == 4 & plot_distribution
							figure(234)
							set(gcf, 'Position', get(0, 'Screensize'));
							subplot(3,6,feature_iterator)
							hist_temp = histogram(zscore(global_bin_power_all_bands(:,feature_iterator)), 100, 'FaceColor', band_color(band_counter,:), 'EdgeColor', band_color(band_counter,:));
							title(sprintf('%s %s', animal_id, x_ticks_array(feature_iterator)))
							caxis_range(2) = prctile(global_bin_power_all_bands(:,feature_iterator), 99);
							
							set(gca,'box','off') 
							xticks(caxis_range)
							xticklabels(fix(caxis_range))
							yaxis_range = [min(hist_temp.Values) max(hist_temp.Values)];

							yaxis_range = ceil(yaxis_range/10)*10;
							xaxis_range = ceil(caxis_range/10)*10;

							yticks(yaxis_range)
							yticklabels(yaxis_range)
							xticks(xaxis_range)
							xticklabels(xaxis_range)

								xlim(xaxis_range)
						end % End if trial_iterator

						if feature_iterator == 1 & trial_iterator == 1
							figure(185)
							set(gcf, 'Position', get(0, 'Screensize'));
							subplot(1,4,trial_iterator)
							scatter(projection_csd(:,1),projection_csd(:,2),5 ,'k', 'filled')
							xticks(global_edges_x)
							yticks(global_edges_y)
							grid on
							ax = gca;
							% ax.GridColor = [0 0 0];
							ax.GridLineStyle = '-';
							ax.GridColor = 'k';
							ax.GridAlpha = 1; % maximum line opacity
							set (gca, 'xticklabel' , {[]});
							set (gca, 'yticklabel' , {[]});
							pbaspect([1 1 1])
							xlim([min(projection_csd(:,1)) max(projection_csd(:,1)) ])
							ylim([min(projection_csd(:,2)) max(projection_csd(:,2)) ])
										

					

						end % End if feature_iterator


 

	 				end % End for trial_iterator



	 				if plot_overlay_maps 
		 				figure(123)
						suptitle(sprintf('%s %s ( %d-%d Hz ) Median Power on CSD', animal_id, x_ticks_array(feature_iterator), bands_array(band_counter,1), bands_array(band_counter,2)))
						filename = sprintf('%s_%s_CSD_power_overlay.png', animal_id, x_ticks_array(feature_iterator));;
						fullname = fullfile(csd_overlay_output_folder,filename);
						saveas(gcf,fullname);


						if make_pdf_images
							pdfname = sprintf('%s_%s_CSD_power_overlay.pdf', animal_id, x_ticks_array(feature_iterator));;
							fullname_pdf = fullfile(csd_overlay_output_folder,pdfname);
							figure(123)
							exportgraphics(gcf,fullname_pdf,'ContentType','image','BackgroundColor','none', 'Resolution',300)
						end


						close(figure(123))

					end % End if plot_overlay_maps

				

					

					if band_counter == 6
						band_counter = 1;
					else
						band_counter = band_counter + 1;
					end


				end % End for feature_iterator

				if plot_distribution

					% figure(234)
					% suptitle(sprintf('%s Distribution of Band Power CSD', animal_id ))
					% filename = sprintf('%s_CSD_power_distribution.png', animal_id);
					% fullname = fullfile(csd_overlay_output_folder,filename);

					% saveas(gcf,fullname);
					% % exportgraphics(gcf,fullname,'Resolution',300)

					% if make_pdf_images
					% 	pdfname = sprintf('%s_CSD_power_distribution.pdf', animal_id);
					% 	fullname_pdf = fullfile(csd_overlay_output_folder,pdfname)
					% 	figure(234)
					% 	exportgraphics(gcf,fullname_pdf,'ContentType','image','BackgroundColor','none', 'Resolution',300)
					% end
				end % End if plot_distribution
					
				output_global_filename = fullfile(csd_overlay_output_folder,'global_median_power_csd.mat');
				save(output_global_filename,'global_median_power_in_trials')



			end % End if feature_overlay_bins_csd


		
		% -------------------------------------------------------------------
		% CSD Overlay
		% -------------------------------------------------------------------		
			if feature_overlay_csd 	

				csd_overlay_output_folder = fullfile(global_output_folder,'overlay_CSD');
				mkdir(csd_overlay_output_folder)

				for trial_iterator = 1 : total_trials
					pv_range = []; 

					range_start = global_states_count_csd(trial_iterator) + 1 ;
					range_end = global_states_count_csd(trial_iterator + 1) ;

					pv_range = range_start : range_end; 

					x = projection_csd(pv_range,1);
					y = projection_csd(pv_range,2);			

					current_trial = char(trial_folders(trial_iterator));

					for layer_iterator = 1 : total_layers

						current_layer = layer_strings(layer_iterator);

						for band_iterator = 1:total_bands

							band_color = []; current_trial_band = [];

							band_index = (layer_iterator - 1) * total_bands + band_iterator;	
						
							cmin = min(global_csd_states(:,band_index));
							cmax = prctile(global_csd_states(:,band_index), 95);

							caxis_range = [cmin, cmax];
							current_trial_band = [];
							current_trial_band = global_csd_states(pv_range, band_index);

							[band_cmap, band_color] = assign_colors(current_trial_band, caxis_range, 10000, 'plasma');

							% % save median power for summary plots
							global_median_power_in_bands_across_trials_csd(band_index, trial_iterator) = median(current_trial_band);			

							% band_color = csd_color_matrix(pv_range,:,band_index);
							csd_overlay_fig_id = 1000 + layer_iterator + trial_iterator;
							plot_map = 0;
							if plot_map
								figure(csd_overlay_fig_id)
								set(gcf, 'Position', get(0, 'Screensize'));
								subplot(3,2,band_iterator)						
								scatter(projection_csd(:,1),projection_csd(:,2),5 ,color_gray, 'filled')
								hold on
								scatter(x,y,5, band_color ,'filled');
				
								colormap(band_cmap)
								caxis(caxis_range)
								h = colorbar();					
								ylabel(h, 'Median Binned Power')
								
								title(sprintf('%0.1f - %0.1f Hz', bands_array(band_iterator,1), bands_array(band_iterator,2)))
								pbaspect([1 1 1]) 
							end % End if plot_map
							% return
						end	% End for band_iterator

						if plot_map
							suptitle(sprintf('%s Trial %s Layer: %s CSD State Space', animal_id , current_trial, current_layer ))
							filename = sprintf('%s_%s_%s_csd_overlay.png', animal_id, current_trial, current_layer);
							fullname = fullfile(csd_overlay_output_folder,filename);
							saveas(csd_overlay_fig_id,fullname);

							close(csd_overlay_fig_id)
						end % End if plot_map

					end % End for layer_iterator

				end % End for trial_iterator

				% Save mean power in each bands across trials
				output_global_filename = fullfile(csd_overlay_output_folder,'global_median_power_in_bands_across_trials.mat');
				save(output_global_filename,'global_median_power_in_bands_across_trials_csd')

			end % End if feature_overlay CSD	


		% -------------------------------------------------------------------
		% LFP Overlay
		% -------------------------------------------------------------------		
			if feature_overlay_lfp 	

				lfp_overlay_output_folder = fullfile(global_output_folder,'overlay_lfp');
				mkdir(lfp_overlay_output_folder)

				for trial_iterator = 1 : total_trials

					range_start = global_states_count_lfp(trial_iterator) + 1 ;
					range_end = global_states_count_lfp(trial_iterator + 1) ;

					pv_range = range_start : range_end; 

					x = projection_lfp(pv_range,1);
					y = projection_lfp(pv_range,2);			

					current_trial = char(trial_folders(trial_iterator));

					for layer_iterator = 1 : total_layers

						current_layer = layer_strings(layer_iterator);

						for band_iterator = 1:total_bands

							band_color = [];

							band_index = (layer_iterator - 1) * total_bands + band_iterator;	
						
							cmin = min(global_lfp_states(:,band_index));
							cmax = prctile(global_lfp_states(:,band_index), 95);

							caxis_range = [cmin, cmax];
							current_trial_band = [];
							current_trial_band = global_lfp_states(pv_range, band_index);

							[band_cmap, band_color] = assign_colors(current_trial_band, caxis_range, 10000, 'viridis');

							% % save median power for summary plots
							global_median_power_in_bands_across_trials_lfp(band_index, trial_iterator) = median(current_trial_band);			

							% band_color = lfp_color_matrix(pv_range,:,band_index);
							plot_map = 0;
							if plot_map
								lfp_overlay_fig_id = 1000 + layer_iterator + trial_iterator;
								figure(lfp_overlay_fig_id)
								set(gcf, 'Position', get(0, 'Screensize'));
								subplot(3,2,band_iterator)						
								scatter(projection_lfp(:,1),projection_lfp(:,2),5 ,color_gray, 'filled')
								hold on
								scatter(x,y,5, band_color ,'filled');
				
								colormap(band_cmap)
								caxis(caxis_range)
								h = colorbar();					
								ylabel(h, 'Median Binned Power')
								
								title(sprintf('%0.1f - %0.1f Hz', bands_array(band_iterator,1), bands_array(band_iterator,2)))
								pbaspect([1 1 1]) 
								% 
							end % End if plot_map	
						end	% End for band_iterator

						if plot_map
							suptitle(sprintf('%s Trial %s Layer: %s LFP State Space', animal_id , current_trial, current_layer ))
							filename = sprintf('%s_%s_%s_lfp_overlay.png', animal_id, current_trial, current_layer);
							fullname = fullfile(lfp_overlay_output_folder,filename);
							saveas(lfp_overlay_fig_id,fullname);

							close(lfp_overlay_fig_id)
						end % End if plot_map	
					end % End for layer_iterator

				end % End for trial_iterator

				% Save mean power in each bands across trials
				output_global_filename = fullfile(lfp_overlay_output_folder,'global_median_power_in_bands_across_trials.mat');
				save(output_global_filename,'global_median_power_in_bands_across_trials_lfp')

			end % End if feature_overlay lfp


		% -------------------------------------------------------------------
		%  LFP Across Trials
		% -------------------------------------------------------------------		
			if lfp_overlay_across_trials 

				power_overlay_folder = fullfile(global_output_folder,'lfp_across_trials');
				mkdir(power_overlay_folder)

				for layer_iterator = 1 : total_layers

					current_layer = layer_strings(layer_iterator);

					for trial_iterator = 1 : total_trials

						range_start = global_states_count_lfp(trial_iterator) + 1;
						range_end = global_states_count_lfp(trial_iterator + 1);

						pv_range = range_start : range_end; 

						x = projection_lfp(pv_range,1);
						y = projection_lfp(pv_range,2);
						
						current_trial = char(trial_folders(trial_iterator));

						% return
						for band_iterator = 1:total_bands

							current_trial_band = []; band_color = [];

							band_index = (layer_iterator - 1) * total_bands + band_iterator;	

							cmin = min(global_lfp_states(:,band_index));
							cmax = prctile(global_lfp_states(:,band_index), 95);

							caxis_range = [cmin, cmax];

							current_trial_band = global_lfp_states(pv_range, band_index);

							[band_cmap, band_color] = assign_colors(current_trial_band, caxis_range, 10000, 'viridis');


							% New fig for each band
							fig_id = 1525 + band_iterator;
							figure(fig_id);
							set(gcf, 'Position', get(0, 'Screensize'));
							subplot(2,2,trial_iterator)						
							scatter(projection_lfp(:,1), projection_lfp(:,2), 5 ,color_gray, 'filled')
							hold on
							scatter(x,y,5, band_color ,'filled');
							colormap(band_cmap)
							colorbar
							caxis(caxis_range)
							h = colorbar;
							if do_rescale_global
								ylabel(h, 'Norm. Bin Power')
							else
								ylabel(h, 'Median Binned Power (uV^2)')
							end
							title(sprintf('%s', current_trial  ))
							pbaspect([1 1 1]) 


							if trial_iterator == total_trials
								suptitle(sprintf('%s %s - %0.1f - %0.1f Hz', animal_id, current_layer, bands_array(band_iterator,1), bands_array(band_iterator,2)))
								filename = sprintf('%s_%s_all_trials_band_%d_overlay.png', animal_id, current_layer, band_iterator );
								fullname = fullfile(power_overlay_folder,filename);
								saveas(figure(fig_id),fullname);
							end

						end	% End for band_iterator

					end % End for trial_iterator			

				end % End for layer_iterator

			end % End if lfp_overlay_across_trials 	


		% -------------------------------------------------------------------
		%  CSD Across Trials
		% -------------------------------------------------------------------	
			if csd_overlay_across_trials 

				csd_overlay_folder = fullfile(global_output_folder,'csd_across_trials');
				mkdir(csd_overlay_folder)

				for layer_iterator = 1 : total_layers

					current_layer = layer_strings(layer_iterator);

					for trial_iterator = 1 : total_trials

						
						range_start = global_states_count_csd(trial_iterator) + 1;
						range_end = global_states_count_csd(trial_iterator + 1);

						pv_range = range_start : range_end; 

						x = projection_csd(pv_range,1);
						y = projection_csd(pv_range,2);

						current_trial = char(trial_folders(trial_iterator));

						% return
						for band_iterator = 1:total_bands

							band_color = []; current_trial_band = [];

							band_index = (layer_iterator - 1) * total_bands + band_iterator;	
						
							cmin = min(global_csd_states(:,band_index));
							cmax = prctile(global_csd_states(:,band_index), 95);

							caxis_range = [cmin, cmax];

							current_trial_band = global_csd_states(pv_range, band_index);

							[band_cmap, band_color] = assign_colors(current_trial_band, caxis_range, 10000, 'plasma');

							% New fig for each band
							fig_id = 1525 + band_iterator;
							figure(fig_id);
							set(gcf, 'Position', get(0, 'Screensize'));
							subplot(2,2,trial_iterator)
							scatter(projection_csd(:,1), projection_csd(:,2), 5, color_gray, 'filled')
							hold on
							scatter(x,y,5, band_color ,'filled');			
							colormap(band_cmap)
							caxis(caxis_range)
							h = colorbar;
							ylabel(h, 'Median Binned Power')
							title(sprintf('%s', current_trial  ))
							pbaspect([1 1 1]) 

							if trial_iterator == total_trials
								suptitle(sprintf('%s %s - %0.1f - %0.1f Hz', animal_id, current_layer, bands_array(band_iterator,1), bands_array(band_iterator,2)))
								filename = sprintf('%s_%s_all_trials_band_%d_overlay.png', animal_id, current_layer, band_iterator);
								fullname = fullfile(csd_overlay_folder,filename);
								saveas(figure(fig_id),fullname);
							end

						end	% End for band_iterator

					end % End for trial_iterator

				end % End for layer_iterator

			end % End if csd_overlay_across_trials	
	

		% -------------------------------------------------------------------
		% Correlation of CSD bands
		% -------------------------------------------------------------------	
			if correlate_csd 

				correlate_csd_folder = fullfile(global_output_folder,'correlate_csd');
				mkdir(correlate_csd_folder)

				corr_vector_across_trials_csd = [];

				% x_ticks_array = [];
				total_freq_bands = length(bands_strings);
				
				global_labels = [];

				sorted_mat_across_trials = {};

				for trial_iterator = 1:total_trials

					% -------------------------------------------------------------------
					% Trial Specific Constants
					% -------------------------------------------------------------------				
						range_start = global_states_count_csd(trial_iterator) + 1 ;
						range_end = global_states_count_csd(trial_iterator + 1) ;

						pv_range = range_start : range_end; 
					
						trial_states_csd = [];
						trial_states_csd = global_csd_states(pv_range,:);

						current_trial = char(trial_folders(trial_iterator));

						corr_mat = corrcoef(trial_states_csd);						

						sorted_indices = sorted_band_indices;
						
						labels_sorted = x_ticks_array(sorted_indices);

						sorted_mat = corr_mat(sorted_indices, sorted_indices);

						sorted_mat_across_trials{1,trial_iterator} = sorted_mat;

						figure(654)
						set(gcf, 'Position', get(0, 'Screensize'));
						subplot(2,2, trial_iterator)
						imagesc(sorted_mat); 
						colormap(brewermap([],'rdylgn'));
						pbaspect([1 1 1])
						colorbar
						xticks(1:18)
						xticklabels(labels_sorted)
						yticks(1:18)
						yticklabels(labels_sorted)
						caxis([-1 1 ])
						title(current_trial)
						xtickangle(90)
						ax = gca;
						ax.XAxis.FontSize = 10 ;

					
					end % End for trial_iterator	

				suptitle(sprintf('%s CSD correlation matrix', animal_id ))

				filename = sprintf('%s_CSD_correlation_matrix.png', animal_id);
				fullname = fullfile(correlate_csd_folder,filename);
				saveas(gcf,fullname);

				% if make_eps_images

				% 	epsname = sprintf('%s_CSD_correlation_matrix.eps', animal_id);
				% 	fullname_eps = fullfile(correlate_csd_folder,epsname)
				% 	figure(654)	
				% 	set(gcf,'renderer','Painters')
				% 	saveas(gcf,fullname_eps,'epsc')

				% end % End if make_eps_images


				global_labels_indices = sorted_band_indices;

				
				output_filename = fullfile(correlate_csd_folder, 'corr_vector_across_trials_csd.mat');
				save(output_filename, 'x_ticks_array', 'global_labels_indices','sorted_mat_across_trials');

				
		
			end % End if correlate_csd


		% -------------------------------------------------------------------
		% Correlation of LFP bands
		% -------------------------------------------------------------------	
			if correlate_lfp 

				correlate_lfp_folder = fullfile(global_output_folder,'correlate_lfp');
				mkdir(correlate_lfp_folder)

				corr_vector_across_trials_lfp = [];

				% x_ticks_array = [];
				total_freq_bands = length(bands_strings);
				
				global_labels = [];

				sorted_mat_across_trials = {};

				for trial_iterator = 1:total_trials

					% -------------------------------------------------------------------
					% Trial Specific Constants
					% -------------------------------------------------------------------				
						range_start = global_states_count_lfp(trial_iterator) + 1 ;
						range_end = global_states_count_lfp(trial_iterator + 1) ;

						pv_range = range_start : range_end; 
					
						trial_states_lfp = [];
						trial_states_lfp = global_lfp_states(pv_range,:);

						current_trial = char(trial_folders(trial_iterator));

						corr_mat = corrcoef(trial_states_lfp);						

						sorted_indices = sorted_band_indices;
						
						labels_sorted = x_ticks_array(sorted_indices);

						sorted_mat = corr_mat(sorted_indices, sorted_indices);

						sorted_mat_across_trials{1,trial_iterator} = sorted_mat;

						figure(654)
						set(gcf, 'Position', get(0, 'Screensize'));
						subplot(2,2, trial_iterator)
						imagesc(sorted_mat); 
						colormap(brewermap([],'brbg'));
						pbaspect([1 1 1])
						colorbar
						xticks(1:18)
						xticklabels(labels_sorted)
						yticks(1:18)
						yticklabels(labels_sorted)
						caxis([-1 1 ])
						title(current_trial)
						xtickangle(90)
						ax = gca;
						ax.XAxis.FontSize = 10 ;

						total_features = 18;

						row_grids = [0 3 6 9 12 15 18];					
						hold on;
						for row = row_grids
						  line([0, total_features+0.5], [row+0.5, row+0.5], 'Color', 'k','lineWidth',2);
						end

						col_grids = [0 3 6 9 12 15 18];					
						hold on;
						for col = col_grids
						  line([col+0.5, col+0.5], [0, total_features+0.5],  'Color', 'k','lineWidth',2);
						end

					
					end % End for trial_iterator	

				suptitle(sprintf('%s LFP correlation matrix', animal_id ))

				filename = sprintf('%s_LFP_correlation_matrix.png', animal_id);
				fullname = fullfile(correlate_lfp_folder,filename);
				saveas(gcf,fullname);

				if make_eps_images

					epsname = sprintf('%s_LFP_correlation_matrix.eps', animal_id);
					fullname_eps = fullfile(correlate_lfp_folder,epsname)
					figure(654)	
					set(gcf,'renderer','Painters')
					saveas(gcf,fullname_eps,'epsc')

				end % End if make_eps_images


				global_labels_indices = sorted_band_indices;

				
				output_filename = fullfile(correlate_lfp_folder, 'corr_vector_across_trials_lfp.mat');
				save(output_filename, 'x_ticks_array', 'global_labels_indices','sorted_mat_across_trials');

				
		
			end % End if correlate_lfp


		% -------------------------------------------------------------------
		% Correlation of LFP bands in bins
		% ------------------------------------------------------------------
			if correlate_lfp_in_bins
				fprintf('Correlate LFP Bands in Bins.yewa..\n')

				corr_lfp_bins_folder = fullfile(global_output_folder,'correlate_lfp_in_bins');
				mkdir(corr_lfp_bins_folder)

				lfp_overlay_output_folder = fullfile(global_output_folder,'overlay_lfp_power_bins');
				load(fullfile(lfp_overlay_output_folder,'is_rem_nrem_lfp.mat'))


			
				no_xbins = 3;
				no_ybins = 3;
				total_bins = no_xbins * no_ybins;

				global_x = projection_lfp(:,1);
				global_y = projection_lfp(:,2);

				% Divide entire umap output into bins
				global_edges_x  = linspace( min(global_x), max(global_x), no_xbins+1 );
				global_edges_y  = linspace( min(global_y), max(global_y), no_ybins+1 );

				global_rem_ratio = rescale(global_lfp_states(:,theta_index) ./ global_lfp_states(:,delta_index));
				global_nonrem_ratio = rescale(global_lfp_states(:,delta_index) .* global_lfp_states(:,spindle_index));

				noisfree_awake_bin_ids = awake_bin_id_lfp{1,animal_iterator};
				noisefree_sleep_bin_ids = sleep_bin_id_lfp{1,animal_iterator};

				global_bin_sorted_matrix = [];
				global_bin_median_rem = [];
				global_bin_median_nonrem = [];
				global_bin_sleep_status = [];
				global_isrem_status = [];
				global_isnonrem_status = [];

				global_counter = 1;

				sorted_indices = sorted_band_indices;
				labels_sorted = x_ticks_array2(sorted_band_indices);

				% Compute density for all trials
				for trial_iterator = 1 : total_trials			

					range_start = global_states_count_lfp(trial_iterator) + 1 ;
					range_end = global_states_count_lfp(trial_iterator + 1) ;

					pv_range = range_start : range_end; 

					x = []; y = []; trial_rem = []; trial_nonrem = []; trial_lfp_states = []; 

					x = projection_lfp(pv_range,1);
					y = projection_lfp(pv_range,2);
					trial_lfp_states = global_lfp_states(pv_range,:);

					trial_isrem = is_rem(pv_range);
					trial_isnonrem = is_nonrem(pv_range);

					trial_rem = global_rem_ratio(pv_range);
					trial_nonrem = global_nonrem_ratio(pv_range);
					current_trial = char(trial_folders(trial_iterator));

					

					% Get xindices 
					[xcounts, x_edges, xindices] = histcounts(x, 'BinEdges', global_edges_x);
					[ycounts, y_edges, yindices] = histcounts(y, 'BinEdges', global_edges_y);

					

					for xbin_iterator = 1:no_xbins
						tempx = []; 
						tempx =	find(xindices == xbin_iterator);

						for ybin_iterator = 1:no_ybins
							tempy = [];
							bin_id = (xbin_iterator - 1) * no_xbins + ybin_iterator;
							tempy =	find(yindices == ybin_iterator);

							bin_states_indices = intersect(tempx,tempy);
							bsi = bin_states_indices;	
							half_size = 0.1 * length(bsi) ;

							if numel(bsi) == 0
								global_isrem_status(global_counter) = NaN;
								global_isnonrem_status(global_counter) = NaN;

								if ismember(trial_iterator,[2,3])
									sleep_status = 0;	
									global_bin_sleep_status(global_counter) = sleep_status;
								else
									sleep_status = 1;
									global_bin_sleep_status(global_counter) = sleep_status;
								end
								global_counter = global_counter + 1;
								continue;
							end		
							
							% Assign sleep status and ignore noisy bins in awake trials
							if ismember(trial_iterator,[2,3])
								sleep_status = 0;	

								if ~ismember(bin_id, noisfree_awake_bin_ids)
									global_counter = global_counter + 1;
									continue;
								end
								global_isrem_status(global_counter) = 0;
								global_isnonrem_status(global_counter) = 0;

									
							else 
								if ~ismember(bin_id, noisefree_sleep_bin_ids)
									sleep_status = 1;
									global_bin_sleep_status(global_counter) = sleep_status;
									global_counter = global_counter + 1;
									continue;
								end
								sleep_status = 1;

								if sum(trial_isrem(bsi)) > half_size
									global_isrem_status(global_counter) = 1;
								else
									global_isrem_status(global_counter) = 0;
								end

								if sum(trial_isnonrem(bsi)) > half_size
									global_isnonrem_status(global_counter) = 1;
								else
									global_isnonrem_status(global_counter) = 0;
								end
							
							end % End if ismember

							

							


							bin_lfp_states = []; bin_corr_matrix = []; bin_sorted_matrix = [];

							bin_lfp_states = trial_lfp_states(bsi,:);

							bin_corr_matrix = corrcoef(bin_lfp_states);

							bin_sorted_matrix = bin_corr_matrix(sorted_indices,sorted_indices);

							global_bin_sorted_matrix(:,:,global_counter) = bin_sorted_matrix;
							global_bin_median_rem(global_counter) = median(trial_rem(bsi));
							global_bin_median_nonrem(global_counter) = median(trial_nonrem(bsi));



							global_bin_sleep_status(global_counter) = sleep_status;

							global_counter = global_counter + 1;

							plot_binframes = 0;
							
							plot_corr_matrix = 0;
							
							if plot_binframes 

								figure_binsquare = 12321;
								figure(figure_binsquare)
								set(gcf, 'Position', get(0, 'Screensize'));
								subplot(no_xbins,no_ybins, bin_id)
								scatter(projection_lfp(:,1),projection_lfp(:,2),5 ,color_gray, 'filled')
								hold on
								scatter(x(bsi), y(bsi), 5 ,'k', 'filled')

								% [xmin, xmax] = bounds(x(bsi)); 
								% [ymin, ymax] = bounds(y(bsi)); 


								xmin = global_edges_x(xbin_iterator);
								xmax = global_edges_x(xbin_iterator+1) ; 
								ymin = global_edges_y(ybin_iterator);
								ymax = global_edges_y(ybin_iterator+1) ; 
								
								width_rectangle = abs(xmax - xmin);
								height_rectangle = abs(ymax - ymin);
								rectangle('Position',[xmin, ymin, width_rectangle, height_rectangle],'EdgeColor','k','LineWidth',1.5)

								title(sprintf('%s %s Bin: %d N-%d R-%d', animal_id,char(trial_folders(trial_iterator)), bin_id, global_isnonrem_status(global_counter-1), global_isrem_status(global_counter-1)))
								pbaspect([1 1 1])
								xticks(x_edges)
								yticks(y_edges)
								grid on
								ax = gca;
								ax.GridLineStyle = '-';
								ax.GridColor = 'k';
								ax.GridAlpha = 1; % maximum line opacity
								set (gca, 'xticklabel' , {[]});
								set (gca, 'yticklabel' , {[]});

							end % End if plot_binframes

							if plot_corr_matrix 

								figure_corr_matrix = 654;
								figure(figure_corr_matrix)
								set(gcf, 'Position', get(0, 'Screensize'));
								subplot(no_xbins, no_ybins, bin_id)
								imagesc(bin_sorted_matrix); 
								colormap(brewermap([],'brbg'));
								pbaspect([1 1 1])
								colorbar
								xticks(1:18)
								xticklabels(labels_sorted)
								yticks(1:18)
								yticklabels(labels_sorted)
								caxis([-1 1 ])
								title(sprintf('%s %s Bin: %d ', animal_id,char(trial_folders(trial_iterator)), bin_id))
								xtickangle(90)
								ax = gca;
								ax.XAxis.FontSize = 10 ;

								total_features = 18;

								row_grids = [0 3 6 9 12 15 18];					
								hold on;
								for row = row_grids
								  line([0, total_features+0.5], [row+0.5, row+0.5], 'Color', 'k','lineWidth',1.2);
								end

								col_grids = [0 3 6 9 12 15 18];					
								hold on;
								for col = col_grids
								  line([col+0.5, col+0.5], [0, total_features+0.5],  'Color', 'k','lineWidth',1.2);
								end

							end % End if plot_corr_matrix

						end % End for ybin_iterator

					end % End for xbin_iterator	


					if plot_binframes	
						figure(figure_binsquare)
						suptitle(sprintf('%s %s Bins ', animal_id, current_trial))
						filename = sprintf('%s_%s_lfp_bins.png', animal_id, current_trial);
						fullname = fullfile(corr_lfp_bins_folder,filename);
						saveas(gcf,fullname);
						
						if make_eps_images 
									
							epsname = sprintf('%s_%s_lfp_bins.eps', animal_id, current_trial);
							fullname_eps = fullfile(corr_lfp_bins_folder,epsname)
							figure(figure_binsquare)	
							set(gcf,'renderer','Painters')
							saveas(gcf,fullname_eps,'epsc')

						end % End if make_eps_images		

					end % End if plot_binframes
					

					if plot_corr_matrix	
						figure(figure_corr_matrix)
						suptitle(sprintf('%s %s Corr Matrix in Bins', animal_id, current_trial))
						filename = sprintf('%s_%s_corr_matrix_bins.png', animal_id, current_trial);
						fullname = fullfile(corr_lfp_bins_folder,filename);
						saveas(gcf,fullname);

						if make_eps_images 
							epsname = sprintf('%s_%s_corr_matrix_bins.eps', animal_id, current_trial);
							fullname_eps = fullfile(corr_lfp_bins_folder,epsname)
							figure(figure_corr_matrix)	
							set(gcf,'renderer','Painters')
							saveas(gcf,fullname_eps,'epsc')
						end % End if make_eps_images

					end % End if plot_corr_matrix	

					close all;


				end % End for trial_iterator

				output_filename = fullfile(corr_lfp_bins_folder, 'corr_vector_across_trials_lfp.mat');
				save(output_filename, 'x_ticks_array2', 'sorted_band_indices','global_bin_sorted_matrix',...
				 'global_bin_median_nonrem', 'global_bin_median_rem','global_bin_sleep_status','global_isrem_status','global_isnonrem_status');

			end % End if correlate_lfp_in_bins


		% -------------------------------------------------------------------
		% Correlation of CSD bands in bins
		% ------------------------------------------------------------------
			if correlate_csd_in_bins
				fprintf('Correlate CSD Bands in Bins...\n')
				corr_csd_bins_folder = fullfile(global_output_folder,'correlate_csd_in_bins');
				mkdir(corr_csd_bins_folder)

				no_xbins = 3;
				no_ybins = 3;
				total_bins = no_xbins * no_ybins;

				global_x = projection_csd(:,1);
				global_y = projection_csd(:,2);

				% Divide entire umap output into bins
				global_edges_x  = linspace( min(global_x), max(global_x), no_xbins+1 );
				global_edges_y  = linspace( min(global_y), max(global_y), no_ybins+1 );

				global_rem_ratio = rescale(global_lfp_for_csd(:,3) ./ global_lfp_for_csd(:,1));
				global_nonrem_ratio = rescale(global_lfp_for_csd(:,1) .* global_lfp_for_csd(:,2));

				noisfree_awake_bin_ids = awake_bin_id_csd{1,animal_iterator};
				noisefree_sleep_bin_ids = sleep_bin_id_csd{1,animal_iterator};

				global_bin_sorted_matrix = [];
				global_bin_median_rem = [];
				global_bin_median_nonrem = [];
				global_bin_sleep_status = [];

				global_counter = 1;

				sorted_indices = sorted_band_indices;
				labels_sorted = x_ticks_array2(sorted_band_indices);

				% Compute density for all trials
				for trial_iterator = 1 : total_trials			

					range_start = global_states_count_csd(trial_iterator) + 1 ;
					range_end = global_states_count_csd(trial_iterator + 1) ;

					pv_range = range_start : range_end; 

					x = []; y = []; trial_rem = []; trial_nonrem = []; trial_csd_states = []; 

					x = projection_csd(pv_range,1);
					y = projection_csd(pv_range,2);
					trial_csd_states = global_csd_states(pv_range,:);

					trial_rem = global_rem_ratio(pv_range);
					trial_nonrem = global_nonrem_ratio(pv_range);
					current_trial = char(trial_folders(trial_iterator));
					

					% Get xindices 
					[xcounts, x_edges, xindices] = histcounts(x, 'BinEdges', global_edges_x);
					[ycounts, y_edges, yindices] = histcounts(y, 'BinEdges', global_edges_y);

					

					for xbin_iterator = 1:no_xbins
						tempx = []; 
						tempx =	find(xindices == xbin_iterator);

						for ybin_iterator = 1:no_ybins
							tempy = [];
							bin_id = (xbin_iterator - 1) * no_xbins + ybin_iterator;
							tempy =	find(yindices == ybin_iterator);

							bin_states_indices = intersect(tempx,tempy);
							bsi = bin_states_indices;	

							if numel(bsi) == 10
								continue;
							end		
							
							% Assign sleep status and ignore noisy bins in awake trials
							if ismember(trial_iterator,[2,3])
								sleep_status = 0;	

								if ~ismember(bin_id, noisfree_awake_bin_ids)
									continue;
								end
							else 
								if ~ismember(bin_id, noisefree_sleep_bin_ids)
									continue;
								end
								sleep_status = 1;
							end % End if ismember


							bin_csd_states = []; bin_corr_matrix = []; bin_sorted_matrix = [];

							bin_csd_states = trial_csd_states(bsi,:);

							bin_corr_matrix = corrcoef(bin_csd_states);

							bin_sorted_matrix = bin_corr_matrix(sorted_indices,sorted_indices);

							global_bin_sorted_matrix(:,:,global_counter) = bin_sorted_matrix;
							global_bin_median_rem(global_counter) = median(trial_rem(bsi));
							global_bin_median_nonrem(global_counter) = median(trial_nonrem(bsi));


							global_bin_sleep_status(global_counter) = sleep_status;

							global_counter = global_counter + 1;

							plot_binframes = 1;
							if plot_binframes 

								figure_binsquare = 12321;
								figure(figure_binsquare)
								set(gcf, 'Position', get(0, 'Screensize'));
								subplot(no_xbins,no_ybins, bin_id)
								scatter(projection_csd(:,1),projection_csd(:,2),5 ,color_gray, 'filled')
								hold on
								scatter(x(bsi), y(bsi), 5 ,color_pink, 'filled')

								% [xmin, xmax] = bounds(x(bsi)); 
								% [ymin, ymax] = bounds(y(bsi)); 


								xmin = global_edges_x(xbin_iterator);
								xmax = global_edges_x(xbin_iterator+1) ; 
								ymin = global_edges_y(ybin_iterator);
								ymax = global_edges_y(ybin_iterator+1) ; 
								
								width_rectangle = abs(xmax - xmin);
								height_rectangle = abs(ymax - ymin);
								rectangle('Position',[xmin, ymin, width_rectangle, height_rectangle],'EdgeColor','k','LineWidth',3)

								title(sprintf('%s %s Bin: %d ', animal_id,char(trial_folders(trial_iterator)), bin_id))
								pbaspect([1 1 1])
								xticks(x_edges)
								yticks(y_edges)
								grid on
								ax = gca;
								ax.GridLineStyle = '-';
								ax.GridColor = 'k';
								ax.GridAlpha = 1; % maximum line opacity
								set (gca, 'xticklabel' , {[]});
								set (gca, 'yticklabel' , {[]});

							end % End if plot_binframes

							plot_corr_matrix = 1;
							if plot_corr_matrix 

								figure_corr_matrix = 654;
								figure(figure_corr_matrix)
								set(gcf, 'Position', get(0, 'Screensize'));
								subplot(no_xbins, no_ybins, bin_id)
								imagesc(bin_sorted_matrix); 
								colormap(brewermap([],'rdylgn'));
								pbaspect([1 1 1])
								colorbar
								xticks(1:18)
								xticklabels(labels_sorted)
								yticks(1:18)
								yticklabels(labels_sorted)
								caxis([-1 1 ])
								title(sprintf('%s %s Bin: %d ', animal_id,char(trial_folders(trial_iterator)), bin_id))
								xtickangle(90)
								ax = gca;
								ax.XAxis.FontSize = 10 ;

							end % End if plot_corr_matrix

						end % End for ybin_iterator

					end % End for xbin_iterator	

					
					if plot_binframes	
						figure(figure_binsquare)
						suptitle(sprintf('%s %s Bins ', animal_id, current_trial))
						filename = sprintf('%s_%s_csd_bins.png', animal_id, current_trial);
						fullname = fullfile(corr_csd_bins_folder,filename);
						saveas(gcf,fullname);
					end % End if plot_binframes
					

					if plot_corr_matrix	
						figure(figure_corr_matrix)
						suptitle(sprintf('%s %s Corr Matrix in Bins', animal_id, current_trial))
						filename = sprintf('%s_%s_corr_matrix_bins.png', animal_id, current_trial);
						fullname = fullfile(corr_csd_bins_folder,filename);
						saveas(gcf,fullname);
					end % End if plot_corr_matrix	

					close all;

				end % End for trial_iterator

				output_filename = fullfile(corr_csd_bins_folder, 'corr_vector_across_trials_csd.mat');
				save(output_filename, 'x_ticks_array2', 'sorted_band_indices','global_bin_sorted_matrix', 'global_bin_median_nonrem', 'global_bin_median_rem','global_bin_sleep_status');

			end % End if correlate_csd_in_bins


		% -------------------------------------------------------------------
		% Correlation of LFP bands over time windows 
		% -------------------------------------------------------------------	
			if correlate_lfp_time 	
				fprintf('Correlation Evolution Over Timeframes LFP...\n');
				foldername = strcat('correlate_lfp_time_',num2str(time_interval_for_correlation_evolution_sec))
				correlate_lfp_time_folder = fullfile(global_output_folder,foldername);
				mkdir(correlate_lfp_time_folder)

				lfp_raw_traces_folder = fullfile(global_output_folder,'lfp_raw_traces');
				mkdir(lfp_raw_traces_folder)

				sorted_indices = sorted_band_indices;

				

				corr_vector_across_trials_lfp = [];

				% x_ticks_array = [];
				total_freq_bands = length(bands_strings);
				
				global_labels = [];

				sorted_mat_across_trials = {};
				global_median_rem = {};
				global_median_nonrem = {};
				global_acceleration_variance = {};
				global_timeframe_indices = {};
				global_rem_peak_locations = {};
				global_timeframe_variance = {};
				global_timeframe_counts = [0];
				global_time_axis = {};

				global_rem_ratio = rescale(global_lfp_states(:,theta_index) ./ global_lfp_states(:,delta_index));
				global_nonrem_ratio = rescale(global_lfp_states(:,delta_index) .* global_lfp_states(:,spindle_index));

				for trial_iterator = 1:total_trials

					% -------------------------------------------------------------------
					% Trial Specific Constants
					% -------------------------------------------------------------------	
						range_start = global_states_count_lfp(trial_iterator) + 1 ;
						range_end = global_states_count_lfp(trial_iterator + 1) ;

						pv_range = [];
						pv_range = range_start : range_end; 
					
						trial_states_lfp = []; 
						trial_states_lfp = global_lfp_states(pv_range,:);

						trial_rawlfp_indices = global_raw_bin_indices_lfp(pv_range,:);

						current_trial = char(trial_folders(trial_iterator));

						trial_acc = global_acceleration_data_lfp(pv_range);

						% Binary Vector 1-Awake, 0-Sleep
						trial_accelerometer_binary = get_accelerometer_variance(trial_acc, acc_var_threshold);

						total_time_bins = length(pv_range);

						timepoints_for_corr_time = fix(1 : time_interval_for_correlation_evolution : total_time_bins);

						total_windows = length(timepoints_for_corr_time) - 1;


						% For plotting activity in specified time on state space
						plots_per_figure = 9;
						subplot_id = 1;
						figure_counter = 1;
						subplot_id_rt = 1;
						figure_counter_rt = 1;
						subplot_id_tfpb = 1;
						figure_counter_tfpb  = 1;

						
						plot_raw_traces = 1;
						plot_timeframes = 0;
						plot_timeframe_power_bands = 0;

						if plot_raw_traces
							fprintf('Loading Raw Traces\n')
							lfp_filepath = fullfile(root_directory, animal_id, csd_root_dirname, current_trial, lfp_filename);
							load(lfp_filepath)

							pyr_lfp = Raw_LFP(lfp_layer_indices(1),:);

							clear Raw_LFP

							raw_traces_colors = viridis(plots_per_figure);

						end % End if plot_raw_traces


						% Compute median REM and NonREM Ratio 
						rem_ratio = []; nonrem_ratio = [];
						rem_ratio = global_rem_ratio(pv_range);
						nonrem_ratio = global_nonrem_ratio(pv_range);

						temp_rem = [];
						temp_nonrem = [];
						temp_variance = [];
						sorted_matrix_array = [];
						temp_sleep_status = [];
						temp_timeframe_indices = [];
						temp_rem_peak_location = [];
						temp_time_axis = [];


						

						for timeframe_iterator = 1 : total_windows

							sorted_mat = []; corr_mat_temp = []; timeframe_variance = [];
							timeframe_sleep_status = [];

							time_start = timepoints_for_corr_time(timeframe_iterator);

							time_end = timepoints_for_corr_time(timeframe_iterator+1);		

							temp_time_axis(timeframe_iterator) = (time_start  + time_interval_for_correlation_evolution / 2 ) * bin_size_sec;				

							bands_in_time_range = trial_states_lfp(time_start:time_end,:);

							corr_mat_temp = corrcoef(bands_in_time_range);

							sorted_mat = corr_mat_temp(sorted_indices, sorted_indices);

							sorted_matrix_array(:,:,timeframe_iterator) = sorted_mat;

							% % Accelerometer data in Time frame
							timeframe_variance = var(trial_acc(time_start:time_end));
							
							if ismember(trial_iterator,[2,3])
								timeframe_sleep_status = 0;
								status_string = 'Awake';
							else
								timeframe_sleep_status = 1;
								status_string = 'Sleep';
							end

							temp_rem(timeframe_iterator)	 = median(rem_ratio(time_start:time_end));
							temp_nonrem(timeframe_iterator)	 = median(nonrem_ratio(time_start:time_end));
							temp_sleep_status(timeframe_iterator) = timeframe_sleep_status;
							temp_timeframe_indices(timeframe_iterator,:) = time_start:time_end;
							temp_variance(timeframe_iterator) = timeframe_variance;



							labels_sorted = x_ticks_array2(sorted_indices);			
							
						
							
							if plot_timeframes
								figure_id_temp = 234 ;
								figure(figure_id_temp)
								set(gcf, 'Position', get(0, 'Screensize'));	
								subplot(2,2,subplot_id)
								imagesc(sorted_mat); 
								colormap(brewermap([],'brbg'));
								pbaspect([1 1 1])
								colorbar
								xticks(1:18)
								xticklabels(labels_sorted)
								yticks(1:18)
								yticklabels(labels_sorted)
								caxis([-1 1 ])
								title(current_trial)
								xtickangle(90)
								ax = gca;
								ax.XAxis.FontSize = 10 ;

								title(sprintf('%d - %d sec' , fix(time_start*bin_size_sec),  fix(time_end * bin_size_sec)))

								pbaspect([1 1 1])


								subplot_id = subplot_id +1 ;

								if subplot_id > plots_per_figure | timeframe_iterator == total_windows
									subplot_id = 1;
									suptitle(sprintf('%s %s', animal_id ,current_trial))
									filename = sprintf('%s_%s_%d_.png', animal_id, current_trial, figure_counter);
									fullname = fullfile(correlate_lfp_time_folder, filename);
									saveas(figure_id_temp,fullname);
									figure_counter = figure_counter + 1;

								end 

							end % End if plot_timeframes

							

							if plot_raw_traces
								timeframe_raw_indices_edges = []; timeframe_lfp_indices = [];
								temp_lfp = []; jump = [];
								timeframe_raw_indices_edges = trial_rawlfp_indices(time_start:time_end,:);
								timeframe_lfp_indices = get_indices_from_edges(timeframe_raw_indices_edges);
							
								jump = diff(timeframe_lfp_indices) ;
								jump(jump>1 ) = 10000;
								jump = jump - 800;


								xjump = 1:numel(jump);
							
								temp_lfp = pyr_lfp(timeframe_lfp_indices);
								temp_lfp = temp_lfp - mean(temp_lfp);

								fig_ig_rt = 1212321321;
							
								figure(fig_ig_rt)
								set(gcf, 'Position', get(0, 'Screensize'));	
								subplot(10,1,subplot_id_rt)
								% plot(temp_lfp, 'Color', raw_traces_colors(subplot_id_rt,:),'lineWidth',1.5);
								plot(temp_lfp, 'Color', 'k','lineWidth',1.5);
								hold on
								plot(xjump, jump, 'Color','k','lineWidth', 3 );
								xlim([0 length(temp_lfp)]);
								ylim([-800 800])
								set(gca,'XColor', 'none','YColor','none')									
							

								title(sprintf('%d - %d sec', fix(time_start*bin_size_sec),  fix(time_end*bin_size_sec) ))
							

								if subplot_id_rt == plots_per_figure
									figure(fig_ig_rt)
									subplot(10,1,subplot_id_rt+1)
									plot([0; 0], [500; -500], '-k',  [0; 1000], [-500; -500], '-k', 'LineWidth', 2)
									xlim([0 length(temp_lfp)]);
									ylim([-800 800])
									set(gca, 'Visible', 'off')
								end

								subplot_id_rt = subplot_id_rt +1 ;

								if subplot_id_rt == plots_per_figure+1
									subplot_id_rt = 1;
									suptitle(sprintf('%s %s', animal_id ,current_trial))
									filename = sprintf('%s_%s_%d_.png', animal_id, current_trial, figure_counter_rt);
									fullname = fullfile(lfp_raw_traces_folder, filename);
									saveas(fig_ig_rt,fullname);
									figure_counter_rt = figure_counter_rt + 1;		


									if make_eps_images & figure_counter_rt - 1 == 22
									
										epsname = sprintf('%s_%s_%d_.eps', animal_id, current_trial, figure_counter_rt-1);
										fullname_eps = fullfile(lfp_raw_traces_folder,epsname)
										figure(fig_ig_rt)	
										set(gcf,'renderer','Painters')
										saveas(gcf,fullname_eps,'epsc')

									end % End if make_eps_images		

								end 

								

								

							end % End if plot_raw_traces


							if plot_timeframe_power_bands
								
								total_features = size(bands_in_time_range,2);
								total_timeframe_samples = size(bands_in_time_range,1);
								zs_bands_tr = zscore(bands_in_time_range)';

								[minzs maxzs] = bounds(zs_bands_tr(:));
								tfpb_fig_id = 23123;
								figure(tfpb_fig_id);
								set(gcf, 'Position', get(0, 'Screensize'));
								subplot(2,2,subplot_id_tfpb)
								imagesc(zs_bands_tr(sorted_band_indices,:));
								colormap('viridis')
								title(sprintf('%d - %d sec' , fix(time_start*bin_size_sec),  fix(time_end * bin_size_sec)))
								yticks(1:18)
								yticklabels(x_ticks_array2(sorted_band_indices))
								xticks(1:100:total_timeframe_samples)

								xl = (time_start:time_end)*bin_size_sec;
								xl = fix(xl(1:100:total_timeframe_samples));
								% xl = fix(time_start*bin_size_sec) : fix(time_end * bin_size_sec);
								% xl = fix(xl);
								xlabel('Time (sec)')
								ylabel('Bands')
								xticklabels(xl)
								cb2 = colorbar;
								ylabel(cb2, 'Zscored Binned Power')
								caxis([minzs maxzs])
								pbaspect([1 1 1])
								subplot_id_tfpb = subplot_id_tfpb +1 ;

								if subplot_id_tfpb > plots_per_figure | timeframe_iterator == total_windows
									subplot_id_tfpb = 1;
									suptitle(sprintf('%s %s Raw Bands', animal_id ,current_trial))
									filename = sprintf('%s_%s_%d_raw_bands.png', animal_id, current_trial, figure_counter_tfpb);
									fullname = fullfile(correlate_lfp_time_folder, filename);
									saveas(tfpb_fig_id,fullname);
									figure_counter_tfpb = figure_counter_tfpb + 1;
									
								end 

							end % End if plot_timeframe_power_bands	
							
						end % End for timeframe_iterator

						temp_rem = smooth_gaussian(temp_rem,50,200);

						temp_colors = viridis(3);
						figure(trial_iterator)
						subplot(3,1,1)
						plot(temp_time_axis, temp_rem,'Color',temp_colors(1,:), 'LineWidth',2)
						title(sprintf('%s REM ratio',current_trial))
						ylim([0 0.3])
						subplot(3,1,2)
						plot(temp_time_axis,temp_nonrem,'Color',temp_colors(2,:), 'LineWidth',2)
						title(sprintf('%s Non-REM ratio',current_trial))
						ylim([0 0.3])
						subplot(3,1,3)
						plot(temp_time_axis,temp_variance ,'Color',temp_colors(3,:), 'LineWidth',2)
						% plot( movmean(temp_sleep_status,3) ,'Color',temp_colors(3,:), 'LineWidth',2)
						title(sprintf('%s Accelerometer Variance',current_trial))
						ylim([0 0.3])

						suptitle(sprintf('%s LFP Time Frames', animal_id ))

						filename = sprintf('%s_rem_vs_acceleration_%s.png', animal_id, current_trial);
						fullname = fullfile(correlate_lfp_time_folder, filename);
						saveas(trial_iterator,fullname);
						close(trial_iterator)



						
						[temp_rem_peaks, temp_rem_peak_location ,w,p] = findpeaks(temp_rem, 'MinPeakHeight', 0.08, 'MinPeakProminence', 0.1);
						

						sorted_mat_across_trials{1,trial_iterator} = sorted_matrix_array;
						global_median_rem{1,trial_iterator} = temp_rem;
						global_median_nonrem{1,trial_iterator} = temp_nonrem;
						global_sleep_status{1,trial_iterator} = temp_sleep_status;
						global_rem_peak_locations{1,trial_iterator} = temp_rem_peak_location;
						global_timeframe_indices{1,trial_iterator} = temp_timeframe_indices;
						global_timeframe_variance{1,trial_iterator} = temp_variance;
						global_timeframe_counts = [global_timeframe_counts ; total_windows];
						global_time_axis{1,trial_iterator} = temp_time_axis;

						close all



						return
					end % End for trial_iterator	

					global_timeframe_counts = cumsum(global_timeframe_counts);
					
				output_filename = fullfile(correlate_lfp_time_folder, 'corr_vector_across_trials_lfp_over_time.mat');
				save(output_filename, 'x_ticks_array', 'sorted_band_indices','global_rem_peak_locations',...
					'sorted_mat_across_trials','global_sleep_status','global_median_rem',...
					'global_median_nonrem', 'global_timeframe_indices', 'time_interval_for_correlation_evolution',...
					'global_timeframe_variance','global_timeframe_counts','global_time_axis');
			

			end % End if correlate_lfp_time	

			
		% -------------------------------------------------------------------
		% Correlation of CSD bands over time windows 
		% -------------------------------------------------------------------	
			if correlate_csd_time 	
				fprintf('Correlation Evolution Over Timeframes CSD...\n');

				correlate_csd_time_folder = fullfile(global_output_folder,'correlate_csd_time');
				mkdir(correlate_csd_time_folder)

				csd_raw_traces_folder = fullfile(global_output_folder,'csd_raw_traces');
				mkdir(csd_raw_traces_folder)

				load(fullfile(root_directory,'AH2',output_folder_string,'correlate_csd', 'corr_vector_across_trials_csd.mat'));
				sorted_indices = global_labels_indices;

				

				corr_vector_across_trials_csd = [];

				% x_ticks_array = [];
				total_freq_bands = length(bands_strings);
				
				global_labels = [];

				sorted_mat_across_trials = {};
				global_median_rem = {};
				global_median_nonrem = {};
				global_acceleration_variance = {};
				global_timeframe_indices = {};
				global_rem_peak_locations = {};
				global_timeframe_variance = {};
				global_timeframe_counts = [0];
				global_time_axis = {};

				global_rem_ratio = rescale(global_lfp_for_csd(:,3) ./ global_lfp_for_csd(:,1));
				global_nonrem_ratio = rescale(global_lfp_for_csd(:,1) .* global_lfp_for_csd(:,2));

				for trial_iterator = 1:total_trials

					% -------------------------------------------------------------------
					% Trial Specific Constants
					% -------------------------------------------------------------------	
						range_start = global_states_count_csd(trial_iterator) + 1 ;
						range_end = global_states_count_csd(trial_iterator + 1) ;

						pv_range = [];
						pv_range = range_start : range_end; 
					
						trial_states_csd = []; 
						trial_states_csd = global_csd_states(pv_range,:);

						trial_rawcsd_indices = global_raw_bin_indices_csd{1,trial_iterator};

						current_trial = char(trial_folders(trial_iterator));

						trial_acc = global_acceleration_data_csd(pv_range);

						% Binary Vector 1-Awake, 0-Sleep
						trial_accelerometer_binary = get_accelerometer_variance(trial_acc, acc_var_threshold);

						total_time_bins = length(pv_range);

						timepoints_for_corr_time = fix(1 : time_interval_for_correlation_evolution : total_time_bins);

						total_windows = length(timepoints_for_corr_time) - 1;


						% For plotting activity in specified time on state space
						plots_per_figure = 9;
						subplot_id = 1;
						figure_counter = 1;
						subplot_id_rt = 1;
						figure_counter_rt = 1;
						subplot_id_tfpb = 1;
						figure_counter_tfpb  = 1;

						
						plot_raw_traces = 1;
						plot_timeframes = 1;
						plot_timeframe_power_bands = 1;

						if plot_raw_traces
							fprintf('Loading Raw Traces\n')
							csd_filepath = fullfile(root_directory, animal_id, csd_root_dirname, current_trial, csd_filename);
							load(csd_filepath)

							pyr_csd = Raw_CSD(csd_layer_indices(1),:);

							clear Raw_CSD

							raw_traces_colors = viridis(plots_per_figure);

						end % End if plot_raw_traces


						% Compute median REM and NonREM Ratio 
						rem_ratio = []; nonrem_ratio = [];
						rem_ratio = global_rem_ratio(pv_range);
						nonrem_ratio = global_nonrem_ratio(pv_range);

						temp_rem = [];
						temp_nonrem = [];
						temp_variance = [];
						sorted_matrix_array = [];
						temp_sleep_status = [];
						temp_timeframe_indices = [];
						temp_rem_peak_location = [];
						temp_time_axis = [];


						

						for timeframe_iterator = 1 : total_windows

							sorted_mat = []; corr_mat_temp = []; timeframe_variance = [];
							timeframe_sleep_status = [];

							time_start = timepoints_for_corr_time(timeframe_iterator);

							time_end = timepoints_for_corr_time(timeframe_iterator+1);		

							temp_time_axis(timeframe_iterator) = (time_start  + time_interval_for_correlation_evolution / 2 ) * bin_size_sec;				

							bands_in_time_range = trial_states_csd(time_start:time_end,:);

							corr_mat_temp = corrcoef(bands_in_time_range);

							sorted_mat = corr_mat_temp(sorted_indices, sorted_indices);

							sorted_matrix_array(:,:,timeframe_iterator) = sorted_mat;

							% % Accelerometer data in Time frame
							timeframe_variance = var(trial_acc(time_start:time_end));
							
							if ismember(trial_iterator,[2,3])
								timeframe_sleep_status = 0;
								status_string = 'Awake';
							else
								timeframe_sleep_status = 1;
								status_string = 'Sleep';
							end

							temp_rem(timeframe_iterator)	 = median(rem_ratio(time_start:time_end));
							temp_nonrem(timeframe_iterator)	 = median(nonrem_ratio(time_start:time_end));
							temp_sleep_status(timeframe_iterator) = timeframe_sleep_status;
							temp_timeframe_indices(timeframe_iterator,:) = time_start:time_end;
							temp_variance(timeframe_iterator) = timeframe_variance;



							labels_sorted = x_ticks_array2(sorted_indices);			
							
						
							
							if plot_timeframes
								figure_id_temp = 234 ;
								figure(figure_id_temp)
								set(gcf, 'Position', get(0, 'Screensize'));	
								subplot(3,3,subplot_id)
								imagesc(sorted_mat); 
								colormap(brewermap([],'brbg'));
								pbaspect([1 1 1])
								colorbar
								xticks(1:18)
								xticklabels(labels_sorted)
								yticks(1:18)
								yticklabels(labels_sorted)
								caxis([-1 1 ])
								title(current_trial)
								xtickangle(90)
								ax = gca;
								ax.XAxis.FontSize = 10 ;

								title(sprintf('%d - %d sec' , fix(time_start*bin_size_sec),  fix(time_end * bin_size_sec)))

								pbaspect([1 1 1])


								subplot_id = subplot_id +1 ;

								if subplot_id > plots_per_figure | timeframe_iterator == total_windows
									subplot_id = 1;
									suptitle(sprintf('%s %s CSD' , animal_id ,current_trial))
									filename = sprintf('%s_%s_%d_csd.png', animal_id, current_trial, figure_counter);
									fullname = fullfile(correlate_csd_time_folder, filename);
									saveas(figure_id_temp,fullname);
									figure_counter = figure_counter + 1;

								end 

							end % End if plot_timeframes

							

							if plot_raw_traces
								timeframe_raw_indices_edges = []; timeframe_csd_indices = [];
								temp_csd = []; jump = [];
								timeframe_raw_indices_edges = trial_rawcsd_indices(time_start:time_end,:);
								timeframe_csd_indices = get_indices_from_edges(timeframe_raw_indices_edges);
							
								jump = diff(timeframe_csd_indices) ;
								jump(jump>1 ) = 10000;
								jump = jump - 500;


								xjump = 1:numel(jump);
							
								temp_csd = pyr_csd(timeframe_csd_indices);
								temp_csd = temp_csd - mean(temp_csd);

								fig_ig_rt = 1212321321;
							
								figure(fig_ig_rt)
								set(gcf, 'Position', get(0, 'Screensize'));	
								subplot(10,1,subplot_id_rt)
								plot(temp_csd, 'Color', raw_traces_colors(subplot_id_rt,:),'lineWidth',1.5);
								hold on
								plot(xjump, jump, 'Color','k','lineWidth', 3 );
								xlim([0 length(temp_csd)]);
								ylim([-500 500])
								set(gca,'XColor', 'none','YColor','none')									
							

								title(sprintf('%d - %d sec', fix(time_start*bin_size_sec),  fix(time_end*bin_size_sec) ))
							

								if subplot_id_rt == plots_per_figure
									figure(fig_ig_rt)
									subplot(10,1,subplot_id_rt+1)
									plot([0; 0], [500; -500], '-k',  [0; 1000], [-500; -500], '-k', 'LineWidth', 2)
									xlim([0 length(temp_csd)]);
									ylim([-500 500])
									set(gca, 'Visible', 'off')
								end

								subplot_id_rt = subplot_id_rt +1 ;

								if subplot_id_rt == plots_per_figure+1
									subplot_id_rt = 1;
									suptitle(sprintf('%s %s CSD', animal_id ,current_trial))
									filename = sprintf('%s_%s_%d_csd.png', animal_id, current_trial, figure_counter_rt);
									fullname = fullfile(csd_raw_traces_folder, filename);
									saveas(fig_ig_rt,fullname);
									figure_counter_rt = figure_counter_rt + 1;									
								end 


								% return

							end % End if plot_raw_traces


							if plot_timeframe_power_bands
								
								total_features = size(bands_in_time_range,2);
								total_timeframe_samples = size(bands_in_time_range,1);
								zs_bands_tr = zscore(bands_in_time_range)';

								[minzs maxzs] = bounds(zs_bands_tr(:));
								tfpb_fig_id = 23123;
								figure(tfpb_fig_id);
								set(gcf, 'Position', get(0, 'Screensize'));
								subplot(3,3,subplot_id_tfpb)
								imagesc(zs_bands_tr);
								colormap('viridis')
								title(sprintf('%d - %d sec' , fix(time_start*bin_size_sec),  fix(time_end * bin_size_sec)))
								yticks(1:18)
								yticklabels(x_ticks_array2)
								xticks(1:20:total_timeframe_samples)

								xl = (time_start:time_end)*bin_size_sec;
								xl = fix(xl(1:20:total_timeframe_samples));
								% xl = fix(time_start*bin_size_sec) : fix(time_end * bin_size_sec);
								% xl = fix(xl);
								xlabel('Time (sec)')
								ylabel('Bands')
								xticklabels(xl)
								cb2 = colorbar;
								ylabel(cb2, 'Zscored Binned Power')
								caxis([minzs maxzs])
								pbaspect([1 1 1])
								subplot_id_tfpb = subplot_id_tfpb +1 ;

								if subplot_id_tfpb > plots_per_figure | timeframe_iterator == total_windows
									subplot_id_tfpb = 1;
									suptitle(sprintf('%s %s CSD Raw Bands', animal_id ,current_trial))
									filename = sprintf('%s_%s_%d_csd_raw_bands.png', animal_id, current_trial, figure_counter_tfpb);
									fullname = fullfile(correlate_csd_time_folder, filename);
									saveas(tfpb_fig_id,fullname);
									figure_counter_tfpb = figure_counter_tfpb + 1;
									
								end 

							end % End if plot_timeframe_power_bands	

						end % End for timeframe_iterator


						temp_colors = viridis(3);
						figure(trial_iterator)
						subplot(3,1,1)
						plot(temp_time_axis, temp_rem,'Color',temp_colors(1,:), 'LineWidth',2)
						title(sprintf('%s REM ratio',current_trial))
						ylim([0 0.3])
						subplot(3,1,2)
						plot(temp_time_axis,temp_nonrem,'Color',temp_colors(2,:), 'LineWidth',2)
						title(sprintf('%s Non-REM ratio',current_trial))
						ylim([0 0.3])
						subplot(3,1,3)
						plot(temp_time_axis,temp_variance ,'Color',temp_colors(3,:), 'LineWidth',2)
						% plot( movmean(temp_sleep_status,3) ,'Color',temp_colors(3,:), 'LineWidth',2)
						title(sprintf('%s Accelerometer Variance',current_trial))
						ylim([0 0.3])

						suptitle(sprintf('%s CSD Time Frames', animal_id ))

						filename = sprintf('%s_rem_vs_acceleration_%s.png', animal_id, current_trial);
						fullname = fullfile(correlate_csd_time_folder, filename);
						saveas(trial_iterator,fullname);
						close(trial_iterator)


						
						[temp_rem_peaks, temp_rem_peak_location ,w,p] = findpeaks(temp_rem, 'MinPeakHeight', 0.08, 'MinPeakProminence', 0.1);
						

						sorted_mat_across_trials{1,trial_iterator} = sorted_matrix_array;
						global_median_rem{1,trial_iterator} = temp_rem;
						global_median_nonrem{1,trial_iterator} = temp_nonrem;
						global_sleep_status{1,trial_iterator} = temp_sleep_status;
						global_rem_peak_locations{1,trial_iterator} = temp_rem_peak_location;
						global_timeframe_indices{1,trial_iterator} = temp_timeframe_indices;
						global_timeframe_variance{1,trial_iterator} = temp_variance;
						global_timeframe_counts = [global_timeframe_counts ; total_windows];
						global_time_axis{1,trial_iterator} = temp_time_axis;

						close all
					end % End for trial_iterator	

					global_timeframe_counts = cumsum(global_timeframe_counts);
					
				output_filename = fullfile(correlate_csd_time_folder, 'corr_vector_across_trials_csd_over_time.mat');
				save(output_filename, 'x_ticks_array', 'global_labels_indices','global_rem_peak_locations',...
					'sorted_mat_across_trials','global_sleep_status','global_median_rem',...
					'global_median_nonrem', 'global_timeframe_indices', 'time_interval_for_correlation_evolution',...
					'global_timeframe_variance','global_timeframe_counts','global_time_axis');
			

			end % End if correlate_csd_time	


		% -------------------------------------------------------------------
		% Make trajectory video CSD
		% -------------------------------------------------------------------	
			if make_video_csd

				trajectory_video_folder = fullfile(global_output_folder,'trajectory_videos');
				mkdir(trajectory_video_folder)

				for trial_iterator = 1 : total_trials

					current_trial = char(trial_folders(trial_iterator));

					

					range_start = global_states_count_csd(trial_iterator) + 1;
					range_end = global_states_count_csd(trial_iterator + 1);

					pv_range = range_start : range_end; 

					x = projection_csd(pv_range,1);
					y = projection_csd(pv_range,2);
					

					video_filename = sprintf('%s_dataset_csd_trial_%s.',animal_id, current_trial);
					video_savepath = fullfile(trajectory_video_folder,video_filename);

					
					PlotComet_2D(x(1:end), y(1:end), 'timeStep', bin_size_sec, 'outputPath', video_savepath, 'plotSpeed', plot_speed_vector(trial_iterator))
				


				end % End trial_iterator

			end % End if make video


		% -------------------------------------------------------------------
		% Make trajectory video LFP
		% -------------------------------------------------------------------	
			if make_video_lfp

				trajectory_video_folder = fullfile(global_output_folder,'trajectory_videos');
				mkdir(trajectory_video_folder)

				for trial_iterator = 1 : total_trials

					current_trial = char(trial_folders(trial_iterator));

					

					range_start = global_states_count_lfp(trial_iterator) + 1;
					range_end = global_states_count_lfp(trial_iterator + 1);

					pv_range = range_start : range_end; 

					x = projection_lfp(pv_range,1);
					y = projection_lfp(pv_range,2);
					

					video_filename = sprintf('%s_dataset_lfp_trial_%s.',animal_id, current_trial);
					video_savepath = fullfile(trajectory_video_folder,video_filename);

					
					PlotComet_2D(x(1:end), y(1:end), 'timeStep', bin_size_sec, 'outputPath', video_savepath, 'plotSpeed', plot_speed_vector(trial_iterator))
					


				end % End trial_iterator

			end % End if make video	
	

		% % -------------------------------------------------------------------
		% % REM overlay on CSD
		% % -------------------------------------------------------------------	
		% 	if rem_overlay_csd
		% 		fprintf('--REM on CSD...\n')
				

		% 		all_theta = global_lfp_for_csd(:,3);
		% 		all_delta = global_lfp_for_csd(:,1);

		% 		global_rem_ratio = rescale( all_theta ./ all_delta  );

		% 		lfp_rem_ratio = rescale(global_lfp_states(:,3) ./ global_lfp_states(:,1));
				
		% 		caxis_range = [0 prctile(lfp_rem_ratio,95)];
		% 		rem_colors = [];
		% 		[rem_cmap, rem_colors] = assign_colors(global_rem_ratio, caxis_range, 10, 'magma' );

		% 		rem_output_folder = fullfile(global_output_folder,'bands_ratio_overlay');
		% 		mkdir(rem_output_folder)

		% 		for trial_iterator = 1 : total_trials

		% 			current_trial = char(trial_folders(trial_iterator));
		% 			range_start = global_states_count_csd(trial_iterator) + 1;
		% 			range_end = global_states_count_csd(trial_iterator + 1);
		% 			pv_range = range_start : range_end; 

		% 			x = projection_csd(pv_range,1);
		% 			y = projection_csd(pv_range,2);

		% 			band_ratio_color = [];
										
		% 			band_ratio_color = rem_colors(pv_range,:);

		% 			figure(1001552)
		% 			set(gcf, 'Position', get(0, 'Screensize'));
		% 			subplot(2,2,trial_iterator)
				
		% 			scatter(projection_csd(:,1), projection_csd(:,2),5 ,color_gray, 'filled')
		% 			hold on
		% 			scatter(x,y,5, band_ratio_color ,'filled');		
		% 			colormap(rem_cmap)
		% 			caxis(caxis_range)
		% 			h = colorbar;
		% 			ylabel(h, 'Norm. Theta/Delta')
		% 			title(sprintf('%s',current_trial))
		% 			pbaspect([1 1 1]) 

		% 			if trial_iterator == total_trials	
		% 				suptitle(sprintf('%s Distribution of theta/delta ratio on manifold', animal_id ))
		% 				filename = sprintf('%s_rem_on_csd.png', animal_id );
		% 				fullname = fullfile(rem_output_folder,filename);
		% 				saveas(gcf,fullname);
		% 			end


		% 		end % End for trial_iterator	

		% 	end % End if rem_overlay_csd.


		% % -------------------------------------------------------------------
		% % Non-REM overlay on CSD
		% % -------------------------------------------------------------------	
		% 	if nonrem_overlay_csd
		% 		fprintf('--Non-REM on CSD...\n')
				
				
		% 		all_delta = global_lfp_for_csd(:,1);
		% 		all_spindle = global_lfp_for_csd(:,2);

		% 		lfp_nonrem_ratio = rescale(global_lfp_states(:,1) .* global_lfp_states(:,2));

		% 		global_nonrem_ratio = rescale(all_delta .* all_spindle);
				
		% 		caxis_range = [0 prctile(lfp_nonrem_ratio,95)];
		% 		nonrem_colors = [];
		% 		[nonrem_cmap, nonrem_colors] = assign_colors(global_nonrem_ratio, caxis_range, 100, 'pink' );

		% 		for trial_iterator = 1 : total_trials

		% 			current_trial = char(trial_folders(trial_iterator));
		% 			range_start = global_states_count_csd(trial_iterator) + 1;
		% 			range_end = global_states_count_csd(trial_iterator + 1);
		% 			pv_range = range_start : range_end; 
		% 			x = projection_csd(pv_range,1);
		% 			y = projection_csd(pv_range,2);
					

		% 			band_ratio_color = [];

		% 			current_trial = char(trial_folders(trial_iterator));

		% 			band_ratio_color = nonrem_colors(pv_range,:);


		% 			figure(10015521)
		% 			set(gcf, 'Position', get(0, 'Screensize'));
		% 			subplot(2,2,trial_iterator)
		% 			scatter(projection_csd(:,1), projection_csd(:,2),5 ,color_gray, 'filled')
		% 			hold on
		% 			scatter(x,y,5, band_ratio_color ,'filled');
		% 			colormap(nonrem_cmap)
		% 			colorbar
		% 			caxis(caxis_range)
		% 			h = colorbar;
		% 			ylabel(h, 'Norm. delta x spindle')
		% 			title(sprintf('%s',current_trial))
		% 			pbaspect([1 1 1]) 
					

		% 			if trial_iterator == total_trials	
		% 				suptitle(sprintf('%s Distribution of Delta x Spindle on CSD manifold', animal_id));
		% 				filename = sprintf('%s_nonrem_on_csd.png',animal_id);
		% 				fullname = fullfile(rem_output_folder,filename);
		% 				saveas(gcf,fullname);
		% 			end


		% 		end % End for trial_iterator

		% 	end % End if nonrem_overlay_csd


		% % -------------------------------------------------------------------
		% % REM overlay on LFP
		% % -------------------------------------------------------------------	
		% 	if rem_overlay_lfp
		% 		fprintf('--REM on LFP...\n')

		% 		rem_color = [227, 14, 85] / 255;	

		% 		all_theta = global_lfp_states(:,theta_index);
		% 		all_delta = global_lfp_states(:,delta_index);

		% 		global_rem_ratio = zscore( all_theta ./ all_delta  );
				
		% 		rem_threshold_std = prctile(global_rem_ratio,80)

		% 		is_rem = global_rem_ratio > rem_threshold_std;

		% 		caxis_range = [0 prctile(global_rem_ratio,95)];
		% 		rem_colors = [];
		% 		[rem_cmap, rem_colors] = assign_colors(global_rem_ratio, caxis_range, 10, 'magma' );

		% 		rem_state_indices = find(is_rem == 1);
		% 		remaining_states = find(is_rem == 0);

		% 		is_rem_colors(rem_state_indices,:) = repmat(rem_color, length(rem_state_indices), 1);
		% 		is_rem_colors(remaining_states,:) = repmat(color_gray, length(remaining_states), 1);

		% 		rem_output_folder = fullfile(global_output_folder,'bands_ratio_overlay');
		% 		mkdir(rem_output_folder)

		% 		all_pv = 1:global_states_count_lfp(5);

		% 		for trial_iterator = 1 : total_trials

		% 			current_trial = char(trial_folders(trial_iterator));
		% 			range_start = global_states_count_lfp(trial_iterator) + 1;
		% 			range_end = global_states_count_lfp(trial_iterator + 1);
		% 			pv_range = range_start : range_end; 

					
		% 			non_pv_range = setdiff(all_pv, pv_range);

		% 			x = projection_lfp(pv_range,1);
		% 			y = projection_lfp(pv_range,2);
					

		% 			band_ratio_color = [];
										
		% 			band_ratio_color = rem_colors(pv_range,:);

		% 			figure(1001552)
		% 			set(gcf, 'Position', get(0, 'Screensize'));
		% 			subplot(1,4,trial_iterator)			
		% 			scatter(projection_lfp(non_pv_range,1), projection_lfp(non_pv_range,2),5 ,color_gray, 'filled')
		% 			hold on
		% 			scatter(x,y,5, band_ratio_color ,'filled');
		% 			colormap(rem_cmap)
		% 			caxis(caxis_range)
		% 			h = colorbar;
		% 			ylabel(h, 'Theta/Delta')
		% 			title(sprintf('%s',current_trial))
		% 			pbaspect([1 1 1]) 



		% 			figure(1001512)
		% 			set(gcf, 'Position', get(0, 'Screensize'));
		% 			subplot(1,4,trial_iterator)			
		% 			scatter(projection_lfp(non_pv_range,1), projection_lfp(non_pv_range,2),5 ,color_gray, 'filled')
		% 			hold on
		% 			scatter(x,y,5, is_rem_colors(pv_range,:) ,'filled');
		% 			% colormap(is_rem_cmap)
		% 			% caxis(caxis_range)
		% 			% h = colorbar;
		% 			% ylabel(h, 'Theta/Delta')
		% 			title(sprintf('%s',current_trial))
		% 			pbaspect([1 1 1]) 

		% 			if trial_iterator == total_trials
		% 				figure(1001552)	
		% 				suptitle(sprintf('%s Distribution of theta/delta ratio on LFP manifold', animal_id))
		% 				filename = sprintf('%s_rem_on_LFP.png', animal_id );
		% 				fullname = fullfile(rem_output_folder,filename);
		% 				saveas(gcf,fullname);

		% 				figure(1001512)
		% 				suptitle(sprintf('%s Is REM LFP', animal_id))
		% 				filename = sprintf('%s_isrem_on_LFP.png', animal_id );
		% 				fullname = fullfile(rem_output_folder,filename);
		% 				saveas(gcf,fullname);

		% 				if make_eps_images

		% 					epsname = sprintf('%s_is_rem_overlay_lfp.eps', animal_id);
		% 					fullname_eps = fullfile(rem_output_folder,epsname)
		% 					figure(1001512)	
		% 					set(gcf,'renderer','Painters')
		% 					saveas(gcf,fullname_eps,'epsc')

		% 				end % End if make_eps_images

		% 			end




		% 		end % End for trial_iterator	

		% 	end % End if rem_overlay_lfp.

		% % -------------------------------------------------------------------
		% % Non-REM overlay on LFP
		% % -------------------------------------------------------------------	
		% 	if nonrem_overlay_lfp
		% 		fprintf('--Non-REM on LFP...\n')

		% 		nonrem_color = [26, 130, 150] / 255;	

		% 		all_spindle = global_lfp_states(:,spindle_index);
		% 		all_delta = global_lfp_states(:,delta_index);

		% 		global_nonrem_ratio = zscore( all_spindle .* all_delta  );
				
		% 		nonrem_threshold_std = prctile(global_nonrem_ratio,60)

		% 		is_nonrem = global_nonrem_ratio > nonrem_threshold_std;

		% 		caxis_range = [0 prctile(global_nonrem_ratio,95)];
		% 		nonrem_colors = [];
		% 		[nonrem_cmap, nonrem_colors] = assign_colors(global_nonrem_ratio, caxis_range, 10, 'pink' );

		% 		nonrem_state_indices = find(is_nonrem == 1);
		% 		nonremaining_states = find(is_nonrem == 0);

		% 		is_nonrem_colors(nonrem_state_indices,:) = repmat(nonrem_color, length(nonrem_state_indices), 1);
		% 		is_nonrem_colors(nonremaining_states,:) = repmat(color_gray, length(nonremaining_states), 1);

		% 		nonrem_output_folder = fullfile(global_output_folder,'bands_ratio_overlay');
		% 		mkdir(nonrem_output_folder)

		% 		all_pv = 1:global_states_count_lfp(5);

		% 		for trial_iterator = 1 : total_trials

		% 			current_trial = char(trial_folders(trial_iterator));
		% 			range_start = global_states_count_lfp(trial_iterator) + 1;
		% 			range_end = global_states_count_lfp(trial_iterator + 1);
		% 			pv_range = range_start : range_end; 

					
		% 			non_pv_range = setdiff(all_pv, pv_range);

		% 			x = projection_lfp(pv_range,1);
		% 			y = projection_lfp(pv_range,2);
					

		% 			band_ratio_color = [];
										
		% 			band_ratio_color = nonrem_colors(pv_range,:);

		% 			figure(1001552)
		% 			set(gcf, 'Position', get(0, 'Screensize'));
		% 			subplot(1,4,trial_iterator)			
		% 			scatter(projection_lfp(non_pv_range,1), projection_lfp(non_pv_range,2),5 ,color_gray, 'filled')
		% 			hold on
		% 			scatter(x,y,5, band_ratio_color ,'filled');
		% 			colormap(nonrem_cmap)
		% 			caxis(caxis_range)
		% 			h = colorbar;
		% 			ylabel(h, 'Theta/Delta')
		% 			title(sprintf('%s',current_trial))
		% 			pbaspect([1 1 1]) 



		% 			figure(1001512)
		% 			set(gcf, 'Position', get(0, 'Screensize'));
		% 			subplot(1,4,trial_iterator)			
		% 			scatter(projection_lfp(non_pv_range,1), projection_lfp(non_pv_range,2),5 ,color_gray, 'filled')
		% 			hold on
		% 			scatter(x,y,5, is_nonrem_colors(pv_range,:) ,'filled');
		% 			% colormap(is_nonrem_cmap)
		% 			% caxis(caxis_range)
		% 			% h = colorbar;
		% 			% ylabel(h, 'Theta/Delta')
		% 			title(sprintf('%s',current_trial))
		% 			pbaspect([1 1 1]) 

		% 			if trial_iterator == total_trials
		% 				figure(1001552)	
		% 				suptitle(sprintf('%s Distribution of spindle * delta ratio on LFP ', animal_id))
		% 				filename = sprintf('%s_nonrem_on_LFP.png', animal_id );
		% 				fullname = fullfile(nonrem_output_folder,filename);
		% 				saveas(gcf,fullname);

		% 				figure(1001512)
		% 				suptitle(sprintf('%s Is Non-REM LFP', animal_id))
		% 				filename = sprintf('%s_isnonrem_on_LFP.png', animal_id );
		% 				fullname = fullfile(nonrem_output_folder,filename);
		% 				saveas(gcf,fullname);

		% 				if make_eps_images

		% 					epsname = sprintf('%s_is_nonrem_overlay_lfp.eps', animal_id);
		% 					fullname_eps = fullfile(nonrem_output_folder,epsname)
		% 					figure(1001512)	
		% 					set(gcf,'renderer','Painters')
		% 					saveas(gcf,fullname_eps,'epsc')

		% 				end % End if make_eps_images

		% 			end




		% 		end % End for trial_iterator	

		% 	end % End if nonrem_overlay_lfp.
	
		% -------------------------------------------------------------------
		% Area Coverage, Speed and Acceleration on state space -- LFP
		% -------------------------------------------------------------------
			if area_coverage_lfp
				fprintf('Area Coverage LFP...\n')

				area_output_folder = fullfile(global_output_folder,'area_coverage_lfp');
				mkdir(area_output_folder);
			

				global_fraction_of_area = {};
				global_speed_of_coverage = {};
				global_acceleration_of_coverage = {};
				global_time_axis = {};
				global_median_rem = {};
				global_median_nonrem = {};
				global_median_acceleration_animal = {};
				global_rem_peak_locations = {};
				global_median_power = {};
				global_power_variance = {};


	
				% global_acceleration_data_lfp

				all_theta = global_lfp_states_og(:,theta_index);
				all_delta = global_lfp_states_og(:,delta_index);
				all_spindle = global_lfp_states_og(:,spindle_index);

				global_rem_ratio = rescale(all_theta ./ all_delta);
				global_nonrem_ratio = rescale(all_spindle .* all_delta);


				% Time interval between two measurements of area coverage
				% in no. of bins;
				% time_interval_for_area = 200;

				for trial_iterator = 1 : total_trials

					current_trial = char(trial_folders(trial_iterator));

					% Save snapshots of activity for each trial
					trial_area_folder  = fullfile(area_output_folder, current_trial);
					mkdir(trial_area_folder);

					current_trial = char(trial_folders(trial_iterator));

					

					range_start = global_states_count_lfp(trial_iterator) + 1;
					range_end = global_states_count_lfp(trial_iterator + 1);

					pv_range = range_start : range_end; 

					trial_lfp_states = global_lfp_states(pv_range,:);

					x = projection_lfp(pv_range,1);
					y = projection_lfp(pv_range,2);
					
					no_xbins = 50;
					no_ybins = 50;

					area_x_bins = linspace(min(projection_lfp(:,1)), max(projection_lfp(:,1)), no_xbins);
					area_y_bins = linspace(min(projection_lfp(:,2)), max(projection_lfp(:,2)), no_ybins);

					area_edges = {area_x_bins,  area_y_bins};		

					total_time_bins = length(x);

					fraction_of_area = [];
					inst_speed_coverage = [];
					temp_time_axis = [];
					temp_rem = [];
					temp_nonrem = [];
					temp_power_matrix = [];
					temp_power_variance_matrix = [];


					global_xy = [projection_lfp(:,1), projection_lfp(:,2)];

					xy_umap = [x y];

					% Check bins that are forever zero
					all_occupancy = hist3(global_xy, area_edges);

					forever_zeros_bins = length(find(all_occupancy(:) == 0));

					% Subtract bins that are forever zeros to get the acutal bins of state space
					total_area_bins = no_xbins * no_ybins - forever_zeros_bins; 

					timepoints_for_area = fix(1:time_interval_for_area:total_time_bins);


					% For plotting activity in specified time on state space
					plots_per_figure = 9;
					subplot_id = 1;
					figure_counter = 1;

					% Compute median REM and NonREM Ratio 
					rem_ratio = []; nonrem_ratio = [];
					rem_ratio = global_rem_ratio(pv_range);
					nonrem_ratio = global_nonrem_ratio(pv_range);

					temp_rem_peak_location = [];

					total_windows = length(timepoints_for_area) - 1;

					% Measure area covered at each of these timepoints
					for time_iterator = 1 : total_windows

						time_start = timepoints_for_area(time_iterator);

						time_end = timepoints_for_area(time_iterator+1);

						temp_time_axis(time_iterator) = (time_start  + time_interval_for_area / 2 ) * bin_size_sec;

						% state space occupancy for computing speed of coverage
						ss_occupancy = [];

						ss_occupancy = hist3(xy_umap(time_start:time_end,:) , 'Edges', area_edges);

						ss_occupancy(ss_occupancy  > 0 ) = 1;

						inst_speed_coverage(time_iterator) = sum(ss_occupancy(:)) ;



						% state space occupancy for computing speed of coverage
						ss_occupancy = [];

						ss_occupancy = hist3(xy_umap(1:time_end,:) , 'Edges', area_edges);

						ss_occupancy(ss_occupancy  > 0 ) = 1;

						fraction_of_area(time_iterator) = sum(ss_occupancy(:)) / total_area_bins;

						temp_rem(time_iterator)	 = median(rem_ratio(time_start:time_end));
						temp_nonrem(time_iterator)	 = median(nonrem_ratio(time_start:time_end));

						plot_timeframes = 0;

						


						if plot_timeframes
							plot_points = time_start : time_end;
							tf_colors = magma(length(plot_points));
							figure_id_temp = 234 ;
							figure(figure_id_temp)
							set(gcf, 'Position', get(0, 'Screensize'));	
							subplot(3,3,subplot_id)
							scatter(x,y, 4, color_gray, 'filled')
							hold on;
							scatter(x(plot_points), y(plot_points), 10, tf_colors, 'filled' );
							title(sprintf('%d - %d ', fix(time_start*bin_size_sec),  fix(time_end * bin_size_sec) ))
							xlim([-15 10])
							ylim([-10 10])
					

							pbaspect([1 1 1])

							subplot_id = subplot_id +1 ;

							if subplot_id > plots_per_figure | time_iterator == total_windows
								subplot_id = 1;
								suptitle(sprintf('%s', current_trial))
								filename = sprintf('%s_%s_%d_.png', animal_id, current_trial, figure_counter);
								fullname = fullfile(trial_area_folder,filename);
								saveas(figure_id_temp,fullname);
								figure_counter = figure_counter + 1;
								

								if make_eps_images & ismember(figure_counter - 1, [14,16])
									
									epsname = sprintf('%s_%s_%d_.eps', animal_id, current_trial, figure_counter-1);
									fullname_eps = fullfile(trial_area_folder,epsname);
									set(gcf,'renderer','Painters')
									saveas(gcf,fullname_eps,'epsc')

								end %End if make_eps_images
								
								close(figure_id_temp)

							end % End if subplot_id > plots_per_figure

						end % End if plot_timeframes

						timeframe_bands = [];

						timeframe_bands = trial_lfp_states(time_start:time_end,:);

						temp_power_matrix(:,time_iterator) = median(timeframe_bands);

						temp_power_variance_matrix(:,time_iterator) = var(timeframe_bands);



					end % End for time_iterator

			


					[temp_rem_peaks, temp_rem_peak_location ,w,p] = findpeaks(temp_rem, 'MinPeakHeight', 0.08, 'MinPeakProminence', 0.1);
	
					global_fraction_of_area{trial_iterator} = fraction_of_area;
					speed_coverage = inst_speed_coverage;
					acceleration_coverage = [0 diff(speed_coverage)] ;

					global_speed_of_coverage{trial_iterator} = speed_coverage;
					global_acceleration_of_coverage{trial_iterator} = acceleration_coverage;
					global_time_axis{trial_iterator} =  temp_time_axis;

					global_median_rem{trial_iterator} = temp_rem;
					global_median_nonrem{trial_iterator} = temp_nonrem;
					global_rem_peak_locations{1,trial_iterator} = temp_rem_peak_location;

					global_median_power{1,trial_iterator} = temp_power_matrix;
					global_power_variance{1,trial_iterator}  = temp_power_variance_matrix;
					
					temp_power_matrix = temp_power_matrix(sorted_band_indices,:);
					temp_power_matrix = [speed_coverage; temp_rem;  temp_nonrem; temp_power_matrix];	
					temp_power_matrix = zscore(temp_power_matrix, 0, 2)	;


					temp_power_variance_matrix = temp_power_variance_matrix(sorted_band_indices,:);
					temp_power_variance_matrix = [speed_coverage; temp_rem;  temp_nonrem; temp_power_variance_matrix];	
					temp_power_variance_matrix = zscore(temp_power_variance_matrix, 0, 2)	;


					total_features = size(trial_lfp_states,2);	

					zscored_speed = temp_power_matrix(1,:);
					corr_power_vs_speed = [];
					corr_var_vs_speed = [];

					% Compute correlation between coverage speed and median power
					for band_iterator = 1:total_features

						corr_power_vs_speed(band_iterator) = xcorr(zscored_speed, temp_power_matrix(band_iterator+3,:), 0,'coeff');

						corr_var_vs_speed(band_iterator) = xcorr(zscored_speed, temp_power_variance_matrix(band_iterator+3,:), 0,'coeff');

					end % End for band_iterator

					global_speed_vs_power{1,trial_iterator} = corr_power_vs_speed;
					global_speed_vs_variance{1,trial_iterator} = corr_var_vs_speed;


					temp_colors2 = cividis(2);
					fig_speed_vs_power = 852;
					figure(fig_speed_vs_power)
					set(gcf, 'Position', get(0, 'Screensize'));	
					subplot(2,2,trial_iterator)
					plot(corr_power_vs_speed,'Color',temp_colors2(1,:),'lineWidth',2 )
					hold on 
					plot(corr_var_vs_speed,'Color',temp_colors2(2,:),'lineWidth',2 )
					xticks(1:total_features);
					xticklabels(x_ticks_array(sorted_band_indices))
					xtickangle(90)
					ylabel('Correlation')
					ylim([-0.8 0.8])
					grid on
					legend('Median Power', 'Power Variance')
					pbaspect([1 1 1])
					title(current_trial)





					temp_labels = ["Coverage Speed", "REM Ratio", "Non-REM Ratio", x_ticks_array(sorted_band_indices)];


					figure_power_matrix = 980;
					figure(figure_power_matrix)
					set(gcf, 'Position', get(0, 'Screensize'));	
					subplot(2,2,trial_iterator)
					imagesc(temp_power_matrix)
					% colormap('cividis')
					colormap(brewermap([],'*rdgy'))
					maxval = max(abs(temp_power_matrix(:)));
					caxis([-maxval maxval])
					cbh = colorbar;
					ylabel(cbh, 'Zscored Rows')
					title(sprintf('%s %s Coverage Speed vs Median Power in Time Frames', animal_id, current_trial))
					yticks(1:total_features+3);
					yticklabels(temp_labels)



					xticks_indices = 1:50:total_windows;
					xlabels_temp = fix(temp_time_axis(xticks_indices)) -5;
					xticks(xticks_indices)
					xticklabels(xlabels_temp)
					xlabel('Time (seconds)')
					
					row_grids = [4:3:4+total_features];					
					hold on;
					for row = row_grids
					  line([1, total_windows], [row-0.5, row-0.5], 'Color', 'k','lineWidth',2);

					 
					end
					
			
					
					figure_power_variance_matrix = 981;
					figure(figure_power_variance_matrix)
					set(gcf, 'Position', get(0, 'Screensize'));	
					subplot(2,2,trial_iterator)
					imagesc(temp_power_variance_matrix)
					colormap(brewermap([],'prgn'))
					maxval = max(abs(temp_power_variance_matrix(:)));
					caxis([-maxval maxval])
					% colormap('magma')
					cbh = colorbar;
					ylabel(cbh,'Zscored Rows')
					title(sprintf('%s %s Coverage Speed vs Power Variance in Time Frames', animal_id, current_trial))
					yticks(1:size(temp_power_matrix,1));
					yticklabels(temp_labels)
					xlabel('Time (seconds)')
					yticks(1:total_features+3);
					yticklabels(temp_labels)

					xticks_indices = 1:fix(total_windows/5):total_windows;
					xticks(xticks_indices)
					xticklabels(fix(temp_time_axis(xticks_indices)))
					xlabel('Time (seconds)')
					
					row_grids = [4:3:4+total_features];					
					hold on;
					for row = row_grids
					  line([1, total_windows], [row-0.5, row-0.5], 'Color', 'k','lineWidth',2);
					end
					
			

					temp_colors = viridis(3);
					figure_binned_state_space = 12321;
					figure(figure_binned_state_space)
					subplot(2,2,trial_iterator)
					set(gcf, 'Position', get(0, 'Screensize'));	
					scatter(x,y, 5, color_gray, 'filled')
					title(sprintf('%s', char(trial_folders(trial_iterator))))
					xticks(area_x_bins)
					yticks(area_y_bins)
					grid on
					ax = gca;
					% ax.GridColor = [0 0 0];
					ax.GridLineStyle = '-';
					ax.GridColor = 'k';
					ax.GridAlpha = 0.5; % maximum line opacity
					set (gca, 'xticklabel' , {[]});
					set (gca, 'yticklabel' , {[]});
					pbaspect([1 1 1])
					xlim([min(area_x_bins) max(area_x_bins)])
					ylim([min(area_y_bins) max(area_y_bins)])


						
					fig_id_fraction = 123 ;
					figure(fig_id_fraction)
					subplot(2,2, trial_iterator)
					plot(temp_time_axis,fraction_of_area ,'lineWidth',2)
					title('Fraction of area swept');
					xlabel('Time (seconds)')
					ylabel('Fraction')
					ylim([0 1 ])
					grid on
					pbaspect([1 1 1])

					figure_inst_speed = 1123213;
					figure(figure_inst_speed)
					set(gcf, 'Position', get(0, 'Screensize'));	
					subplot(3,1, 1)
					plot(temp_time_axis, speed_coverage, 'Color',temp_colors(1,:), 'lineWidth', 2);
					title('Inst. Speed of Coverage');
					grid on
					xlabel('Time (seconds)')
					ylabel(sprintf('No. of bins swept / %d sec', fix(time_interval_for_area) * bin_size_sec ))

				
					subplot(3,1,2)
					plot(temp_time_axis, temp_rem,'Color',temp_colors(2,:), 'LineWidth',2)
					title(sprintf('%s REM ratio',current_trial))
					ylim([0 0.3])
					subplot(3,1,3)
					plot(temp_time_axis,temp_nonrem,'Color',temp_colors(3,:), 'LineWidth',2)
					title(sprintf('%s Non-REM ratio',current_trial))
					ylim([0 0.3])

		
					suptitle(sprintf('%s', current_trial))
					filename = sprintf('%s_%s_speed_rem_nonrem_ratios.png', animal_id, current_trial);
					fullname = fullfile(area_output_folder,filename);
					saveas(figure(figure_inst_speed),fullname);
					close(figure(figure_inst_speed))

				
					if trial_iterator == total_trials
						figure(fig_id_fraction)
						suptitle(sprintf('%s Fraction of Area Swept', animal_id ))
						filename = sprintf('%s_fraction_area.png', animal_id);
						fullname = fullfile(area_output_folder,filename);
						saveas(figure(fig_id_fraction),fullname);

						figure(figure_binned_state_space)
						suptitle(sprintf('%s Binned State Space', animal_id ))
						filename = sprintf('%s_binned_state_space.png', animal_id);
						fullname = fullfile(area_output_folder,filename);
						saveas(figure(figure_binned_state_space),fullname);

						if make_eps_images 
									
							epsname = sprintf('%s_binned_state_space.eps', animal_id);
							fullname_eps = fullfile(area_output_folder,epsname);
							set(gcf,'renderer','Painters')
							saveas(gcf,fullname_eps,'epsc')

						end %End if make_eps_images


						figure(figure_power_matrix)
						suptitle(sprintf('%s Coverage Speed vs Median Power in Time Frames', animal_id ))
						filename = sprintf('%s_coverage_speed_vs_median_power.png', animal_id);
						fullname = fullfile(area_output_folder,filename);
						saveas(figure(figure_power_matrix),fullname);

						if make_eps_images 
									
							epsname = sprintf('%s_coverage_speed_vs_median_power.eps', animal_id);
							fullname_eps = fullfile(area_output_folder,epsname);
							set(gcf,'renderer','Painters')
							saveas(gcf,fullname_eps,'epsc')

						end %End if make_eps_images


						figure(figure_power_variance_matrix)
						suptitle(sprintf('%s Coverage Speed vs Median Power Variance in Time Frames', animal_id ))
						filename = sprintf('%s_coverage_speed_vs_power_variance.png', animal_id);
						fullname = fullfile(area_output_folder,filename);
						saveas(figure(figure_power_variance_matrix),fullname);

						figure(fig_speed_vs_power)
						suptitle(sprintf('%s Coverage Speed vs Power Variance and Median Power ', animal_id))
						filename = sprintf('%s_corr_coverage_speed_vs_power_and_variance.png', animal_id);
						fullname = fullfile(area_output_folder,filename);
						saveas(figure(fig_speed_vs_power),fullname);


					end
					
				
				end % End for trial_iterator


				output_global_filename = fullfile(area_output_folder,'global_coverage_variables_lfp.mat');
				save(output_global_filename,'global_fraction_of_area','global_speed_of_coverage','global_rem_peak_locations','global_acceleration_of_coverage',...
				 'global_time_axis', 'global_median_rem','global_median_nonrem','global_median_power','global_power_variance','global_speed_vs_power',...
				 'global_speed_vs_variance' ,'sorted_band_indices','x_ticks_array','time_interval_for_area_sec','bin_size_sec');		
				
				close all
			end % End if area_coverage_lfp
			

		% -------------------------------------------------------------------
		% Area Coverage, Speed and Acceleration on state space -- CSD
		% -------------------------------------------------------------------
			if area_coverage_csd

				area_output_folder = fullfile(global_output_folder,'area_coverage_csd');
				mkdir(area_output_folder);
			

				global_fraction_of_area = {};
				global_speed_of_coverage = {};
				global_acceleration_of_coverage = {};
				global_time_axis = {};
				global_median_rem = {};
				global_median_nonrem = {};
				global_median_acceleration_animal = {};
				global_rem_peak_locations = {};


	
				% global_acceleration_data_csd

				global_rem_ratio = rescale(global_lfp_for_csd(:,3) ./ global_lfp_for_csd(:,1));
				global_nonrem_ratio = rescale(global_lfp_for_csd(:,1) .* global_lfp_for_csd(:,2));


				% Time interval between two measurements of area coverage
				% in no. of bins;
				% time_interval_for_area = 200;

				for trial_iterator = 1 : total_trials

					current_trial = char(trial_folders(trial_iterator));

					% Save snapshots of activity for each trial
					trial_area_folder  = fullfile(area_output_folder, current_trial);
					mkdir(trial_area_folder);

					current_trial = char(trial_folders(trial_iterator));

					

					range_start = global_states_count_csd(trial_iterator) + 1;
					range_end = global_states_count_csd(trial_iterator + 1);

					pv_range = range_start : range_end; 

					trial_csd_states = global_csd_states(pv_range,:);

					x = projection_csd(pv_range,1);
					y = projection_csd(pv_range,2);
					
					no_xbins = 50;
					no_ybins = 50;

					area_x_bins = linspace(min(projection_csd(:,1)), max(projection_csd(:,1)), no_xbins);
					area_y_bins = linspace(min(projection_csd(:,2)), max(projection_csd(:,2)), no_ybins);

					area_edges = {area_x_bins,  area_y_bins};		

					total_time_bins = length(x);

					fraction_of_area = [];
					inst_speed_coverage = [];
					temp_time_axis = [];
					temp_rem = [];
					temp_nonrem = [];


					global_xy = [projection_csd(:,1), projection_csd(:,2)];

					xy_umap = [x y];

					% Check bins that are forever zero
					all_occupancy = hist3(global_xy, area_edges);

					forever_zeros_bins = length(find(all_occupancy(:) == 0));

					% Subtract bins that are forever zeros to get the acutal bins of state space
					total_area_bins = no_xbins * no_ybins - forever_zeros_bins; 

					timepoints_for_area = fix(1:time_interval_for_area:total_time_bins);


					% For plotting activity in specified time on state space
					plots_per_figure = 9;
					subplot_id = 1;
					figure_counter = 1;

					% Compute median REM and NonREM Ratio 
					rem_ratio = []; nonrem_ratio = [];
					rem_ratio = global_rem_ratio(pv_range);
					nonrem_ratio = global_nonrem_ratio(pv_range);

					temp_rem_peak_location = [];



					% Measure area covered at each of these timepoints
					for time_iterator = 1 : length(timepoints_for_area) - 1 

						time_start = timepoints_for_area(time_iterator);

						time_end = timepoints_for_area(time_iterator+1);

						temp_time_axis(time_iterator) = (time_start  + time_interval_for_area / 2 ) * bin_size_sec;

						% state space occupancy for computing speed of coverage
						ss_occupancy = [];

						ss_occupancy = hist3(xy_umap(time_start:time_end,:) , 'Edges', area_edges);

						ss_occupancy(ss_occupancy  > 0 ) = 1;

						inst_speed_coverage(time_iterator) = sum(ss_occupancy(:)) ;



						% state space occupancy for computing speed of coverage
						ss_occupancy = [];

						ss_occupancy = hist3(xy_umap(1:time_end,:) , 'Edges', area_edges);

						ss_occupancy(ss_occupancy  > 0 ) = 1;

						fraction_of_area(time_iterator) = sum(ss_occupancy(:)) / total_area_bins;

						temp_rem(time_iterator)	 = median(rem_ratio(time_start:time_end));
						temp_nonrem(time_iterator)	 = median(nonrem_ratio(time_start:time_end));

						plot_timeframes = 0;

						if plot_timeframes
							plot_points = time_start : time_end;
							figure_id_temp = 234 ;
							figure(figure_id_temp)
							set(gcf, 'Position', get(0, 'Screensize'));	
							subplot(3,3,subplot_id)
							scatter(x,y, 4, color_gray, 'filled')
							hold on;
							scatter(x(plot_points), y(plot_points), 10, color_pink, 'filled' );
							title(sprintf('%d - %d ', fix(time_start*bin_size_sec),  fix(time_end * bin_size_sec) ))

							pbaspect([1 1 1])


							subplot_id = subplot_id +1 ;

							if subplot_id > plots_per_figure
								subplot_id = 1;
								suptitle(sprintf('%s', current_trial))
								filename = sprintf('%s_%s_%d_.png', animal_id, current_trial, figure_counter);
								fullname = fullfile(trial_area_folder,filename);
								saveas(figure_id_temp,fullname);
								figure_counter = figure_counter + 1;
								close(figure_id_temp)
							end 

						end % End if plot_timeframes

					end % End for time_iterator


					[temp_rem_peaks, temp_rem_peak_location ,w,p] = findpeaks(temp_rem, 'MinPeakHeight', 0.08, 'MinPeakProminence', 0.1);
	
					global_fraction_of_area{trial_iterator} = fraction_of_area;
					speed_coverage = inst_speed_coverage;
					acceleration_coverage = [0 diff(speed_coverage)] ;

					global_speed_of_coverage{trial_iterator} = speed_coverage;
					global_acceleration_of_coverage{trial_iterator} = acceleration_coverage;
					global_time_axis{trial_iterator} =  temp_time_axis;

					global_median_rem{trial_iterator} = temp_rem;
					global_median_nonrem{trial_iterator} = temp_nonrem;
					global_rem_peak_locations{1,trial_iterator} = temp_rem_peak_location;


					temp_colors = viridis(3);

					figure(12321)
					subplot(2,2,trial_iterator)
					set(gcf, 'Position', get(0, 'Screensize'));	
					scatter(x,y, 5, 'k', 'filled')
					title(sprintf('%s', char(trial_folders(trial_iterator))))
					xticks(area_x_bins)
					yticks(area_y_bins)
					grid on
					ax = gca;
					% ax.GridColor = [0 0 0];
					ax.GridLineStyle = '-';
					ax.GridColor = 'k';
					ax.GridAlpha = 1; % maximum line opacity
					set (gca, 'xticklabel' , {[]});
					set (gca, 'yticklabel' , {[]});
					pbaspect([1 1 1])


					
					fig_id_fraction = 123 ;
					figure(fig_id_fraction)
					subplot(2,2, trial_iterator)
					plot(temp_time_axis,fraction_of_area ,'lineWidth',2)
					title('Fraction of area swept');
					xlabel('Time (seconds)')
					ylabel('Fraction')
					ylim([0 1 ])
					grid on
					pbaspect([1 1 1])


					figure(1123213)
					set(gcf, 'Position', get(0, 'Screensize'));	
					subplot(3,1, 1)
					plot(temp_time_axis, speed_coverage, 'Color',temp_colors(1,:), 'lineWidth', 2);
					title('Inst. Speed of Coverage');
					grid on
					xlabel('Time (seconds)')
					ylabel(sprintf('No. of bins swept / %d sec', fix(time_interval_for_area) * bin_size_sec ))

				
					subplot(3,1,2)
					plot(temp_time_axis, temp_rem,'Color',temp_colors(2,:), 'LineWidth',2)
					title(sprintf('%s REM ratio',current_trial))
					ylim([0 0.3])
					subplot(3,1,3)
					plot(temp_time_axis,temp_nonrem,'Color',temp_colors(3,:), 'LineWidth',2)
					title(sprintf('%s Non-REM ratio',current_trial))
					ylim([0 0.3])


					% subplot(3,1, 3)
					% plot(temp_time_axis, acceleration_coverage, 'lineWidth', 2);
					% grid on
					% title('Inst. Acceleration of Coverage');
					% xlabel('Time (seconds)')
		
					suptitle(sprintf('%s', current_trial))
					filename = sprintf('%s_%s_speed_and_acceleration.png', animal_id, current_trial);
					fullname = fullfile(area_output_folder,filename);
					saveas(figure(1123213),fullname);
					close(figure(1123213))

				
					if trial_iterator == total_trials
						suptitle(sprintf('%s Fraction of Area Swept', animal_id ))
						filename = sprintf('%s_fraction_area.png', animal_id);
						fullname = fullfile(area_output_folder,filename);
						saveas(figure(fig_id_fraction),fullname);


						suptitle(sprintf('%s Binned State Space', animal_id ))
						filename = sprintf('%s_binned_state_space.png', animal_id);
						fullname = fullfile(area_output_folder,filename);
						saveas(figure(12321),fullname);
					end
					
				
				end % End for trial_iterator

				output_global_filename = fullfile(area_output_folder,'global_coverage_variables_csd.mat');
				save(output_global_filename,'global_fraction_of_area','global_speed_of_coverage','global_rem_peak_locations','global_acceleration_of_coverage', 'global_time_axis', 'global_median_rem','global_median_nonrem' );		
				
				close all
			end % End if area_coverage_csd
			
		% -------------------------------------------------------------------
		%  Quantify movements LFP
		% -------------------------------------------------------------------	
			if quantify_movements_lfp

				quantify_movements_folder = fullfile(global_output_folder,'quantify_movements_awake_lfp');
				mkdir(quantify_movements_folder)
		

				for trial_iterator = 1 : total_trials

					current_trial = char(trial_folders(trial_iterator));
			

					

					range_start = global_states_count(trial_iterator) + 1;
					range_end = global_states_count(trial_iterator + 1);

					pv_range = range_start : range_end; 

					population_vector = global_lfp_states(pv_range,:);
					total_pop_vec = length(pv_range);

					x = projection_lfp(pv_range,1);
					y = projection_lfp(pv_range,2);
					
				
					time_color = magma(length(x));

					total_time_regime = 25;

					time_regime = floor(linspace(1,total_pop_vec, total_time_regime+1));

					regime_median_power = [];

					temp_time_axis = [];

					for regime_iterator = 1 : total_time_regime

						% Lower and upper limit
						regime_ll = time_regime(regime_iterator);
						regime_ul = time_regime(regime_iterator+1);

						regime_indices = regime_ll : regime_ul;

						regime_median_power(regime_iterator,:) = median(population_vector(regime_indices,:));	


						time_start = floor(regime_ll * bin_size_sec - bin_size_sec); 
						time_end = floor(regime_ul * bin_size_sec);

						temp_time_axis(regime_iterator) = fix(time_start  + (time_end-time_start)/2);

						fig_id = 5311 + trial_iterator;
						figure(fig_id)
						set(gcf, 'Position', get(0, 'Screensize'))
						subplot(5,5,regime_iterator)
						scatter(x,y,3,color_gray, 'filled')
						hold on
						scatter(x(regime_indices), y(regime_indices) ,3, time_color(regime_indices,:) ,'filled');
						title(sprintf('%d - %d sec', time_start, time_end));
						pbaspect([1 1 1]) 


						if regime_iterator == total_time_regime
							suptitle(sprintf('%s %s Time overlay on state space', animal_id, char(trial_folders(trial_iterator)) ));
							filename = sprintf('%s_%s_time_overlay_regimes.png', animal_id, char(trial_folders(trial_iterator)) );
							fullname = fullfile(quantify_movements_folder,filename);
							saveas(fig_id,fullname);
						end 


						regime_ll = [];
						regime_ul = [];
						regime_indices = [];

					end % End for regime_iterator
					close(fig_id)


					fig_id = 221 + trial_iterator;

					% Plot each band median power in regime
					for band_iterator = 1:total_bands
						% smooth median power
						regime_median_power(:,band_iterator) = zscore(smooth_gaussian(regime_median_power(:,band_iterator), 1, 2));

						figure(fig_id)
						subplot(3,2, band_iterator)
						plot(temp_time_axis, regime_median_power(:,band_iterator), 'lineWidth', 2);
						title(char(bands_strings(band_iterator)))
						grid on



					end % End for band_iterator

					suptitle(sprintf('%s %s Median power in time regimes', animal_id,  char(trial_folders(trial_iterator)) ));

					filename = sprintf('%s_trial_%s_median_power_in_time.png',animal_id, char(trial_folders(trial_iterator)) );
					fullname = fullfile(quantify_movements_folder,filename);
					saveas(gcf,fullname);

					global_regime_median_power{trial_iterator} = regime_median_power;

					time_regime = [];


				end % End for trial_iterator



				output_global_filename = fullfile(quantify_movements_folder,'global_regime_median_power.mat');
				save(output_global_filename,'global_regime_median_power')		
				
				global_regime_median_power = [];
			end % End if quantify_movements_lfp		


		% -------------------------------------------------------------------
		%  Quantify movements csd
		% -------------------------------------------------------------------	
			if quantify_movements_csd

				quantify_movements_folder = fullfile(global_output_folder,'quantify_movements_awake_csd');
				mkdir(quantify_movements_folder)
		

				for trial_iterator = 1 : total_trials

					current_trial = char(trial_folders(trial_iterator));
			

					

					range_start = global_states_count(trial_iterator) + 1;
					range_end = global_states_count(trial_iterator + 1);

					pv_range = range_start : range_end; 

					population_vector = global_csd_states(pv_range,:);
					total_pop_vec = length(pv_range);

					x = projection_csd(pv_range,1);
					y = projection_csd(pv_range,2);
					
				
					time_color = magma(length(x));

					total_time_regime = 25;

					time_regime = floor(linspace(1,total_pop_vec, total_time_regime+1));

					regime_median_power = [];

					for regime_iterator = 1 : total_time_regime

						% Lower and upper limit
						regime_ll = time_regime(regime_iterator);
						regime_ul = time_regime(regime_iterator+1);

						regime_indices = regime_ll : regime_ul;

						regime_median_power(regime_iterator,:) = median(population_vector(regime_indices,:));	


						time_start = floor(regime_ll * bin_size_sec - bin_size_sec); 
						time_end = floor(regime_ul * bin_size_sec);

						temp_time_axis(regime_iterator) = fix(time_start  + (time_end-time_start)/2);

						fig_id = 5311 + trial_iterator;
						figure(fig_id)
						set(gcf, 'Position', get(0, 'Screensize'))
						subplot(5,5,regime_iterator)
						scatter(x,y,3,color_gray, 'filled')
						hold on
						scatter(x(regime_indices), y(regime_indices) ,3, time_color(regime_indices,:) ,'filled');
						title(sprintf('%d - %d sec', time_start, time_end));
						pbaspect([1 1 1]) 


						if regime_iterator == total_time_regime
							suptitle(sprintf('%s %s Time overlay on state space', animal_id, char(trial_folders(trial_iterator)) ));
							filename = sprintf('%s_%s_time_overlay_regimes.png', animal_id, char(trial_folders(trial_iterator)) );
							fullname = fullfile(quantify_movements_folder,filename);
							saveas(fig_id,fullname);
						end 


						regime_ll = [];
						regime_ul = [];
						regime_indices = [];

					end % End for regime_iterator
					close(fig_id)


					fig_id = 221 + trial_iterator;

					% Plot each band median power in regime
					for band_iterator = 1:total_bands
						% smooth median power
						regime_median_power(:,band_iterator) = zscore(smooth_gaussian(regime_median_power(:,band_iterator), 1, 2));

						figure(fig_id)
						subplot(3,2, band_iterator)
						plot(temp_time_axis, regime_median_power(:,band_iterator), 'lineWidth', 2);
						title(char(bands_strings(band_iterator)))
						grid on



					end % End for band_iterator

					suptitle(sprintf('%s %s Median power in time regimes', animal_id,  char(trial_folders(trial_iterator)) ));

					filename = sprintf('%s_trial_%s_median_power_in_time.png',animal_id, char(trial_folders(trial_iterator)) );
					fullname = fullfile(quantify_movements_folder,filename);
					saveas(gcf,fullname);

					global_regime_median_power{trial_iterator} = regime_median_power;

					time_regime = [];


				end % End for trial_iterator


				output_global_filename = fullfile(quantify_movements_folder,'global_regime_median_power.mat');
				save(output_global_filename,'global_regime_median_power')		
				
				global_regime_median_power = [];
			end % End if quantify_movements_csd		
		
		

			return
	end % End for animal_iterator

toc;
