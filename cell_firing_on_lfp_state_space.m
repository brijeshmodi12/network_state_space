% ===================================================================
	% Overlay cell firing rate on LFP state space
	% Plot trajectory dependent firing of each cell

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
fprintf('Cell Firing on LFP State Space\n')
fprintf('------------------------------------------------\n')

	
	load_precomputed_mean_firing = 01;


	plot_cells = 0;
	eps_cells = [8, 27, 39, 44, 51, 80, 81 ];

	plot_cells_awake_local_colors = 0;

	plot_cells_awake_global_colors = 0;

	plot_cells_sleep_global_colors = 0;

	

	project_signature_space_combined_umap = 0;		

	overlay_cell_params_on_signature_space = 0;

	sparsity_based_grouping = 	1;


	make_eps_images = 01;

	% Oudated

	classify_cell_signatures_combined_clustering = 0;

	compute_average_group_signatures_combined = 0;

	compute_average_group_signatures_combined_global_scale = 0;

	classify_cell_signatures_awake = 0;

	compute_average_group_signatures_awake = 0;

% -------------------------------------------------------------------
% General Constants
% -------------------------------------------------------------------
	% Paths and Filenames
		root_directory = 'C:\Users\Matteo\Desktop\Brijesh\Open Field Dataset\';
		root_directory = 'D:\matteo_datasets\Open Field Datasets\';

		animal_folders_list = ["AH2","AH3","AH4","AH5"];

		lfp_root_state_space = 'network_state_space_bin200';
		lfp_root_raw_data = 'CSD';


		output_foldername_lfp = 'cell_firing_overlay_lfp'; 
		
		state_space_filename = 'state_space_data.mat';

		spikedata_filename = 'spike_times.mat';

		placemaps_data_filename = 'placemaps_data.mat';
		pyr_indices_filename = 'idx_pyramidal.mat';
		int_indices_filename = 'idx_interneuron.mat';

		cell_classification_foldername = 'cell_classification_information';
		placemaps_data_foldername = 'place_maps_speed_threshold_0';



		% Output files / Pre computed data file names
		binned_spiketrains_filename = 'binned_spiketrains.mat';
		mean_firing_bins_all_cells_filename = 'mean_firing_bins_all_cells_trials.mat';

		cell_indices_sorted_filename = 'lfp_cell_signature_sorted_indices.mat';
		umap_data_filename = 'lfp_umap_data_all_trials.mat';
		umap_data_cellsign_combined_filename = 'lfp_umap_data_cell_signatures_combined.mat';
		lfp_signature_combined_group_filename = 'lfp_signature_combined_group.mat';
		lfp_signature_combined_group_kmeans_filename = 'lfp_signature_combined_group_kmeans.mat';
		lfp_signature_placecells_combined_group_filename = 'lfp_signature_placecells_combined_group.mat';

		group_cell_data_filename = 'group_cell_data.mat';
		sparsity_group_cell_data_filename = 'sparsity_group_cell_data.mat';

		lfp_signature_awake_group_filename = 'lfp_signature_awake_group.mat';
		awake_group_cell_data_filename = 'awake_group_cell_data.mat';


		total_folders = length(animal_folders_list);

		color_gray = [134, 138, 145 ] ./ 255;

		% Combine all the bands in an single array so you can apply loops 
		bands_strings = ["Delta", "Spindle", "Theta", "Slow gamma", "Medium gamma", "Fast gamma"];

		bin_size_sec = 0.2;

		lfp_filename = 'C_Raw_LFP.mat';


		downsampling_frequency = 1000;


		% in number of lfp samples
		bin_size_lfp = bin_size_sec * downsampling_frequency;

		% Set it to high number to include all cells
		firing_threshold = 1000; % In Hz


		% Smoothing binned spiketrains
		smoothing_std = 3;
		smoothing_width = smoothing_std * 6;

		smoothing_paramters.std = smoothing_std;
		smoothing_paramters.width = smoothing_width;

		% For binning on state space ; to compute firing rate
		no_xbins = 15;
		no_ybins = 15;

		% For umap
		n_components = 2;

		awake_trial_list = [2,3];
		sleep_trial_list = [1,4];


% -------------------------------------------------------------------
% Core Loop
% ------------------------------------------------------------------
	for animal_iterator = 1:total_folders

		animal_id = char(animal_folders_list(animal_iterator))
		full_filepath = fullfile(root_directory, animal_id, lfp_root_state_space, state_space_filename);

		% Load state space data
		try
			load(full_filepath);
		catch
			fprintf('State space data not found\n');
			return
		end


		% -------------------------------------------------------------------
		% Animal specific constants
		% -------------------------------------------------------------------	
			total_trials = length(global_states_count_lfp) - 1;

			switch animal_id
				case 'AH2'
					trial_folders = ["preSleep", "1Rectangle", "2Circle", "postSleep"]; % AH2 
					output_folder_string = 'network_state_space_bin200';

					% trial_folders = [ "1Rectangle", "2Circle"]; % AH2 
					% output_folder_string = 'all_trials_power_manifold_awake';
					% only_awake = 1;

					% trial_folders = ["preSleep", "postSleep"]; % AH2 
					% output_folder_string = 	'all_trials_power_manifold_sleep';
					% only_sleep = 1;

					


				case 'AH3'
					trial_folders = ["preSleep", "1Square", "2Circle", "postSleep"]; % AH3 
					output_folder_string = 	'network_state_space_bin200';
					
					% trial_folders = [ "1Square", "2Circle" ]; % AH3
					% output_folder_string = 	'all_trials_power_manifold_awake';
					% only_awake = 1;

					% trial_folders = ["preSleep", "postSleep"]; % AH3  
					% output_folder_string = 	'all_trials_power_manifold_sleep';
					% only_sleep = 1;


				case 'AH4'	
					trial_folders = ["preSleep", "1Circle", "2Square", "postSleep"]; % AH4 
					output_folder_string = 'network_state_space_bin200';

					% trial_folders = [ "1Circle", "2Square" ]; % AH4 
					% output_folder_string = 'all_trials_power_manifold_awake';
					% only_awake = 1;

					% trial_folders = ["preSleep", "postSleep"]; % AH4 
					% output_folder_string = 	'all_trials_power_manifold_sleep';
					% only_sleep = 1;

				case 'AH5'	
					trial_folders = ["preSleep", "1Circle", "2Square", "postSleep"]; % AH5 
					output_folder_string = 'network_state_space_bin200';

					% trial_folders = ["1Circle", "2Square"]; % AH5 
					% output_folder_string = 'all_trials_power_manifold_awake';
					% only_awake = 1;

					% trial_folders = ["preSleep", "postSleep"]; % AH5
					% output_folder_string = 	'all_trials_power_manifold_sleep';
					% only_sleep = 1;
					

				otherwise
					fprintf('Incorrect animal id\n')
					return
			end % End switch animal_id

			lfp_output_folder = fullfile(root_directory,animal_id, lfp_root_state_space, output_foldername_lfp);
			mkdir(lfp_output_folder);

			cell_signature_output_folder = fullfile(root_directory,animal_id,lfp_root_state_space,'cell_signatures_lfp');
			mkdir(cell_signature_output_folder)

			cell_classification_folderpath = fullfile(root_directory,animal_id,cell_classification_foldername);
			

			% global_trial_duration = [];
		

		% -------------------------------------------------------------------
		% Load Precomputed binned spiketrains
		% -------------------------------------------------------------------	
			% binned_spiketrains_filepath = fullfile(root_directory, animal_id, lfp_root_state_space, binned_spiketrains_filename);

			% if load_precomputed_binned_spikes 
			% 	try
			% 		load(binned_spiketrains_filepath)
			% 		fprintf('Binned Spiketrains Found\n')
			% 	catch
			% 		fprintf('Binned Spiketrains Missing\n')
			% 		return
			% 	end
			% end % End if load_precomputed_binned_spikes


		% -------------------------------------------------------------------
		% Compute binned spiketrains
		% -------------------------------------------------------------------	
			% if ~load_precomputed_binned_spikes

			% 	% -------------------------------------------------------------------
			% 	% Process each trial for binning spike times 
			% 	% -------------------------------------------------------------------	
			% 	for trial_iterator = 1:total_trials

			% 		current_trial = char(trial_folders(trial_iterator));

			% 		% -------------------------------------------------------------------
			% 		% Load Trial specific lfp files and Spike data
			% 		% -------------------------------------------------------------------
			% 			try
			% 				lfp_filepath = fullfile(root_directory, animal_id, lfp_root_raw_data, current_trial, lfp_filename);
							
			% 				load(lfp_filepath)

			% 			catch
			% 				fprintf('Raw LFP not found\n')
			% 				return
			% 			end

			% 			try
			% 				spikedata_filepath = fullfile(root_directory, animal_id, current_trial, spikedata_filename);
							
			% 				load(spikedata_filepath)

			% 			catch
			% 				fprintf('Spike data not found\n')
			% 				return
			% 			end


			% 		% -------------------------------------------------------------------
			% 		% Trial Specific Constants
			% 		% -------------------------------------------------------------------
			% 			total_lfp_samples = size(Raw_LFP,2);

			% 			clear Raw_LFP

			% 			total_cells = length(spike_times);

			% 			% total bins for spike train binning
			% 			total_bins = ceil(total_lfp_samples / bin_size_lfp);

			% 			% Edges for binning spikes
			% 			bin_edges = linspace(1 , total_lfp_samples, total_bins);

			% 			bin_start = bin_edges(1:end-1);
						
			% 			bin_end = bin_edges(2:end);

			% 			bin_edges = [bin_start' bin_end'];

			% 			% in seconds
			% 			trial_duration = total_lfp_samples / downsampling_frequency;
		
			% 			global_trial_duration(trial_iterator) = trial_duration;

						
			% 		% -------------------------------------------------------------------
			% 		% Pre process spike timings
			% 		% -------------------------------------------------------------------
			
			% 			cells_list = spike_times;

			% 			% Adjusting occurence of spike time on lfp samples

			% 			for cell_iterator = 1:total_cells
			% 				cells_list{cell_iterator} = floor(cells_list{cell_iterator} * downsampling_frequency);
			% 			end


			% 		% -------------------------------------------------------------------
			% 		% Bin spike times for each cell
			% 		% -------------------------------------------------------------------
			% 			binned_spiketimes_cells = [];



			% 			binned_spiketimes_cells = create_pv_spiketrain(cells_list, bin_edges, trial_duration, firing_threshold, smoothing_paramters);

			% 			global_binned_spikes_lfp = [global_binned_spikes_lfp ; binned_spiketimes_cells];
					

						

			% 	end % End for trial_iterator

			% 	% Save computed binned spiketrains
			% 	binned_spiketrains_filepath = fullfile(root_directory, animal_id, lfp_root_state_space, binned_spiketrains_filename);
			% 	save(binned_spiketrains_filepath, 'global_binned_spikes_lfp', 'global_states_count_lfp', 'global_trial_duration')		
			% 	fprintf('Spiketrains binnned and saved.\n')
				
			% 	return
			% end % End if ~load_precomputed_binned_spikes
	
			

		% -------------------------------------------------------------------
		% Compute mean firing in bin across trials - lfp
		% -------------------------------------------------------------------	
			total_cells = size(global_binned_spikes_lfp,2);
			total_pv = size(global_binned_spikes_lfp,1);

			mean_firing_bins_all_cells = [];

			cell_signatures_across_trials = {[],[],[],[]};

			% Counts no. of bins in each trials
			global_trial_bin_count = [0];

			% Store indices of states in each bins for all trials
			global_trial_bsi = {};

			
			if ~load_precomputed_mean_firing	

				for cell_iterator = 1:total_cells 

					fprintf('Processing Cell : %d\n', cell_iterator);
					
					current_cell = [];
					current_cell = global_binned_spikes_lfp(:,cell_iterator);
					current_cell = smooth_gaussian(current_cell, smoothing_std, smoothing_width);

					% for storing mean firing rates in each bins across trials
					global_mean_firing_bins = [];

					
					for trial_iterator = 1:total_trials 

						current_trial = char(trial_folders(trial_iterator));
						
						pv_range = [];
						x = [];  y = [];

						cell_trial_signature = [];

						range_start = global_states_count_lfp(trial_iterator) + 1 ;
						range_end = global_states_count_lfp(trial_iterator + 1) ;
						pv_range = range_start : range_end; 

						% Trial specific binned spikes
						spiketimes_binned = current_cell(pv_range);

						x = projection_lfp(pv_range,1);
						y = projection_lfp(pv_range,2);
										
						current_trial_duration = length(pv_range) * bin_size_sec;
					
						
						total_spikes = sum(spiketimes_binned);
					
						mean_firing_rate = total_spikes / current_trial_duration;

						mean_firing_matrix(cell_iterator,trial_iterator) = mean_firing_rate;


						% -------------------------------------------------------------------
						% Compute mean firing and cell signatures in each bin
						% -------------------------------------------------------------------	
							% X and Y axis reference points for umap
							x_edges = linspace(min(projection_lfp(:,1)), max(projection_lfp(:,1)), no_xbins);
							y_edges = linspace(min(projection_lfp(:,2)), max(projection_lfp(:,2)), no_ybins);


							% Get xindices
							[xcounts, x_edges, xindices] = histcounts(x, 'BinEdges', x_edges);
							[ycounts, y_edges, yindices] = histcounts(y, 'BinEdges', y_edges);

							total_bins = no_xbins * no_ybins;

							mean_firing_in_bins = [];			

							bin_counter = 1;
														
							for xbin_iterator = 1:no_xbins

								tempx =	find(xindices == xbin_iterator);

								for ybin_iterator = 1:no_ybins

									bin_id = (xbin_iterator - 1) * no_xbins + ybin_iterator;

									tempy =	find(yindices == ybin_iterator);

									bin_states_indices = intersect(tempx,tempy);
									bsi = bin_states_indices;								

									if length(bsi) < 5 
										mean_firing_temp = NaN;
										mean_firing_in_bins(bsi) = mean_firing_temp;
										global_trial_bsi{trial_iterator,bin_counter} = bsi;
										cell_trial_signature(bin_counter) = mean_firing_temp;

										bin_counter = bin_counter + 1;

										continue;
									end

									total_spikes_in_bins = sum(spiketimes_binned(bsi));
									time_in_bins = length(bsi) * bin_size_sec;

									mean_firing_temp = total_spikes_in_bins / time_in_bins ;												


									% in Hz
									mean_firing_in_bins(bsi) = mean_firing_temp;

									cell_trial_signature(bin_counter) = mean_firing_temp;

									if cell_iterator == 1;
										global_trial_bsi{trial_iterator,bin_counter} = bsi;
									end

									bin_counter = bin_counter + 1;
								
								end % End for ybin_iterator

							end % End for xbin_iterator	

							% Save cell signature into global collection
							cell_sig_temp = cell_signatures_across_trials{1,trial_iterator};

							% Append current cell signature to list
							cell_sig_temp = [cell_sig_temp  ; cell_trial_signature];

							% Update global collection
							cell_signatures_across_trials{1,trial_iterator} = cell_sig_temp;

							% Save mean firing in bins for each trial
							global_mean_firing_bins = [global_mean_firing_bins  mean_firing_in_bins];

							if cell_iterator == 1 
								global_trial_bin_count = [global_trial_bin_count bin_counter-1];
							end
								
					end % End for trial_iterator		

					mean_firing_bins_all_cells(:, cell_iterator) = reshape(global_mean_firing_bins, [], 1);

				end % End for cell_iterator	
				
				global_trial_bin_count = cumsum(global_trial_bin_count);
				fullname = fullfile(lfp_output_folder, mean_firing_bins_all_cells_filename);
				save(fullname, 'mean_firing_bins_all_cells', 'global_states_count_lfp', 'mean_firing_matrix','cell_signatures_across_trials','global_trial_bin_count','global_trial_bsi');

			end % End if load_precomputed_mean_firing


		% -------------------------------------------------------------------
		% Plot each cells 
		% -------------------------------------------------------------------	
			if plot_cells
				try 
					load(fullfile(lfp_output_folder, mean_firing_bins_all_cells_filename));
				catch
					fprintf('Mean Firing in Bins All Cells Not Found\n')
					return
				end

				total_cells =  size(mean_firing_bins_all_cells,2);
				total_pv = size(mean_firing_bins_all_cells,1);

				for cell_iterator = 1:total_cells
					fprintf('Plotting cell %d...\n',cell_iterator)

					cell_firing_across_trials = [];
					cell_firing_across_trials = mean_firing_bins_all_cells(: ,cell_iterator);

					cell_colors = []; sort_ind = [];

					% Assign global colors to cell firing across trials
					% [s, sort_ind] = sort(cell_firing_across_trials);
					% cell_colors(sort_ind,:) = plasma(length(cell_firing_across_trials));

					cell_firing_range = [min(cell_firing_across_trials), prctile(cell_firing_across_trials,95)];


					[cell_cmap cell_colors] = assign_colors(cell_firing_across_trials, cell_firing_range, 100, 'plasma');


					nan_index = find(isnan(cell_firing_across_trials));

					if nan_index
						cell_colors(nan_index,:) = repmat(color_gray, length(nan_index),1);
					end

					all_pv = 1:global_states_count_lfp(5);
					

					for trial_iterator = 1:total_trials	

						pv_range = []; x = []; y = [];	trial_colors = [];

						current_trial = char(trial_folders(trial_iterator));

						range_start = global_states_count_lfp(trial_iterator) + 1 ;
						range_end = global_states_count_lfp(trial_iterator + 1) ;

						pv_range = range_start : range_end; 
						non_pv_range = setdiff(all_pv, pv_range);

						x = projection_lfp(pv_range,1);
						y = projection_lfp(pv_range,2);
						trial_colors = cell_colors(pv_range,:);

						figure(185)
						set(gcf, 'Position', get(0, 'Screensize'));
						subplot(1,4,trial_iterator)
						scatter(projection_lfp(non_pv_range,1), projection_lfp(non_pv_range,2), 5, color_gray, 'filled')
						hold on 
						scatter(x,y,5, trial_colors ,'filled')
						title(sprintf('%s (%0.1f Hz)' ,  current_trial, mean_firing_matrix(cell_iterator,trial_iterator) ) );
						colormap(cell_cmap)
						caxis(cell_firing_range)
						h = colorbar;
						ylabel(h, 'Frequency (Hz)')
						h.FontSize = 12;
						pbaspect([1 1 1]) 
						% return

					end % End for trial_iterator

					suptitle(sprintf('%s Cell firing on LFP ; Cell %d ',animal_id, cell_iterator))

					cell_filename = sprintf('%s_lfp_cell_%d.png', animal_id, cell_iterator);
					cell_image_filepath = fullfile(lfp_output_folder, cell_filename);
					saveas(gcf, cell_image_filepath);	


					if make_eps_images & ismember(cell_iterator , eps_cells)
									
						epsname = sprintf('%s_lfp_cell_%d.eps', animal_id, cell_iterator);
						fullname_eps = fullfile(lfp_output_folder,epsname);
						set(gcf,'renderer','Painters')
						saveas(gcf,fullname_eps,'epsc')
						

					end % End if make_eps_images	

				
				end % End for cell_iterator	

			end % End plot_cells


	

		% -------------------------------------------------------------------
		% Classify cell signatures combined clustering
		% -------------------------------------------------------------------
			if classify_cell_signatures_combined_clustering
				fprintf('Clustering Cell Signatures Combined\n')
				try 
					load(fullfile(lfp_output_folder, mean_firing_bins_all_cells_filename));
					fprintf('Cell Signatures Found\n')
				catch
					fprintf('Cell Signatures Not Found\n')
					return
				end

				total_cells = size(cell_signatures_across_trials{1,1}, 1);

				global_labels_indices = [];

				x_ticks_array = 1:total_cells;

				lfp_cell_signatures_combined = [];
				lfp_cell_signatures_combined_og = [];

				% for trial_iterator = 1:total_trials
				for trial_iterator = 1:total_trials	

					lfp_cell_signatures_combined_og = [lfp_cell_signatures_combined_og cell_signatures_across_trials{1, trial_iterator} ] ;
					if ismember(trial_iterator, [2,3])
						continue;
					end
		
					lfp_cell_signatures_combined = [lfp_cell_signatures_combined cell_signatures_across_trials{1, trial_iterator} ] ;

				end % End for trial_iterator



				nan_columns = find(isnan(sum(lfp_cell_signatures_combined )) );
				lfp_cell_signatures_combined(:,nan_columns) = [];

			

				
				% % % across cells, normalizes each column
				% lfp_cell_signatures_combined = zscore(lfp_cell_signatures_combined, 0, 1);

				% across trials, normalizes each row
				lfp_cell_signatures_combined = zscore(lfp_cell_signatures_combined, 0, 2);

				% % across cells, normalizes each column
				% lfp_cell_signatures_combined = zscore(lfp_cell_signatures_combined, 0, 1);
				

				% pairwise_distance = pdist(lfp_cell_signatures_combined,'cosine');
			
				% % Just put the results in matrix form and add some small noise to avoid coincident points
				% pairwise_distance = squareform(pairwise_distance + rand(size(pairwise_distance))*0.0001 );

				

				correlation_matrix = corrcoef(lfp_cell_signatures_combined');

				core_mat = correlation_matrix;

				%Re-order your matrices 
				% no_of_groups = Number of groups to sort your matrix into
				no_of_groups = 6;
				
				%# Remove diagonal elements
				corrMat = core_mat; %MATRICES WILL BE RE-ORDERED ON THE BASE OF THIS ONE 
				corrMat = corrMat - eye(size(corrMat));

				%# and convert to a vector (as pdist)
				dissimilarity = 1 - corrMat';

				%# decide on a cutoff
				%# remember that 0.4 corresponds to corr of 0.6!
				% cutoff = 0.1; 
				%# perform complete linkage clustering
				Z = linkage(dissimilarity,'complete');

				%# group the data into clusters
				groups = cluster(Z,'maxclust',no_of_groups);

				[~, sorted_indices] = sort(groups);
				
				labels_sorted = x_ticks_array(sorted_indices);

				sorted_corr_matrix = core_mat(sorted_indices, sorted_indices);

				figure(654)
				set(gcf, 'Position', get(0, 'Screensize'));
				imagesc(sorted_corr_matrix); 
				ax = gca;
				colormap(brewermap([],'*spectral'));
				pbaspect([1 1 1])
				h = colorbar;
				ylabel(h,'Correlation Coefficient')
				% xticks(1:total_cells)
				% xticklabels(labels_sorted)
				% yticks(1:total_cells)
				% yticklabels(labels_sorted)
				xlabel('Cells')
				ylabel('Cells')
				caxis([-1 1 ])
				title('Sorted Correlation Matrix')
				xtickangle(90)
				
				ax.XAxis.FontSize = 12;
				ax.YAxis.FontSize = 12;

				suptitle(sprintf('%s LFP Cell Signatures Combined Trials (Groups = %d)', animal_id, no_of_groups))
			
				cell_filename = sprintf('%s_lfp_cell_signatures_combined.png', animal_id );
				cell_image_filepath = fullfile(cell_signature_output_folder, cell_filename);
				saveas(gcf, cell_image_filepath);
				
				lfp_cell_signatures_combined = lfp_cell_signatures_combined_og;


				cell_indices_sorted_filepath = fullfile(cell_signature_output_folder, lfp_signature_combined_group_filename);
				save(cell_indices_sorted_filepath,  'groups', 'lfp_cell_signatures_combined')

							
			end % End if classify_cell_signatures_combined_clustering


		



		% -------------------------------------------------------------------
		% Project cell signatures combined with UMAP (Combined)
		% -------------------------------------------------------------------
			if project_signature_space_combined_umap 

				fprintf('Projecting Signature Space using UMAP\n')

				try 
					load(fullfile(lfp_output_folder, mean_firing_bins_all_cells_filename));
					fprintf('Cell Signatures Found\n')
				catch
					fprintf('Cell Signatures Not Found\n')
					return
				end

				total_cells = size(cell_signatures_across_trials{1,1}, 1);

				global_labels_indices = [];

				lfp_cell_signatures_combined_raw = [];
				lfp_cell_signatures_combined_nanfree = [];
				lfp_sleep_signatures_raw = [];

				for trial_iterator = 1:total_trials

					if ismember(trial_iterator, [1 4])
						lfp_sleep_signatures_raw = [lfp_sleep_signatures_raw cell_signatures_across_trials{1, trial_iterator} ];
						
					end

					current_trial = char(trial_folders(trial_iterator));
		
					lfp_cell_signatures_combined_raw = [lfp_cell_signatures_combined_raw cell_signatures_across_trials{1, trial_iterator} ] ;


				end % End for trial_iterator

				nan_columns = find(isnan(sum(lfp_sleep_signatures_raw)) );

				lfp_cell_signatures_combined_nanfree = lfp_sleep_signatures_raw;

				lfp_cell_signatures_combined_nanfree(:,nan_columns) = [];

			
				% across trials, normalizes each row
				% lfp_cell_signatures_combined = zscore(lfp_cell_signatures_combined, 0, 2);

				% mean(lfp_cell_signatures_combined(1,:))

				% across cells, normalizes each column
				lfp_cell_signatures_combined_nanfree = zscore(lfp_cell_signatures_combined_nanfree, 0, 1);

			

									
					% ------------------------
					% Run UMAP 
					% ------------------------
						[umap_output_cell_signatures_combined, umap_params] = run_umap(lfp_cell_signatures_combined_nanfree, 'n_components', n_components, 'n_neighbors', 10 ,'min_dist',0.5, 'metric', 'cosine' );
 

						x = umap_output_cell_signatures_combined(:,1);
						y = umap_output_cell_signatures_combined(:,2);
						

						figure(4)
						set(gcf, 'Position', get(0, 'Screensize'));
						scatter(x,y, 50,  color_gray ,'filled');	
						ylim([-10 10])
						xlim([-5 15])
							

						pbaspect([1 1 1]) 
						title(sprintf('%dD ',size(lfp_cell_signatures_combined_nanfree,2) ))
					
				

				suptitle(sprintf('%s LFP Signatures Combined UMAP (%d cells)', animal_id, total_cells )); 

				
				cell_filename = sprintf('%s_lfp_cell_signatures_combined_UMAP.png', animal_id);
				cell_image_filepath = fullfile(cell_signature_output_folder, cell_filename);
				saveas(gcf, cell_image_filepath);

				if make_eps_images 
									
					epsname = sprintf('%s_lfp_cell_signatures_combined_UMAP.eps', animal_id);
					fullname_eps = fullfile(cell_signature_output_folder,epsname);
					set(gcf,'renderer','Painters')
					saveas(gcf,fullname_eps,'epsc')
					

				end % End if make_eps_images	



				lfp_cell_signatures_combined = lfp_cell_signatures_combined_raw;

				
				umap_data_filepath = fullfile(cell_signature_output_folder, umap_data_cellsign_combined_filename);
				save(umap_data_filepath, 'umap_output_cell_signatures_combined', 'lfp_cell_signatures_combined')

			end % End if project_signature_space_combined_umap

	
		% -------------------------------------------------------------------
		% Average signatures of groups (Combined)
		% -------------------------------------------------------------------	
			if compute_average_group_signatures_combined

				fprintf('Compute Average signatures of groups\n')

				try 
					load(fullfile(cell_signature_output_folder, lfp_signature_combined_group_filename));
					load(fullfile(lfp_output_folder, mean_firing_bins_all_cells_filename));
					fprintf('Cell Signatures, Group Detials and Bin States Info Found\n')
				catch
					fprintf('Cell Signatures, Group Detials and Bin States Info Not Found\n')
					return
				end

				total_groups = length(unique(groups));

				group_average_matrix = [];

				group_list = unique(groups);

				

				group_color = magma(total_groups);

				group_cell_data = {};

				nan_columns = find(isnan(sum(lfp_cell_signatures_combined )) );
				non_nan_columms = find(~isnan(sum(lfp_cell_signatures_combined )) );

				for group_iterator = 1:total_groups

					group_cells_indices = []; group_average = []; sort_ind = [];
					group_bin_colors = []; group_cells = [];

					group_cells_indices = find(groups == group_list(group_iterator));

					group_cell_data{1, group_iterator} = group_cells_indices;

					temp_signatures_without_nan = lfp_cell_signatures_combined(:,non_nan_columms);

					group_cells = zscore(temp_signatures_without_nan(group_cells_indices,:)');
					group_cells = group_cells';

					temp_signatures_with_nan = lfp_cell_signatures_combined(group_cells_indices,:);
					temp_signatures_with_nan(:,non_nan_columms) = group_cells;

					if length(group_cells_indices) > 1
				
						group_average = nanmean(temp_signatures_with_nan);
					else

						group_average = temp_signatures_with_nan;
					end


					group_average_matrix = [group_average_matrix; group_average];

					group_average_without_nan = group_average(find(~isnan(group_average)));

					[s,sort_ind] = sort(group_average_without_nan);

					% Colors for all bins across trials
					group_bin_colors(sort_ind,:) = plasma(length(group_average_without_nan));

					final_group_colors(non_nan_columms,:) = group_bin_colors;

					final_group_colors(nan_columns,:) = repmat(color_gray,length(nan_columns),1);



					% Assign colors to each states in a bin
					for trial_iterator = 1 : total_trials

						current_trial = char(trial_folders(trial_iterator));

						trial_group_colors = []; trial_bin_range = [];
						trial_state_colors = [];

						trial_bin_start = global_trial_bin_count(trial_iterator) + 1;
						trial_bin_end =  global_trial_bin_count(trial_iterator + 1);

						trial_bin_range = trial_bin_start:trial_bin_end;

						trial_group_colors = final_group_colors(trial_bin_range,:);

						total_bins = length(trial_bin_range);

						for bin_iterator = 1 : total_bins

							bin_states_indices = global_trial_bsi{trial_iterator, bin_iterator};
							
							trial_state_colors(bin_states_indices,:) = repmat(trial_group_colors(bin_iterator,:), length(bin_states_indices), 1);

						end % End for bin_iterator

						

						% Plot trial average for this group
						range_start = global_states_count_lfp(trial_iterator) + 1 ;
						range_end = global_states_count_lfp(trial_iterator + 1) ;

						pv_range = range_start : range_end; 

						x = projection_lfp(pv_range,1);
						y = projection_lfp(pv_range,2);
						cmin = min(group_average);
						cmax = max(group_average);

						figure(185)
						set(gcf, 'Position', get(0, 'Screensize'));
						subplot(2,2,trial_iterator)
						scatter(projection_lfp(:,1), projection_lfp(:,2), 5, color_gray, 'filled')
						hold on 
						scatter(x,y,5, trial_state_colors ,'filled')
						title(sprintf('%s' ,  current_trial ));
						colormap('plasma')
						% colormap(brewermap([],'pubugn'))
						caxis( [ cmin cmax ] )
						h = colorbar;
						ylabel(h, 'Norm. Frequency')
						pbaspect([1 1 1]) 


					end % End for trial_iterator

					suptitle(sprintf('%s Group %d LFP Avg. Signature (n = %d cells)',animal_id, group_list(group_iterator), length(group_cells_indices) ));

					cell_filename = sprintf('%s_avg_signature_group_%d.png', animal_id, group_list(group_iterator) );
					cell_image_filepath = fullfile(cell_signature_output_folder, cell_filename);
					saveas(gcf, cell_image_filepath);

					

					figure(198)
					h = plot(group_average,'lineWidth',1);
				    set(h, {'color'}, num2cell(group_colors(group_iterator,:), 2));
				    cbh = colorbar;
				    set(cbh,'YTick',[1:total_groups])
				    ylabel(cbh,'Group ID')
				    colormap('magma')
				    caxis([1 total_groups]);
					hold on;
					xlabel('Dimensions')
					ylabel('Values of Avg Group Signature')
					
					group_average = smooth_gaussian(group_average, 5, 20);

					figure(195)
					h = plot(group_average,'lineWidth',1);
				    set(h, {'color'}, num2cell(group_colors(group_iterator,:), 2));
				    cbh = colorbar;
				    set(cbh,'YTick',[1:total_groups])
				    ylabel(cbh,'Group ID')
				    colormap('magma')
				    caxis([1 total_groups]);
					hold on;
					xlabel('Dimensions')
					ylabel('Values of Avg Group Signature')


					% return

			
				end % End for group iterator


				suptitle(sprintf('%s LFP Signature in High Dimensional Space (%d Groups)  ',animal_id , total_groups));

				cell_filename = sprintf('%s_avg_signature_group_high_dim.png', animal_id  );
				cell_image_filepath = fullfile(cell_signature_output_folder, cell_filename);
				saveas(figure(195), cell_image_filepath);

				group_cell_data_filepath = fullfile(cell_signature_output_folder, group_cell_data_filename);
				save(group_cell_data_filepath, 'group_cell_data')
				
			end % End if compute_average_group_signatures_combined
			

		% -------------------------------------------------------------------
		% Average signatures of groups (Combined)
		% -------------------------------------------------------------------	
			if compute_average_group_signatures_combined_global_scale	

				fprintf('Compute Average signatures of groups\n')

				try 
					load(fullfile(cell_signature_output_folder, lfp_signature_combined_group_filename));
					load(fullfile(lfp_output_folder, mean_firing_bins_all_cells_filename));
					fprintf('Cell Signatures, Group Detials and Bin States Info Found\n')
				catch
					fprintf('Cell Signatures, Group Detials and Bin States Info Not Found\n')
					return
				end

				total_groups = length(unique(groups));

				group_average_matrix = [];

				group_list = unique(groups);

				group_colors = magma(total_groups);

				group_cell_data = {};

				nan_columns = find(isnan(sum(lfp_cell_signatures_combined )) );
				non_nan_columms = find(~isnan(sum(lfp_cell_signatures_combined )) );

				for group_iterator = 1:total_groups

					group_cells_indices = []; group_average = []; sort_ind = [];
					group_bin_colors = []; group_cells = [];

					group_cells_indices = find(groups == group_list(group_iterator));

					group_cell_data{1, group_iterator} = group_cells_indices;

					temp_signatures_without_nan = lfp_cell_signatures_combined(:,non_nan_columms);

					group_cells = zscore(temp_signatures_without_nan(group_cells_indices,:)');
					group_cells = group_cells';

					temp_signatures_with_nan = lfp_cell_signatures_combined(group_cells_indices,:);
					temp_signatures_with_nan(:,non_nan_columms) = group_cells;

					if length(group_cells_indices) > 1
						group_average = nanmean(temp_signatures_with_nan);
					else
						group_average = temp_signatures_with_nan;
					end


					group_average_matrix = [group_average_matrix; group_average];

					group_average_without_nan = group_average(find(~isnan(group_average)));

					[s,sort_ind] = sort(group_average_without_nan);

					% Colors for all bins across trials
					group_bin_colors(sort_ind,:) = plasma(length(group_average_without_nan));

					final_group_colors(non_nan_columms,:) = group_bin_colors;

					final_group_colors(nan_columns,:) = repmat(color_gray,length(nan_columns),1);



					% Assign colors to each states in a bin
					for trial_iterator = 1 : total_trials

						current_trial = char(trial_folders(trial_iterator));

						trial_group_colors = []; trial_bin_range = [];
						trial_state_colors = [];

						trial_bin_start = global_trial_bin_count(trial_iterator) + 1;
						trial_bin_end =  global_trial_bin_count(trial_iterator + 1);

						trial_bin_range = trial_bin_start:trial_bin_end;

						trial_group_colors = final_group_colors(trial_bin_range,:);

						total_bins = length(trial_bin_range);

						for bin_iterator = 1 : total_bins

							bin_states_indices = global_trial_bsi{trial_iterator, bin_iterator};
							
							trial_state_colors(bin_states_indices,:) = repmat(trial_group_colors(bin_iterator,:), length(bin_states_indices), 1);

						end % End for bin_iterator

						

						% Plot trial average for this group
						range_start = global_states_count_lfp(trial_iterator) + 1 ;
						range_end = global_states_count_lfp(trial_iterator + 1) ;

						pv_range = range_start : range_end; 

						x = projection_lfp(pv_range,1);
						y = projection_lfp(pv_range,2);
						cmin = min(group_average);
						cmax = max(group_average);

						figure(185)
						set(gcf, 'Position', get(0, 'Screensize'));
						subplot(2,2,trial_iterator)
						scatter(projection_lfp(:,1), projection_lfp(:,2), 5, color_gray, 'filled')
						hold on 
						scatter(x,y,5, trial_state_colors ,'filled')
						title(sprintf('%s' ,  current_trial ));
						colormap('plasma')
						% colormap(brewermap([],'pubugn'))
						caxis( [ cmin cmax ] )
						h = colorbar;
						ylabel(h, 'Norm. Frequency')
						pbaspect([1 1 1]) 


					end % End for trial_iterator

					suptitle(sprintf('%s Group %d LFP Avg. Signature (n = %d cells)',animal_id, group_list(group_iterator), length(group_cells_indices) ));

					cell_filename = sprintf('%s_avg_signature_group_%d.png', animal_id, group_list(group_iterator) );
					cell_image_filepath = fullfile(cell_signature_output_folder, cell_filename);
					saveas(gcf, cell_image_filepath);

					

					figure(198)
					h = plot(group_average,'lineWidth',1);
				    set(h, {'color'}, num2cell(group_colors(group_iterator,:), 2));
				    cbh = colorbar;
				    set(cbh,'YTick',[1:total_groups])
				    ylabel(cbh,'Group ID')
				    colormap('magma')
				    caxis([1 total_groups]);
					hold on;
					xlabel('Dimensions')
					ylabel('Values of Avg Group Signature')
					
					group_average = smooth_gaussian(group_average, 5, 20);

					figure(195)
					h = plot(group_average,'lineWidth',1);
				    set(h, {'color'}, num2cell(group_colors(group_iterator,:), 2));
				    cbh = colorbar;
				    set(cbh,'YTick',[1:total_groups])
				    ylabel(cbh,'Group ID')
				    colormap('magma')
				    caxis([1 total_groups]);
					hold on;
					xlabel('Dimensions')
					ylabel('Values of Avg Group Signature')


					% return

			
				end % End for group iterator


				suptitle(sprintf('%s LFP Signature in High Dimensional Space (%d Groups)  ',animal_id , total_groups));

				cell_filename = sprintf('%s_avg_signature_group_high_dim.png', animal_id  );
				cell_image_filepath = fullfile(cell_signature_output_folder, cell_filename);
				saveas(figure(195), cell_image_filepath);

				group_cell_data_filepath = fullfile(cell_signature_output_folder, group_cell_data_filename);
				save(group_cell_data_filepath, 'group_cell_data')
				
			end % End if compute_average_group_signatures_combined_global_scale


		% -------------------------------------------------------------------
		% Classify cell signatures awake 
		% -------------------------------------------------------------------
			if classify_cell_signatures_awake
				fprintf('Classifying Cells in Awake Trials\n')

				try 
					load(fullfile(lfp_output_folder, mean_firing_bins_all_cells_filename));
					fprintf('Cell Signatures Found\n')
				catch
					fprintf('Cell Signatures Not Found\n')
					return
				end

				total_cells = size(cell_signatures_across_trials{1,1}, 1);

				global_labels_indices = [];

				x_ticks_array = 1:total_cells;

				lfp_cell_signatures_awake = [];

				% for trial_iterator = 1:total_trials
				for trial_iterator = awake_trial_list	

					current_trial = char(trial_folders(trial_iterator));
		
					lfp_cell_signatures_awake = [lfp_cell_signatures_awake cell_signatures_across_trials{1, trial_iterator} ] ;

				end % End for trial_iterator

				lfp_cell_signatures_awake = zscore(lfp_cell_signatures_awake, 0, 2);

				lfp_cell_signatures_awake = zscore(lfp_cell_signatures_awake, 0, 1);

				% pairwise_distance = pdist(lfp_cell_signatures_awake,'cosine');
			
				% % Just put the results in matrix form and add some small noise to avoid coincident points
				% pairwise_distance = squareform(pairwise_distance + rand(size(pairwise_distance))*0.0001 );

				correlation_matrix = corrcoef(lfp_cell_signatures_awake');

				core_mat = correlation_matrix;

				%Re-order your matrices 
				no_of_groups = 5;
				
				%# Remove diagonal elements
				corrMat = core_mat; %MATRICES WILL BE RE-ORDERED ON THE BASE OF THIS ONE 
				corrMat = corrMat - eye(size(corrMat));

				%# and convert to a vector (as pdist)
				dissimilarity = 1 - corrMat';

				%# decide on a cutoff
				%# remember that 0.4 corresponds to corr of 0.6!
				% cutoff = 0.1; 
				%# perform complete linkage clustering
				Z = linkage(dissimilarity,'complete');

				%# group the data into clusters
				groups = cluster(Z,'maxclust',no_of_groups);

				[~, sorted_indices] = sort(groups);
				
				labels_sorted = x_ticks_array(sorted_indices);

				sorted_corr_matrix = correlation_matrix(sorted_indices, sorted_indices);

				figure(654)
				set(gcf, 'Position', get(0, 'Screensize'));
				imagesc(sorted_corr_matrix); 
				ax = gca;
				colormap(brewermap([],'*PuOr'));
				pbaspect([1 1 1])
				h = colorbar;
				ylabel(h,'Correlation Coefficient')
				% xticks(1:total_cells)
				% xticklabels(labels_sorted)
				% yticks(1:total_cells)
				% yticklabels(labels_sorted)
				xlabel('Cells')
				ylabel('Cells')
				caxis([-1 1 ])
				title('Sorted Correlation Matrix Awake')
				xtickangle(90)
				
				ax.XAxis.FontSize = 12;
				ax.YAxis.FontSize = 12;

				suptitle(sprintf('%s LFP Cell Signatures Awake Trials (Groups = %d)', animal_id, no_of_groups))
			
				cell_filename = sprintf('%s_lfp_cell_signatures_awake.png', animal_id );
				cell_image_filepath = fullfile(cell_signature_output_folder, cell_filename);
				saveas(gcf, cell_image_filepath);
				
				cell_indices_sorted_filepath = fullfile(cell_signature_output_folder, lfp_signature_awake_group_filename);
				save(cell_indices_sorted_filepath, 'global_labels_indices', 'groups', 'lfp_cell_signatures_awake')

							
			end % End if classify_cell_signatures_awake



		% -------------------------------------------------------------------
		% Average signatures of groups (Awake)
		% -------------------------------------------------------------------	
			if compute_average_group_signatures_awake

				fprintf('Compute Average signatures of groups -- Awake Trails\n')

				try 
					load(fullfile(cell_signature_output_folder, lfp_signature_awake_group_filename));
					load(fullfile(lfp_output_folder, mean_firing_bins_all_cells_filename));
					fprintf('Cell Signatures, Group Detials and Bin States Info Found\n')
				catch
					fprintf('Cell Signatures, Group Detials and Bin States Info Not Found\n')
					return
				end

				total_groups = length(unique(groups));

				group_average_matrix = [];

				group_list = unique(groups);

				group_colors = magma(total_groups);

				group_cell_data = {};

				for group_iterator = 1:total_groups

					group_cells_indices = []; group_average = []; sort_ind = [];
					group_bin_colors = [];

					group_cells_indices = find(groups == group_list(group_iterator));

					awake_group_cell_data{1, group_iterator} = group_cells_indices;

					if length(group_cells_indices) > 1
						group_average = mean(lfp_cell_signatures_awake(group_cells_indices,:));
					else
						group_average = lfp_cell_signatures_awake(group_cells_indices,:);
					end



					group_average_matrix = [group_average_matrix; group_average];

					[s,sort_ind] = sort(group_average);

					% Colors for all bins across trials
					group_bin_colors(sort_ind,:) = plasma(length(group_average));

					% Assign colors to each states in a bin

					for trial_iterator = awake_trial_list

						current_trial = char(trial_folders(trial_iterator));

						trial_group_colors = []; trial_bin_range = [];
						trial_state_colors = [];

						trial_bin_start = global_trial_bin_count(trial_iterator) + 1 - global_trial_bin_count(2);
						trial_bin_end =  global_trial_bin_count(trial_iterator + 1) - global_trial_bin_count(2);

						trial_bin_range = trial_bin_start:trial_bin_end;

						
						trial_group_colors = group_bin_colors(trial_bin_range,:);

						total_bins = length(trial_bin_range);

						for bin_iterator = 1 : total_bins

							bin_states_indices = global_trial_bsi{trial_iterator, bin_iterator};

							trial_state_colors(bin_states_indices,:) = repmat(trial_group_colors(bin_iterator,:), length(bin_states_indices), 1);

						end % End for bin_iterator

						% Plot trial average for this group
						range_start = global_states_count_lfp(trial_iterator) + 1 ;
						range_end = global_states_count_lfp(trial_iterator + 1) ;

						pv_range = range_start : range_end; 

						x = projection_lfp(pv_range,1);
						y = projection_lfp(pv_range,2);
						cmin = min(group_average);
						cmax = max(group_average);

						figure(185)
						set(gcf, 'Position', get(0, 'Screensize'));
						subplot(1,2,trial_iterator-1)
						scatter(projection_lfp(:,1), projection_lfp(:,2), 10, color_gray, 'filled')
						hold on 
						scatter(x,y,10, trial_state_colors ,'filled')
						title(sprintf('%s' ,  current_trial ));
						colormap('plasma')
						% colormap(brewermap([],'pubugn'))
						caxis( [ cmin cmax ] )
						h = colorbar;
						ylabel(h, 'Norm. Frequency')
						pbaspect([1 1 1]) 


					end % End for trial_iterator

					suptitle(sprintf('%s Awake Group %d LFP Avg. Signature (n = %d cells)',animal_id, group_list(group_iterator), length(group_cells_indices) ));

					cell_filename = sprintf('%s_awake_avg_signature_group_%d.png', animal_id, group_list(group_iterator) );
					cell_image_filepath = fullfile(cell_signature_output_folder, cell_filename);
					saveas(gcf, cell_image_filepath);

					
					group_average = smooth_gaussian(group_average,5,20);

					figure(195)
					h = plot(group_average,'lineWidth',1);
				    set(h, {'color'}, num2cell(group_colors(group_iterator,:), 2));
				    cbh = colorbar;
				    set(cbh,'YTick',[1:total_groups])
				    ylabel(cbh,'Group ID')
				    colormap('magma')
				    caxis([1 total_groups]);
					hold on;
					xlabel('Dimensions')
					ylabel('Values of Avg Group Signature')



			
				end % End for group iterator


				suptitle(sprintf('%s Awake LFP Signature in High Dimensional Space (%d Groups)  ',animal_id , total_groups));

				cell_filename = sprintf('%s_awake_avg_signature_group_high_dim.png', animal_id  );
				cell_image_filepath = fullfile(cell_signature_output_folder, cell_filename);
				saveas(figure(195), cell_image_filepath);

				group_cell_data_filepath = fullfile(cell_signature_output_folder, awake_group_cell_data_filename);
				save(group_cell_data_filepath, 'awake_group_cell_data')


				
			end % End if compute_average_group_signatures_awake




		% -------------------------------------------------------------------
		% Overlay cell parameters on signature space
		% -------------------------------------------------------------------
			if overlay_cell_params_on_signature_space

				% Import parameters (pyramidal or interneuron, sparsity, si, mean firing)

				fprintf('\nOverlay Cell parameters on signature space\n')

				try 
					
					load(fullfile(cell_signature_output_folder, umap_data_cellsign_combined_filename));
					fprintf('UMAP Projection of Signature Space Found\n')
					% load(fullfile(root_directory, animal_id, cell_classification_foldername, pyr_indices_filename));
					% load(fullfile(root_directory, animal_id, cell_classification_foldername, int_indices_filename));
					% fprintf('Pyramidal, Interneuron Indices Found\n')
					load(fullfile(root_directory, animal_id, placemaps_data_foldername, placemaps_data_filename ));
					fprintf('Placemaps Data Found\n')
					load(fullfile(lfp_output_folder, mean_firing_bins_all_cells_filename));
					fprintf('Mean Firing Found\n')
					% load(fullfile(cell_signature_output_folder, lfp_signature_combined_group_filename));
					% fprintf('Cell Groups Found\n')

				catch
					fprintf('UMAP Projection of Signature Space or Parameters Not Found\n')
					return
				end
				
				total_cells = size(mean_firing_matrix, 1);
				x = umap_output_cell_signatures_combined(:,1);
				y = umap_output_cell_signatures_combined(:,2);


				% -----------------------------------
				% Overlay cell mean firing
				% ----------------------------------- 	
					total_trials = size(mean_firing_matrix, 2);
					
					% for trial_iterator = 1 : total_trials

						% current_trial = char(trial_folders(trial_iterator));

						current_mean_firing = mean(mean_firing_matrix');

						[cmin, cmax] = bounds(current_mean_firing);

						sort_ind = []; mean_firing_colors = [];
						% Assign Colors
						[firing_cmap, mean_firing_colors] = assign_colors(current_mean_firing, [cmin, cmax], 10, 'viridis' );
						

						

						figure(212)
						% subplot(2,2,trial_iterator)
						scatter(x,y,50,mean_firing_colors,'filled')
						colorbar
						colormap(firing_cmap)
						caxis([cmin, cmax])
						% title(sprintf('%s', current_trial))
						pbaspect([1 1 1]) 
						ylim([-10 10])
						xlim([-5 15])


					% end % End for trial_iterator

					suptitle(sprintf('%s Mean Firing Overlay on LFP Signature Space', animal_id ))

					cell_filename = sprintf('%s_mean_firing_on_lfp_signature_space.png', animal_id);
					cell_image_filepath = fullfile(cell_signature_output_folder, cell_filename);
					saveas(gcf, cell_image_filepath);

					if make_eps_images 
									
						epsname = sprintf('%s_mean_firing_on_lfp_signature_space.eps', animal_id);
						fullname_eps = fullfile(cell_signature_output_folder,epsname);
						set(gcf,'renderer','Painters')
						saveas(gcf,fullname_eps,'epsc')
						

					end % End if make_eps_images	
						

				% -----------------------------------
				% Overlay spatial information
				% -----------------------------------	
					% total_trials = size(si_matrix, 2);
					
					% for trial_iterator = 1 : total_trials

					% 	current_trial = char(trial_folders(trial_iterator+1));

					% 	current_si = si_matrix(:,trial_iterator);

					% 	sort_ind = []; si_colors = [];
					% 	% Assign Colors
					% 	[s, sort_ind] = sort(current_si);
					% 	si_colors(sort_ind,:) = viridis(total_cells);

					% 	[cmin, cmax] = bounds(current_si);

					% 	figure(2122)
					% 	subplot(1,2,trial_iterator)
					% 	scatter(x,y,20,si_colors)
					% 	colorbar
					% 	colormap('viridis')
					% 	caxis([cmin, cmax])
					% 	title(sprintf('%s', current_trial))
					% 	pbaspect([1 1 1]) 


					% end % End for trial_iterator

					% suptitle(sprintf('%s Spatial Information Overlay on LFP Signature Space', animal_id ))

					% cell_filename = sprintf('%s_spatial_information_on_lfp_signature_space.png', animal_id);
					% cell_image_filepath = fullfile(cell_signature_output_folder, cell_filename);
					% saveas(gcf, cell_image_filepath);



				% -----------------------------------
				% Overlay sparsity
				% -----------------------------------	
					total_trials = size(sparsity_matrix, 2);
					
					% for trial_iterator = 1 : total_trials
% 
						% current_trial = char(trial_folders(trial_iterator+1));

						mean_sparsity = mean(sparsity_matrix');

						[sparsity_cmap , sparsity_colors ] = assign_colors(mean_sparsity, [0 1], 100, 'cividis');
				

						figure(21221)
						scatter(x,y,50,sparsity_colors,'filled')
						colorbar
						colormap(sparsity_cmap)
						caxis([0 1])
						ylim([-10 10])
						xlim([-5 15])
						% title(sprintf('%s', current_trial))
						pbaspect([1 1 1]) 


					% end % End for trial_iterator

					suptitle(sprintf('%s Sparsity Overlay on LFP Signature Space', animal_id ))

					cell_filename = sprintf('%s_sparsity_on_lfp_signature_space.png', animal_id);
					cell_image_filepath = fullfile(cell_signature_output_folder, cell_filename);
					saveas(gcf, cell_image_filepath);

					if make_eps_images 
									
						epsname = sprintf('%s_sparsity_on_lfp_signature_space.eps', animal_id);
						fullname_eps = fullfile(cell_signature_output_folder,epsname);
						set(gcf,'renderer','Painters')
						saveas(gcf,fullname_eps,'epsc')
						

					end % End if make_eps_images	

				% -----------------------------------
				% Overlay cell type (pyramidal, interneuron, unclassified)
				% -----------------------------------	
					% cell_type_vector = zeros(total_cells,1);
					
					% cell_type_vector(find(idx_pyramidal)) = 1;
					% cell_type_vector(find(idx_interneuron)) = 2;


					% idx_unclassified = find(cell_type_vector == 0);
			

					% sort_ind = []; celltype_colors = [];
					
					% % Assign Colors
					% % [s, sort_ind] = sort(cell_type_vector);

					% celltype_colors(find(idx_pyramidal),:) = repmat([1 0 0], sum(idx_pyramidal), 1);
					% celltype_colors(find(idx_interneuron),:) = repmat([0 1 0], sum(idx_interneuron), 1);
					% celltype_colors(idx_unclassified,:) = repmat([0 0 1], sum(cell_type_vector == 0), 1);

					% % [cmin, cmax] = bounds(cell_type_vector);

					% figure(21212)
					% scatter(x,y,20,celltype_colors,'filled')
					% % colorbar
					% % colormap('tab20b')
					% % caxis([cmin, cmax])
					% pbaspect([1 1 1]) 


					% suptitle(sprintf('%s Cell Type Overlay on LFP Signature Space', animal_id ))

					% cell_filename = sprintf('%s_cell_type_on_lfp_signature_space.png', animal_id);
					% cell_image_filepath = fullfile(cell_signature_output_folder, cell_filename);
					% saveas(gcf, cell_image_filepath);


				% -----------------------------------
				% Overlay Group ID
				% -----------------------------------	
					% total_groups = length(unique(groups));

					% group_list = unique(groups);

					% group_colors = zeros(total_cells,3);

					% gc_panel = inferno(total_groups);

					% for group_iterator = 1:total_groups
					% 	group_indices = [];

					% 	group_indices = find(groups == group_list(group_iterator));

					% 	group_colors(group_indices,:) = repmat(gc_panel(group_iterator,:), length(group_indices), 1);

					% end % End


					% figure(21215)
					% scatter(x,y,50,group_colors,'filled')
					% pbaspect([1 1 1]) 
					% cbh = colorbar;
					% colormap('inferno')
					% title(sprintf('Group ID overlay on LFP signatures'))
					% set(cbh,'YTick',[1:total_groups])
				 %    ylabel(cbh,'Group ID')
				 %    caxis([1 total_groups]);

				 %    cell_filename = sprintf('%s_group_id_on_lfp_signature_space.png', animal_id);
					% cell_image_filepath = fullfile(cell_signature_output_folder, cell_filename);
					% saveas(gcf, cell_image_filepath);


			end % End if overlay_cell_params_on_signature_space


		% -------------------------------------------------------------------
		% Sparsity Based Grouping
		% -------------------------------------------------------------------	
			if sparsity_based_grouping

				fprintf('\nSparsity Based Grouping\n')

				try 
					% load(fullfile(cell_signature_output_folder, lfp_signature_combined_group_filename));
					load(fullfile(lfp_output_folder, mean_firing_bins_all_cells_filename));
					fprintf('Cell Signatures, Group Detials and Bin States Info Found\n')
					load(fullfile(root_directory, animal_id, placemaps_data_foldername, placemaps_data_filename ));
					fprintf('Placemaps Data Found\n')
					load(fullfile(cell_signature_output_folder, umap_data_cellsign_combined_filename));
					fprintf('UMAP Projection of Signature Space Found\n')
				catch
					fprintf('Cell Signatures, Group Detials and Bin States Info Not Found\n')
					return
				end
		
				total_trials = length(trial_folders);

				% sparsity_groups = [0 0.3; 0.3 0.5; 0.5 0.7; 0.7 0.9; 0.9 1];
				sparsity_groups = [0 0.5; 0.5 0.8 ; 0.8 1];

				total_groups = size(sparsity_groups,1);

				groups = zeros(total_cells,1);

				minimal_sparsity = mean(sparsity_matrix');

				for group_iterator = 1:total_groups

					sparsity_ll = []; sparsity_ul = []; gi_temp = [];
					sparsity_ll = sparsity_groups(group_iterator,1);
					sparsity_ul = sparsity_groups(group_iterator,2);

					gi_temp = find(minimal_sparsity > sparsity_ll & minimal_sparsity <= sparsity_ul);

					groups(gi_temp) = group_iterator;

				end % End for group_iterator


				group_average_matrix = [];

				group_list = unique(groups);

				total_groups = length(group_list);

				group_colors = magma(total_groups);

				group_cell_data = {};

				nan_columns = find(isnan(sum(lfp_cell_signatures_combined )) );
				awake_columns = [15*15+1 : 3*15*15];

				nan_columns = union(nan_columns, awake_columns);



				non_nan_columms = setdiff(1:4*15*15, nan_columns);
				

				for group_iterator = 1:total_groups

					group_cells_indices = []; group_average = []; sort_ind = [];
					group_bin_colors = []; group_cells = []; temp_signatures_without_nan = [];

					group_cells_indices = find(groups == group_list(group_iterator));

					group_cell_data{1, group_iterator} = group_cells_indices;

					temp_signatures_without_nan = lfp_cell_signatures_combined(:,non_nan_columms);

					% group_cells = zscore(temp_signatures_without_nan(group_cells_indices,:)');
					% group_cells = group_cells';
					
					group_cells = temp_signatures_without_nan(group_cells_indices,:);

					temp_signatures_with_nan = lfp_cell_signatures_combined(group_cells_indices,:);
					temp_signatures_with_nan(:,non_nan_columms) = group_cells;

					if length(group_cells_indices) > 1
						group_average = nanmean(temp_signatures_with_nan);
					else
						group_average = temp_signatures_with_nan;
					end


					group_average_matrix = [group_average_matrix; group_average];

					% group_average_without_nan = zscore(group_average(non_nan_columms));
					% caxis_range = [prctile(group_average_without_nan(:),1) prctile(group_average_without_nan(:),95)];

					group_average_without_nan = zscore(group_average(non_nan_columms));


					% caxis_temp = group_average_without_nan;
					% caxis_temp = lfp_cell_signatures_combined(:,non_nan_columms);
					% caxis_range = [prctile(caxis_temp(:),1) prctile(caxis_temp(:),95)];

					caxis_range = [-3 3];
					
					[group_cmap,  group_bin_colors] = assign_colors(group_average_without_nan,  caxis_range , 100, 'brbg');

					
					final_group_colors(non_nan_columms,:) = group_bin_colors;

					% Gray color for nan columns and awake columns
					final_group_colors(nan_columns,:) = repmat(color_gray,length(nan_columns),1);


					plot_id = 1;
					% Assign colors to each states in a bin
					for trial_iterator = 1 : total_trials

						if ismember(trial_iterator,awake_trial_list)

							continue
						end

						current_trial = char(trial_folders(trial_iterator));

						trial_group_colors = []; trial_bin_range = [];
						trial_state_colors = [];

						trial_bin_start = global_trial_bin_count(trial_iterator) + 1;
						trial_bin_end =  global_trial_bin_count(trial_iterator + 1);

						trial_bin_range = trial_bin_start:trial_bin_end;

						trial_group_colors = final_group_colors(trial_bin_range,:);

						total_bins = length(trial_bin_range);

						for bin_iterator = 1 : total_bins

							bin_states_indices = global_trial_bsi{trial_iterator, bin_iterator};
							
							trial_state_colors(bin_states_indices,:) = repmat(trial_group_colors(bin_iterator,:), length(bin_states_indices), 1);

						end % End for bin_iterator
						

						% Plot trial average for this group
						range_start = global_states_count_lfp(trial_iterator) + 1 ;
						range_end = global_states_count_lfp(trial_iterator + 1) ;

						pv_range = range_start : range_end; 

						x = projection_lfp(pv_range,1);
						y = projection_lfp(pv_range,2);
						% cmin = min(group_average);
						% cmax = max(group_average);

						figure(185)
						set(gcf, 'Position', get(0, 'Screensize'));
						subplot(1,2,plot_id )
						scatter(projection_lfp(:,1), projection_lfp(:,2), 10, color_gray, 'filled')
						hold on 
						scatter(x,y,10, trial_state_colors ,'filled')
						title(sprintf('%s' ,  current_trial ));
						colormap(group_cmap)
						% colormap(brewermap([],'pubugn'))
						caxis( caxis_range )
						h = colorbar;
						ylabel(h, 'Stand. Frequency')
						pbaspect([1 1 1]) 

						plot_id = plot_id + 1;


					end % End for trial_iterator


					suptitle(sprintf('%s Sparsity Group %d LFP Avg. Signature RAW SCALE(n = %d cells) ',animal_id, group_list(group_iterator), length(group_cells_indices) ));

					cell_filename = sprintf('%s_sparsity_avg_signature_group_%d_raw.png', animal_id, group_list(group_iterator) );
					cell_image_filepath = fullfile(cell_signature_output_folder, cell_filename);
					saveas(gcf, cell_image_filepath);


					if make_eps_images

						fprintf('1');

						epsname = sprintf('%s_sparsity_avg_signature_group_%d_raw.eps', animal_id, group_iterator );			
						fullname_eps = fullfile(cell_signature_output_folder,epsname);
						set(gcf,'renderer','Painters')
						saveas(gcf,fullname_eps,'epsc')
						

					end % End if make_eps_images	
					

					figure(198)
					h = plot(group_average,'lineWidth',1);
				    set(h, {'color'}, num2cell(group_colors(group_iterator,:), 2));
				    cbh = colorbar;
				    set(cbh,'YTick',[1:total_groups])
				    ylabel(cbh,'Group ID')
				    colormap('magma')
				    caxis([1 total_groups]);
					hold on;
					xlabel('Dimensions')
					ylabel('Values of Avg Group Signature')
					
					group_average = smooth_gaussian(group_average, 5, 20);

					figure(195)
					h = plot(group_average,'lineWidth',1);
				    set(h, {'color'}, num2cell(group_colors(group_iterator,:), 2));
				    cbh = colorbar;
				    set(cbh,'YTick',[1:total_groups])
				    ylabel(cbh,'Group ID')
				    colormap('magma')
				    caxis([1 total_groups]);
					hold on;
					xlabel('Dimensions')
					ylabel('Values of Avg Group Signature')



			
				end % End for group iterator


				suptitle(sprintf('%s LFP Signature in High Dimensional Space (%d Groups)  ',animal_id , total_groups));

				cell_filename = sprintf('%s_avg_signature_group_high_dim.png', animal_id  );
				cell_image_filepath = fullfile(cell_signature_output_folder, cell_filename);
				saveas(figure(195), cell_image_filepath);

				group_cell_data_filepath = fullfile(cell_signature_output_folder, sparsity_group_cell_data_filename);
				save(group_cell_data_filepath, 'group_cell_data')


				% -----------------------------------
				% Overlay Group ID
				% -----------------------------------	
					total_groups = length(unique(groups));

					group_list = unique(groups);

					group_colors = zeros(total_cells,3);

					% gc_panel = inferno(total_groups);


					gc_panel = [63, 133, 143;
								147, 171, 67;
								207, 103, 117;
								];

					gc_panel = gc_panel / 255;


					for group_iterator = 1:total_groups
						group_indices = [];

						group_indices = find(groups == group_list(group_iterator));

						group_colors(group_indices,:) = repmat(gc_panel(group_iterator,:), length(group_indices), 1);

					end % End

					total_cells = size(mean_firing_matrix, 1);
					x = umap_output_cell_signatures_combined(:,1);
					y = umap_output_cell_signatures_combined(:,2);


					figure(21215)
					scatter(x,y,70,group_colors,'filled')
					pbaspect([1 1 1]) 
					cbh = colorbar;
					colormap('inferno')
					title(sprintf('Sparsity Group ID overlay on LFP signatures'))
					set(cbh,'YTick',[1:total_groups])
				    ylabel(cbh,'Group ID')
				    caxis([1 total_groups]);
				   ylim([-10 10])
						xlim([-5 15])



				    cell_filename = sprintf('%s_sparsity_group_id_on_lfp_signature_space.png', animal_id);
					cell_image_filepath = fullfile(cell_signature_output_folder, cell_filename);
					saveas(gcf, cell_image_filepath);


					if make_eps_images

						fprintf('1');

						epsname = sprintf('%s_sparsity_group_id_on_lfp_signature_space.eps', animal_id);			
						fullname_eps = fullfile(cell_signature_output_folder,epsname);
						set(gcf,'renderer','Painters')
						saveas(gcf,fullname_eps,'epsc')
						

					end % End if make_eps_images	
				
			end % End if sparsity_based_grouping
		


	






			return
	end % End for animal_iterator

toc;

