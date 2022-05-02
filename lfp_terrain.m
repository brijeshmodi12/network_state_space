% ===================================================================
	% Compute and Plot the terrain of state space LFP
	% 

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

	% Date Created: 10 Nov 2020
% ===================================================================

clear all; close all; clc; tic;

% -------------------------------------------------------------------
% General Constants
% -------------------------------------------------------------------
	% Paths and Filenames
		% root_directory = 'C:\Users\Matteo\Desktop\Brijesh\Open Field Dataset\';
		root_directory = 'D:\matteo_datasets\Open Field Datasets\';

		animal_folders_list = ["AH2","AH3","AH4","AH5"];

		% lfp_root = 'csd_manifold_awake_sleep_pca_global_rescale';
		lfp_root = 'network_state_space_bin200';

		output_foldername = 'lfp_terrain_video'; 

		file_to_get = 'state_space_data.mat';

		total_folders = length(animal_folders_list);

		color_gray = [134, 138, 145 ] ./ 255;
		color_pink = [189, 11, 73] ./ 255;
		color_hist =  [134, 138, 145 ] ./ 255;	

		% Combine all the bands in an single array so you can apply loops 
		bands_strings = ["Delta", "Theta", "Spindle", "Slow gamma", "Medium gamma", "Fast gamma"];

		bin_size_sec = 0.2;

		% 1 for processing
		compute_eccentricity_and_plot = 01;
		
		compute_transition_matrix = 0;

		% For shuffled sleep trials 
		shuffle_transition_trajectories = 0;

		plot_difference_matrix = 0;

		plot_transition_probability = 0;

		compute_power_transitions = 0;

		make_eps_images = 00;
	

		do_sleep_trials = 1;
		do_all_trials = 0;

		% outdated
		overlay_eccentricity =  0;

		next_states_jump = 5;

		layer_letters =  ["P", "R", "S"];
		total_layers  = length(layer_letters);
		total_freq_bands = 6;

		% Make strings for layers
		for layer_iterator = 1:total_layers

			for band_iterator = 1:total_freq_bands

				ind =  (layer_iterator - 1)*total_freq_bands + band_iterator;

				% x_axis_string = strcat(layer_letters(layer_iterator),'-',bands_strings2(band_iterator));

				% x_ticks_array(ind) = x_axis_string;

				x_axis_string2 = strcat(layer_letters(layer_iterator),num2str(band_iterator));

				x_ticks_array2(ind) = x_axis_string2;

			end % End for band_iterator

		end % End for layer_iterator

		% Sorts all delta bands together, all theta bands together and so on.
		sorted_band_indices = reshape(1:18, total_freq_bands, []);
		sorted_band_indices = reshape(sorted_band_indices', [], 1 ) ;

		

% -------------------------------------------------------------------
% Core Loop
% -------------------------------------------------------------------
	% Core Loop
	for animal_iterator = 1:total_folders


		animal_id = char(animal_folders_list(animal_iterator))
		full_filepath = fullfile(root_directory, animal_id, lfp_root, file_to_get);

		try
			load(full_filepath);
		catch
			fprintf('File not found\n');
			return
		end

		
		% -------------------------------------------------------------------
		% Animal specific Constants
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

					% lfp_filename = '135_CH42.continuous'; % AH2 

					


				case 'AH3'
					trial_folders = ["preSleep", "1Square", "2Circle", "postSleep"]; % AH3 
					output_folder_string = 	'network_state_space_bin200';
					
					% trial_folders = [ "1Square", "2Circle" ]; % AH3
					% output_folder_string = 	'all_trials_power_manifold_awake';
					% only_awake = 1;

					% trial_folders = ["preSleep", "postSleep"]; % AH3  
					% output_folder_string = 	'all_trials_power_manifold_sleep';
					% only_sleep = 1;

					% lfp_filename = '159_CH46.continuous'; % AH3 

					
					

				case 'AH4'	
					trial_folders = ["preSleep", "1Circle", "2Square", "postSleep"]; % AH4 
					output_folder_string = 'network_state_space_bin200';

					% trial_folders = [ "1Circle", "2Square" ]; % AH4 
					% output_folder_string = 'all_trials_power_manifold_awake';
					% only_awake = 1;

					% trial_folders = ["preSleep", "postSleep"]; % AH4 
					% output_folder_string = 	'all_trials_power_manifold_sleep';
					% only_sleep = 1;

					% lfp_filename = '126_CH46.continuous'; % AH4

					


				case 'AH5'	
					trial_folders = ["preSleep", "1Circle", "2Square", "postSleep"]; % AH5 
					output_folder_string = 'network_state_space_bin200';

					% trial_folders = ["1Circle", "2Square"]; % AH5 
					% output_folder_string = 'all_trials_power_manifold_awake';
					% only_awake = 1;

					% trial_folders = ["preSleep", "postSleep"]; % AH5
					% output_folder_string = 	'all_trials_power_manifold_sleep';
					% only_sleep = 1;
					
					% lfp_filename = '101_CH42.continuous'; % AH5

					% lfp_cut_off_vector
					

				otherwise
					fprintf('Incorrect animal id\n')
					return
			end % End switch animal_id
		
			trial_colors = magma(total_trials+5);
			trial_colors = trial_colors(end-3:end,:);

			incoming_tm_across_trials = [];
			outgoing_tm_across_trials = [];

			% Bin-id of states for all trials
			global_bin_id_of_states = {};


			new_tm_across_trials = [];

		
			trial_list = [];

			if do_sleep_trials
				trial_list = [1,4];
				% For Trajectories
				no_xbins = 10;
				no_ybins = 10;	
				eps_bins = [4 9];

				trial_string = 'Sleep';
				increment_grid = 1;

				fprintf('Processing Sleep Trials\n')
			else
				trial_list = [2,3];
				% For Trajectories
				no_xbins = 10;
				no_ybins = 10;
				eps_bins = [49 ,68];
				increment_grid = 10;
				trial_string = 'Awake';
				fprintf('Processing Awake Trials\n')

			end

			% if do_all_trials 
			% 	trial_list = 1:total_trials;
			% end 
			total_bins = no_xbins * no_ybins;

			bin_number = no_xbins;
			bin_folder_name = sprintf('bin_%d_%s', bin_number, trial_string);

			output_root_folder = fullfile(root_directory,animal_id,lfp_root,output_foldername);
			mkdir(output_root_folder);

			bin_output_folder = fullfile(root_directory,animal_id,lfp_root,output_foldername, bin_folder_name);

			global_output_folder = bin_output_folder;
			mkdir(global_output_folder);

			temp_counter = 1;


		% -------------------------------------------------------------------
		% Process each trial
		% -------------------------------------------------------------------	
			for trial_iterator = trial_list

				current_trial = char(trial_folders(trial_iterator));

				% -------------------------------------------------------------------
				% Trial Specific Constants
				% -------------------------------------------------------------------	
					n_components = 2;

					range_start = global_states_count_lfp(trial_iterator) + 1 ;
					range_end = global_states_count_lfp(trial_iterator + 1) ;

					pv_range = range_start : range_end; 
					x = [];  y = []; z = [];

					x = projection_lfp(pv_range,1);
					y = projection_lfp(pv_range,2);
					% if n_components > 2
					% 	z = projection_lfp(pv_range,3);
					% end

					% % umap output for this trial
					% Y = [x, y];

					population_vector = [];
					population_vector = global_lfp_states(pv_range,:);

					total_pop_vec = size(population_vector,1);
					total_bands = size(population_vector,2);

					incoming_eccentricities = zeros(1,total_pop_vec);
					outgoing_eccentricities = zeros(1,total_pop_vec);

						

				% -------------------------------------------------------------------
				% Compute Eccentricity and Plot Bin Trajectories
				% -------------------------------------------------------------------	
					if compute_eccentricity_and_plot


						% X and Y axis reference points for umap
						% X and Y axis reference points for umap
						x_edges = linspace(min(projection_lfp(:,1)), max(projection_lfp(:,1)), no_xbins+1);
						y_edges = linspace(min(projection_lfp(:,2)), max(projection_lfp(:,2)), no_ybins+1);


						% Get xindices
						[xcounts, x_edges, xindices] = histcounts(x, 'BinEdges', x_edges);
						[ycounts, y_edges, yindices] = histcounts(y, 'BinEdges', y_edges);

						total_bins = no_xbins * no_ybins;

						plots_per_fig = 9;	
						subplot_id = 1;
						fig_id = 1;

						next_states_jump = 5;
						bin_id = 1;
					
						for xbin_iterator = 1:no_xbins

							tempx =	find(xindices == xbin_iterator);

							for ybin_iterator = 1:no_ybins

								bin_id = (xbin_iterator - 1) * no_xbins + ybin_iterator;

								tempy =	find(yindices == ybin_iterator);

								bin_states_indices = intersect(tempx,tempy);
								bsi = bin_states_indices;								

								if length(bsi) == 0 
									continue;
								end
								
								next_states = bsi + 1;
								previous_states = bsi - 1;

								next_states_2 = bsi + next_states_jump;

								next_states(next_states > total_pop_vec) = total_pop_vec;
								previous_states(previous_states<1) = 1;

								next_states_2(next_states_2 > total_pop_vec) = total_pop_vec;

				
								total_trajectories = length(bin_states_indices);

								incoming_colors = viridis(total_trajectories);
								outgoing_colors = magma(total_trajectories);

								trajectory_start = bsi - next_states_jump;
								trajectory_end = bsi + next_states_jump;

								trajectory_length = 2*next_states_jump + 1;
								trajectory_length = 6;
								incoming_t_colors = viridis(trajectory_length+5);
								incoming_t_colors  = incoming_t_colors(6:end,:);

								outgoing_t_colors = magma(trajectory_length+5);
								outgoing_t_colors  = outgoing_t_colors(6:end,:);

								trajectory_fig_id = 100 + bin_id; 


								% ----------------------------------------	
								% Incoming trajectories
								% ----------------------------------------

									figure(trajectory_fig_id)
									subplot(1,2,1)
									set(gcf, 'Position', get(0, 'Screensize'));
									scatter(projection_lfp(pv_range,1),projection_lfp(pv_range,2), 20, color_gray, 'filled')
									hold on

									bin_incoming_indices = [];

									for trajectory_iterator = 1 : total_trajectories

										trajectory_indices = trajectory_start(trajectory_iterator) : bsi(trajectory_iterator);
										ti = trajectory_indices;

										ti(ti<1) = 1;
										ti(ti>total_pop_vec) = total_pop_vec;
										bin_incoming_indices = [bin_incoming_indices ti];
									
										plot(x(ti),y(ti),  'Color', incoming_colors(trajectory_iterator,:),'lineWidth', 1)
										hold on
										% scatter(x(ti), y(ti), 20, incoming_t_colors, 'filled')
										% hold on
										% scatter(x(bsi(trajectory_iterator)), y(bsi(trajectory_iterator)), 50, 'k', 'filled')
										% hold on
									
									end % End for trajectory_iterator
									xmin = x_edges(xbin_iterator);
									xmax = x_edges(xbin_iterator+1) ; 
									ymin = y_edges(ybin_iterator);
									ymax = y_edges(ybin_iterator+1) ; 
									
									width_rectangle = abs(xmax - xmin);
									height_rectangle = abs(ymax - ymin);
									rectangle('Position',[xmin, ymin, width_rectangle, height_rectangle],'EdgeColor','k','LineWidth',3)
									hold on

									% scatter(median(x(bsi)), median(y(bsi)), 50, 'k', 'filled');

									title(sprintf('%s Bin: %d Incoming', char(trial_folders(trial_iterator)), bin_id))
									xticks(x_edges)
									yticks(y_edges)
									xlim([min(x_edges) max(x_edges)])
									ylim([min(y_edges) max(y_edges)])
									grid on
									ax = gca;
									ax.GridLineStyle = '-';
									ax.GridColor = 'k';
									ax.GridAlpha = 1; % maximum line opacity
									set (gca, 'xticklabel' , {[]});
									set (gca, 'yticklabel' , {[]});
									pbaspect([1 1 1]) 

							

								
													

								% ----------------------------------------	
								% Outgoing trajectories
								% ----------------------------------------
									figure(trajectory_fig_id)
									subplot(1,2,2)
									scatter(projection_lfp(pv_range,1),projection_lfp(pv_range,2), 20, color_gray, 'filled')
									hold on

									bin_outgoing_indices = [];

									for trajectory_iterator = 1 : total_trajectories

										trajectory_indices =  bsi(trajectory_iterator) : trajectory_end(trajectory_iterator);
										ti = trajectory_indices;

										ti(ti<1) = 1;
										ti(ti>total_pop_vec) = total_pop_vec;
										bin_outgoing_indices = [bin_outgoing_indices ti];
									
										plot(x(ti),y(ti),  'Color', outgoing_colors(trajectory_iterator,:),'lineWidth', 1)
										hold on
										% scatter(x(ti), y(ti), 20, outgoing_t_colors, 'filled')
										% hold on
										% scatter(x(bsi(trajectory_iterator)), y(bsi(trajectory_iterator)), 50, 'k', 'filled')
										% hold on

									end % End for trajectory_iterator

									xmin = x_edges(xbin_iterator);
									xmax = x_edges(xbin_iterator+1) ; 
									ymin = y_edges(ybin_iterator);
									ymax = y_edges(ybin_iterator+1) ; 
									
									width_rectangle = abs(xmax - xmin);
									height_rectangle = abs(ymax - ymin);
									rectangle('Position',[xmin, ymin, width_rectangle, height_rectangle],'EdgeColor','k','LineWidth',3)
									hold on 
									% scatter(median(x(bsi)), median(y(bsi)), 50, 'k', 'filled');

									title(sprintf('%s Bin: %d Outgoing', char(trial_folders(trial_iterator)), bin_id))
									xticks(x_edges)
									yticks(y_edges)
									grid on
									ax = gca;
									ax.GridLineStyle = '-';
									ax.GridColor = 'k';
									ax.GridAlpha = 1; % maximum line opacity
									set (gca, 'xticklabel' , {[]});
									set (gca, 'yticklabel' , {[]});
									pbaspect([1 1 1]) 
									xlim([min(x_edges) max(x_edges)])
									ylim([min(y_edges) max(y_edges)])	

									suptitle(sprintf('%s %s Bin: %d',animal_id, char(trial_folders(trial_iterator)), bin_id))


									

									filename = sprintf('%s_%s_bin_%d.png', animal_id, char(trial_folders(trial_iterator)), bin_id );
									fullname = fullfile(global_output_folder,filename);
									saveas(figure(trajectory_fig_id),fullname);
									

									if make_eps_images & ismember(bin_id , eps_bins)
									
										epsname = sprintf('%s_%s_bin_%d.eps', animal_id, char(trial_folders(trial_iterator)), bin_id );
										fullname_eps = fullfile(global_output_folder,epsname);
										set(gcf,'renderer','Painters')
										saveas(gcf,fullname_eps,'epsc')

									end % End if make_eps_images		
									close all;		

							end % End for ybin_iterator

							% return

						end % End for xbin_iterator	


						% filename = sprintf('%s_Eccenctricity_%s.mat', animal_id, char(trial_folders(trial_iterator)) );
						% fullname = fullfile(global_output_folder,filename);
						% save(fullname, 'incoming_eccentricities', 'outgoing_eccentricities');

					end % End if compute_eccentricity_and_plot



				% -------------------------------------------------------------------
				% Compute Transition Matrix
				% -------------------------------------------------------------------
					if compute_transition_matrix

						incoming_transition_matrix = [];
						outgoing_transition_matrix = [];

						fprintf('Computing Transition Matrix\n')

						% X and Y axis reference points for umap
						x_edges = linspace(min(projection_lfp(:,1)), max(projection_lfp(:,1)), no_xbins+1);
						y_edges = linspace(min(projection_lfp(:,2)), max(projection_lfp(:,2)), no_ybins+1);


						% Get xindices
						[xcounts, x_edges, xindices] = histcounts(x, 'BinEdges', x_edges);
						[ycounts, y_edges, yindices] = histcounts(y, 'BinEdges', y_edges);

						% next_states_jump = 15;
						
						
						% Vector to store bin id for each state
						bin_id_of_states = [];

						% Bin states indices for each bin
						bin_bsi_array = cell(1, no_xbins *  no_ybins);

						% Collect Bin ID for each states in bin
							for xbin_iterator = 1:no_xbins

								tempx = [];
								tempx =	find(xindices == xbin_iterator);

								for ybin_iterator = 1:no_ybins

									bin_id = (xbin_iterator - 1) * no_xbins + ybin_iterator ;

									tempy = [];
									tempy =	find(yindices == ybin_iterator);

									bin_states_indices = intersect(tempx,tempy);
									bsi = bin_states_indices;
					
									if length(bsi) == 0 
										continue;
									end

									bin_id_of_states(bsi) = bin_id;

									bin_bsi_array{1,bin_id} = bsi;

								end % End for ybin_iterator
							
							end % End for xbin_iterator

						global_bin_id_of_states{1,temp_counter} = bin_id_of_states;

						
						total_histbins = size(bin_bsi_array,2);

						hist_edges_start = 1:total_histbins;
						hist_edges_start = hist_edges_start - 0.5;
						
						hist_edges_end = hist_edges_start + 1;

						edges = [hist_edges_start' hist_edges_end'];

						trajectory_start1 = bin_id_of_states(1:end - next_states_jump);
						trajectory_end1 = bin_id_of_states(next_states_jump + 1 : end );

						% Compuate new transition matrix
							transition_matrix =  make_transition_matrix(trajectory_start1, trajectory_end1, edges,total_bins);

							new_tm_across_trials(:,:,temp_counter) = transition_matrix;						


						% Plot transition matrix
							figure(123)
							set(gcf, 'Position', get(0, 'Screensize'));
							subplot(1,2,temp_counter)
							imagesc(transition_matrix)
							pbaspect([1 1 1])
							colormap('viridis')
							cbh = colorbar;
							ylabel(cbh,'Probability')
							caxis([0 1]);
							title(sprintf('%s',current_trial))
							% Add Grids
							grid_edges = 0.5:increment_grid:total_bins - 0.5;
							xticks(grid_edges)
							yticks(grid_edges)
							grid on
							ax = gca;
							ax.GridLineStyle = '-';
							ax.GridColor = 'k';
							ax.GridAlpha = 1; % maximum line opacity
							set (gca, 'xticklabel' , {[]});
							set (gca, 'yticklabel' , {[]});


							temp_counter = temp_counter + 1;	


							
							

						% % Compute transition probability for each bin
							% 	for bin_id = 1:total_histbins

							% 		current_bsi = bin_bsi_array{1,bin_id};

							% 		if length(current_bsi) < 10
							% 			incoming_transition_matrix = [incoming_transition_matrix ; NaN * zeros(1,total_histbins)];
							% 			outgoing_transition_matrix = [outgoing_transition_matrix ; NaN * zeros(1,total_histbins)];
							% 			continue;
							% 		end

							% 		total_bsi = length(current_bsi);

							% 		incoming_bsi = []; incoming_bin_id = [];
							% 		outgoing_bsi = []; outgoing_bin_id = [];
							% 		incoming_bin_colors = []; outgoing_bin_colors = [];
							% 		incoming_states_colors = []; outgoing_states_colors = []; 

							% 		total_trajectories = total_bsi;
							% 		trajectory_start_incoming = current_bsi - next_states_jump;
							% 		trajectory_end_outgoing = current_bsi + next_states_jump;


							% 		for trajectory_iterator = 1 : total_trajectories
							% 			ti = [];trajectory_indices = []; 
							% 			% trajectory_indices = trajectory_start_incoming(trajectory_iterator) : current_bsi(trajectory_iterator);
							% 			trajectory_indices = [trajectory_start_incoming(trajectory_iterator) ] ;
							% 			ti = trajectory_indices;

							% 			ti(ti<1) = 1;
							% 			ti(ti>total_pop_vec) = total_pop_vec;
							% 			incoming_bsi = [incoming_bsi ti];
										

							% 			ti = []; trajectory_indices = [];
							% 			% trajectory_indices = current_bsi(trajectory_iterator): trajectory_end_outgoing(trajectory_iterator) ;
							% 			trajectory_indices = [ trajectory_end_outgoing(trajectory_iterator) ];
							% 			ti = trajectory_indices;

							% 			ti(ti<1) = 1;
							% 			ti(ti>total_pop_vec) = total_pop_vec;
							% 			outgoing_bsi = [outgoing_bsi ti];
										
							% 		end % End for trajectory_iterator

							% 		incoming_bsi(incoming_bsi<1) = 1;
							% 		outgoing_bsi(outgoing_bsi > total_pop_vec) = total_pop_vec;

							% 		% Id of incoming bins
							% 		incoming_bin_id = bin_id_of_states(incoming_bsi);

							% 		% Id of outgoing bins
							% 		outgoing_bin_id = bin_id_of_states(outgoing_bsi);				
									
							% 		% Compute probability of transition to all bins
							% 		incoming_transition_prob = custom_histcounts(incoming_bin_id, edges);
							% 		incoming_transition_prob = incoming_transition_prob / sum(incoming_transition_prob);

							% 		outgoing_transition_prob = custom_histcounts(outgoing_bin_id, edges);
							% 		outgoing_transition_prob = outgoing_transition_prob / sum(outgoing_transition_prob);

							% 		incoming_transition_matrix = [incoming_transition_matrix ; incoming_transition_prob'];
							% 		outgoing_transition_matrix = [outgoing_transition_matrix ; outgoing_transition_prob'];

							% 		% Replace zeros with NaN for coloring
							% 			zero_itp_ind = find(incoming_transition_prob == 0);
							% 			zero_otp_ind = find(outgoing_transition_prob == 0);

							% 			incoming_transition_prob(zero_itp_ind) = NaN;
							% 			outgoing_transition_prob(zero_itp_ind) = NaN;
								

							% 		% Assign Colors -- Incoming
							% 			[incoming_colormap, incoming_bin_colors] = assign_colors(incoming_transition_prob, [0 1], 100, 'viridis');
										
							% 			for bin_iterator = 1:total_histbins

							% 				temp_states_ind = [];
							% 				temp_states_ind = bin_bsi_array{1,bin_iterator};
							% 				if isempty(temp_states_ind)
							% 					continue;
							% 				end

							% 				incoming_states_colors(temp_states_ind,:) = repmat(incoming_bin_colors(bin_iterator,:) , length(temp_states_ind),1);

							% 			end % End for bin_iterator

										
							% 		% Plot Incoming
							% 			if plot_transition_probability
							% 				tp_fig_id = 821;
							% 				figure(tp_fig_id)
							% 				set(gcf, 'Position', get(0, 'Screensize'));
							% 				subplot(1,2,1)
							% 				scatter(projection_lfp(:,1), projection_lfp(:,2), 10, color_gray, 'filled')
							% 				hold on
							% 				scatter(x,y, 10, incoming_states_colors, 'filled')
							% 				pbaspect([1 1 1])
							% 				% colormap(brewermap('*ylgnbu'))
							% 				colormap(gca, incoming_colormap)
							% 				cbh = colorbar;
							% 				ylabel(cbh,'Probability')
							% 				% [cmin, cmax] = bounds(incoming_transition_prob);
							% 				% set(cbh, 'ylim', [cmin cmax])
							% 				caxis([0 1]);					
										
							% 				% Plot square
							% 				[xmin, xmax] = bounds(x(current_bsi)); 
							% 				[ymin, ymax] = bounds(y(current_bsi)); 

							% 				width_rectangle = abs(xmax - xmin);
							% 				height_rectangle = abs(ymax - ymin);
							% 				rectangle('Position',[xmin, ymin, width_rectangle, height_rectangle],'EdgeColor','k','LineWidth',3)

							% 				title(sprintf('%s Bin: %d Incoming', char(trial_folders(trial_iterator)), bin_id))
							% 				xticks(x_edges)
							% 				yticks(y_edges)
							% 				grid on
							% 				ax = gca;
							% 				ax.GridLineStyle = '-';
							% 				ax.GridColor = 'k';
							% 				ax.GridAlpha = 1; % maximum line opacity
							% 				set (gca, 'xticklabel' , {[]});
							% 				set (gca, 'yticklabel' , {[]});

											
							% 			end % End if plot_transition_probability


							% 		% Assign Colors -- Outgoing
							% 			[outgoing_colormap, outgoing_bin_colors] = assign_colors(outgoing_transition_prob, [0 1], 100, 'magma');

							% 			for bin_iterator = 1:total_histbins

							% 				temp_states_ind = [];
							% 				temp_states_ind = bin_bsi_array{1,bin_iterator};
							% 				if isempty(temp_states_ind)
							% 					continue;
							% 				end

							% 				outgoing_states_colors(temp_states_ind,:) = repmat(outgoing_bin_colors(bin_iterator,:) , length(temp_states_ind),1);

							% 			end % End for bin_iterator


							% 		% Plot Outgoing
							% 			if plot_transition_probability
							% 				tp_fig_id = 821;
							% 				figure(tp_fig_id)
							% 				subplot(1,2,2)
							% 				scatter(projection_lfp(:,1), projection_lfp(:,2), 10, color_gray, 'filled')
							% 				hold on
							% 				scatter(x,y, 10, outgoing_states_colors, 'filled')
							% 				pbaspect([1 1 1])
							% 				% colormap(brewermap('*ylgnbu'))
							% 				colormap(gca, outgoing_colormap)
							% 				cbh = colorbar;
							% 				ylabel(cbh,'Probability')
							% 				% [cmin, cmax] = bounds(outgoing_transition_prob);
							% 				% set(cbh, 'ylim', [cmin cmax])
							% 				caxis([0 1]);
														
							% 				% Plot square
							% 					[xmin, xmax] = bounds(x(current_bsi)); 
							% 					[ymin, ymax] = bounds(y(current_bsi)); 

							% 					width_rectangle = abs(xmax - xmin);
							% 					height_rectangle = abs(ymax - ymin);
							% 					rectangle('Position',[xmin, ymin, width_rectangle, height_rectangle],'EdgeColor','k','LineWidth',3)

							% 				title(sprintf('%s Bin: %d Outgoing', char(trial_folders(trial_iterator)), bin_id))
							% 				xticks(x_edges)
							% 				yticks(y_edges)
							% 				grid on
							% 				ax = gca;
							% 				ax.GridLineStyle = '-';
							% 				ax.GridColor = 'k';
							% 				ax.GridAlpha = 1; % maximum line opacity
							% 				set (gca, 'xticklabel' , {[]});
							% 				set (gca, 'yticklabel' , {[]});

										

							% 				% Save Plot
							% 				suptitle(sprintf('%s %s Bin %d Transition Probability LFP', animal_id, current_trial, bin_id ))
							% 				filename = sprintf('%s_%s_bin_%d_transition_map.png', animal_id, current_trial ,bin_id);
							% 				fullname = fullfile(global_output_folder,filename);
							% 				saveas(figure(tp_fig_id),fullname);

							% 				close(figure(tp_fig_id))
							% 			end % End if plot_transition_probability



										
							% 	end % End for bin_id 


							% incoming_transition_matrix(incoming_transition_matrix > 0.5) = 0;
							% 	itm_copy = incoming_transition_matrix;
							% 	nan_rows = find(isnan(sum(incoming_transition_matrix,2)));
							% 	incoming_transition_matrix(nan_rows,:) = [];
							% 	incoming_transition_matrix(:, nan_rows) = [];


							% % Plot Incoming Transition Matrix
							% 	tm_fig_id = 123;
							% 	figure(tm_fig_id)
							% 	% set(gcf, 'Position', get(0, 'Screensize'));
							% 	subplot(1,2,1)
							% 	imagesc(incoming_transition_matrix)
							% 	title(sprintf('%s Incoming',current_trial))
							% 	pbaspect([1 1 1])
							% 	colormap(gca,viridis)
							% 	cbh = colorbar;
							% 	ylabel(cbh,'Probability')
							% 	[cmin, cmax] = bounds(incoming_transition_matrix(:));
							% 	caxis([0 1]);
							% 	cbh.FontSize = 12;

							% otm_copy = outgoing_transition_matrix;
							% outgoing_transition_matrix(nan_rows,:) = [];
							% outgoing_transition_matrix(:, nan_rows) = [];

							% % Plot Outgoing Transition Matrix
							% 	tm_fig_id = 123;
							% 	figure(tm_fig_id)
							% 	subplot(1,2,2)
							% 	imagesc(outgoing_transition_matrix)
							% 	title(sprintf('%s Outgoing',current_trial))
							% 	pbaspect([1 1 1])
							% 	colormap(gca,magma)
							% 	cbh = colorbar;
							% 	ylabel(cbh,'Probability')
							% 	caxis([0 1]);
							% 	cbh.FontSize = 12;

							% 	% Save Plot
							% 	suptitle(sprintf('%s %s Transition Matrix LFP', animal_id, current_trial ))
							% 	filename = sprintf('%s_%s_transition_matrix.png', animal_id, current_trial);
							% 	fullname = fullfile(global_output_folder,filename);
							% 	saveas(figure(tm_fig_id),fullname);
								

							% 	if make_eps_images

							% 		epsname = sprintf('%s_%s_transition_matrix.eps', animal_id, current_trial);
							% 		fullname_eps = fullfile(global_output_folder,epsname);
							% 		set(gcf,'renderer','Painters')
							% 		saveas(gcf,fullname_eps,'epsc')

							% 	end % End if make_eps_images

							% 	close(figure(tm_fig_id))




							% incoming_tm_across_trials(:,:,trial_iterator)  = itm_copy;
							% outgoing_tm_across_trials(:,:,trial_iterator)  = otm_copy;


						if trial_iterator >= 3
							if trial_iterator == 3
								str = 'awake';
							else
								str = 'sleep';
							end 

							filename = sprintf('%s_global_transition_data_%s.mat', animal_id, str );
							fullname = fullfile(global_output_folder,filename);
							save(fullname, 'global_bin_id_of_states','total_histbins', 'new_tm_across_trials');
							fprintf('Data Saved\n')


							figure(123)
							filename = sprintf('%s_transition_matrix_%s.png', animal_id, trial_string);
							fullname = fullfile(global_output_folder,filename);
							saveas(gcf,fullname)


							if make_eps_images

								epsname = sprintf('%s_transition_matrix_%s.eps', animal_id, trial_string);
								fullname_eps = fullfile(global_output_folder,epsname);
								set(gcf,'renderer','Painters')
								saveas(gcf,fullname_eps,'epsc')

							end % End if make_eps_images	


						end % End if trial_iterator




					end % End if compute_transition_matrix



				% -------------------------------------------------------------------
				% Compute Trajectory and Power Transitions
				% -------------------------------------------------------------------
					if compute_power_transitions

						fprintf('Computing Power Transitions\n')

						% X and Y axis reference points for umap
						x_edges = linspace(min(projection_lfp(:,1)), max(projection_lfp(:,1)), no_xbins+1);
						y_edges = linspace(min(projection_lfp(:,2)), max(projection_lfp(:,2)), no_ybins+1);


						% Get xindices
						[xcounts, x_edges, xindices] = histcounts(x, 'BinEdges', x_edges);
						[ycounts, y_edges, yindices] = histcounts(y, 'BinEdges', y_edges);

						% next_states_jump = 15;
						
						
						% Vector to store bin id for each state
						bin_id_of_states = [];

						% Bin states indices for each bin
						bin_bsi_array = {};

						% Collect Bin ID for each states in bin
							for xbin_iterator = 1:no_xbins

								tempx = [];
								tempx =	find(xindices == xbin_iterator);

								for ybin_iterator = 1:no_ybins

									bin_id = (xbin_iterator - 1) * no_xbins + ybin_iterator ;

									tempy = [];
									tempy =	find(yindices == ybin_iterator);

									bin_states_indices = intersect(tempx,tempy);
									bsi = bin_states_indices;
					
									if length(bsi) == 0 
										continue;
									end

									bin_id_of_states(bsi) = bin_id;

									bin_bsi_array{1,bin_id} = bsi;



								end % End for ybin_iterator
							
							end % End for xbin_iterator


						total_histbins = size(bin_bsi_array,2);

		
						% Compute transition probability for each bin
							for bin_id = 1:total_histbins
								fprintf('Processing Bin %d...\n', bin_id)

								current_bsi = bin_bsi_array{1,bin_id};

								if numel(current_bsi)  < 1
									continue;
								end

								total_bsi = length(current_bsi);
							
								outgoing_bsi = []; outgoing_bin_id = [];

								outgoing_bsi = current_bsi + next_states_jump;	

								outgoing_bsi(outgoing_bsi > total_pop_vec) = total_pop_vec;

								% Id of outgoing bins
								outgoing_bin_id = bin_id_of_states(outgoing_bsi);

								all_output_bins = unique(outgoing_bin_id);

								total_output_bins = numel(all_output_bins);

								% Compute power transitionitions for each output bins
								for output_bin_iterator = 1:total_output_bins

									current_output_bin_id =  all_output_bins(output_bin_iterator);

									% Skip tranisiton within the same bin
									if current_output_bin_id == bin_id 
										continue;
									end

									% indices where trajectory ends in current output bin
									temp_ind = find(outgoing_bin_id == current_output_bin_id);

									% use these indices to get trajectories
									trajectories_start_indices = current_bsi(temp_ind);
									trajectories_end_indices = outgoing_bsi(temp_ind);

									total_trajectories = length(trajectories_start_indices);
									trajectory_colors = magma(total_trajectories);

									trajectory_power = [];average_power = []; zscored_avg_power = [];

									% Plot trajectory to confirm 
									trajectory_fig_id = 987;
									figure(trajectory_fig_id)
									subplot(3,3,output_bin_iterator)
									set(gcf, 'Position', get(0, 'Screensize'));
									scatter(projection_lfp(:,1),projection_lfp(:,2), 20, color_gray, 'filled')
									hold on

									tpower_index = 1;
									
									for trajectory_iterator = 1:total_trajectories

										ti = trajectories_start_indices(trajectory_iterator) : trajectories_end_indices(trajectory_iterator);
										if numel(ti) <= next_states_jump
											continue
										end
										plot(x(ti),y(ti),  'Color', trajectory_colors(trajectory_iterator,:),'lineWidth', 1)
										hold on

										trajectory_power(:,:,tpower_index) = population_vector(ti,:);
										tpower_index = tpower_index + 1;

									end % End for trajectory_iterator

									title(sprintf('Output ID %d', current_output_bin_id))
									
									[xmin, xmax] = bounds(x(current_bsi)); 
									[ymin, ymax] = bounds(y(current_bsi)); 

									
									width_rectangle = abs(xmax - xmin);
									height_rectangle = abs(ymax - ymin);
									rectangle('Position',[xmin, ymin, width_rectangle, height_rectangle],'EdgeColor','k','LineWidth',3)
									
									output_bin_bsi = bin_bsi_array{1,current_output_bin_id};
									[xmin, xmax] = bounds(x(output_bin_bsi)); 
									[ymin, ymax] = bounds(y(output_bin_bsi)); 

									
									width_rectangle = abs(xmax - xmin);
									height_rectangle = abs(ymax - ymin);
									rectangle('Position',[xmin, ymin, width_rectangle, height_rectangle],'EdgeColor',color_pink,'LineWidth',3)
									pbaspect([1 1 1])

									% Add Grids
									xticks(x_edges)
									yticks(y_edges)
									grid on
									ax = gca;
									ax.GridLineStyle = '-';
									ax.GridColor = 'k';
									ax.GridAlpha = 1; % maximum line opacity
									set (gca, 'xticklabel' , {[]});
									set (gca, 'yticklabel' , {[]});
									
						

									% Plot Average Trajectory Power
									average_power = mean(trajectory_power,3);

									zscored_avg_power = zscore(average_power);

									zscored_avg_power = zscored_avg_power(:,sorted_band_indices);
									xlabels =  x_ticks_array2(sorted_band_indices);

									power_figure = 12312321;
									figure(power_figure)
									subplot(3,3,output_bin_iterator)
									set(gcf, 'Position', get(0, 'Screensize'));
									imagesc(zscored_avg_power')
									colormap('viridis')
									pbaspect([1 1 1])
									cbh = colorbar;
									ylabel(cbh,'Zscored Power')
									cmin = min(zscored_avg_power(:));
									cmax = max(zscored_avg_power(:));
									caxis([cmin cmax])
									yticks(1:18)
									yticklabels(xlabels)
									xticks(1:1+next_states_jump)
									xticklabels([0:0.2:1]);
									xlabel('Time (seconds)')
									title(sprintf('Output ID %d', current_output_bin_id))
									total_features = 18;

									row_grids = [ 3 6 9 12 15 ];					
									hold on;
									for row = row_grids
									  line([0, total_features+0.5], [row+0.5, row+0.5], 'Color', 'k','lineWidth',2);
									end



						

								end % End output_bin_iterator

								figure(power_figure)
								suptitle(sprintf('%s %s, Average Power during Transition from Bin: %d ',animal_id, current_trial, bin_id));
								filename = sprintf('%s_%s_bin_%d_average_power_transition.png', animal_id, current_trial ,bin_id);
								fullname = fullfile(global_output_folder,filename);
								saveas(figure(power_figure),fullname);


								figure(trajectory_fig_id)
								suptitle(sprintf('%s %s, Trasition Trajectories from Bin %d ',animal_id, current_trial, bin_id));								
								filename = sprintf('%s_%s_bin_%d_transition_trajectories.png', animal_id, current_trial ,bin_id);
								fullname = fullfile(global_output_folder,filename);
								saveas(figure(trajectory_fig_id),fullname);								
							
							
								close all
									
							
									
							end % End for bin_id 

						

					end % End if compute_power_transitions

					
			end % End for trial_iterator	





			% -------------------------------------------------------------------
			% Compute Shuffled Sleep and Awake Trials
			% -------------------------------------------------------------------
				if shuffle_transition_trajectories
					
					fprintf('Shuffling Data \n')

					if do_sleep_trials
						sleep_filename = sprintf('%s_global_transition_data_sleep.mat', animal_id );
						shuffled_filename = fullfile(global_output_folder,sleep_filename);
					else
						awake_filename = sprintf('%s_global_transition_data_awake.mat', animal_id );
						shuffled_filename = fullfile(global_output_folder,awake_filename);
					end
					
					try
						load(shuffled_filename)
						fprintf('Shuffled Raw Data Found\n')
					catch
						fprintf('Shuffled Raw Data Not Found\n')
					end

					
					temp_diff = new_tm_across_trials(:,:,2) - new_tm_across_trials(:,:,1);

					mean_abs_real_data = mean(abs(temp_diff(:)))

					max_val = max(abs(temp_diff(:)));

					[real_rem2rem, real_nrem2nrem] = get_rem_nonrem_transition(abs(temp_diff), animal_id);

				

					figure(124)
					set(gcf, 'Position', get(0, 'Screensize'));
					imagesc(temp_diff)
					pbaspect([1 1 1])
					colormap(gca,brewermap([],'spectral'))
					cbh = colorbar;
					ylabel(cbh,'Probability')
					caxis([ -max_val max_val ]);
					title(sprintf('Difference Matrix %d',1))	
					% Add Grids
					grid_edges = 0.5:increment_grid:total_bins - 0.5;
					xticks(grid_edges)
					yticks(grid_edges)
					grid on
					ax = gca;
					ax.GridLineStyle = '-';
					ax.GridColor = 'k';
					ax.GridAlpha = 1; % maximum line opacity
					set (gca, 'xticklabel' , {[]});
					set (gca, 'yticklabel' , {[]});
					set(gca,'box','off') 

					filename = sprintf('%s_diff_matrix.png', animal_id );
					fullname = fullfile(global_output_folder,filename);
					saveas(gcf,fullname)


					if make_eps_images

						epsname = sprintf('%s_%s_diff_matrix.eps', animal_id, trial_string);
						fullname_eps = fullfile(global_output_folder,epsname);
						set(gcf,'renderer','Painters')
						saveas(gcf,fullname_eps,'epsc')

					end % End if make_eps_images	


					total_iterations = 1000; 

					all_states_bin_id = [global_bin_id_of_states{1,1}, global_bin_id_of_states{1,2} ];



					total_states = length(all_states_bin_id(1:end - next_states_jump));

					all_states_indices = 1:total_states;

					% Make Edges for Histogram (probablility of transitions)
					hist_edges_start = 1:total_histbins;
					hist_edges_start = hist_edges_start - 0.5;
					hist_edges_end = hist_edges_start + 1;
					edges = [hist_edges_start' hist_edges_end'];

					avg_abs_change = zeros(1,total_iterations);

					global_surrogate_rem2rem = [];
					global_surrogate_nonrem2nonrem = [];

					fprintf('Processing Random Iteration\n')
					for random_iterator = 1:total_iterations
						if ~mod(random_iterator, 100) 
							fprintf('%d...',random_iterator)
						end

						surr_rem2rem = []; surr_nrem2nrem = [];
					
						random_indices1 = randi([1 total_states],1, round(total_states/2));

						random_indices2 = setdiff(all_states_indices, random_indices1);

						shuffle_1_start = all_states_bin_id(random_indices1);

						shuffle_2_start = all_states_bin_id(random_indices2);

						shuffle_1_end = all_states_bin_id(random_indices1 + next_states_jump);

						transition_matrix1 = make_transition_matrix(shuffle_1_start, shuffle_1_end, edges, total_bins);	

						shuffle_2_end = all_states_bin_id(random_indices2 + next_states_jump);

						transition_matrix2 = make_transition_matrix(shuffle_2_start, shuffle_2_end, edges, total_bins);	

						diff_matt = transition_matrix2 - transition_matrix1;

						[surr_rem2rem, surr_nrem2nrem] = get_rem_nonrem_transition(abs(diff_matt), animal_id);

						global_surrogate_rem2rem = [global_surrogate_rem2rem; mean(surr_rem2rem)];

						global_surrogate_nonrem2nonrem = [global_surrogate_nonrem2nonrem; mean(surr_nrem2nrem)];

						avg_abs_change(random_iterator) = mean(abs(diff_matt(:)));

					end % End for random_iterator	
					fprintf('\nShuffling Completed\n')



					pval = 95;
					pval_data = prctile(avg_abs_change,pval)
					figure
					histogram([avg_abs_change mean_abs_real_data] ,100,'FaceColor', color_hist, 'EdgeColor', color_hist);
					hold on
					xline(pval_data,'--r','lineWidth',2)
					hold on 
					xline(mean_abs_real_data,'--','Color','#77AC30','lineWidth',2)
					% hold on
					% xline()

					title(sprintf('Real Data = %0.4f, Shuffled 95 percentile = %0.4f', mean_abs_real_data, pval_data ));
					set(gca,'box','off') 
					temp = [avg_abs_change mean_abs_real_data];
					max_val = max(temp)
					max_val = max_val + max_val/2;
					xlim([0 max_val])
					filename = sprintf('%s_shuffle_hist.png', animal_id );
					fullname = fullfile(global_output_folder,filename);
					saveas(gcf,fullname)

					if make_eps_images

						epsname = sprintf('%s_%s_shuffle_hist.eps', animal_id, trial_string);
						fullname_eps = fullfile(global_output_folder,epsname);
						set(gcf,'renderer','Painters')
						saveas(gcf,fullname_eps,'epsc')

					end % End if make_eps_images	


					% Plot surrogate rem2rem and nonrem2nonrem
					pval = 95;
					pval_data = prctile(global_surrogate_rem2rem, pval)
					figure(12321)
					histogram([global_surrogate_rem2rem ; mean(real_rem2rem)] ,100,'FaceColor', color_hist, 'EdgeColor', color_hist);
					hold on
					xline(pval_data,'--r','lineWidth',2)
					hold on 
					xline(mean(real_rem2rem),'--','Color','#77AC30','lineWidth',2)
					title(sprintf('REM Real Data = %0.4f, Shuffled 95 percentile = %0.4f', mean(real_rem2rem), pval_data ));
					real_colors = magma(length(real_rem2rem));

					% for i = 1:length(real_rem2rem)
						
						% hold on
					% end

					
					

					if make_eps_images

						epsname = sprintf('%s_%s_shuffle_hist_rem2rem.eps', animal_id, trial_string);
						fullname_eps = fullfile(global_output_folder,epsname);
						set(gcf,'renderer','Painters')
						saveas(gcf,fullname_eps,'epsc')

					end % End if make_eps_images	



					pval = 95;
					pval_data = prctile(global_surrogate_nonrem2nonrem, pval)
					figure(7898542)
					histogram([global_surrogate_nonrem2nonrem ; mean(real_nrem2nrem)] ,100,'FaceColor', color_hist, 'EdgeColor', color_hist);
					hold on
					xline(pval_data,'--r','lineWidth',2)
					hold on 
					xline(mean(real_nrem2nrem),'--','Color','#77AC30','lineWidth',2)
					title(sprintf('NONREM Real Data = %0.4f, Shuffled 95 percentile = %0.4f', mean(real_nrem2nrem), pval_data ));

					if make_eps_images

						epsname = sprintf('%s_%s_shuffle_hist_nonrem2nonrem.eps', animal_id, trial_string);
						fullname_eps = fullfile(global_output_folder,epsname);
						set(gcf,'renderer','Painters')
						saveas(gcf,fullname_eps,'epsc')

					end % End if make_eps_images	

					return




					


					return

				end % End shuffle_transition_trajectories
						

			% -------------------------------------------------------------------
			% PLot differencein Pre and Post sleep transition Matrix
			% -------------------------------------------------------------------
				if plot_difference_matrix
					ind2 = trial_list(2) ;
					ind1 = trial_list(1);
					incoming_diff = incoming_tm_across_trials(:,:,ind2) - incoming_tm_across_trials(:,:,ind1);
					outgoing_diff = outgoing_tm_across_trials(:,:,ind2) - outgoing_tm_across_trials(:,:,ind1);

					nan_rows = find(isnan(sum(incoming_diff,2)));
					incoming_diff(nan_rows,:) = [];
					incoming_diff(:, nan_rows) = [];

					outgoing_diff(nan_rows,:) = [];
					outgoing_diff(:, nan_rows) = [];

					resultant_incoming_change = sum(abs(incoming_diff(:))) / size(incoming_diff,1);
					resultant_outgoing_change = sum(abs(outgoing_diff(:))) / size(outgoing_diff,1);

					allvalues = [ incoming_diff(:), outgoing_diff(:)];
					
					cmax = max(abs(allvalues(:)));				
					cmin = -cmax;

					tm_fig_id = 1213;
					figure(tm_fig_id)
					% set(gcf, 'Position', get(0, 'Screensize'));

					subplot(1,2,1)
					imagesc(incoming_diff)
					title(sprintf('%s Incoming Diff',current_trial))
					pbaspect([1 1 1])
					colormap(gca,brewermap([],'spectral'))
					% colormap( [0 0 0; brewermap([],'spectral')] )
					cbh = colorbar;
					ylabel(cbh,'Difference in Probability')
					caxis([cmin cmax]);
					cbh.FontSize = 12;


						

					figure(tm_fig_id)
					subplot(1,2,2)
					imagesc(outgoing_diff)
					title(sprintf('%s Outgoing',current_trial))
					pbaspect([1 1 1])
					colormap(gca,brewermap([],'spectral'))
					% colormap( [0 0 0; brewermap([],'spectral')] )
					cbh = colorbar;
					ylabel(cbh,'Difference in Probability')
					caxis([cmin cmax]);
					cbh.FontSize = 12;

					% Save Plot
					suptitle(sprintf('%s %s  Diff Matrix  LFP', animal_id, current_trial ))
					filename = sprintf('%s_%s_diff_matrix.png', animal_id, current_trial);
					fullname = fullfile(global_output_folder,filename);
					saveas(figure(tm_fig_id),fullname);

					if make_eps_images

						epsname = sprintf('%s_%s_diff_matrix.eps', animal_id, current_trial);
						fullname_eps = fullfile(global_output_folder,epsname);
						set(gcf,'renderer','Painters')
						saveas(gcf,fullname_eps,'epsc')

					end % End if make_eps_images	



					filename = sprintf('%s_resultant_diff.mat', animal_id );
					fullname = fullfile(global_output_folder,filename);
					save(fullname, 'resultant_outgoing_change', 'resultant_incoming_change');

					return
					close(figure(tm_fig_id))

				end % End if plot_difference_matrix

				
			% -------------------------------------------------------------------
			% Overlay eccentricity on state space
			% -------------------------------------------------------------------
				if overlay_eccentricity	
					global_incoming = [];
					global_outgoing = [];

					for trial_iterator = 1:total_trials
					
						filename = sprintf('%s_Eccenctricity_%s.mat', animal_id, char(trial_folders(trial_iterator)) );
						fullname = fullfile(global_output_folder,filename);
						load(fullname)
					
						global_incoming = [global_incoming incoming_eccentricities ];
						global_outgoing = [global_outgoing outgoing_eccentricities ];

					end % End trial_iterator

					% Assign global colors to eccentrcity
					% [s, sort_ind] = sort(global_incoming);
					% in_ecc_color(sort_ind,:) = viridis(length(global_incoming));

					[in_ecc_min, in_ecc_max] = bounds(global_incoming);
					[out_ecc_min, out_ecc_max] = bounds(global_outgoing);

					in_ecc_bounds = [in_ecc_min, in_ecc_max];
					out_ecc_bounds = [out_ecc_min, out_ecc_max];

					[in_ecc_cmap in_ecc_color] = assign_colors(global_incoming, in_ecc_bounds, [], 'viridis');
					[out_ecc_cmap out_ecc_color] = assign_colors(global_outgoing, out_ecc_bounds, [], 'magma');

					% s = []; sort_ind = []; 
					% [s, sort_ind] = sort(global_outgoing);
					% out_ecc_color(sort_ind,:) = magma(length(global_outgoing));


					for trial_iterator = 1:total_trials		

						range_start = global_states_count_lfp(trial_iterator) + 1 ;
						range_end = global_states_count_lfp(trial_iterator + 1) ;

						pv_range = range_start : range_end; 

						x = projection_lfp(pv_range,1);
						y = projection_lfp(pv_range,2);
						temp_inc_colors = in_ecc_color(pv_range,:);
						temp_out_colors = out_ecc_color(pv_range,:);

						figure(185)
						subplot(2,2,trial_iterator)
						scatter(projection_lfp(:,1), projection_lfp(:,2),5,color_gray,'filled')
						hold on 
						scatter(x,y,5, temp_inc_colors,'filled')
						title(char(trial_folders(trial_iterator)))
						colormap(in_ecc_cmap)
						cbh = colorbar;
						pbaspect([1 1 1]) 
						caxis(gca,in_ecc_bounds);


						figure(186)
						subplot(2,2,trial_iterator)
						scatter(projection_lfp(:,1), projection_lfp(:,2),5, color_gray, 'filled')
						hold on 
						scatter(x,y,5, temp_out_colors,'filled')
						title(char(trial_folders(trial_iterator)))
						colormap(out_ecc_cmap)
						cbh = colorbar;
						pbaspect([1 1 1]) 
						caxis(gca, out_ecc_bounds);

						temp_out_colors = []; x = []; y = [];
						temp_inc_colors = [];

					end % End trial_iterator

					figure(185)
					suptitle('Incoming Eccentrcity LFP')
					filename = sprintf('Incoming_Eccenctricity_overlay_lfp_%s.png', animal_id );
					fullname = fullfile(global_output_folder,filename);
					saveas(figure(185),fullname);


					figure(186)
					suptitle('Outgoing Eccentrcity LFP')
					filename = sprintf('Outgoing_Eccenctricity_overlay_lfp_%s.png', animal_id );
					fullname = fullfile(global_output_folder,filename);
					saveas(figure(186),fullname);
				

				end % End overlay_eccentricity
 
			return
	end % End for animal_iterator

toc;
