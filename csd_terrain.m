% ===================================================================
	% Compute and Plot the terrain of state space CSD
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
		root_directory = 'C:\Users\Matteo\Desktop\Brijesh\Open Field Dataset\';
		root_directory = 'D:\matteo_datasets\Open Field Datasets\';

		animal_folders_list = ["AH2","AH3","AH4","AH5"];

		% csd_root = 'csd_manifold_awake_sleep_pca_global_rescale';
		csd_root = 'network_state_space_bin200';

		output_foldername = 'csd_terrain'; 

		file_to_get = 'state_space_data.mat';

		total_folders = length(animal_folders_list);

		color_gray = [134, 138, 145 ] ./ 255;

		% Combine all the bands in an single array so you can apply loops 
		bands_strings = ["Delta", "Spindle", "Theta", "Slow gamma", "Medium gamma", "Fast gamma"];

		bin_size_sec = 0.2;

		% 1 for processing
		compute_eccentricity_and_plot = 1;

		overlay_eccentricity =  1;
	

		compute_transition_matrix = 1;

		plot_transition_probability = 1;	


		% For Trajectories
		no_xbins = 2;
		no_ybins = 2;	

		next_states_jump = 5;
		

% -------------------------------------------------------------------
% Core Loop
% -------------------------------------------------------------------
	% Core Loop
	for animal_iterator = 1:total_folders


		animal_id = char(animal_folders_list(animal_iterator))
		full_filepath = fullfile(root_directory, animal_id, csd_root, file_to_get);

		try
			load(full_filepath);
		catch
			fprintf('File not found\n');
			return
		end

		
		% -------------------------------------------------------------------
		% Animal specific Constants
		% -------------------------------------------------------------------		
			total_trials = length(global_states_count_csd) - 1;

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

					% csd_filename = '135_CH42.continuous'; % AH2 

					


				case 'AH3'
					trial_folders = ["preSleep", "1Square", "2Circle", "postSleep"]; % AH3 
					output_folder_string = 	'network_state_space_bin200';
					
					% trial_folders = [ "1Square", "2Circle" ]; % AH3
					% output_folder_string = 	'all_trials_power_manifold_awake';
					% only_awake = 1;

					% trial_folders = ["preSleep", "postSleep"]; % AH3  
					% output_folder_string = 	'all_trials_power_manifold_sleep';
					% only_sleep = 1;

					% csd_filename = '159_CH46.continuous'; % AH3 

					
					

				case 'AH4'	
					trial_folders = ["preSleep", "1Circle", "2Square", "postSleep"]; % AH4 
					output_folder_string = 'network_state_space_bin200';

					% trial_folders = [ "1Circle", "2Square" ]; % AH4 
					% output_folder_string = 'all_trials_power_manifold_awake';
					% only_awake = 1;

					% trial_folders = ["preSleep", "postSleep"]; % AH4 
					% output_folder_string = 	'all_trials_power_manifold_sleep';
					% only_sleep = 1;

					% csd_filename = '126_CH46.continuous'; % AH4

					


				case 'AH5'	
					trial_folders = ["preSleep", "1Circle", "2Square", "postSleep"]; % AH5 
					output_folder_string = 'network_state_space_bin200';

					% trial_folders = ["1Circle", "2Square"]; % AH5 
					% output_folder_string = 'all_trials_power_manifold_awake';
					% only_awake = 1;

					% trial_folders = ["preSleep", "postSleep"]; % AH5
					% output_folder_string = 	'all_trials_power_manifold_sleep';
					% only_sleep = 1;
					
					% csd_filename = '101_CH42.continuous'; % AH5

					csd_cut_off_vector
					

				otherwise
					fprintf('Incorrect animal id\n')
					return
			end % End switch animal_id

			bin_number = no_xbins
			bin_folder_name = strcat('bin_', num2str(bin_number));

			output_root_folder = fullfile(root_directory,animal_id,csd_root,output_foldername);
			mkdir(output_root_folder);

			bin_output_folder = fullfile(root_directory,animal_id,csd_root,output_foldername, bin_folder_name);

			global_output_folder = bin_output_folder;
			mkdir(global_output_folder);

		 

		
			trial_colors = magma(total_trials+5);
			trial_colors = trial_colors(end-3:end,:);

			


		% -------------------------------------------------------------------
		% Process each trial
		% -------------------------------------------------------------------	
			for trial_iterator = 1:total_trials

				current_trial = char(trial_folders(trial_iterator));

				% -------------------------------------------------------------------
				% Trial Specific Constants
				% -------------------------------------------------------------------	
					n_components = 2;

					range_start = global_states_count_csd(trial_iterator) + 1 ;
					range_end = global_states_count_csd(trial_iterator + 1) ;

					pv_range = range_start : range_end; 
					x = [];  y = []; z = [];

					x = projection_csd(pv_range,1);
					y = projection_csd(pv_range,2);
					% if n_components > 2
					% 	z = projection_csd(pv_range,3);
					% end

					% % umap output for this trial
					% Y = [x, y];

					population_vector = [];
					population_vector = global_csd_states(pv_range,:);

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
						x_edges = linspace(min(projection_csd(:,1)), max(projection_csd(:,1)), no_xbins+1);
						y_edges = linspace(min(projection_csd(:,2)), max(projection_csd(:,2)), no_ybins+1);


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
									subplot(2,2,1)
									set(gcf, 'Position', get(0, 'Screensize'));
									scatter(projection_csd(:,1),projection_csd(:,2), 20, color_gray, 'filled')
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
									[xmin, xmax] = bounds(x(bsi)); 
									[ymin, ymax] = bounds(y(bsi)); 

									width_rectangle = abs(xmax - xmin);
									height_rectangle = abs(ymax - ymin);
									rectangle('Position',[xmin, ymin, width_rectangle, height_rectangle],'EdgeColor','k','LineWidth',3)
									hold on

									scatter(median(x(bsi)), median(y(bsi)), 50, 'k', 'filled');

									title(sprintf('%s Bin: %d Incoming', char(trial_folders(trial_iterator)), bin_id))
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

							

									% Compute flow for incoming
										bin_center = [median(x(bsi)), median(y(bsi)) ];

										x_shifted = x - bin_center(1);
										y_shifted = y - bin_center(2);

										[theta, rho] = cart2pol(x_shifted, y_shifted);

										nbins_theta = 30;
																		
										figure(trajectory_fig_id)
										subplot(2,2,3)
										incoming_probability = polarhistogram(theta(bin_incoming_indices), nbins_theta, 'Normalization', 'probability','FaceColor','green','FaceAlpha',.3);
										hold on

										prob_values = incoming_probability.Values;

										theta_bins = incoming_probability.BinEdges;
										theta_bins = theta_bins(1:end-1) + (theta_bins(2)-theta_bins(1))/2;

										[max_val, max_ind] = max(prob_values);

										max_ind_opposite = mod(max_ind + fix(nbins_theta/2) ,nbins_theta) + 1;

										axis_1_ind = [max_ind  max_ind_opposite];
										axis_2_ind = [mod(max_ind + fix(nbins_theta/4)  ,nbins_theta) + 1 mod(max_ind + fix(3*nbins_theta/4), nbins_theta) + 1];

										axis1_length = sum(abs(prob_values(axis_1_ind)));
										axis2_length = sum(abs(prob_values(axis_2_ind)));


										if axis1_length > axis2_length
											majoraxis = axis1_length;
											minoraxis = axis2_length;
										else
											majoraxis = axis2_length;
											minoraxis = axis1_length;
										end

										if minoraxis == 0
											minoraxis = 0.1 * majoraxis;
										end



										% Eccentricity of ellipse
										ecc = axes2ecc(majoraxis,minoraxis);
										title(sprintf('Transition Probability; e : %0.3f', ecc ))
										hold off

										incoming_eccentricities(bin_states_indices) = ecc;

										ecc = [];
										majoraxis = [];  minoraxis = [];

													

								% ----------------------------------------	
								% Outgoing trajectories
								% ----------------------------------------
									figure(trajectory_fig_id)
									subplot(2,2,2)
									scatter(projection_csd(:,1),projection_csd(:,2), 20, color_gray, 'filled')
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

									[xmin, xmax] = bounds(x(bsi)); 
									[ymin, ymax] = bounds(y(bsi)); 

									width_rectangle = abs(xmax - xmin);
									height_rectangle = abs(ymax - ymin);
									rectangle('Position',[xmin, ymin, width_rectangle, height_rectangle],'EdgeColor','k','LineWidth',3)
									hold on 
									scatter(median(x(bsi)), median(y(bsi)), 50, 'k', 'filled');

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

									suptitle(sprintf('%s %s Bin: %d',animal_id, char(trial_folders(trial_iterator)), bin_id))


									figure(trajectory_fig_id)
									subplot(2,2,4)
									outgoing_probability = polarhistogram(theta(bin_outgoing_indices), nbins_theta, 'Normalization', 'probability', 'FaceColor','red','FaceAlpha',.3);
									hold on

									prob_values = outgoing_probability.Values;

									theta_bins = outgoing_probability.BinEdges;
									theta_bins = theta_bins(1:end-1) + (theta_bins(2)-theta_bins(1))/2;

									[max_val, max_ind] = max(prob_values);

									max_ind_opposite = mod(max_ind + fix(nbins_theta/2) ,nbins_theta) + 1;

									axis_1_ind = [max_ind  max_ind_opposite];
									axis_2_ind = [mod(max_ind + fix(nbins_theta/4)  ,nbins_theta) + 1 mod(max_ind + fix(3*nbins_theta/4), nbins_theta) + 1];


									axis1_length = sum(abs(prob_values(axis_1_ind)));
									axis2_length = sum(abs(prob_values(axis_2_ind)));

									theta_bins = outgoing_probability.BinEdges;
									theta_bins = theta_bins(1:end-1) + (theta_bins(2)-theta_bins(1))/2;

								
									axis1_length = sum(abs(outgoing_probability.Values(axis_1_ind)));
									axis2_length = sum(abs(outgoing_probability.Values(axis_2_ind)));


									if axis1_length > axis2_length
										majoraxis = axis1_length;
										minoraxis = axis2_length;
									else
										majoraxis = axis2_length;
										minoraxis = axis1_length;
									end

									if minoraxis == 0
										minoraxis = 0.1 * majoraxis;
									end


									ecc = axes2ecc(majoraxis,minoraxis);
									title(sprintf('Transition Probability; e : %0.3f', ecc ))
									hold off

									outgoing_eccentricities(bin_states_indices) = ecc;

									ecc = [];
									majoraxis = [];  minoraxis = [];
									incoming_probability = [];
									outgoing_probability = [];


									filename = sprintf('%s_%s_bin_%d.png', animal_id, char(trial_folders(trial_iterator)), bin_id );
									fullname = fullfile(global_output_folder,filename);
									saveas(figure(trajectory_fig_id),fullname);
									
									close all;		

							end % End for ybin_iterator

							% return

						end % End for xbin_iterator	


						filename = sprintf('%s_Eccenctricity_%s.mat', animal_id, char(trial_folders(trial_iterator)) );
						fullname = fullfile(global_output_folder,filename);
						save(fullname, 'incoming_eccentricities', 'outgoing_eccentricities');

					end % End if compute_eccentricity_and_plot



				% -------------------------------------------------------------------
				% Compute Transition Matrix and Plot Transition Probabilities
				% -------------------------------------------------------------------
					if compute_transition_matrix

						incoming_transition_matrix = [];
						outgoing_transition_matrix = [];

						fprintf('Computing Transition Matrix\n')

						% X and Y axis reference points for umap
						x_edges = linspace(min(projection_csd(:,1)), max(projection_csd(:,1)), no_xbins+1);
						y_edges = linspace(min(projection_csd(:,2)), max(projection_csd(:,2)), no_ybins+1);


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

						hist_edges_start = 1:total_histbins;
						hist_edges_start = hist_edges_start - 0.5;
						
						hist_edges_end = hist_edges_start + 1;

						edges = [hist_edges_start' hist_edges_end'];


						% Compute transition probability for each bin
							for bin_id = 1:total_histbins

								current_bsi = bin_bsi_array{1,bin_id};

								if isempty(current_bsi)
									incoming_transition_matrix = [incoming_transition_matrix ; zeros(1,total_histbins)];
									outgoing_transition_matrix = [outgoing_transition_matrix ; zeros(1,total_histbins)];
									continue;
								end

								total_bsi = length(current_bsi);

								incoming_bsi = []; incoming_bin_id = [];
								outgoing_bsi = []; outgoing_bin_id = [];
								incoming_bin_colors = []; outgoing_bin_colors = [];
								incoming_states_colors = []; outgoing_states_colors = []; 

								total_trajectories = total_bsi;
								trajectory_start_incoming = current_bsi - next_states_jump;
								trajectory_end_outgoing = current_bsi + next_states_jump;


								for trajectory_iterator = 1 : total_trajectories
									ti = [];trajectory_indices = []; 
									trajectory_indices = trajectory_start_incoming(trajectory_iterator) : current_bsi(trajectory_iterator);
									ti = trajectory_indices;

									ti(ti<1) = 1;
									ti(ti>total_pop_vec) = total_pop_vec;
									incoming_bsi = [incoming_bsi ti];
									

									ti = []; trajectory_indices = [];
									trajectory_indices = current_bsi(trajectory_iterator): trajectory_end_outgoing(trajectory_iterator) ;
									ti = trajectory_indices;

									ti(ti<1) = 1;
									ti(ti>total_pop_vec) = total_pop_vec;
									outgoing_bsi = [outgoing_bsi ti];
									
								end % End for trajectory_iterator

								incoming_bsi(incoming_bsi<1) = 1;
								outgoing_bsi(outgoing_bsi > total_pop_vec) = total_pop_vec;

								% Id of incoming bins
								incoming_bin_id = bin_id_of_states(incoming_bsi);
								% Id of outgoing bins
								outgoing_bin_id = bin_id_of_states(outgoing_bsi);				
								
								% Compute probability of transition to all bins
								incoming_transition_prob = custom_histcounts(incoming_bin_id, edges);
								incoming_transition_prob = incoming_transition_prob / sum(incoming_transition_prob);

								outgoing_transition_prob = custom_histcounts(outgoing_bin_id, edges);
								outgoing_transition_prob = outgoing_transition_prob / sum(outgoing_transition_prob);

								incoming_transition_matrix = [incoming_transition_matrix ; incoming_transition_prob'];
								outgoing_transition_matrix = [outgoing_transition_matrix ; outgoing_transition_prob'];

								% Replace zeros with NaN for coloring
									zero_itp_ind = find(incoming_transition_prob == 0);
									zero_otp_ind = find(outgoing_transition_prob == 0);

									incoming_transition_prob(zero_itp_ind) = NaN;
									outgoing_transition_prob(zero_itp_ind) = NaN;
							
								% Assign Colors -- Incoming
									[incoming_colormap, incoming_bin_colors] = assign_colors(incoming_transition_prob, [0 1], 100, 'viridis');
									
									for bin_iterator = 1:total_histbins

										temp_states_ind = [];
										temp_states_ind = bin_bsi_array{1,bin_iterator};
										if isempty(temp_states_ind)
											continue;
										end

										incoming_states_colors(temp_states_ind,:) = repmat(incoming_bin_colors(bin_iterator,:) , length(temp_states_ind),1);

									end % End for bin_iterator

									
								% Plot Incoming
									if plot_transition_probability
										tp_fig_id = 821;
										figure(tp_fig_id)
										set(gcf, 'Position', get(0, 'Screensize'));
										subplot(1,2,1)
										scatter(projection_csd(:,1), projection_csd(:,2), 10, color_gray, 'filled')
										hold on
										scatter(x,y, 10, incoming_states_colors, 'filled')
										pbaspect([1 1 1])
										% colormap(brewermap('*ylgnbu'))
										colormap(gca, incoming_colormap)
										cbh = colorbar;
										ylabel(cbh,'Probability')
										% [cmin, cmax] = bounds(incoming_transition_prob);
										% set(cbh, 'ylim', [cmin cmax])
										caxis([0 1]);					
									
										% Plot square
										[xmin, xmax] = bounds(x(current_bsi)); 
										[ymin, ymax] = bounds(y(current_bsi)); 

										width_rectangle = abs(xmax - xmin);
										height_rectangle = abs(ymax - ymin);
										rectangle('Position',[xmin, ymin, width_rectangle, height_rectangle],'EdgeColor','k','LineWidth',3)

										title(sprintf('%s Bin: %d Incoming', char(trial_folders(trial_iterator)), bin_id))
										xticks(x_edges)
										yticks(y_edges)
										grid on
										ax = gca;
										ax.GridLineStyle = '-';
										ax.GridColor = 'k';
										ax.GridAlpha = 1; % maximum line opacity
										set (gca, 'xticklabel' , {[]});
										set (gca, 'yticklabel' , {[]});

										
									end % End if plot_transition_probability


								% Assign Colors -- Outgoing
									[outgoing_colormap, outgoing_bin_colors] = assign_colors(outgoing_transition_prob, [0 1], 100, 'magma');

									for bin_iterator = 1:total_histbins

										temp_states_ind = [];
										temp_states_ind = bin_bsi_array{1,bin_iterator};
										if isempty(temp_states_ind)
											continue;
										end

										outgoing_states_colors(temp_states_ind,:) = repmat(outgoing_bin_colors(bin_iterator,:) , length(temp_states_ind),1);

									end % End for bin_iterator


								% Plot Outgoing
									if plot_transition_probability
										tp_fig_id = 821;
										figure(tp_fig_id)
										subplot(1,2,2)
										scatter(projection_csd(:,1), projection_csd(:,2), 10, color_gray, 'filled')
										hold on
										scatter(x,y, 10, outgoing_states_colors, 'filled')
										pbaspect([1 1 1])
										% colormap(brewermap('*ylgnbu'))
										colormap(gca, outgoing_colormap)
										cbh = colorbar;
										ylabel(cbh,'Probability')
										% [cmin, cmax] = bounds(outgoing_transition_prob);
										% set(cbh, 'ylim', [cmin cmax])
										caxis([0 1]);
													
										% Plot square
											[xmin, xmax] = bounds(x(current_bsi)); 
											[ymin, ymax] = bounds(y(current_bsi)); 

											width_rectangle = abs(xmax - xmin);
											height_rectangle = abs(ymax - ymin);
											rectangle('Position',[xmin, ymin, width_rectangle, height_rectangle],'EdgeColor','k','LineWidth',3)

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

									end % End if plot_transition_probability


								% Save Plot
									suptitle(sprintf('%s %s Bin %d Transition Probability CSD', animal_id, current_trial, bin_id ))
									filename = sprintf('%s_%s_bin_%d_transition_map.png', animal_id, current_trial ,bin_id);
									fullname = fullfile(global_output_folder,filename);
									saveas(figure(tp_fig_id),fullname);

									close(figure(tp_fig_id))
							
									
							end % End for bin_id 


						% incoming_transition_matrix(incoming_transition_matrix > 0.5) = 0;
					


						% Plot Incoming Transition Matrix
							tm_fig_id = 123;
							figure(tm_fig_id)
							set(gcf, 'Position', get(0, 'Screensize'));
							subplot(1,2,1)
							imagesc(incoming_transition_matrix)
							title(sprintf('%s Incoming',current_trial))
							pbaspect([1 1 1])
							colormap(gca,viridis)
							cbh = colorbar;
							ylabel(cbh,'Probability')
							[cmin, cmax] = bounds(incoming_transition_matrix(:));
							caxis([0 1]);
							cbh.FontSize = 12;



						% Plot Outgoing Transition Matrix
							tm_fig_id = 123;
							figure(tm_fig_id)
							subplot(1,2,2)
							imagesc(outgoing_transition_matrix)
							title(sprintf('%s Outgoing',current_trial))
							pbaspect([1 1 1])
							colormap(gca,magma)
							cbh = colorbar;
							ylabel(cbh,'Probability')
							caxis([0 1]);
							cbh.FontSize = 12;

							% Save Plot
							suptitle(sprintf('%s %s Transition Matrix CSD', animal_id, current_trial ))
							filename = sprintf('%s_%s_transition_matrix.png', animal_id, current_trial);
							fullname = fullfile(global_output_folder,filename);
							saveas(figure(tm_fig_id),fullname);
							
							close(figure(tm_fig_id))
							

					end % End if compute_transition_matrix

					
			end % End for trial_iterator	


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

						range_start = global_states_count_csd(trial_iterator) + 1 ;
						range_end = global_states_count_csd(trial_iterator + 1) ;

						pv_range = range_start : range_end; 

						x = projection_csd(pv_range,1);
						y = projection_csd(pv_range,2);
						temp_inc_colors = in_ecc_color(pv_range,:);
						temp_out_colors = out_ecc_color(pv_range,:);

						figure(185)
						subplot(2,2,trial_iterator)
						scatter(projection_csd(:,1), projection_csd(:,2),5,color_gray,'filled')
						hold on 
						scatter(x,y,5, temp_inc_colors,'filled')
						title(char(trial_folders(trial_iterator)))
						colormap(in_ecc_cmap)
						cbh = colorbar;
						pbaspect([1 1 1]) 
						caxis(gca,in_ecc_bounds);


						figure(186)
						subplot(2,2,trial_iterator)
						scatter(projection_csd(:,1), projection_csd(:,2),5, color_gray, 'filled')
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
					suptitle('Incoming Eccentrcity CSD')
					filename = sprintf('Incoming_Eccenctricity_overlay_csd_%s.png', animal_id );
					fullname = fullfile(global_output_folder,filename);
					saveas(figure(185),fullname);


					figure(186)
					suptitle('Outgoing Eccentrcity CSD')
					filename = sprintf('Outgoing_Eccenctricity_overlay_csd_%s.png', animal_id );
					fullname = fullfile(global_output_folder,filename);
					saveas(figure(186),fullname);
				

				end % End overlay_eccentricity
 
			return
	end % End for animal_iterator

toc;
