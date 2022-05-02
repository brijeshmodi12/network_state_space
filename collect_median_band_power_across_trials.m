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

	overlay_lfp_folder = 'overlay_lfp_power_bins';
	overlay_csd_folder = 'overlay_csd_power_bins';

	global_csd_filename = 'global_median_power_csd.mat';
	global_lfp_filename = 'global_median_power_lfp.mat';

	total_folders = length(animal_folders_list);

	lfp_trials_across_animals = {};
	csd_trials_across_animals = {};


	animal_colors = plasma(total_folders);

	color_matrix = [];

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


	output_folder = fullfile(root_directory,'median_band_power_across_trials')
	mkdir(output_folder);

% Core Loop
	for animal_iterator = 1:total_folders

		animal_id = char(animal_folders_list(animal_iterator))
		full_filepath = fullfile(root_directory, animal_id, sub_directory_name1, overlay_lfp_folder, global_lfp_filename);

		try
			load(full_filepath);
		catch
			fprintf('File not found\n');
			return
		end

		total_bands = size(global_median_power_in_trials,1);

		for band_iterator = 1 : total_bands 

			 if animal_iterator <= 1 ;
			 	lfp_trials_across_animals{1,band_iterator} = [];

			 end

			 temp = [];

			 temp = lfp_trials_across_animals{1,band_iterator};

			 row = []; average_row = [];
			 
			 row = global_median_power_in_trials(band_iterator,:);

			 average_row = [mean(row([1,4])) mean(row([2,3])) ];

			 diff1 = diff(average_row) / sum(average_row);

			 average_row = [average_row diff1];

			 temp = [temp; average_row];

			 

			 lfp_trials_across_animals{1,band_iterator} = temp;

		end % End for band_iterator

	
		
		full_filepath = fullfile(root_directory, animal_id, sub_directory_name1, overlay_csd_folder, global_csd_filename);

		try
			load(full_filepath);
		catch
			fprintf('File not found\n');
			return
		end
		

		for band_iterator = 1 : total_bands 

			 if animal_iterator <= 1 ;
			 	csd_trials_across_animals{1,band_iterator} = [];
			 end

			 temp = [];

			 temp = csd_trials_across_animals{1,band_iterator};

			 row = []; average_row = [];
			 
			 row = global_median_power_in_trials(band_iterator,:);

			 average_row = [mean(row([1,4])) mean(row([2,3])) ];

			 temp = [temp; average_row];

			 csd_trials_across_animals{1,band_iterator} = temp;

		end % End for band_iterator
	

	end % End for animal_iterator

	output_global_filename = fullfile(output_folder,'median_power_across_animals.mat');
			save(output_global_filename,'csd_trials_across_animals','lfp_trials_across_animals')

	return


	% -------------------------------------------------------------------
	% Plot 
	% -------------------------------------------------------------------

		total_bands = length(bands_strings);

		for layer_iterator = 1 : total_layers

			current_layer = layer_strings(layer_iterator);

			for band_iterator = 1:total_bands

				band_color = [];

				band_index = (layer_iterator - 1) * total_bands + band_iterator;

				current_plot_matrix = [];
				current_plot_matrix = lfp_trials_across_animals{1,band_index};

				lfp_fig_id = 1000 + layer_iterator;
				figure(lfp_fig_id)
				set(gcf, 'Position', get(0, 'Screensize'));
				subplot(3,2,band_iterator)
				h = plot(current_plot_matrix','lineWidth',3);
				set(h, {'color'}, num2cell(magma(4),2) );
				cbh = colorbar;
				colormap('magma')
				caxis([1 total_folders])
				set(cbh,'YTick',1:total_folders);
				% set(h, 'XTickLabel',{'PreSleep','Arena1','Arena2','PostSleep'});



				ylabel(cbh,'Animals')
				title(sprintf('%0.1f - %0.1f Hz', bands_array(band_iterator,1), bands_array(band_iterator,2)))
				pbaspect([1 1 1]) 
				xlim([0 5])

			

				current_plot_matrix = [];
				current_plot_matrix = csd_trials_across_animals{1,band_index};

				csd_fig_id = 2000 + layer_iterator;
				figure(csd_fig_id)
				set(gcf, 'Position', get(0, 'Screensize'));
				subplot(3,2,band_iterator)
				h = plot(current_plot_matrix','lineWidth',3);
				set(h, {'color'}, num2cell(magma(4),2) );
				cbh = colorbar;
				colormap('magma')
				caxis([1 total_folders])
				set(cbh,'YTick',1:total_folders);

				ylabel(cbh,'Animals')
				title(sprintf('%0.1f - %0.1f Hz', bands_array(band_iterator,1), bands_array(band_iterator,2)))
				pbaspect([1 1 1]) 
				xlim([0 5])
			
				% sleep_data = []; awake_data = [];
				% sleep_data = reshape(current_plot_matrix(:,[1,4]), [], 1 );
				% awake_data = reshape(current_plot_matrix(:,[2,3]), [], 1 );

				% violinplot_matrix = [sleep_data awake_data];
				% vp_names = {'Sleep', 'Awake'}
				% vp_plot_id = 1100 + layer_iterator;
				% figure(vp_plot_id)

				% violinplot(violinplot_matrix, vp_names, 'ViolinColor', num2cell(magma(2),2) )
				% return
				
				% return
			end	% End for band_iterator


			figure(lfp_fig_id);
			suptitle(sprintf('LFP Layer: %s ',  current_layer ))
			filename = sprintf('lfp_%s.png',  current_layer);
			fullname = fullfile(output_folder,filename);
			saveas(lfp_fig_id,fullname);

			close(lfp_fig_id)


			figure(csd_fig_id);

			suptitle(sprintf('CSD Layer: %s ',  current_layer ))
			filename = sprintf('csd_%s.png',  current_layer);
			fullname = fullfile(output_folder,filename);
			saveas(csd_fig_id,fullname);

			close(csd_fig_id)

			% close(csd_overlay_fig_id)

		end % End for layer_iterator

		% fname = sprintf('lfp_corr_mat_across_trials_across_animals.png', animal_id);
		% lfp_filepath = fullfile(output_folder, fname);
		% saveas(gcf, lfp_filepath);

		% Save mean power in each bands across trials
			



	% -------------------------------------------------------------------
	% CSD
	% -------------------------------------------------------------------
