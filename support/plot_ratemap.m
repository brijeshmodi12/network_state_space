%% plot_ratemap: Plots rate map using scaled color, and white for unvisited bins

function plot_ratemap(frequency_map, duration_map, title_string, colormap_string)

	% Test arguments
		% frequency_map = magic(5);

		% smoothed_duration = magic(5);

		% smoothed_duration(1,1) = 0;

		% duration_map = smoothed_duration;

	% Handle astronomical values in frequencies 	
		frequency_map(frequency_map > 1000 ) = 0;

	% Get original dimension of frequency map
		[nrows, ncols] = size(frequency_map);

	% No of color levels between 0 and maximum firing rate.
		total_levels = 30;

	% Get color RGB values for colors
		[color_levels, brewermap_flag] = get_color(colormap_string, total_levels);

	% Color for unvisited bins --  White
		color_unvisited = [1 1 1];

	% Create a copy 
		rate_map = frequency_map;
		temp_duration_map = duration_map(:);

	% Convert matrix into a column vector for easier operations
		rate_map = rate_map(:);

	% Maximum firing rate
		max_rate = max(rate_map);

	% Normalize rate from 0-1;
		norm_rate_map = rate_map / max_rate;

	% Assign a 'level' to normalized rate map
	% We then use this level to assign a color for this corresponding firing rate
		level_map = round(total_levels *  norm_rate_map) + 1;
		level_map(level_map > total_levels) = total_levels;

	% Assign colors
		color_map = color_levels(level_map,:);

	% Assign white colors to unvisited bins
				
	% Find unvisited bins
		unvisited_indices = find(temp_duration_map == 0);

	% Overwrite the values	
		color_map(unvisited_indices,:) = repmat(color_unvisited, [length(unvisited_indices), 1]);

		% Reshape each column of color_map to 2D matrix form.
		% Size of that matrix corresponds to size of original frequency map
		for rgb_iterator = 1:3

			colored_matrix(:,:,rgb_iterator) = reshape(color_map(:,rgb_iterator) , nrows, ncols);

		end % End	

		image(colored_matrix);
		if ~brewermap_flag
			colormap(colormap_string)
		else
			colormap(brewermap(colormap_string))
		end
		colorbar
	   	caxis([min(norm_rate_map) max(norm_rate_map)]);
   	 	caxis([0 max_rate]);
	   	% axis off
		set(gca,'YDir','reverse');
		axis image;
		title(title_string);




	
	
	% return




	% outputs = ;
%% get_color: function description
function [color_levels, brewermap_flag] = get_color(colormap_string, n)
	brewermap_flag = 0;

	switch colormap_string
		case 'jet'
			color_levels = jet(n);

		case 'viridis'
			color_levels = viridis(n);

		case 'inferno'
			color_levels = inferno(n);	

		case 'parula'
			color_levels = parula(n);	

		case 'magma'
			color_levels = magma(n);

		otherwise
			color_levels = brewermap(n, colormap_string);
			brewermap_flag = 1;
	end

