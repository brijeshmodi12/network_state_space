%% assign_colors: function description
function [colormap_handle, color_values ] = assign_colors(input_values, caxis_range, resolution, colormap_string)
	% Input values are subset of range

	color_gray = [134, 138, 145 ] ./ 255;

	nan_color = color_gray;
	% nan_color = [1 0 0];

	nan_indices = find(isnan(input_values));

	non_nan_indices = find(~isnan(input_values));

	total_input = numel(input_values);

	[min_range max_range] = bounds(caxis_range);

	if isempty(resolution)
		resolution = max((max_range-min_range), 1000);
	end

	colormap_handle = get_colormap(resolution, colormap_string);

	binned_range_values = linspace(min_range, max_range, resolution);

	% Rounds off input values to nearest values in the range at a given resolution
	[~ , nearest_value_index] = get_nearest_values(binned_range_values, input_values(non_nan_indices));

	[nvim ,nvimax ] = bounds(nearest_value_index);

	color_values(non_nan_indices, :) = colormap_handle(nearest_value_index,:);

	if ~isempty(nan_indices)
		% fprintf('NaN values found while assigning colors...\n')
		color_values(nan_indices,:) = repmat(nan_color, length(nan_indices), 1);
	end


function [colormap_handle] = get_colormap(resolution, colormap_string)

	switch colormap_string

		case 'cividis'
			colormap_handle = cividis(resolution);

		case 'viridis'
			colormap_handle = viridis(resolution);	

		case 'magma'
			colormap_handle = magma(resolution);

		case 'plasma'
			colormap_handle = plasma(resolution);

		case 'inferno'
			colormap_handle = inferno(resolution);

		case 'twilight'
			colormap_handle = twilight(resolution);		

		case 'pink'
			colormap_handle = pink(resolution);		

		case 'bone'
			colormap_handle = bone(resolution);		

		case 'winter'
			colormap_handle = winter(resolution);		

		case 'summer'
			colormap_handle = summer(resolution);
			
		otherwise
			colormap_handle = brewermap(resolution, colormap_string);

	end



% Rounds up values to preset levels in vector
% kinda like analog to digital conversion
function [nearest_value_vector, index_vector] = get_nearest_values(ref_vector, values_vector)
	% ref_vector = 1:100;
	% values_vector = [5 60];
	total_vals = length(values_vector);

	% Initialize vector to store indices of nearest values
	index_vector = zeros(size(values_vector));

	for val_iterator = 1 : total_vals

		current_value = values_vector(val_iterator);

		[~,  ind] = min(abs(ref_vector -  current_value));

		index_vector(val_iterator) = ind;

		nearest_value_vector(val_iterator) = ref_vector(ind);

	end % End for val_iterator