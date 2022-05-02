%% get_nearest_timestamp: Returns the index of nearest value in the vector for a given val (scalar or vector)
% Input;
	% timestamps_vector  - possible timestamps 
	% value 
function [index_vector] = get_nearest_bin_center(vector, values_vector)

	total_vals = length(values_vector);

	% Initialize vector to store indices of nearest values
	index_vector = zeros(size(values_vector));

	for val_iterator = 1 : total_vals

		current_value = values_vector(val_iterator);

		[~,  ind] = min(abs(vector -  current_value));

		index_vector(val_iterator) = ind;

	end % End for val_iterator


	
