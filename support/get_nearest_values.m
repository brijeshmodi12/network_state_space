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


	
