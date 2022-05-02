%% find_start_end: function description.
% Inputs 
	% - Vector containing ones and zeros.
% Outputs
	% - start and end indices of each continous patterns of ones
function [start_indices, end_indices] = find_start_end(vector)

	if max(vector) >  1
		fprintf('Input vector should only contain 1s and 0s')
		start_indices = 0;
		end_indices = 0;
		return
	end

	% make column vector
	vector = reshape(vector, [length(vector), 1]);

	% Find all ones
	ones_indices = find(vector == 1);

	% Subtract adjacent elements
	ones_indices_diff = diff(ones_indices);


	
	temp_indices = find(ones_indices_diff >  1) + 1 ;

	start_indices = [ ones_indices(1); ones_indices(temp_indices) ] ;

	
	temp_indices = temp_indices - 1;

	end_indices =  [ ones_indices(temp_indices) ; ones_indices(end)] ;


	outputs = [start_indices end_indices];
