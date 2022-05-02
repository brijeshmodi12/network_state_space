%% custom_binning_mean: function description
% Computes mean of all elements in user defined bin edges for given input vector 
% (bins can be overlapping) 
% Inputs 
% - base_vector  	: timestamps for quantity vector
% - edges   		: N X 2 matrix with first column as start value and second column as end values of bins;
% - quantity_vector	: Quantity vector to take mean	

% Output
% - bin_mean				: N X 1 vector containing total elements in each bins
% - bin_elements_indices	: N X 1 cells containig indices of elements in each bin.

function [bin_means, bin_elements_indices] = custom_binning_mean(base_vector, base_edges, quantity_vector)

	total_bins = size(base_edges,1);

	bin_means = zeros(total_bins,1);

	bin_elements_indices = {};
	
	for bin_iterator = 1 : total_bins

		bin_start = base_edges(bin_iterator, 1);
		bin_end = base_edges(bin_iterator, 2);
					
		bin_elements_indices{bin_iterator} = find( base_vector >= bin_start & base_vector <= bin_end);

		bin_means(bin_iterator) = mean(quantity_vector(bin_elements_indices{bin_iterator}));

	end % End bin_iterator



