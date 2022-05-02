%% custom_histcounts: function description
% Computes total elements in user defined bin edges for given input vector 
% (bins can be overlapping) 
% Inputs 
% - vector  : row/column vector data
% - edges   : N X 2 matrix with first column as start value and second column as end values of bins;

% Output
% - bin_counts				: N X 1 vector containing total elements in each bins
% - bin_elements_indices	: N X 1 cells containig indices of elements in each bin.

function [bin_counts, bin_elements_indices] = custom_histcounts(vector, edges)

	total_bins = size(edges,1);

	bin_counts = zeros(total_bins,1);

	bin_elements_indices = {};
	
	for bin_iterator = 1 : total_bins

		bin_start = edges(bin_iterator, 1);
		bin_end = edges(bin_iterator, 2);
					
		bin_elements_indices{bin_iterator} = find( vector >= bin_start & vector < bin_end);

		bin_counts(bin_iterator) = numel(bin_elements_indices{bin_iterator});

	end % End bin_iterator






