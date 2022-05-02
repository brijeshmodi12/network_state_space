%% create_pv_spiketrain: Create population vector given with given list of cells using bin_edges;
% Inputs;
%  cells_list		: Cell array containing spike timings of each spike 
%  bin_edges		: Bin edges to be used for binning spike trains. (can be overlapping)
%  trial_duration	: Trial duration in seconds 
%  firing_threshold	: Firing threshold to skip high firing neurons

% Ouput vector
% population_vector : matrix containing cells in columns, data points in rows.


function [population_vector] = create_pv_spiketrain_raw(cells_list, bin_edges)

	total_cells = numel(cells_list);

	% Initalize 
	population_vector = [];

	% Iterate over each cells to make population vector
	for cell_iterator = 1:total_cells

		current_cell_spiketimes = cells_list{cell_iterator};
			

		% Initialize
		spiketrain_binned = [];

		[spiketrain_binned, ~] = custom_histcounts(current_cell_spiketimes, bin_edges);

		% Convert more than one spikes in a bin to 1..
		% spiketrain_binned(spiketrain_binned >= 1 ) = 1;	

		% spiketrain_binned = smooth_gaussian(spiketrain_binned, smoothing_std, smoothing_width);

		% Make sure the spiketrains are column vector
		spiketrain_binned =  reshape(spiketrain_binned,[length(spiketrain_binned), 1 ]);
		
		population_vector = [population_vector   spiketrain_binned];

		
		
	end % End for cell_iterator

