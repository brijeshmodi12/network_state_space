

%% make_transition_matrix: function description
function [transition_matrix] = make_transition_matrix(trajectory_start_bin_id, trajectory_end_bin_id, edges, total_bins)

	% initialize output 
	transition_matrix = [];

	if numel(trajectory_start_bin_id) ~=  numel(trajectory_end_bin_id)
		fprintf('Input Mismatch\n')
		return
	end

	bins_list = unique(trajectory_start_bin_id);

	% total_bins = length(bins_list);

	for bin_iterator = 1:total_bins

		bin_iterator;

		bin_states_indices = find(trajectory_start_bin_id == bin_iterator);

		% Consider with with atleast 50 states (10 seconds of data)
		if numel(bin_states_indices) < 50
			transition_matrix(bin_iterator,:) = zeros(1,total_bins);
			continue
		end

		outgoing_bin_id = trajectory_end_bin_id(bin_states_indices);

		outgoing_transition_prob = custom_histcounts(outgoing_bin_id, edges);

		outgoing_transition_prob = outgoing_transition_prob / sum(outgoing_transition_prob);

		transition_matrix(bin_iterator,:) = outgoing_transition_prob;

	end

