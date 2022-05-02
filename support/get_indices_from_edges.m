%% get_indices_from_edges: function description
% Returns all the indices given the start and end points.
%
% Input : 
% -	edges n x 2 matrix containing starting and end point of each edge
%
% Output
% Row vector containing all the numbers between the edges
%
% Example Input : [1  11 ;  15 20]
% Example Output : [1:11 15:20]
function [indices] = get_indices_from_edges(edges)

	

	indices = [];
	total_edges = size(edges,1);

	for edge_iterator = 1:total_edges
		
		indices = [indices edges(edge_iterator,1) : edges(edge_iterator,2)];

	end % End edge_iterator

	% indices = unique(indices);


