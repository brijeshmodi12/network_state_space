%% linkage_clustering: function description
function [sorted_indices] = linkage_clustering(input_matrix,  no_of_groups)
				
	%# Remove diagonal elements
	input_matrix = input_matrix - eye(size(input_matrix));

	%# and convert to a vector (as pdist)
	dissimilarity = 1 - input_matrix';

	%# perform complete linkage clustering
	Z = linkage(dissimilarity,'complete');
	
	%# group the data into clusters
	groups = cluster(Z,'maxclust', no_of_groups);

	[~, sorted_indices] = sort(groups);

	

