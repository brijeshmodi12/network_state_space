%% compute_difference: Computes distance between all points in a vector and a reference point
% Input: 
% - vector - Each point is a row, columns (dimension of points)
% - reference_point - 1 row, same dimension as vector
function [distance_vector] = compute_difference(vector, reference_point)
	% clc; clear all
	% vector = [-1, 0; 0,1; 2,0; 0,3]
	% ref = [1,0]

	if size(vector,2) ~= size(reference_point, 2)
		fprintf('Reference point and vector dimension mismatch\n')
		return
	end

	% Make a matrix by repeating the reference point.
	ref_point_matrix = repmat(reference_point, size(vector,1), 1);

	squared_diff_between_coordinates = (ref_point_matrix - vector) .^2;

	distance_vector = sqrt(sum(squared_diff_between_coordinates,2));




	
