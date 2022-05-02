% Upsamples the data from 30 Hz to 30000 Hz
% Inputs
	% trajectory_obj - trajectory matrix N x 5 
	% path - path to save upsampled file.


function upsampledata(trajectory_obj, path)
	
	upsampling_factor = 1000;
	trajectory_obj1 = [];


	% Apply Interpolation on position
	trajectory_obj1(:,1) = interp(trajectory_obj(:,1), upsampling_factor);
	trajectory_obj1(:,2) = interp(trajectory_obj(:,2), upsampling_factor);

	% Upsample Occupancy, if present 
	if size(trajectory_obj, 2)  > 2 

		trajectory_obj1(:,3) = interp(trajectory_obj(:,3), upsampling_factor);
		trajectory_obj1(:,4) = interp(trajectory_obj(:,4), upsampling_factor);
		trajectory_obj1(:,5) = interp(trajectory_obj(:,5), upsampling_factor);

		% remap
		for temp_iterator = 3:5
			temp_column = trajectory_obj1(:,temp_iterator);
			temp_column(temp_column < 0.5) = 0;
			temp_column(temp_column  > 0.5) = 1;

			trajectory_obj1(:,temp_iterator) = temp_column;	
		end % temp_iterator

		
	end % End if 

	filepath = fullfile(path,'trajectory_obj_upsample.mat');

	trajectory_obj = trajectory_obj1;

	save(filepath, 'trajectory_obj');
	fprintf('File Upsampled\n')

	outputs = 1;


