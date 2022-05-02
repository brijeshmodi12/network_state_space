	clear all ; close all; clc; tic

% Root Path - 
	root_directory = 'C:\Users\Matteo\Desktop\Brijesh\animal_4\KnierimDay1\';

% Folders to be considered for analysis
	trial_folders = ["Hab1", "Hab2","Hab3","Hab4","Hab5", "Hab6"];


	total_trials = length(trial_folders);

	% Process each trial tracking to extract time spent in each spatial bin
	for trial_iterator = 1:total_trials

		current_trial_path = fullfile(root_directory, char(trial_folders(trial_iterator)));

		% % Load trajectory and object information
		fprintf('Loading trial %d data... \n', trial_iterator)
		xlfile = xlsread(fullfile(current_trial_path, 'Hab.csv'));

		upsampledata(xlfile, current_trial_path)

		return

	end % End trial_iterator


	toc;