%% get_accelerometer_variance: function description
% trial_acc_variance :: binary vector (for sleep/awake)
function [trial_acc_variance] = get_accelerometer_variance(trial_binned_acceleration, threshold)

	trial_acc_variance = movvar(trial_binned_acceleration, 10);

	trial_acc_variance(trial_acc_variance > threshold) = 1;
	trial_acc_variance(trial_acc_variance < threshold) = 0;
