%% zscore_mod: function description
% Computes Z score with custom std value
% Used for normalizing lfp bands to its origin source LFP
function [output_signal] = zscore_mod(input_signal, std_val)

	if std_val == 0
		output_signal = zscore(input_signal);
		return
	end

	output_signal = input_signal - nanmean(input_signal);
	output_signal = input_signal ./ std_val;
end


