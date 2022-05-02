%% detect_swr: function description
function [swr_event_indices, swr_presence_vector] = detect_swr(lfp_swr_band, lfp_pyramidal, plot_swr)

		downsampling_frequency = 1000;
		% plot_swr = 0;
	% Parameters For Detecting Ripples
	% (Adapted from Jadhav et al, 2012)
	% https://science.sciencemag.org/content/sci/suppl/2012/05/02/science.1217230.DC1/Jadhav.SM.pdf
		% Standard deviation of gaussian kernel(in milliseconds)
		std_gaussian_swr = 4;

		% in sample numbers
		std_gaussian_swr = (std_gaussian_swr / 1000) * downsampling_frequency;

		gaussian_width = std_gaussian_swr * 6;


		% Detection threshold (std.)	
		swr_detection_threshold = 4;

		% Event duration threshold 
		event_duration_threshold_ms = 30;

		event_duration_threshold = event_duration_threshold_ms / 1000;


		% SWR event padding 
		event_padding_ms = 100;

		event_padding_sec = event_padding_ms / 1000;

		event_padding_samples = event_padding_sec * downsampling_frequency;	

		total_lfp_samples = length(lfp_swr_band);

		lfp_pyramidal = zscore(lfp_pyramidal);

	% Core	
		swr_indices = [];

		x = reshape(lfp_swr_band, [], 1);

		swr_presence_vector = zeros(numel(x),1);

		
		swr_envelope = abs(hilbert(x)) .^ 2;

		% Smooth the envelope with gaussian kernel and zscore it.
		smoothed_envelope = zscore(smooth_gaussian(swr_envelope, std_gaussian_swr, gaussian_width));


		% Binary vector containing 1 where envelope goes above threshold
		swr_indices = smoothed_envelope > swr_detection_threshold;

		sum(swr_indices)

		% Find start and end of events
		[swr_start_indices, swr_end_indices] = find_start_end(swr_indices);

		event_duration = [swr_end_indices - swr_start_indices] / downsampling_frequency; 

		% Find valid events
		valid_events = find(event_duration  > event_duration_threshold);

		swr_start_indices = swr_start_indices(valid_events);
		swr_end_indices = swr_end_indices(valid_events);

		total_swr_events = length(swr_start_indices);

		if total_swr_events == 0
			fprintf('No sharp wave ripples found \n');
		else
			swr_event_indices = [swr_start_indices swr_end_indices];
		end % 

		swr_presence_indices = get_indices_from_edges(swr_event_indices);

		% Binary vector indicating presence of SWR in LFP
		swr_presence_vector(swr_presence_indices) = 1;


		

		if plot_swr
			% Collect indices of each swr event.
			for event_iterator = 1 : min(total_swr_events,100)

				event_start = max([swr_start_indices(event_iterator) - event_padding_samples, 1]);
				event_end = min([swr_end_indices(event_iterator) + event_padding_samples, total_lfp_samples]);

				edges = [event_start, event_end];

				event_indices = get_indices_from_edges(edges);	

				swr_indices_array{event_iterator} = event_indices; 


				figure(2)
				subplot(10,10,event_iterator)
				plot(lfp_pyramidal(event_indices))
				axis('off')
				
			end % End for event_iterator
			return
		end % End if plot_swr