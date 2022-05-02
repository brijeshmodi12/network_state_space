%% gaussian_smooth for 2D maps:
%% gaussian2d: function description
function [smoothed_output] = gaussian2d(x, mu, sigma)

	smoothed_output =  exp(-0.5 *((x-mu) / sigma).^2) ./ ( sigma*sqrt(2*pi));;


