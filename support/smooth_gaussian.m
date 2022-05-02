function [smoothed_data] = smooth_gaussian(data, gaussian_std, gaussian_width)
	% Gaussian smoothing function 
	% Input :
	% -	data 						Row/Column vector
	% - gaussian_std				Standard deviation for gaussian kernel (Units : No. of samples)
	% - gaussian_width				Width for gaussian kernel (Units : No. of samples)

	% Output :
	% - smoothed_data				Row/Column vector smoothed version of data.

	% Create basis of gaussian kernel with user-defined width.
	basis = round(-gaussian_width / 2)  : round(gaussian_width / 2);

	% Create Gaussian kernel on that basis with user-defined standard deviation.
	kernel = normpdf(basis, 0 , gaussian_std);
	% figure
	% plot(kernel)
	% Returns smooth data of same length as original data
	smoothed_data = conv(data, kernel, 'same');

