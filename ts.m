
clear all; close all ; clc
meta_data_root_folder = 'D:\Silvia Datasets\singol planes';
tif_root = 'D:\Silvia Datasets\singol planes\images';


ctrl_string = 'CTR'; 
plx_string = 'PLX';
syn_string = 'SYN';
dre_string = 'DRE';
	
% tif_filename = 'MAX_C2-Composite_SCSR_PLX001_big.tif 5P.tif';

toprocess = 'SYN';

if toprocess == syn_string
	image_prestring = '';
	image_poststring = '';
	output_folder = fullfile(meta_data_root_folder,'results_syn');
else
	image_prestring = '';
	image_poststring = '';
	output_folder = fullfile(meta_data_root_folder,'results_dre');
end 

mkdir(output_folder)

xcol = 6; 
ycol = 7;
voxelcol = 2;

binsize = [20,20];

% total_random_points = 500;
% total_random_iteration = 9;

% lower_tri_matrix = tril(magic(total_random_points),-1);
% lower_tri_matrix = lower_tri_matrix > 0;

color_gray = [134, 138, 145 ] ./ 255;
color_pink = [189, 11, 73] ./ 255;
color_blue = [12, 235, 227] ./ 255;

listings = dir(meta_data_root_folder);
filenames = {listings.name};
total_files = length(filenames);
% return



global_matrix_storage = [];
global_type_storage = [];
global_density_distribution = [];
global_type_density = [];
global_density = [];

global_peak_density = [];



for file_iterator = 1:total_files

	current_file = char(filenames{1,file_iterator});

	density_distribution = [];
	current_type = [];

	[~,cfname,~] = fileparts(current_file);

	if length(current_file) < 3;
		continue;
	end

	
	if current_file(1:3) ~= toprocess
		% fprintf('File Skipped\n')
		continue;
	end

	fprintf('Prcessing file: %s\n',current_file)
	stringcheck = strfind(current_file, plx_string);
	if isempty(stringcheck)
		file_type = 'CTR';
		current_type = 0;
		temp_string = sprintf('%s',cfname);
		tif_filename = sprintf('%s%s%s.tif',image_prestring,  temp_string, image_poststring)
		
	else
		file_type = 'PLX'
		current_type = 1;
		temp_string = sprintf('%s',cfname);
		tif_filename = sprintf('%s%s%s.tif',image_prestring,  temp_string, image_poststring)
	end


	meta_data_filepath = fullfile(meta_data_root_folder, current_file)

	meta_data_table = readtable(meta_data_filepath,'Format','auto');



	xval = cell2mat(table2cell(meta_data_table(:,xcol)));
	yval = cell2mat(table2cell(meta_data_table(:,ycol)));
	voxelval = cell2mat(table2cell(meta_data_table(:,voxelcol)));

	% xval = randi([1 1000], 5000, 1);
	% yval = randi([1 1000], 5000, 1);

	points = [xval yval];
	% caxis_range = [min(voxelval) prctile(voxelval,99)];
	% [points_cmap, points_colors] = assign_colors(voxelval, caxis_range, 10000, 'brbg'); 

	outliers_ind = find(voxelval>prctile(voxelval,95));
	normal_ind = find(voxelval<=prctile(voxelval,95));

	points_colors = [];
	points_colors(normal_ind,:) = repmat(color_gray, length(normal_ind),1);
	points_colors(outliers_ind,:) = repmat(color_pink, length(outliers_ind),1);

	points_size = rescale(voxelval, 20, 100);
	% return


	total_points = size(points,1);

	x_edges = linspace(min(xval),max(xval)+0.5, binsize(1)+1);
	y_edges = linspace(min(yval),max(yval)+0.5, binsize(1)+1);
	edges = {x_edges, y_edges};

	hist3d_values = hist3(points, 'Edges', edges);

	hist3copy = hist3d_values;


	hist3d_values(binsize(1)+1,:) = [];
	hist3d_values(:,binsize(1)+1) = [];

	edges1d = linspace(0,50,30);


	density_distribution = histogram(hist3d_values(:), 'BinEdges', edges1d);

	density_distribution = density_distribution.Values;
	figure(98)
	set(gcf, 'Position', get(0, 'Screensize'));
	subplot(1,2,1)
	scatter(points(:,1), points(:,2), points_size, points_colors,'filled',  'MarkerFaceAlpha',.7)
	pbaspect([1 1 1])
	% colorbar
	% caxis(caxis_range);
	% return
	% return
	
	xticks(x_edges)
	yticks(y_edges)
	xlim([min(xval) max(xval)])
	ylim([min(yval) max(yval)])
	grid on
	ax = gca;
	ax.GridLineStyle = '-';
	ax.GridColor = color_blue;
	ax.GridAlpha = 1; % maximum line opacity
	set (gca, 'xticklabel' , {[]});
	set (gca, 'yticklabel' , {[]});
	set(gca,'ydir','reverse')
	% return

	subplot(1,2,2)
	tif_path = fullfile(tif_root,tif_filename)
	a = Tiff(tif_path,'r');
	ab =read(a);
	imshow(ab)

	% binsize_pixels = 100;

	total_image_rows = size(ab,1);
	total_image_cols = size(ab,2);

	row_edges = linspace(1,total_image_rows,binsize(1)+1);
	col_edges = linspace(1,total_image_cols,binsize(2)+1);

	hold on;
	for row = row_edges
	  line([1, total_image_cols], [row, row], 'Color', 'w');
	end
	for col = col_edges
	  line([col, col], [1, total_image_rows], 'Color', 'w');
	end

	% set(gca,'xdir','reverse')
	

	% plot_matrix = imgaussfilt(flipud(hist3d_values'), 'FilterSize', 5);
	% % figure
	% subplot(1,4,3)
	% imagesc(plot_matrix)
	% title('Norm. Density Map')
	% colormap('viridis');
	% pbaspect([1 1 1])
	% colorbar
	% caxis([min(plot_matrix(:)) max(plot_matrix(:))])

	% subplot(1,4,4)
	% histogram(hist3d_values(:), 'BinEdges', edges1d);
	% xlim([0 50])
	% xlabel('Local Density (synapsin / unit space)')
	% ylabel('Counts')
	% % ylim([0 100])
	% pbaspect([1 1 1])
	% title(sprintf('%s',cfname));

	suptitle(sprintf('%s',cfname))

	
	fname = sprintf('%s_local_density.png', cfname);
	filepath = fullfile(output_folder, fname);
	saveas(gcf, filepath);

	% return
	close(figure(98))
	% return

	smoothed_density_distribution = reshape(smooth(density_distribution,5), 1, []);

	global_density_distribution = [global_density_distribution ; smoothed_density_distribution];
	global_type_density = [global_type_density; current_type];
	global_density = [global_density; total_points];
	global_peak_density = [global_peak_density; max(hist3d_values(:))];
	
	% return


	% for random_iterator = 1:total_random_iteration
	% 	fprintf('%d...',random_iterator)


	% 	randomset = []; pairwise_distance = []; temp_matrix = []; distance_matrix = [];
		
	% 	randomset = sort(randperm(total_points,total_random_points));

	% 	subset = points(randomset,:);

	% 	pairwise_distance = pdist(subset);

	% 	distance_matrix = squareform(pairwise_distance + rand(size(pairwise_distance))*0.00001 );


	% 	temp_matrix = distance_matrix(lower_tri_matrix);

	% 	temp_matrix = reshape(temp_matrix, 1, [] ) ;

	% 	global_matrix_storage = [global_matrix_storage; temp_matrix] ;

	% 	global_type_storage = [global_type_storage; current_type];

	% 	plotfigure = 0;

	% 	if plotfigure
	% 		figure(1)
	% 		subplot(3,3,random_iterator)
	% 		imagesc(distance_matrix)
	% 		colormap('viridis')
	% 		title(sprintf('RI %d',random_iterator))
	% 		pbaspect([1 1 1])

	% 		figure(2)
	% 		subplot(3,3,random_iterator)
	% 		scatter(points(:,1), points(:,2), 5, color_gray,'filled')
	% 		hold on 
	% 		scatter(points(randomset,1), points(randomset,2), 10, color_pink,'filled')
	% 		pbaspect([1 1 1])
	% 		title(sprintf('RI %d',random_iterator))

	% 	end

	% end % End for random_iterator	

	% fprintf('\n')

	% if plotfigure
	% 	figure(1)
	% 	set(gcf, 'Position', get(0, 'Screensize'));
	% 	suptitle(sprintf('%s',cfname))
	% 	fname = sprintf('%s_random_points.png', cfname);
	% 	filepath = fullfile(output_folder, fname);
	% 	saveas(gcf, filepath);


	% 	figure(2)
	% 	set(gcf, 'Position', get(0, 'Screensize'));
	% 	suptitle(sprintf('%s',cfname))
	% 	fname = sprintf('%s_distance_matrix.png', cfname);
	% 	filepath = fullfile(output_folder, fname);
	% 	saveas(gcf, filepath);

	% 	close all;
	% end

end % End if file_iterator



output_data_filename = sprintf('global_results.mat');
output_data_path = fullfile(output_folder, output_data_filename);
save(output_data_path, 'global_density_distribution','edges1d' ,'global_density', 'global_type_density','ctrl_string','plx_string')

return
[projection, ~] = run_umap(global_density_distribution,  'n_components', 2, 'n_neighbors', 10 ,'min_dist', 0.1, 'metric', 'cosine' );

[projection_cmap, projection_colors] = assign_colors(global_type_density, [0 1], 2, 'viridis');
figure(33)
scatter(projection(:,1), projection(:,2), 50, projection_colors,'filled')
xlim([-5 20])
ylim([-20 0])
colormap(projection_cmap)



rescaled_density_distribution = rescale(global_density_distribution');
rescaled_density_distribution = rescaled_density_distribution';

ctrl_distribution_avg = mean(rescaled_density_distribution(1:9,:));
plx_distribution_avg = mean(rescaled_density_distribution(10:18,:));
figure
plot(ctrl_distribution_avg,'r');
hold on 
plot(plx_distribution_avg, 'b')