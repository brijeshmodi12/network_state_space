function varargout = PlotComet_2D(x_data, y_data,  varargin)

%function creates a moving comet plot in 3 dimensions using the specificed
%x, y, and z line data
%
%PlotComet_3D(x_data, y_data, z_data)
%
%Optional Input Arguements:
%cFigure, handle to the current figure or axes to create the plot in (will
%   add the line data to the plot regardless of whether hold is on
%
%Frequency, playback frequency, if not provided assumed to be 10 Hz
%
%blockSize, number of points in the comet tail, assumed to be 1/20 of
%   length(x_data) if not provided
%
%PlotComet_3D(x_data, y_data, z_data, 'cFigure',FigureHandle,...
%     'Frequency',freq, 'blockSize',num_points)
%
%tailFormat, a structured variable containing the fomrating requirements of
%   the tail (if different than the default red solid line)
%   ex: Tail = struct('LineWidth',2,'Color','g','LineStyle','*');
%   PlotComet_3D(x_data,y_data,z_data,'tailFormat',Tail);
%
%headFormat, a structued vairable containing the formating requirements of
%   the head (if different than the default blue circle)
%   ex: Head = struct('LineStyle','none','MarkerSize',8,'Marker','^','Color','k');
%   PlotComet_3D(x_data,y_data,z_data,'headFormat',head);
%   Note: 'LineSytle','none' is required since the head is plotted as a
%   single point!
%
%Optional Output Arguements
%hFigure = figure handle to the comet plot
%
%Example: Create a surface plot and construct a moving line plot through it
%with connected cyan stars and dashed line on a 50 points long tail
%and large green square for a head.
%
%%create some data for the surface
%[X, Y] = meshgrid(linspace(-1, 1, 25), linspace(-1, 1, 25));
%Z = X .* exp(-X.^2 - Y.^2);
%%open a new figure and plot the surface
%surf(X, Y, Z);
%fig = gcf;
%%create the line data
%t = -pi:pi/500:pi;
%x = sin(5*t);
%y = cos(3*t);
%z = t;
%%define the tail formated structure
%Tail = struct('LineStyle','--','Marker','*','Color','c',...
%   'LineWidth',2,'MarkerSize',4);
%%define the head formated structure
%Head = struct('LineStyle','none','Marker','s','Color','g','MarkerSize',10);
%
%%Invoke the comet plot with 50 Hz playback, and a tail length of 50 pts
%PlotComet_3D(x,y,z,'cFigure',fig,'blockSize',50,'Frequency',50,...
%   'headFormat',Head,'tailFormat',Tail);
%
%Nick Hansford
%09/26/2013

%initialize the inputs

color_gray = [134, 138, 145 ] ./ 255;
head_color = [189, 11, 73] ./ 255;
color_tail = [133, 133, 133] ./ 255;
freq = 100;
blockSize = floor(1/1*length(x_data));
blockSize = 2;
tailFormat = struct('LineWidth',2,'Color',head_color,'LineStyle','-');
headFormat = struct('LineStyle','none','Marker','o','MarkerSize',10,...
    'MarkerFaceColor','black');


% movie_vector = [];

%parse out the inputs
argCount = nargin - 2;

for i = 1:2:argCount
    %         caseVar = varargin{i}
    switch varargin{i}
        case 'cFigure'
            cFigure = varargin{i+1};
            
            %get the original size:
            cAxes = get(cFigure,'CurrentAxes');
            xLim = get(cAxes,'XLim');
            yLim = get(cAxes,'YLim');
            % zLim = get(cAxes,'ZLim');
            
            %resize if the new plot will exceed them
            if xLim(1) > min(x_data)
                xLim(1) = min(x_data);
            end
            
            if xLim(2) < max(x_data)
                xLim(2) = max(x_data);
            end
            
            if yLim(1) > min(y_data)
                yLim(1) = min(y_data);
            end
            
            if yLim(2) < max(y_data)
                yLim(2) = max(y_data);
            end
           
            axis([xLim, yLim]);
            
        case 'blockSize'
            blockSize = varargin{i+1};
        case 'Frequency'
            freq = varargin{i+1};
            
        case 'tailFormat'
            tailFormat = varargin{i+1};
            
        case 'headFormat'
            headFormat = varargin{i+1};

        case 'timeStep'
            timeStep = varargin{i+1};    

        case 'outputPath'    
            outputPath = varargin{i+1}
            
        case 'plotSpeed' 
            plot_speed = varargin{i+1}
    end
end
    
    
%make sure the figure exists, if not, create it
if ~exist('cFigure')
    cFigure = figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    view(2);
    axis([min(x_data), max(x_data),...
        min(y_data), max(y_data)]);
end

    
%user can pass in current axes, if so, get the parent for the figure window
%and use that instead...should make this compatible with GUIs?
if strcmp(get(cFigure,'Type'),'axes')
    cFigure = get(cFigure,'Parent');
end

%activate the figure window
figure(cFigure);
cAxes = get(cFigure,'CurrentAxes');

oldNextPlot = get(cAxes,'NextPlot');
set(cAxes,'NextPlot','add');

pauseTime = 1./freq;

n_start = 1;
n_stop = 1;
%put on the starting point
% plot(x_data(n_start:n_stop), y_data(n_start:n_stop) ,tailFormat);
% plot(x_data(n_stop), y_data(n_stop), headFormat);
scatter(x_data, y_data, 15, color_gray, 'filled');
pbaspect([1 1 1]) 
hold on

frame_counter = 1;
movie_vector(frame_counter) = getframe(gcf, get(0, 'Screensize'));

% points per frame
% plot_speed = 100;

for n = 2:1:length(x_data)
    frame_counter = frame_counter + 1;

    if n <= blockSize
        n_start = 1;
        n_stop = n;
    else
        n_start = n_start + plot_speed;
        n_stop = n_start + plot_speed;
    end

    if n_stop >= length(x_data)
        break
    end

    % p1;
    plot(x_data(n_start:n_stop), y_data(n_start:n_stop) ,tailFormat);
    scatter(x_data(n_stop), y_data(n_stop), 50, head_color, 'filled');

    current_time = n * plot_speed * timeStep ;
    seconds = fix(current_time);

    % seconds = fix(mod(current_time, 60));
    % minutes = floor(seconds / 60);

    title(sprintf('Time = %d sec ',seconds))
    grid on
    % view([-30 56])
    drawnow;

    movie_vector(frame_counter) = getframe(gcf, get(0, 'Screensize'));
    
    %update playback refresh rate
    % b = toc(a);
    % b = 0;
    % pause(pauseTime-b);
    pause(0.000001)


end


set(cAxes,'NextPlot',oldNextPlot);
% fprintf('Finished\n');


% writer_object = VideoWriter(filepath)
writer_object = VideoWriter(outputPath)
% writer_object.FrameRate = 1/timeStep;
writer_object.FrameRate = 30;

open(writer_object);
writeVideo(writer_object, movie_vector);
close(writer_object);

if nargout == 1
    varargout{1} = cFigure;
end



% % writer_object = VideoWriter(filepath)
% writer_object = VideoWriter('output_realtime')
% writer_object.FrameRate = 1/timeStep;
% % writer_object.FrameRate = 30;

% open(writer_object);
% writeVideo(writer_object, movie_vector);
% close(writer_object);

% if nargout == 1
%     varargout{1} = cFigure;
% end