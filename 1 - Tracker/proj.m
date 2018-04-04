%% Script for tracking particles in given file
% The objective of this script is to initiate the detection of the selected
% dataset and subsequently track the positions of the detected objects
% throughout the available frames.
% File Dependencies:
% 1) DIP Image
% 2) dipfiltering.m
% 3) tracker.m
% 4) nearestneighborlinker.m

% Author: Shaurya Gupta
% Date: March 11th, 2018

clear
close all

%% Read in files/data and Detect Objectes in each frame
% This section sets the dimension for the tracking problem. The following
% code calls the detection algorithm to detect the objects in the current
% working directory.
% Note: The detection algorithm requires the DIP Image Library to be initialized

% Dimensions for the current problem
dim = 2;

% Initiate the DIP Image Library
addpath('/Applications/dip/common/dipimage');
dip_initialise

% Call the detection algorithm
% Note: The detection algorithm does not require any input as it operates on
% all the '*.tif' images present in the current working directory.
% The output of the detection algorithm is stroed in a 'dim' dimensional array.
% For the purposes of current application, 'points' is a 2-dimensional array
% with the follwing structure: | x | | y |;
% where 'x' and 'y' are the x and y coordinates of the centroids of the detected
% points in the cartesian coordinate reference frame
addpath('/Users/shauryagupta/Documents/image_processing/2 - Detection/')
points = dipfilteringcontrol();

% Number of frames to track points
n_frames = numel(points);

%% Plot the points
% This section of the script produces two plots. The first plot is produced to
% visualize all the detected points in the available frames of the dataset.
% Here, all the points detected in the 1st frame are labelled as 'x' followed
% by the frame number in which they were detected.
figure(1)
hold on
for frame = 1:1:n_frames
    str = num2str(frame);

    for j_point = 1:1:size(points{frame}, 1)
        pos = points{frame}(j_point, :);
        plot(pos(1), pos(2), 'x','Markersize',6,'MarkerEdgeColor','red')
        text('Position', pos, 'String', str, 'FontSize', 6)
    end
    xlabel('X Position')
    ylabel('Y Position')
    title('Detected Objects in All Frames of Dataset')
end
hold off

% The second plot is produced to visualize the path or 'tracks' of the detected
% object as it moves from the first frame to the last frame. Here, 'x' marks the
% begining of a track, whereas 'o' marks the last position of the object.
% The 'hold on' MATLAB command is used to hold on to the figure in order to plot
% the tracks after they have been determined by the tracking algorithm.
figure(2)
clf
hold on
for frame = 1:(n_frames-1):n_frames
    if frame == 1
      for j_point = 1:1:size(points{frame}, 1)
        pos = points{frame}(j_point, :);
        plot(pos(1), pos(2), '-x','Markersize',10,'MarkerEdgeColor','red',...
        'MarkerFaceColor',[1 .6 .6])
      end
    else
      for j_point = 1:1:size(points{frame}, 1)
        pos = points{frame}(j_point, :);
        plot(pos(1), pos(2), '-o','Markersize',10,'MarkerEdgeColor','red',...
        'MarkerFaceColor',[1 .6 .6])
      end
    end
    xlabel('X Position')
    ylabel('Y Position')
    title('Object Tracks')
end

%% Track points
% This section of the script initiates the tracking algorithm. It is a
% requirement of this script to have the 'tracking.m' in the same working
% directory as the main proj.m file. In order to initialize, the algorithm
% requires the definition of two input variables: 'max_linking_distance' and
% 'max_gap_closing'. Find the definitions of these variables in the
% 'tracker.m' file.
% The output of the tracking algorithm is stored in the form of a cell array,
% named 'adjacency_tracks'
% 'adjacency_tracks' is a n x 1 cell array, where n is the total number of
% tracks. Each cell in the cell array contains the indices for the detected
% points in the order that they were detected by the detection algorithm.
% (See tracker.m for more information)

% Control Case
max_linking_distance = 150;
max_gap_closing = 3;

% Real Data
% max_linking_distance = 50;
% max_gap_closing = 1;

[adjacency_tracks] = tracker(points,...
    'MaxLinkingDistance', max_linking_distance, ...
    'MaxGapClosing', max_gap_closing);

%% Plot tracks
% This section of the script overlays the object tracks on Figure 2 (containing
% the start and end positions of the objects).
% Each seperate track is indicated by a different color.
% In order to simplify the implementation of the plot, the 'adjacency_tracks'
% variable (containing the indices of the detected objects) is used to find the
% object and its coordinate from an array containing all detected points (in the
% order that they were orginally detected).

n_tracks = numel(adjacency_tracks);

% Creating a color palette
colors = hsv(n_tracks);

% Creating a array containing all the coordinates for all the centroids of all
% the objects detected, in succession.
all_points = vertcat(points{:});

% Looping through all the object tracks
for i_track = 1 : n_tracks
    % Picking out individual tracks from adjacency_tracks cell array
    track = adjacency_tracks{i_track};
    % Picking out all the points in a particular track
    track_points = all_points(track, :);

    plot(track_points(:,1), track_points(:, 2), 'Color', colors(i_track, :))
end
