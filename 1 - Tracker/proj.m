%% Script for tracking particles in given file

% Author: Shaurya Gupta
% Date: March 11th, 2018

%% Read in files/data and Detect Objectes in each frame
% This section uses built in a MATLAB function to detect objects in a binary
% image. The function used here is: regionprops

% Dimensions for the current problem
dim = 2;

% Call Daniel's Detection code
addpath('/Users/shauryagupta/Documents/image_processing/2 - Detection/')
points = dipfilteringcontrol();

% Number of frames to track points
n_frames = numel(points);

% Estimate of the number of points per frame
points_per_frame = 4;

%% Plot the points
% We plot a 'x' at each point location, and an index of the frame they are
% in next to the mark.

figure(1)
clf
hold on
for frame = 1:1:n_frames

    str = num2str(frame);
    for j_point = 1:1:size(points{frame}, 1)
        pos = points{frame}(j_point, :);
        plot(pos(1), pos(2), 'x')
        text('Position', pos, 'String', str)

        pause(0.01)
    end

end

%% Track points

max_linking_distance = 150;
max_gap_closing = 3;

[tracks adjacency_tracks] = tracker(points,...
    'MaxLinkingDistance', max_linking_distance, ...
    'MaxGapClosing', max_gap_closing);

%% Plot tracks
% We want to plot each track in a given color. Normally we would have to
% retrieve the points coordinates in the given |points| initiall cell
% arrat, for each point in frame. To skip this, we simple use the
% adjacency_tracks, that can pick points directly in the concatenated
% points array |all_points|.

n_tracks = numel(tracks);
colors = hsv(n_tracks);

all_points = vertcat(points{:});

for i_track = 1 : n_tracks

    % We use the adjacency tracks to retrieve the points coordinates. It
    % saves us a loop.

    track = adjacency_tracks{i_track};
    track_points = all_points(track, :);

    plot(track_points(:,1), track_points(:, 2), 'Color', colors(i_track, :))
    pause(0.2)

end
