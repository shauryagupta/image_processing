%% Script for tracking particles in given file

% Author: Shaurya Gupta
% Date: March 11th, 2018

%% Read in files/data

% % Get files from the Control_data folder
% folderName = '/Users/shauryagupta/Documents/MATLAB/SimpleTracker/Control_data/';
% cd(folderName)
% addpath(pwd)
%
% % Store file names in a variable
% % The files are of type .png. Change value here if the input file type changes
% imdir = dir('*.png');
%
% % Create variable to store image stack
% image = zeros(400,400,numel(imdir));
%
% % Read in image files and store in variable
% for i = 1:1:numel(imdir)
%   imcurr = im2bw(imread(strcat(folderName,imdir(i).name)),0.5);
%   image(:,:,i) = imcrop(imcurr,[136 46 399 399]);
% end
%
% % Go back to parent directory
% cd ..

load('image.mat')

%% Detecting Objectes in each frame
% This section uses built in a MATLAB function to detect objects in a binary
% image. The function used here is: regionprops

% Dimensions for the current problem
dim = 2;

% Number of frames to track points
n_frames = size(image,3);

% Estimate of the number of points per frame
points_per_frame = 35;

% Create variable to store points
points = cell(n_frames,1);

for iframe = 1:1:n_frames
  % Detecting points in each frame
  cc = bwconncomp(image(:,:,iframe));
  stats = regionprops(cc,'Centroid');

  frame_points = zeros(cc.NumObjects,dim);

  for j = 1:1:cc.NumObjects
    frame_points(j,:) = stats(j).Centroid;
  end

  % Storing value in variable
  points{iframe} = frame_points;
end

%% Plot the random points
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

        %pause(0.1)
    end

end

%% Track them
% Finally! A one liner. We add some information to the output, and allow
% gap closing to happen all the way through.

max_linking_distance = 10;
max_gap_closing = 3;

[ tracks adjacency_tracks ] = tracker(points,...
    'MaxLinkingDistance', max_linking_distance, ...
    'MaxGapClosing', max_gap_closing);

%% Plot tracks
% We want to plot eahc track in a given color. Normally we would have to
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
    pause(0.01)

end
