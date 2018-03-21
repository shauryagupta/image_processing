function [tracks adjacency_tracks] = tracker(points, varagrin)
% Tracker: This function is based on the SIMPLETRACKER implementation by
% Jean-Yves Tinevez. This function implements a tracking algorithm to
% link particles between frames (in addition to dealing with gaps)
%
% Dependencies: In additon to the parameters outlined above, this code requires
% the following files to function properly:
% nearestneighborlinker.m: Implementation of a nearest-neighbor algorithm for
% finding the same particles in adjacent frames in oder to perform frame-to-
% frame linking. This algorithm achieves a local optimum for a pair of points.
%
% Input parameters:
% points - represents the matrix of the coordinates of the particles found in
% each frame of the input images (image.mat)
%   - points must be a cell, with one cell dedicated for each frame
%   - each cell has the following structure:
%     |#points x_coordinate y_coordinate|
%
% /*** CHANGE BEFORE SUBMISSION ***/
% 'MaxLinkingDistance' - a positive number, by default Inifity.
% Defines a maximal distance for particle linking. Two particles will not
% be linked (even if they are the remaining closest pair) if their distance
% is larger than this value. By default, it is infinite, not preventing nay
% linking.
%
% 'MaxGapClosing' - a positive integer, by default 3
% Defines a maximal frame distance in gap-closing. Frames further way than
% this value will not be investigated for gap closing. By default, it has
% the value of 3.
% /*** CHANGE BEFORE SUBMISSION ***/
%
% Output parameters:
% /*** CHANGE BEFORE SUBMISSION ***/
% Enter description text
% /*** CHANGE BEFORE SUBMISSION ***/

  % Parse parameters
  p = inputParser;

  defaultMaxGapClosing = 3;
  defaultMaxLinkingDistance = Inf;

  % Update parameters if supplied by user
  p.addParamValue('MaxGapClosing', defaultMaxGapClosing, @isnumeric);
  p.addParamValue('MaxLinkingDistance', defaultMaxLinkingDistance, @isnumeric);

  p.parse(varargin{:});

  max_gap_closing = p.Results.MaxGapClosing;
  max_linking_distance = p.Results.MaxLinkingDistance;

  % Performing frame to frame linking
  disp('Performing frame to frame linking using NearestNeighbor Algorithm');

  num_slices = numel(points);

  % Initializing variables
  curr_slice_i = 0;
  row_i = cell(num_slices,1);
  column_i = cell(num_slices,1);
  unmatched_targets = cell(num_slices,1);
  unmatched_sources = cell(num_slices,1);

  % The number of detected points in each frame
  n_cells = cellfun(@(x) size(x,1),points);

  for i = 1:1:(num_slices-1)
    % clears print message from previous iteration
    fprintf(remapt('\b',1,7));
    fprintf('%03d/%03d',i,n_slices-1);

    source = points{i};
    target = points{i+1};

    % Using the nearest-neighbor algorithm for frame to frame linking
    [target_i,~,unmatched_targets{i+1}] = ...
      nearestneighborlinker(source,target,max_linking_distance);

    % For the nearest-neighbor algorithm, all indices are initially set to -1 and
    % target distaces are set to NaN. These are then changed in the loop based on
    % valid iterations
    unmatched_sources{i} = find(target_i == -1);

    % Initialize variables to hold vaid links generated above
    num_links = sum(target_i ~= -1);
    row_i{i} = NaN(num_links,1);
    column_i{i} = NaN(num_links,1);

    % Implementing an adjacecy matrix: a square matrix to represnt a fitnite
    % graph
    index = 1;
    for j = 1:1:numel(target_i)
      % Skip iteration if proper link to target not found
      if target_i(j) == -1
        continue
      end

      % The source node in the adjacecy matrix
      row_i{i}(index) = curr_slice_i + j;

      % The target node in the adjacency matrix
      column_i{i}(index) = curr_slice_i + n_cells(i) + target_i(j);

      % Update index
      index = index + 1;
    end

    curr_slice_i = curr_slice_i + n_cells(i);
  end

  % Convert cells to matrix form
  row_i = vertcat(row_i{:});
  column_i = vertcat(column_i{:});
  
