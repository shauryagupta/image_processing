function [tracks adjacency_tracks] = tracker(points, varargin)
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

  p.parse( varargin{:} );

  max_gap_closing = p.Results.MaxGapClosing;
  max_linking_distance = p.Results.MaxLinkingDistance;

  %% FRAME LINKING: Performing frame to frame linking
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
    fprintf(repmat('\b',1,7));
    fprintf('%03d/%03d',i,num_slices-1);

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

  % Create variables for the linking step
  link = ones(numel(row_i),1);
  n_cells_total = sum(n_cells);

  fprintf('\nCreating %d links for %d points.\n', numel(link), n_cells_total)

  % Creating the link using the MATLAB sparse function
  A = sparse(row_i,column_i,link,n_cells_total,n_cells_total);

  %% GAP CLOSING: Performing Gap closing
  disp('Performing gap-closing')

  % Resetting the current slice index
  curr_slice_i = 0;

  for i = 1:1:(num_slices-2)

    % The following lines of code find a target in subsequent frames (starting
    % at i+2) and parsing over targets that are not present in the link
    curr_target_slice_i = curr_slice_i + n_cells(i) + n_cells(i+1);

    for j = i+2:min(i+max_gap_closing, num_slices)
      source = points{i}(unmatched_sources{i},:);
      target = points{j}(unmatched_targets{j},:);

      % Check to see if the matrix is empty
      if isempty(source) || isempty(target)
        curr_target_slice_i = curr_target_slice_i + n_cells(j);
        continue
      end

      % Run another iteration of the nearestneighbor algorithm to find the
      % links between previously unmatched targets

      target_i = nearestneighborlinker(source, target, max_linking_distance);

      % Update the adjacecy matrix with the resulsts obtained above
      for k = 1:1:numel(target_i)

        % Look for unmatched targets
        if target_i(k) == -1
          continue
        end

        % Print message to track which point was matched with which other point
        fprintf('Creating a link between point %d of frame %d and point %d of frame %d.\n',...
            unmatched_sources{i}(k),i,unmatched_targets{j}(target_i(k)),j);
        
        % The source line number in the adjacency matrix
        row_i = curr_slice_i + unmatched_sources{i}(k);
        
        % The target line number in the adjacency matrix
        column_i = curr_target_slice_i + unmatched_targets{j}(target_i(k));

        % Update the adjacency matrix
        A(row_i,column_i) = 1;
      end

      % Get indices of matched (round 2) points
      new_target_i = target_i ~= -1;

      % Delete already matched sources so that they are not accidently changed
      unmatched_sources{i}(new_target_i) = [];

      % Delete already matched targets so that they are not accidently changed
      unmatched_targets{j}(target_i(new_target_i)) = [];

      % Update current frame index
      % This line broke the code
      curr_target_slice_i = curr_target_slice_i + n_cells(j);
    end

    curr_slice_i = curr_slice_i + n_cells(i);
  end

    %% BUILD TRACKS: Parse Adjacecy Matrix to build tracks
    disp('Building tracks:')

    % Initializing holding variable
    cells_no_source = [];

    % Finding columns that contain 0: indicates that the cell has no source
    for i = 1:1:size(A,2)
      if length(find(A(:,i))) == 0
        cells_no_source = [cells_no_source; i];
      end
    end

    % Initializing varibales
    n_tracks = numel(cells_no_source);
    adjacency_tracks = cell(n_tracks,1);

    A_t = A';

    for i = 1:1:n_tracks
      tmp = NaN(n_cells_total,1);
      target = cells_no_source(i);

      index = 1;
      while ~isempty(target)
        tmp(index) = target;
        target = find(A_t(:,target), 1, 'first');

        % Update index
        index = index + 1;
      end

      adjacency_tracks{i} = tmp(~isnan(tmp));
    end

    %% Reparse adjacency_tracks to update indices so that it refers to the points
    % in the original array.

    tracks = cell(n_tracks, 1);

    for i = 1:1:n_tracks
      adjacency_track = adjacency_tracks{i};
      track = NaN(num_slices, 1);

      for j = 1:1:numel(adjacency_track)
        cell_i = adjacency_track(j);

        % Determine the corresponding frame for the index j
        tmp_i = cell_i;
        frame_i = 1;

        while tmp_i > 0
          tmp_i = tmp_i - n_cells(frame_i);
          frame_i = frame_i + 1;
        end

        frame_i = frame_i - 1;
        in_frame_cell_i = tmp_i + n_cells(frame_i);

        track(frame_i) = in_frame_cell_i;
      end
      tracks{i} = track;
    end
  end
