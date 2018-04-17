    function [adjacency_tracks] = tracker(points, varargin)
% Tracker: This function implements a tracking algorithm to
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
%     |x_coordinate y_coordinate|
%
% 'MaxLinkingDistance'
% This variable defines a maximum distance between objects to prevent linking.
% Thus, even if the particles are the remaining closest pair, they will not
% be linked if the distance between them is larger than the value specified by
% MaxLinkingDistance.
% By default, the value is set to infinity in order to promote linking.
%
% 'MaxGapClosing'
% This variable defines the maximum number of frames that the object can be
% lost before it is considered as a new object.
% i.e. If an object disapprears in a frame and reapprears in a frame that is
% less than 'MaxGapClosing' frames away from the starting frame, then the
% object in condideration will be linked. Otherwise, the object in the
% subsequent frame will be considered as a new object.
% By default, the value of this variable is set to 3.
%
% Output parameters:
% 'adjacency_tracks'
% This variable is a n x 1 cell array, where n is the total number of
% tracks. Each cell in the cell array contains the indices for the detected
% points in the order that they were detected by the detection algorithm.

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
    [target_i,unmatched_targets{i+1}] = ...
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

        % The source (row) line number in the adjacency matrix
        row_i = curr_slice_i + unmatched_sources{i}(k);

        % The target (column) line number in the adjacency matrix
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
  disp('Building Tracks')

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
  disp('Building Complete')
end

%% SUBROUTINES
% Defininig the Nearest Neighbour Algorithm
function [target_indices unassigned_targets] = nearestneighborlinker(source, target, max_distance)
%NEARESTNEIGHBORLINKER link two lists of points based on nearest neighbor.
%
% target_indices = NEARESTNEIGHBORLINKER(source, target) finds for each
% point in 'source' the closest point in 'target'. These 2 inputs must be
% arrays with one point per row, and have their cartesian coordinates in
% each column (1D, 2D, 3D, ...). Nearest neighbor matching is based on
% euclidean distance. The two arrays might not have the same number of
% points.
%
% The indices of the 'target' points are returned in an array
% 'target_indices', so that each row in 'source' matches the corresponding
% row in 'target(target_indices, :)'.
%
% The linking is exclusive: one source point is linked to at most one
% target point, and conversely. The linking is only locally optimal: the
% two closest points amongst the two sets are sought for first, then the
% second closest pair, excluding the first, etc... This ensures that the
% resulting linking will not depend on the order of the points in each set.
%
% target_indices = NEARESTNEIGHBORLINKER(source, target, max_distance) adds
% a condition on distance. If the nearest neighbor is found to be at a
% distance larger than the given 'max_distance', they are not linked, and
% the 'target_indices' receive the value -1 for this source point. The same
% happens if all target points are exhausted.
%
% [target_indices unmatched_targets]=
%                                   NEARESTNEIGHBORLINKER(source, target)
% additionaly return the indices of the points in 'target' that have not
% been linked.
%
% This is the cheapest (in term of accuracy) algorithm for linking that can
% be made. In particular, it is not guaranteed (and it is generally not the
% case) that the returned linking is an optimum for the sum of distances.
% Each source point is matched regardless of the others, there is no global
% optimization here (the Hungarian algorithm does that). Also, there exists
% refinement to nearest neighbor searches, such as the use of KD-trees;
% this contribution is exempt of such developments.
%
% This code is an adaptation of code written by:
% Jean-Yves Tinevez (2012)

    if nargin < 3
        max_distance = Inf;
    end

    n_source_points = size(source, 1);
    n_target_points = size(target, 1);

    D = NaN(n_source_points, n_target_points);

    % Build distance matrix
    for i = 1 : n_source_points
        % Pick one source point
        current_point = source(i, :);

        % Compute square distance to all target points
        diff_coords = target - repmat(current_point, n_target_points, 1);
        square_dist = sum(diff_coords.^2, 2);

        % Store them
        D(i, :) = square_dist;
    end

    % Deal with maximal linking distance: we simply mark these links as already
    % treated, so that they can never generate a link.
    D ( D > max_distance * max_distance ) = Inf;

    target_indices = -1 * ones(n_source_points, 1);
    target_distances = NaN(n_source_points, 1);

    % Parse distance matrix
    while ~all(isinf(D(:)))

        [ min_D closest_targets ] = min(D, [], 2); % index of the closest target for each source points
        [ ~, sorted_index ] = sort(min_D);

        for i = 1 : numel(sorted_index)
            source_index =  sorted_index(i);
            target_index =  closest_targets ( sorted_index(i) );

            % Did we already assigned this target to a source?
            if any ( target_index == target_indices )
                % Yes, then exit the loop and change the distance matrix to
                % prevent this assignment
                break
            else
                % No, then store this assignment
                target_indices( source_index ) = target_index;
                target_distances ( source_index ) = sqrt ( min_D (  sorted_index(i) ) );

                % And make it impossible to find it again by putting the target
                % point to infinity in the distance matrix
                D(:, target_index) = Inf;
                % And the same for the source line
                D(source_index, :) = Inf;
                if all(isinf(D(:)))
                    break
                end
            end
        end
    end
    unassigned_targets = setdiff ( 1 : n_target_points , target_indices );
end
