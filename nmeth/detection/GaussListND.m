function gaussList = GaussListND(coordList,sigma,center,intNorm)
%GAUSSLISTND calculates the value of a N-D Gaussian at specific pixel/voxel coordinates
%
% SYNOPSIS gaussList = GaussList23D(coordList,sigma,center,intNorm)
%
% INPUT    coordList : m-by-n list of coordinates, where m is the number of
%                      coordinates and n the number of dimensions
%          sigma     : 1-by-n (or scalar): sigma of Gaussian
%          center    : (opt) 1-by-n vector of center of Gaussian.
%                      Default: zeros(1,n)
%          intNorm   : (opt) switch for how the Gaussian should be normed
%                      Default: 0
%                      0 - no norming. Max of Gaussian == 1
%                      1 - normed so that integral of infinite Gaussian = 1
%
% OUTPUT   gaussList : m-by-1 list of intensities. Intensity is the
%                      integral of the Gaussian over the pixel/voxel
%
% REMARKS  The code assumes that a pixel has the edge length 1!
%This file is part of u-track.
%
%    u-track is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%    u-track is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with u-track.  If not, see <http://www.gnu.org/licenses/>.
%
% Copyright: 2/05 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%======================
% TEST INPUT
%======================

% check number of input arguments
nIn = nargin;
% the following doesn't work with Matlab 6.5.0
% error(nargchk(2,4,nIn,'struct'));
if nIn < 2 || nIn > 4
    error('wrong number of input arguments!')
end

% check dimensionality of coordList. 
if isempty(coordList)
    error('you have to supply a list of coordinates for GaussList23D')
else
    [nCoords,nDims] = size(coordList);
end

% sigma
ls = length(sigma);
switch ls
    case nDims
        % make as long as coords
        sigma = repmat(sigma,[nCoords,1]);
    case 1
        sigma = repmat(sigma,[nCoords,nDims]);
    otherwise
        error('sigma has to be a scalar or a 1-by-n vector!')
end

% center
if nIn < 3 || isempty(center)
    center = zeros(nCoords,nDims);
else
    lc = length(center);
    switch lc
        case nDims
            center = repmat(center, [nCoords,1]);
        case 1
            center = repmat(center, [nCoords,3]);
        otherwise
            error('center has to be a scalar or a 1-by-n vector!')
    end
end

% intNorm
if nIn < 4 || isempty(intNorm)
    intNorm = 0;
end

%======================

%======================
% CALC GAUSSLIST
%======================

% 0.5*erfc(-(x+0.5)/sqrt(2))-0.5*erfc(-(x-0.5)/sqrt(2)) gives the integral on the
% pixel at 1 of a Gaussian with mean 0 and sigma 1

% convert coordList to 0/1
coordList = (coordList - center)./sigma;

% double coordList as preparation for erfc
%fixed bug: must divide the 0.5 by sigma - KJ
coordList = cat(3,coordList-0.5./sigma, coordList+0.5./sigma);

% calculate gaussList
%Jonas was missing the minus sign in erfc. I corrected that - KJ
gaussList = diff(0.5 * erfc(-coordList/sqrt(2)),1,3);
gaussList = prod(gaussList,2);

% norm gaussList
switch intNorm
    case 0
        % "un-norm" Gaussian
        gaussList = gaussList*((2*pi)^(0.5*nDims)*prod(sigma(1,:)));
    case 1
        % gaussList is already normed
end