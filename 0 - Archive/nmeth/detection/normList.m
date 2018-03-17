function [listOfNorms,normedVectors]=normList(vectors)
%calculates the norm of a list of vectors
%
%SYNOPSIS [listOfNorms,normedVectors]=normList(vectors)
%
%INPUT list of vectors (nVectorsXdimension)
%
%OUTPUT listOfNorms: list (nX1) containing the norms of the vectors
%       normedVectors: list (nXdim) containing the normed vectors
%
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
%Copyright: 1/03 Jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nVectors,nDims]=size(vectors);

listOfNorms=zeros(nVectors,1);
normedVectors=zeros(nVectors,nDims);

listOfNorms=sqrt(sum(vectors.^2,2));
goodVectors=find(listOfNorms);

normedVectors(goodVectors,:)=vectors(goodVectors,:)./(repmat(listOfNorms(goodVectors),[1,nDims]));
