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
%Copyright Jaqaman 01/2008
%% movie information

movieParam.imageDir = 'C:\Users\Shaurya\Desktop\image_processing\nmeth\example\'; %directory where images are
movieParam.filenameBase = 'crop_071017_37CLNB_'; %image file name base
movieParam.firstImageNum = 1; %number of first image in movie
movieParam.lastImageNum = 40; %number of last image in movie
movieParam.digits4Enum = 4; %number of digits used for frame enumeration (1-4).

%% detection parameters

%Camera bit-depth
detectionParam.bitDepth = 16;

%The standard deviation of the point spread function is defined
%as 0.21*(emission wavelength)/(numerical aperture). If the wavelength is
%given in nanometers, this will be in nanometers. To convert to pixels,
%divide by the pixel side length (which should also be in nanometers).
detectionParam.psfSigma = 1.7;

%Number of frames before and after a frame for time averaging
%For no time averaging, set to 0
detectionParam.integWindow = 1;

%Alpha-value for initial detection of local maxima
detectionParam.alphaLocMax = 0.1;

%Maximum number of iterations for PSF sigma estimation for detected local
%maxima
%To use the input sigma without modification, set to 0
detectionParam.numSigmaIter = 0;

%1 to attempt to fit more than 1 kernel in a local maximum, 0 to fit only 1
%kernel per local maximum
%If psfSigma is < 1 pixel, set doMMF to 0, not 1. There is no point
%in attempting to fit additional kernels in one local maximum under such
%low spatial resolution
detectionParam.doMMF = 1;

%Alpha-values for statistical tests in mixture-model fitting step
detectionParam.testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0);

%1 to visualize detection results, frame by frame, 0 otherwise. Use 1 only
%for small movies. In the resulting images, blue dots indicate local
%maxima, red dots indicate local maxima surviving the mixture-model fitting
%step, pink dots indicate where red dots overlap with blue dots
detectionParam.visual = 0;

%% save results

saveResults.dir = 'C:\Users\Shaurya\Desktop\image_processing\nmeth\example\'; %directory where to save input and output
saveResults.filename = 'testDetection.mat'; %name of file where input and output are saved

%% run the detection function

[movieInfo,exceptions,localMaxima,background,psfSigma] = ...
    detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults);

%% Output variables

%The important output variable is movieInfo, which contains the detected
%particle information

%for a movie with N frames, movieInfo is a structure array with N entries.
%Every entry has the fields xCoord, yCoord, zCoord (if 3D) and amp.
%If there are M features in frame i, each one of these fields in
%moveiInfo(i) will be an Mx2 array, where the first column is the value
%(e.g. x-coordinate in xCoord and intensity in amp) and the second column
%is the standard deviation.

