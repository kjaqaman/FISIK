function getImageStackFromTracksDirect_New(tracksSim,timeInfo,intensityInfo,spaceInfo,saveInfo)
% GETIMAGESTACKFROMTRACKSDIRECT generates images from simulated compTracks 
%
%   SYNPOSIS: getImageStackFromTracksDirect_New(tracksSim,timeInfo,intensityInfo,spaceInfo,saveInfo)
%
%   Input:
%
%       tracksSim     : Compound tracks, e.g. from field in receptorInfoLabeled
%
%       timeInfo      : Structure with fields:
%
%           .simTimeStep    : Simulation time step (s).
%           .sampleStep     : Sampling time step, in same units as
%                             timeStep. Mostly relevant for
%                             simulated data where simulation time step
%                             might be 0.01 s but sampling time step of
%                             interest is e.g. 0.1 s.
%           .cutOffTime     : Total time of simulation (s) to be considered.
%
%       intensityInfo : Structure with fields:
%
%           .bgav           : average background level. In grayscale values,
%                             assuming a 16-bit camera.
%           .bgnoise        : std of backgound noise. In grayscale values,
%                             assuming a 16-bit camera.
%           .scaleFactor    : scale factor for particle intensities to
%                             convert to grayscale units, assuming a 16-bit
%                             camera.
%
%        spaceInfo    : Structure with fields:
%
%           .pixelSize      : pixel size in microns
%           .psfSigma       : width of point-spread function (in pixels),
%                             where sigma = (1/3)*(radius of Airy disc)
%           .imsize         : size of image [sx,sy]
%
%         saveInfo    : Structure with fields:
%
%           .saveVar      : variable that indicates whether or not tif files
%                           are saved to file (1/0)
%           .saveDir      : directory where the images will be saved.
%
%
%
% Modified September 2016, Luciana de Oliveira
% 
% Copyright (C) 2019, Jaqaman Lab - UTSouthwestern
%
% This file is part of FISIK.
%
% FISIK is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% FISIK is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with FISIK.  If not, see <http://www.gnu.org/licenses/>.
%
%

%% Input

% time info
simTimeStep = timeInfo.simTimeStep;
sampleStep = timeInfo.sampleStep;
cutOffTime = timeInfo.cutOffTime;

%intensity info
bgav=intensityInfo.bgav;
bgnoise=intensityInfo.bgnoise;
scaleFactor=intensityInfo.scaleFactor;

%space info
pixelSize =spaceInfo.pixelSize;
psfSigma =spaceInfo.psfSigma;
imsize =spaceInfo.imsize;

%save info
saveVar =saveInfo.saveVar;
saveDir =saveInfo.saveDir;

%% Calculations

%convert tracks structure to matrix
trackedFeatureInfo = convStruct2MatIgnoreMS(tracksSim,0);

%Use a portion of the simulation if desired (truncating) given by the value
% of endIter,i.e, only the first total endIter=cutOffTime*simTimeStep will
% be considered for the analysis
endIter=cutOffTime/simTimeStep;
trackedFeatureInfo = trackedFeatureInfo(:,1:8*endIter);

%Scale up intensities
trackedFeatureInfo(:,4:8:end) = trackedFeatureInfo(:,4:8:end) * scaleFactor;

%rescaling in function of pixel size
trackedFeatureInfo(:,1:8:end) = trackedFeatureInfo(:,1:8:end)/pixelSize;
trackedFeatureInfo(:,2:8:end) = trackedFeatureInfo(:,2:8:end)/pixelSize;

%remove 0 columns
trackedFeatureInfo(:,3:8:end) = [];
trackedFeatureInfo(:,4:7:end) = [];
trackedFeatureInfo(:,4:6:end) = [];
trackedFeatureInfo(:,4:5:end) = [];
trackedFeatureInfo(:,4:4:end) = [];

%Replace zeros with NaNs
trackedFeatureInfo(trackedFeatureInfo == 0) = NaN;

%Subsample whole simulation every convStep iterations
convStep = 3*round(sampleStep/simTimeStep);
xCoord = trackedFeatureInfo(:,1:convStep:end);
yCoord = trackedFeatureInfo(:,2:convStep:end);
amps = trackedFeatureInfo(:,3:convStep:end);
[numTracks,numIters] = size(xCoord);

% Alleviate boundary effects on detection: add 10 pixels*psfSigma to
% particle positions; a few lines down image size will be expanded as
% well
xCoord = xCoord + 10*psfSigma;
yCoord = yCoord + 10*psfSigma;

%reformat track information following function input
trackInfo = reshape([xCoord;yCoord;amps],numTracks,3*numIters);

%Convert image size from micron square to pixels also add 10*psfSigma pixels
%to each image size to alleviate boundary effects on detection

imsize(1:2) = round(imsize/pixelSize) + 20*psfSigma;

% rad is the radius used for the generation of the Airy disc
% for each object, in increments of psfSigma

rad=3*psfSigma;

%call the function that will make the image and save it in saveDir
makeAiryImageFromMPM(trackInfo,bgav,bgnoise,psfSigma,imsize,rad,saveVar,[],saveDir);

end



