function [clusterStats,clusterStatsMotion] = clusterNumbersFromCompTracksFromClustHistory(clustHistoryMerged,infoSpaceTime)
%CLUSTERNUMBERSFROMCOMPTRACKSFROMCLUSTHISTORY calculates number, fraction and density of each cluster size from clustHistory
%
%   SYNPOSIS: clusterStats = clusterNumbersFromCompTracksFromClustHistory(clustHistoryMerged,infoSpaceTime)
%
%   Input:
%       clustHistoryMerged : See description in clusterOnOffRatesAndDensityFromClustHistory
%
%       infoSpaceTime: Structure with fields:
%           .probDim        : Problem dimensionality.
%           .areaSideLen    : Simulation/image side length values,
%                             which can be a single value or a value per
%                             side. In units of interest (e.g. um).
%           .timeStep       : Time between frames/time points. In units of
%                             interest (e.g. s).
%           .sampleStep     : Sampling time step, in same units as
%                             timeStep. Mostly relevant for
%                             simulated data where simulation time step
%                             might be 0.01 s but sampling time step of
%                             interest is e.g. 0.1 s.
%                             Optional. If not input, then sampleStep =
%                             timeStep.
%           .firstLastTP    : Row vector of first and last time points to
%                             use for calculating rates and densities. In
%                             same units as timeStep.
%                             If only one value is input, it is taken as
%                             the last time point.
%                             If no value is input, then all time points
%                             are used.
%
%   Output:
%       clusterStats: a struct with the following fields
%           a) clusterCount (max cluster x num time points): number of
%              clusters per size per time point.
%           b) clusterFrac (max cluster x num time points): fraction
%              of clusters per size per time point.
%           c) clusterDensity (max cluster x num time points): density
%              of clusters per size per time point.
%           d) receptorCount (1 x num time points): number of receptors
%              per time point.
%           e) receptorDensity (1 x num time points): density of receptors
%              per time point.
%           f) largestClusterSize: the largest cluster size.
%           g) infoSpaceTime: the (input) structure infoSpaceTime.
%
%       clusterStatsMotion: Equivalent to clusterStats, but with numbers
%               broken down by motion type (diffusion mode) along the 3rd
%               dimension. This applies to fields clusterCount, clusterFrac
%               and clusterDensity. It has additional field motionFrac
%               that reports fraction of motion types per cluster size.
%
%   Khuloud Jaqaman, May 2015. Based on calcClusterStatsFromCompTracks.
%   Modified by Luciana de Oliveira, Oct 2016. To calculate the clustStats
%   directly from clustHistory.
%   Expanded by KJ in May 2019 to break down numbers by motion type.
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

%get basic space and time information
probDim = infoSpaceTime.probDim;
timeStep = infoSpaceTime.timeStep;

%get sampling and cropping information
if isfield(infoSpaceTime,'sampleStep')
    sampleStep = infoSpaceTime.sampleStep;
else
    sampleStep = timeStep;
end
convStep = round(sampleStep/timeStep);
if isfield(infoSpaceTime,'firstLastTP')
    firstLastTP = infoSpaceTime.firstLastTP;
    if length(firstLastTP)==1
        firstLastTP = [0 firstLastTP];
    end
else
    firstLastTP = [];
end
firstLastFrame = 1 + round(firstLastTP/sampleStep);

%Determine area
if isfield(infoSpaceTime,'areaSideLen')
    areaSideLen = infoSpaceTime.areaSideLen;
    numSideLenVals = length(areaSideLen);
    if (numSideLenVals == 1)
        simArea = areaSideLen ^ probDim;
    else
        simArea = prod(areaSideLen);
    end
else
    simArea = [];
end
if isfield(infoSpaceTime,'area')
    area = infoSpaceTime.area;
    if ~isempty(simArea) && any(area ~= simArea)
        error('clusterNumbersFromCompTracksFromClustHistory: area and areaSideLen in infoSpaceTime must be consistent if both are supplied!');
    end
    simArea = sum(area);
end


%% Calculations

%Get number of tracks and number of frames/iterations
[numOfEvents,numCol] = size(clustHistoryMerged);
numIters = round(max(clustHistoryMerged(:,4)));

%add columns to cluster history if motion info is not supplied
clustHistoryMerged = [clustHistoryMerged NaN(numOfEvents,11-numCol)];

%define default first and last frame if not input
if isempty(firstLastTP)
    firstLastFrame = [1 round(numIters/convStep)];
end

%get maximum cluster size
maxClustSize = max(clustHistoryMerged(:,2));

%get number of diffusion modes
clustHistoryMerged(isnan(clustHistoryMerged(:,10)),10) = 0;
numMode = max(clustHistoryMerged(:,10)) + 1;

%reserve memory using a very large cluster size to be safe
clusterCount = zeros(maxClustSize,numIters);
clusterCountMotion = zeros(maxClustSize,numIters,numMode);

%go over each track and count clusters at each time point
for iTrack =  1 : numOfEvents
    
    %get events from clustHistory
    eventTrack = clustHistoryMerged(iTrack,:);
    
    %get the time span that this track covers
    colStart = eventTrack(3);
    colEnd = eventTrack(4);
    
    %get the motion type
    diffMode = eventTrack(10);
    
    %count clusters and add to overall matrix
    clusterCount(eventTrack(2),colStart:colEnd) = ...
        clusterCount(eventTrack(2),colStart:colEnd) + 1;
    
    %count clusters per motion type and add to overall matrix
    clusterCountMotion(eventTrack(2),colStart:colEnd,numMode-diffMode) = ...
        clusterCountMotion(eventTrack(2),colStart:colEnd,numMode-diffMode) + 1;
    
end

%sub-sample time step and crop time as requested
clusterCount = clusterCount(:,1:convStep:end);
clusterCount = clusterCount(:,firstLastFrame(1):firstLastFrame(2));
numIters = size(clusterCount,2);
clusterCountMotion = clusterCountMotion(:,1:convStep:end,:);
clusterCountMotion = clusterCountMotion(:,firstLastFrame(1):firstLastFrame(2),:);

%calculate cluster fractions at each iteration
clusterFrac = clusterCount ./ repmat(sum(clusterCount),maxClustSize,1);
clusterFracMotion = clusterCountMotion ./ repmat(sum(clusterCountMotion,1),[maxClustSize 1 1]);
motionFracCluster = clusterCountMotion ./ repmat(sum(clusterCountMotion,3),[1 1 numMode]);

%calculate cluster density at each iteration
clusterDensity = clusterCount / simArea;
clusterDensityMotion = clusterCountMotion / simArea;

%calculate receptor density at each iteraction
receptorCount = sum( clusterCount .* repmat((1:maxClustSize)',1,numIters) );
receptorDensity = receptorCount / simArea;

%Save values for return
clusterStats.clusterCount = clusterCount;
clusterStats.clusterFrac = clusterFrac;
clusterStats.clusterDensity = clusterDensity;
clusterStats.receptorCount = receptorCount;
clusterStats.receptorDensity = receptorDensity;

clusterStats.largestClustSize = maxClustSize;
clusterStats.infoSpaceTime = infoSpaceTime;

clusterStatsMotion.clusterCount = clusterCountMotion;
clusterStatsMotion.clusterFrac = clusterFracMotion;
clusterStatsMotion.motionFrac = motionFracCluster;
clusterStatsMotion.clusterDensity = clusterDensity;

%% ~~~ the end ~~~

