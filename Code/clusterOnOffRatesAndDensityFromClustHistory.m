function [rateOnPerClust,rateOffPerClust,densityPerClust,numClustForRateCalc,clustStats] = ...
    clusterOnOffRatesAndDensityFromClustHistory(clustHistoryMerged,infoSpaceTime)
%CLUSTERONOFFRATESANDDENSITYFROMCLUSTHISTORY calculates cluster on and off rates and densities from clusterHistory
%
%   SYNPOSIS: [rateOnPerClust,rateOffPerClust,densityPerClust,numClustForRateCalc,clustStats] = ...
%    clusterOnOffRatesAndDensityFromClustHistory(clustHistoryMerged,infoSpaceTime)
%
% INPUT:
%      clustHistoryMerged: Second output of
%                          clusterHistoryFromCompTracks_aggregState_motion.
%                          This could be the aggregate of multiple (N)
%                          movies/simulations.
%       infoSpaceTime    : Structure with fields:
%           .probDim        : Problem dimensionality.
%           .timeStep       : Time between frames/time points. In units of
%                             interest (e.g. s).
%               Either one of these two fields:
%           .areaSideLen    : Simulation/image side length values,
%                             which can be a single value or a value per
%                             side. In units of interest (e.g. um).
%               OR
%           .area           : Nx1 vector (N = number of movies contributing
%                             to clustHistoryMerged) storing area where
%                             molecules reside, in units of interest (e.g. pixels^2
%                             or um^2). This is needed because masked cell
%                             areas are not rectangular, adn different
%                             between movies. 
%               If structure contains both fields, they must match.
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
%   OUTPUT:
%       rateOnPerClust    :  A 1D array of calculated on rates for clusters
%                            of size 1, 2, 3, etc. Cluster of size 1 gets
%                            NaN, but value kept for ease of reference to
%                            larger clusters. Units: per unit time per (#
%                            molecules/unit area). Area unit = square of
%                            areaSideLen unit.
%       rateOffPerClust   :  A 1D array of calculated off rates for clusters
%                            of size 1, 2, 3, etc. Cluster of size 1 gets
%                            NaN, but value kept for ease of reference to
%                            larger clusters. Units: per unit time.
%       densityPerClust   :  A 1D array of calculated density for clusters
%                            of size 1, 2, 3, etc. Units: # molecules/unit
%                            area. Area unit = square of areaSideLen unit.
%       numClustForRateCalc: First column indicates number of clusters of
%                            each size used to calculate off rate.
%                            Second column indicates mean number of
%                            clusters of each size per iteration,
%                            indirectly used to calculate on rate.
%       clustStats        :  Output of clusterNumbersFromCompTracksFromClustHistory.
%
%   Khuloud Jaqaman, May 2015
%
%   Modified by Luciana de Oliveira, Oct 2016. To calculate the clustStats
%   directly from clustHistory.
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

%get sampling and cropping information
timeStep = infoSpaceTime.timeStep;
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
firstLastFrame = 1 + firstLastTP/timeStep;

%% Densities
%get cluster densities
clustStats = clusterNumbersFromCompTracksFromClustHistory(clustHistoryMerged,infoSpaceTime);

densityPerClust = mean(clustStats.clusterDensity,2);
clustCount = mean(clustStats.clusterCount,2);


%% RATES

% Modification 2016/11/17(LRO): To calculate the rates it is needed to do the
% following steps in clustHistoryMerged.

%KJ 190513: Remove events with unknown starts or ends
clustHistoryMergedChange = clustHistoryMerged(clustHistoryMerged(:,6)~=0 & clustHistoryMerged(:,7)~=0,:);

if ~isempty(firstLastFrame)
    clustHistoryMergedChange = clustHistoryMergedChange(clustHistoryMergedChange(:,3)>=firstLastFrame(1) ...
        & clustHistoryMergedChange(:,4)<=firstLastFrame(2),:);
end

%subsample time
clustHistoryMergedChange(:,3:4) = ceil(clustHistoryMergedChange(:,3:4)/convStep);
clustHistoryMergedChange(:,5) = clustHistoryMergedChange(:,4) - clustHistoryMergedChange(:,3);

%remove clusters which start and end at exactly the same time point -
%these would be not detectable with the sub-sampling time step
clustHistoryMergedChange = clustHistoryMergedChange(clustHistoryMergedChange(:,5)>0,:);

%convert from iterations/frames to real time units
clustHistoryMergedChange(:,3:5) = clustHistoryMergedChange(:,3:5) * sampleStep;

%get maximum cluster size
maxClusterSize = min([max(clustHistoryMergedChange(:,2)) length(densityPerClust)]);
rateOffPerClust = NaN(maxClusterSize,1);
rateOnPerClust = NaN(maxClusterSize,1);
numClustForRateCalc = [NaN(maxClusterSize,1) clustCount(1:maxClusterSize)];

%go over each cluster size > 1 and calculate off rate
for iSize = 2 : maxClusterSize
    
    %get lifetimes of clusters of current size, and whether they ended by
    %association or dissociation
    indxClust = clustHistoryMergedChange(:,2)==iSize&clustHistoryMergedChange(:,6)==2;
    clustLft = clustHistoryMergedChange(indxClust,5);%cluster life time
    clustEndType = clustHistoryMergedChange(indxClust,7);
    
    %calculate dissociation rate
    numClusters = length(clustEndType);
    rateOffPerClust(iSize) = ...
        (length(find(clustEndType==1))/numClusters) / mean(clustLft);
    
    %record number of clusters used for off rate calculation
    numClustForRateCalc(iSize,1) = numClusters;
    
end

%calculate on rates from off rates and densities (assumes steady state)
if maxClusterSize > 1
    rateOnPerClust(2:end) = rateOffPerClust(2:end) .* densityPerClust(2:maxClusterSize) ./ ...
        ( densityPerClust(1:maxClusterSize-1) * densityPerClust(1) );
end


%% ~~~ the end ~~~