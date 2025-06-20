function [rateOnPerClust,rateOffPerClust,densityPerClust,numClustForRateCalc,clustStats,...
    rateOnPerClustM,rateOffPerClustM,densityPerClustM,numClustForRateCalcM,clustStatsM,...
    rateOnPerClustF,rateOffPerClustF,densityPerClustF,numClustForRateCalcF,clustStatsF] = ...
    clusterOnOffRatesDensityFromClustHistAndAggregState(clustHistoryMerged,aggregStateMat,...
    infoSpaceTime,fracFollowerMerged)
%CLUSTERONOFFRATESDENSITYFROMCLUSTHISTANDAGGREGSTATE calculates cluster on and off rates and densities from clusterHistory and track aggregation state
%
%   SYNPOSIS: [rateOnPerClust,rateOffPerClust,densityPerClust,numClustForRateCalc,clustStats,...
%    rateOnPerClustM,rateOffPerClustM,densityPerClustM,numClustForRateCalcM,clustStatsM,...
%    rateOnPerClustF,rateOffPerClustF,densityPerClustF,numClustForRateCalcF,clustStatsF] = ...
%    clusterOnOffRatesDensityFromClustHistAndAggregState(clustHistoryMerged,aggregStateMat,...
%    infoSpaceTime,fracFollowerMerged)
%
% INPUT:
%      clustHistoryMerged: Second output of
%                          clusterHistoryFromCompTracks_aggregState_motion.
%                          This could be the concatenation of multiple (N)
%                          movies/simulations.
%      aggregStateMat    : Matrix of aggregation state of the tracks
%                          corresponding to clustHistoryMerged. Every row
%                          represents one track segment, while number of
%                          columns = number of frames in movie. For any
%                          particular set of tracks, this variable can be
%                          obtained by running "convStruct2MatIgnoreMS" on
%                          the default tracks output of
%                          aggregStateFromCompTracksMIQP.
%      infoSpaceTime     : Structure with fields:
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
%                             areas are not rectangular, and different
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
%       rateOnPerClustM   :  5D array where the first dimension has length 
%                            equal to the largest cluster size, the second
%                            has length one, and the other three have 
%                            length equal to the number of motion modes + 
%                            1. Gives on rates of clusters per cluster size
%                            of the newly formed cluster (along the first 
%                            index), according to the mode of the new 
%                            cluster (third index), dependent on the modes 
%                            of the two mergers, one of size 1 (index 
%                            four), and the other being one size smaller 
%                            than the newly formed cluster (index five).
%                            The second index (along the columns) is kept
%                            at length one for now, because we reserve the
%                            addition of the second column for error values
%                            that can be obtained in the script
%                            summarizeMergeSplitVsMotionStats.
%       rateOffPerClustM  :  Equivalent to rateOffPerClust, per motion
%                            mode. Number of columns = number of motion
%                            modes + 1, where modes are listed from
%                            fastest to slowest, and finally unknown.
%       densityPerClustM  :  Equivalent to densityPerClust, per motion
%                            mode. Number of columns = number of motion
%                            modes + 1, where modes are listed from
%                            fastest to slowest, and finally unknown.
%       numClustForRateCalcM: Equivalent to numClustForRateCalc, per motion
%                            mode (stored in 3rd dimension).
%       clustStatsM       :  Output of clusterNumbersFromCompTracksFromClustHistory.
%
%       rateOnPerClustF   :  All NaNs for now (not calculated).
%       rateOffPerClustF  :  Equivalent to rateOffPerClustM, but division
%                            is based on follower presence: Column 1 with, 
%                            Column 2 without.
%       densityPerClustF  :  Equivalent to densityPerClustM, but division
%                            is based on follower presence: Column 1 with, 
%                            Column 2 without.
%       numClustForRateCalcF: Equivalent to numClustForRateCalcF, but
%                            division is based on follower presence
%                            (stored in 3rd dimension).
%       clustStatsF       :  Output of clusterNumbersFromCompTracksFromClustHistory.
%
% Khuloud Jaqaman, April 2024. 
%
% Based on clusterOnOffRatesAndDensityFromClustHistory, but with added use of
% aggregState directly for counting number of molecules and clusters, to
% properly handle gap closing. The previous code, which uses only cluster
% history, results in over-counting of molecules and clusters, and thus
% over estimation of density.
%
% Copyright (C) 2025, Jaqaman Lab - UTSouthwestern
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

%patch up cluster history for motion info if needed
[numRow,numCol] = size(clustHistoryMerged);
clustHistoryMerged = [clustHistoryMerged NaN(numRow,11-numCol)];

%replace motion mode NaN with 0, get number of modes
clustHistoryMerged(isnan(clustHistoryMerged(:,10)),10) = 0;
numMode = max(clustHistoryMerged(:,10)) + 1;

%master-follower analysis
if nargin < 4 || isempty(fracFollowerMerged)
    fracFollowerMerged = [];
end 

%% Densities

%get cluster densities
[clustStats,clustStatsM,clustStatsF] = clusterNumbersFromClustHistAndAggregState(...
    clustHistoryMerged,aggregStateMat,infoSpaceTime,fracFollowerMerged);

densityPerClust = mean(clustStats.clusterDensity,2);
clustCount = mean(clustStats.clusterCount,2);

densityPerClustM = squeeze(mean(clustStatsM.clusterDensity,2));
clustCountM = mean(clustStatsM.clusterCount,2);

densityPerClustF = squeeze(mean(clustStatsF.clusterDensity,2));
clustCountF = mean(clustStatsF.clusterCount,2);

%% RATES

%KJ 190513: Remove events with unknown starts or ends
indxKeep = find( clustHistoryMerged(:,6)~=0 & clustHistoryMerged(:,7)~=0 );
clustHistoryMergedChange = clustHistoryMerged(indxKeep,:);
if ~isempty(fracFollowerMerged)
    fracFollowerMergedChange = fracFollowerMerged(indxKeep);
end

%keep only clusters within frame range of interest
if ~isempty(firstLastFrame)
    indxKeep = find( clustHistoryMergedChange(:,3)>=firstLastFrame(1) ...
        & clustHistoryMergedChange(:,4)<=firstLastFrame(2) );
    clustHistoryMergedChange = clustHistoryMergedChange(indxKeep,:);
    if ~isempty(fracFollowerMerged)
        fracFollowerMergedChange = fracFollowerMergedChange(indxKeep);
    end
end

%subsample in time if desired
if sampleStep ~= timeStep
    
    %mimic time subsampling
    clustHistoryMergedChange(:,3:4) = ceil(clustHistoryMergedChange(:,3:4)/convStep);
    clustHistoryMergedChange(:,5) = clustHistoryMergedChange(:,4) - clustHistoryMergedChange(:,3) + 1;
    
    %remove clusters which start and end at exactly the same time point -
    %these would be not detectable after sub-sampling
    indxKeep = find( clustHistoryMergedChange(:,5)>1 );
    clustHistoryMergedChange = clustHistoryMergedChange(indxKeep,:);
    if ~isempty(fracFollowerMerged)
        fracFollowerMergedChange = fracFollowerMergedChange(indxKeep);
    end
    
end

%convert from iterations/frames to real time units
clustHistoryMergedChange(:,3:5) = clustHistoryMergedChange(:,3:5) * sampleStep;

%get maximum cluster size
maxClusterSize = min([max(clustHistoryMergedChange(:,2)) length(densityPerClust)]);

%initialize variables
rateOffPerClust = NaN(maxClusterSize,1);
rateOnPerClust = NaN(maxClusterSize,1);
numClustForRateCalc = [NaN(maxClusterSize,1) clustCount(1:maxClusterSize)];
rateOffPerClustM = NaN(maxClusterSize,numMode);
rateOnPerClustM = NaN(maxClusterSize,numMode);
numClustForRateCalcM = [NaN(maxClusterSize,1,numMode) clustCountM(1:maxClusterSize,:,:)];
rateOffPerClustF = NaN(maxClusterSize,2);
rateOnPerClustF = NaN(maxClusterSize,2);
numClustForRateCalcF = [NaN(maxClusterSize,1,2) clustCountF(1:maxClusterSize,:,:)];

%define minimum number of clusters needed for rate calculation
minNumClust = 20;

%go over each cluster size > 1 and calculate off rate
for iSize = 2 : maxClusterSize
    
    %get lifetimes of clusters of current size formed by association, and
    %record whether they ended by association or dissociation
    indxClust = find(clustHistoryMergedChange(:,2)==iSize&clustHistoryMergedChange(:,6)==2);
    clustLft = clustHistoryMergedChange(indxClust,5);%cluster life time
    clustEndType = clustHistoryMergedChange(indxClust,7);
    
    %calculate dissociation rate (only if minimum number of clusters is
    %satisfied)
    numClusters = length(indxClust);
    if numClusters >= minNumClust
        rateOffPerClust(iSize) = ...
            (length(find(clustEndType==1))/numClusters) / mean(clustLft);
    end
    
    %record number of clusters used for off rate calculation
    numClustForRateCalc(iSize,1) = numClusters;
    
    %repeat the above but based on motion mode
    clustMotionMode = clustHistoryMergedChange(indxClust,10);
    for iMode = 1 : numMode
        indxMotion = find(clustMotionMode==(numMode-iMode));
        numClustersM = length(indxMotion);
        if numClustersM >= minNumClust
            rateOffPerClustM(iSize,iMode) = ...
                (length(find(clustEndType(indxMotion)==1))/numClustersM) / mean(clustLft(indxMotion));
        end
        numClustForRateCalcM(iSize,1,iMode) = numClustersM;
    end
    
    %repeat the above but based on follower presence
    if ~isempty(fracFollowerMerged)
        followerPresent = double(fracFollowerMergedChange(indxClust) > 0);
        for iFoll = 1 : 2
            indxFollow = find(followerPresent==(2-iFoll));
            numClustersF = length(indxFollow);
            if numClustersF >= minNumClust
                rateOffPerClustF(iSize,iFoll) = ...
                    (length(find(clustEndType(indxFollow)==1))/numClustersF) / mean(clustLft(indxFollow));
            end
            numClustForRateCalcF(iSize,1,iFoll) = numClustersF;
        end
    end
    
end

%calculate on rates from off rates and densities (assuming steady state)
%note that if number of clusters is not sufficient to calculate off rate,
%then on rate will also be NaN
if maxClusterSize > 1
    
    rateOnPerClust(2:end) = rateOffPerClust(2:end) .* densityPerClust(2:maxClusterSize) ./ ...
        ( densityPerClust(1:maxClusterSize-1) * densityPerClust(1) );
    
    % Compute on rate based on motion mode (ZMalik 2020/03/10)
    rateOnPerClustM = NaN(maxClusterSize,1,numMode,numMode,numMode);
    for iSize = 2:maxClusterSize
        for iNewClustMode = 1:numMode
            for iOldReceptMode = 1:numMode
                for iOldClustMode = 1:numMode
                    rateOnPerClustM(iSize,1,iNewClustMode,...
                        iOldReceptMode,iOldClustMode) = ...
                        (rateOffPerClustM(iSize,iNewClustMode)*densityPerClustM(iSize,iNewClustMode))...
                        /(densityPerClustM(iSize-1,iOldClustMode)*densityPerClustM(1,iOldReceptMode));                   
                end
            end
        end
    end
    
    % 20200706 ZMALIK: Multiply dimer on rates by a factor of two
    rateOnPerClust(2) = rateOnPerClust(2)*2;
    rateOnPerClustM(2,:,:,:,:) = rateOnPerClustM(2,:,:,:,:)*2;
    
end


%% ~~~ the end ~~~