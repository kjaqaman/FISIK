function [rateOnPerClust,rateOffPerClust,densityPerClust,...
    numClustForRateCalc,clustHistory,clustStats] = ...
    clusterOnOffRatesAndDensity(compTracksAggregState,infoSpaceTime)
%CLUSTERONOFFRATESANDDENSITY calculates interaction on and off rates and cluster/oligomer densities from dynamic trajectory data
%
%   SYNPOSIS: [rateOnPerClust,rateOffPerClust,densityPerClust,...
%               numClustForRateCalc,clustHistory,clustStats] = ...
%               clusterOnOffRatesAndDensity(compTracksAggregState,infoSpaceTime)
%
%   INPUT:
%       compTracksAggregState:
%                              Compound tracks as output by
%                              aggregStateFromCompIntensity or
%                              aggregStateFromCompTracksMIQP.
%                              THE INPUT MUST BE IN THE DEFAULT FORMAT.
%
%       infoSpaceTime    : Structure with fields
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
%   OUTPUT:
%
%       clustHistory   :  A 1D cell with rows = number of tracks in
%                         compTracks.
%                         Each entry contains a clusterHistory table for a
%                         track in compTracks. In this version cluster
%                         history is recorded for all clusters in the
%                         obervation time.
%                         clusterHistory is a 2D array with each row
%                         corresponding to an association or a dissociation
%                         event. The 9 colums give the following information:
%                         1) Track segment number in default format.
%                         2) Cluster size.
%                         3) Start time (same units as input timeStep etc).
%                         4) End time (same units as input timeStep etc).
%                         5) Lifetime (same units as input timeStep etc).
%                         6) Event that started the cluster
%                           (1 = dissociation, 2 = association, 0 = unknown/appearance).
%                         7) Event that ended the cluster
%                           (1 = dissociation, 2 = association, 0 = unknown/disappearance).
%                         8) Resulting cluster size.
%                         9) Track segment number in alternative format.
%                        10) Diffusion mode - IRRELEVANT HERE.
%                        11) Diffusion coefficient - IRRELEVANT HERE.
%
%        clustStats    :  A struct with the following fields
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
%
%   Khuloud Jaqaman, May 2015
%   Modified Luciana de Oliveira, April 2018
%   Modified Luciana de Oliveira, June 2018
%   Modified Luciana de Oliveira, October 2018
%   Modified Luciana de Oliveira, November 2018
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

%%
% calculate cluster history from the aggregState in the default format.

[clustHistory,clustHistoryMerged] = ...
    clusterHistoryFromCompTracks_aggregState_motion(compTracksAggregState);

% calculate rates and density

[rateOnPerClust,rateOffPerClust,densityPerClust,numClustForRateCalc,clustStats] = ...
    clusterOnOffRatesAndDensityFromClustHistory(clustHistoryMerged,infoSpaceTime);

%% ~~~ the end ~~~