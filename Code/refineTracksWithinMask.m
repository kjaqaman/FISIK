function cleanedTracksFinal = refineTracksWithinMask(tracksFinal,probArtifactMerge2Split,thresholdTime)
%REFINETRACKSWITHINMASK removes merge and split artifacts from tracks
%
%INPUT
%                   tracksFinal: structure with tracks as output by
%                                trackCloseGapsKalmanSparse.
%
%               probArtifactM2S: Probability of a merge-to-split event to
%                                be an artifact, as a function of its
%                                duration. This is a vector with entry i
%                                for a merge-to-split of duration i frames.
%                                For any merge-to-split events with a
%                                longer duration, the probability of being
%                                an artifact is taken as 0.
%
%thresholdTime = struct('threshSt2M',2,'threshSp2M',2,'threshSp2E',1)  [default values]
%                   .threshSt2M : If start-to-merge time <= threshold,
%                                 merge is discarded. 
%                                 Optional. Default = 2.
%                   .threshSp2M : If split-to-merge time <= threshold,
%                                 split and merge are discarded. 
%                                 Optional. Default = 2.
%                   .threshSp2E : If split-to-end time <= threshold,
%                                 split is discarded. 
%                                 Optional. Default = 1.
%                   Note that default values are not filled in here, but
%                   rather in the function
%                   "removeSplitMergeArtifactsChronological" called below.
%
%OUTPUT
%   cleanedTracksFinal: structure with tracks after removing artifacts.
%
% Khuloud Jaqaman, September 2023
%
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


%%

% Step 1: reform tracks with artifacts coming from simultaneous merges and splits
tracksReformSimul = removeSimultaneousMergeSplit(tracksFinal);

% Step 2: reshape the seqOfEvents when a segment split and merge with the same
% segment again. If the oldest segment merge with the newest this function will
% exchange segment index to maintain one long segment and a short that split and merge. 
tracksReformat = reformatSplitsMergesKeepLongerLivedSegment(tracksReformSimul);

% Step3: calculate aggregState, clustHistory and probability to be an artifact. 
% These inputs are needed to reform merges to splits 

% calculate clustering state

%The proper way, to get true clustering state, is the commented out lines.
% intensityQuantum = [0.002 0.0002]; %this is just an example
% tracksAggregState = aggregStateFromCompTracksMIQP(tracksReformat,intensityQuantum);

%But, for the sake of just getting merge to split times, we do not need the
%true clustering state. We just need the cluster lifetimes. This speeds up
%the code tremendously.
tracksAggregState = aggregStateFromCompIntensity(tracksReformat);

% calculate cluster history
clustHistoryAll = clusterHistoryFromCompTracks_aggregState_motion(tracksAggregState.defaultFormatTracks,tracksAggregState.alternativeFormatTracks);

%KJ 230901: Remnants of Luciana's script. Check, update, modify, use
% calculate probability to be an artifact
%  functionCalcProbArtifactsMerge2SplitTimes(saveRoot,rDDir,aPDir,dRDir,outDirNum,lRDir,infoSpaceTime)

% Step 4:  reform tracks with artifacts coming from finite merges and splits
% KJ 230830: Updated order of input to make it compatible with updated function
%probArtifactMerge2Split = [ones(2,1); zeros(298,1)]; %this is just an example
tracksReformM2S = removeMerge2SplitClustHist(tracksReformat,clustHistoryAll,probArtifactMerge2Split);

% Step 5: reform tracks with artifacts coming from birth to merge, split to end or split to merge. 
%thresholdTime = struct('threshSt2M',2,'threshSp2M',2,'threshSp2E',1); %this is just an example
cleanedTracksFinal = removeSplitMergeArtifactsChronological(tracksReformM2S,thresholdTime);

end
