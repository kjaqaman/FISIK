function  tracksReform = removeMerge2SplitClustHist(tracksIn,clustHistoryAll,probArtifactM2S)
% REMOVEMERGE2SPLITCLUSTHIST removes potentially artifactual merges and splits from single-molecule tracks, 
%
% SYNOPSIS
% tracksReform = removeMerge2SplitClustHist(tracksIn,clustHistoryAll,probArtifactM2S)
%
%
% INPUT
%               tracksIn       : Tracks in the format as output by trackCloseGapsKalman.
%               clustHistoryAll: Cluster history of the compound tracks in
%                                tracksIn, as output by 
%                                clusterHistoryFromCompTracks_aggregState_motion.
%               probArtifactM2S: Probability of a merge-to-split event to
%                                be an artifact, as a function of its
%                                duration. This is a vector with entry i
%                                for a merge-to-split of duration i frames.
%                                For any merge-to-split events with a
%                                longer duration, the probability of being
%                                an artifact is taken as 0.
%
% OUTPUT
%               tracksReform   : tracks after removing artifactual
%                                merge-to-split events.
%
%
% Code uses the information in clustHistory and the probability to be an
% artifact as a function of merge-to-split time.
%
% Luciana de Oliveira, December 2017.
% Modified February 2018.
% Updated and modified, Khuloud Jaqaman, August 2023
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

%load compTracks
compTracks=tracksIn;

%initiate the reform compTrack
tracksReform=tracksIn;

%append zeros to probArtifactM2S to cover the total number of frames
seqOfEventsAll = vertcat(tracksIn.seqOfEvents);
numFrames = max(seqOfEventsAll(:,1));
probVecSize = length(probArtifactM2S);
probArtifactM2S = [probArtifactM2S(:); zeros(max(0,numFrames-probVecSize),1)];


%% check where merges to splits are happening

for iTrack = 1: length(compTracks)
    
    % load seqOfEvents    
    seqOfEvents = compTracks(iTrack).seqOfEvents;
    
    %For experimental data, the first iteration (frame) in seqOfEvents
    %can be some value other than 1.  This occurs because an object
    %can appear or disappear at any time, while in simulation all of
    %the objects are present at frame 1. This affects how aggregState
    %is accessed since the frame # is supposed to correspond to a
    %column in aggregState. So, need to shift frame numbers in
    %seqOfEvents, if necessary, when accessing aggregState. The very
    %first frame # is the shift amount.
    firstFrame= seqOfEvents(1,1);
    
    % load tracksFeatIndxCG and tracksCoordAmpCG
    tracksFeatIndxCG= tracksReform(iTrack).tracksFeatIndxCG;
    tracksCoordAmpCG= tracksReform(iTrack).tracksCoordAmpCG;
    
    % load clustHistory
    clustHistory = clustHistoryAll{iTrack,1};
    
    %calculate length seqOfEvents
    lengthSeqOfEvents= size(seqOfEvents,1);
    
    % if the seqOfEvents has more than one segment
    
    if lengthSeqOfEvents > 2
        
        % identify in which rows merge to split events are happening
        indexMerge2Split=find(clustHistory(:,6)==2 & clustHistory(:,7)==1);
        
        % go over all the merge to split events
        for indexTime=1:length(indexMerge2Split)
            
            % take info from clustHistory
            timeMerge = clustHistory(indexMerge2Split(indexTime),3);
            timeSplit = clustHistory(indexMerge2Split(indexTime),4) + 1; %KJ 230828: bug fix (was missing +1)
            lifetime = clustHistory(indexMerge2Split(indexTime),5);
            segmentContinuous = clustHistory(indexMerge2Split(indexTime),1);
            
            % get random number to decide whether this observation will be considered an artifact
            randNumberArtifacts=rand;
            
            % decide if the observation is an artifact, based on the probability
            if randNumberArtifacts<probArtifactM2S(lifetime)
                
                % identify the rows where the merge is happining
                iMerge= seqOfEvents(:,1)==timeMerge & seqOfEvents(:,2)==2 & seqOfEvents(:,4)==segmentContinuous;
                
                %check where these two segments split
                iSplit = seqOfEvents(:,4)==segmentContinuous & seqOfEvents(:,2)==1 & seqOfEvents(:,1)==timeSplit;
                
                % merge and split segment number
                segmentMerge=seqOfEvents(iMerge,3);
                segmentSplit=seqOfEvents(iSplit,3);
                
                % determine merge and split rows
                rowMerge=find(iMerge);
                rowSplit=find(iSplit);
                rowMergeSplit=[rowMerge,rowSplit];
                
                % copy the positions during the merge to the merging segment
                % also copy intensity (and all other information), although
                % this will be adjusted in the next steps
                tracksFeatIndxCG(segmentMerge,timeMerge-firstFrame+1:timeSplit-firstFrame) = tracksFeatIndxCG(segmentContinuous,timeMerge-firstFrame+1:timeSplit-firstFrame);
                tracksCoordAmpCG(segmentMerge,8*(timeMerge-firstFrame)+1:8*(timeSplit-firstFrame)) = tracksCoordAmpCG(segmentContinuous,8*(timeMerge-firstFrame)+1:8*(timeSplit-firstFrame));
                
                % for the intensity, calculate the fraction of intensity distribution for each
                % segment. For that need to take the intensity before the merge and after the split, but if
                % some event happened before the merge or after the split, need to take only the
                % part of the segments that comes just before the merge event or just after the split event.
                
                timeBefAftEvent  = timeBeforeAfterEvent(seqOfEvents,segmentMerge,segmentSplit,timeMerge,segmentContinuous,timeSplit,firstFrame);
                
                %calculate the intensities before and after the merge (bug
                %fix KJ 230828)
                meanIntSegMerge =  nanmean([tracksCoordAmpCG(segmentMerge,8*(timeBefAftEvent(1)-firstFrame)+4:8:8*(timeMerge-firstFrame)),tracksCoordAmpCG(segmentSplit,8*(timeSplit-firstFrame)+4:8:8*(timeBefAftEvent(2)-firstFrame+1))]);
                meanIntSegCont  =  nanmean(tracksCoordAmpCG(segmentContinuous,[8*(timeBefAftEvent(3)-firstFrame)+4:8:8*(timeMerge-firstFrame),8*(timeSplit-firstFrame)+4:8:8*(timeBefAftEvent(4)-firstFrame+1)]));
                
                % intensity distribution
                intMergeProp = meanIntSegMerge/(meanIntSegMerge+meanIntSegCont);
                intContProp  = meanIntSegCont/(meanIntSegMerge+meanIntSegCont);
                
                % use intensity distribution to adjust intensities after breaking apart the merge
                tracksCoordAmpCG(segmentMerge,8*(timeMerge-firstFrame)+4:8:8*(timeSplit-firstFrame-1)+4) = intMergeProp.*tracksCoordAmpCG(segmentContinuous,8*(timeMerge-firstFrame)+4:8:8*(timeSplit-firstFrame-1)+4);
                tracksCoordAmpCG(segmentContinuous,8*(timeMerge-firstFrame)+4:8:8*(timeSplit-firstFrame-1)+4) = intContProp.*tracksCoordAmpCG(segmentContinuous,8*(timeMerge-firstFrame)+4:8:8*(timeSplit-firstFrame-1)+4);
                
                % paiwise tracks
                cellMergeSplit = {[segmentMerge segmentSplit rowMergeSplit]}; %KJ 230828: Converted to cell for input. Verify.
                [tracksFeatIndxCG,tracksCoordAmpCG] = pairwiseSegTracks(tracksFeatIndxCG,tracksCoordAmpCG,cellMergeSplit);
                
                % reform seqOfEvents and tracks
                [seqOfEvents, clustHistory] = reformSeqOfEventsFromMerge2split(seqOfEvents,cellMergeSplit,clustHistory); %KJ 230828: Converted to cell for input (original was outdated)
                
                % save tracks in the final format
                tracksReform(iTrack).tracksFeatIndxCG = tracksFeatIndxCG;
                tracksReform(iTrack).tracksCoordAmpCG = tracksCoordAmpCG;
                tracksReform(iTrack).seqOfEvents = seqOfEvents;
                
            end
        end
    end
    
end


%% regroup tracks

tracksReform = groupSegmentsComptracksAfterReformation(tracksReform);

end