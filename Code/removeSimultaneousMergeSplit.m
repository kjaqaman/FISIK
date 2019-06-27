function [tracksReform,modificationCount] = removeSimultaneousMergeSplit( tracksIn )
%REMOVESIMULTANEOUSMEGESPLIT resolves merges and splits that occur with the same segment at the same time
%
% SYNOPSIS [tracksReform,modificationCount] = removeSimultaneousMergeSplit( tracksIn )
%
% INPUT
%          tracksIn       : Tracks in the format of the output of
%                           trackCloseGapsKalmanSparse (or
%                           u-track for GUI version).
% OUTPUT
%          tracksReform   : Tracks after removing artifacts.
%          modificationCount :   Number of modified segments.
%
% Luciana de Oliveira, January 2018
% Modification March 2018
% Modification August 2018
% This function is divided in 4 parts:
%
% 1) Check if there are more than one event happening in the same frame;
% 2) Check if these events are happining with the same segments;
% 3) Check if the events are a combination of merges and splits.
% 4) Reform the seqOfEvents and compTracks. In this last step it will
% decouple the segments that are not part of the same compTrack anymore and
% correct starting times for the segments and segment numbers in the
% seqOfEvents.
%
%Khuloud Jaqaman, May 2019: Added code to take care of cases when A merges
%with B which simultaneously merges with C, and equivalently A splits from
%B which simultaneously splits from C. In this case, both A and B merge
%with C, or both A and B split from C. This makes more sense and simplifies
%the seqOfEvents to then remove simultaneous merges and splits.
%While working on code, I also fixed bug on line 445, where the merge
%and split pairs are documented.
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
% Modification LRO 2018/08/15-include counting for
% the number of merges and splits corrected
%in the case of multiple merges and splits, the
%number of segments modified will be the difference
%of the merges and splits
% start counts
modificationCount=0;

% load compTracks
compTracks=tracksIn;


%initiate the reform compTrack
tracksReform=tracksIn;

% go over all the tracks

for iTrack = 1: length(compTracks)
    
    %load sequence of events
    
    seqOfEvents= compTracks(iTrack).seqOfEvents;
    
    %For experimental data, the first iteration (frame) in seqOfEvents
    %can be some value other than 1.  This is occurs becuase an object
    %can appear or disappear at any time while in simulation, all of
    %the objects are present at frame 1. This affects how aggregState
    %is accessed since the frame # is supposed to correspond to a
    %column in aggregState. So, need to shift frame numbers in
    %seqOfEvents, if necessary, when accessing aggregState. The very
    %first frame # is the shift amount.
    firstFrame= seqOfEvents(1,1);
    
    
    %load tracksFeatIndxCG and tracksCoordAmpCG
    
    
    tracksFeatIndxCG= tracksReform(iTrack).tracksFeatIndxCG;
    tracksCoordAmpCG= tracksReform(iTrack).tracksCoordAmpCG;
    
    %calculate length seqOfEvents
    
    lengthSeqOfEvents= size(seqOfEvents,1);
    
    
    % if the seqOfEvents have more than two segments
    
    
    if lengthSeqOfEvents > 4
        
        %find indices for which a start of a new track is not due a birth and
        % an end is not due a death
        
        eventRows = ~isnan(seqOfEvents(:,4));
        
        %% Check if the events are happening in the same frame
        
        % find which are the times that appear more than two times in the same
        % seqOfEvents, because these are candidates to have simultaneous merge and
        % split events
        
        % identify the frames where the events are happening
        timeEvents=seqOfEvents(eventRows,1);
        
        %identify the unique times
        uniqueTimeEvents=unique(timeEvents);
        
        timeRepetition=uniqueTimeEvents(1<histc(timeEvents,uniqueTimeEvents));
        
        % if there are multiple events happening in the same frame
        while timeRepetition
            % check for each time that repeats if the events are simultaneous
            
            eventTime=timeRepetition(1);
            
            %KJ 190507: get all the merges at this time point
            mergeEvents = find( seqOfEvents(:,1)==eventTime & seqOfEvents(:,2)==2 & ~isnan(seqOfEvents(:,4)) );
            
            %KJ 190507: get all splits at this time point
            splitEvents = find( seqOfEvents(:,1)==eventTime & seqOfEvents(:,2)==1 & ~isnan(seqOfEvents(:,4)) );
            
            %KJ 190507: fix the case when a segment is merged with but at the same
            %time it itself merges with another segment
            segmentsMergedWith = seqOfEvents(mergeEvents,4);
            segmentsMergedWith = unique(segmentsMergedWith);
            for iSeg = 1 : length(segmentsMergedWith)
                indxMessyMerge = find(seqOfEvents(mergeEvents,3)==segmentsMergedWith(iSeg));
                if ~isempty(indxMessyMerge)
                    survivingSegment = seqOfEvents(mergeEvents(indxMessyMerge),4);
                    replaceRow = find(seqOfEvents(mergeEvents,4)==segmentsMergedWith(iSeg));
                    seqOfEvents(mergeEvents(replaceRow),4) = survivingSegment; %#ok<FNDSB>
                    replaceRow = find(seqOfEvents(splitEvents,4)==segmentsMergedWith(iSeg));
                    seqOfEvents(splitEvents(replaceRow),4) = survivingSegment; %#ok<FNDSB>
                end
            end
            
            %KJ 190507: fix the case when a segment is split from but at the same time
            %it itself splits from another segment
            segmentsSplitFrom = seqOfEvents(splitEvents,4);
            segmentsSplitFrom = unique(segmentsSplitFrom);
            for iSeg = 1 : length(segmentsSplitFrom)
                indxMessySplit = find(seqOfEvents(splitEvents,3)==segmentsSplitFrom(iSeg));
                if ~isempty(indxMessySplit)
                    originalSegment = seqOfEvents(splitEvents(indxMessySplit),4);
                    replaceRow = find(seqOfEvents(splitEvents,4)==segmentsSplitFrom(iSeg));
                    seqOfEvents(splitEvents(replaceRow),4) = originalSegment; %#ok<FNDSB>
                end
            end
            
            % initialize the cell for info merge and split for each event
            % time
            cellMergeSplitInfo={};
            
            %% check if the events are happening with the same segment
            
            
            % identify the possible segments with artifacts
            
            possibleArtifactsTime=~isnan(seqOfEvents(:,4))&seqOfEvents(:,1)==eventTime;
            
            %find if they are interacting with the same segment
            
            %             % modification LRO 2019/04/25, should consider repetitions for
            %             % 3 and 4 columns
            %
            %             segmentInteraction=seqOfEvents(possibleArtifactsTime,3:4);
            %             segmentInteraction=segmentInteraction(:);
            
            %KJ 190507: Because of the changes above (lines 109-140), we only need
            %to look at the 4th column here
            segmentInteraction = seqOfEvents(possibleArtifactsTime,4);
            
            %determine which ones are the same and how much times it appears
            uniqueSegment=unique(segmentInteraction);
            
            %count the repetitions
            segmentRepetition=histc(segmentInteraction,uniqueSegment);
            
            
            for indexSegment=1:length(uniqueSegment)
                
                
                %recalculate the possible artifacts for each segment
                possibleArtifactsTime=find(~isnan(seqOfEvents(:,4))&seqOfEvents(:,1)==eventTime);
                
                
                
                % if there are multiple events happening with the same segment
                
                if segmentRepetition(indexSegment)>1
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Modification LRO 2018/08/15-include counting for
                    % the number of merges and splits corrected
                    %in the case of multiple merges and splits, the
                    %number of segments modified will be the difference
                    %of the merges and splits
                    modificationCount=modificationCount+1;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % check if it is a merge or a split
                    
                    %chnage LRO 2019/04/18
                    % the seqOfEvents is updating in each step, because of
                    % that I need to calculate locally the mergeIndex
                    
                    %KJ 190507: As above, we only need to look at the 4th column here
                    
                    %                     mergeIndexTmp=find(seqOfEvents(possibleArtifactsTime,2)==2 & seqOfEvents(possibleArtifactsTime,3)==uniqueSegment(indexSegment));
                    %                     mergeIndexTmp1=find(seqOfEvents(possibleArtifactsTime,2)==2 & seqOfEvents(possibleArtifactsTime,4)==uniqueSegment(indexSegment));
                    %                     mergeIndex=[mergeIndexTmp;mergeIndexTmp1];
                    mergeIndex=find(seqOfEvents(possibleArtifactsTime,2)==2 & seqOfEvents(possibleArtifactsTime,4)==uniqueSegment(indexSegment));
                    
                    %                     splitIndexTmp=find(seqOfEvents(possibleArtifactsTime,2)==1 & seqOfEvents(possibleArtifactsTime,3)==uniqueSegment(indexSegment));
                    %                     splitIndexTmp1=find(seqOfEvents(possibleArtifactsTime,2)==1 & seqOfEvents(possibleArtifactsTime,4)==uniqueSegment(indexSegment));
                    %                     splitIndex=[splitIndexTmp;splitIndexTmp1];
                    splitIndex=find(seqOfEvents(possibleArtifactsTime,2)==1 & seqOfEvents(possibleArtifactsTime,4)==uniqueSegment(indexSegment));
                    
                    % pair merge and split, because there are cases where the number of
                    % merges and splits are not equal, so the solution is pairwise all the
                    % merges and splits and leave the "extra" events how they are.
                    
                    %%%%%Modification LRO 2019/04/16%%%%%%%%%%%%%%%%%
                    % For cases where there are repetition of the segment
                    % because of multiple merge or splits, but one of them
                    % is empty, I am adding the condition to verify if we
                    % have at least one merge of split for each condition
                    
                    if ~isempty(mergeIndex) && ~isempty(splitIndex)
                        
                        
                        %% reform seqOfevents- replace segment numbers
                        
                        % identify segment merge and split
                        
                        segmentMerge=seqOfEvents(possibleArtifactsTime(mergeIndex),3);
                        % LRO 2019/04/29, add condition the segment being repeated be in the third or fourth coluns:
                        %take the segment that is interacting with the
                        %                         % segment repeated
                        %                         segmentInter=segmentMerge~=uniqueSegment(indexSegment);
                        %                         segmentMerge=segmentMerge(segmentInter);
                        
                        
                        segmentSplit=seqOfEvents(possibleArtifactsTime(splitIndex),3);
                        % LRO 2019/04/29, add condition the segment being repeated be in the third or fourth coluns:
                        %take the segment that is interacting with the
                        %                         % segment repeated
                        %                         segmentInter=segmentSplit~=uniqueSegment(indexSegment);
                        %                         segmentSplit=segmentSplit(segmentInter);
                        
                        
                        % test if there are multiple merges and split events
                        % happening simultaneously
                        if length(mergeIndex)==1 && length(splitIndex)==1 % if there is only one merge and one split
                            
                            %find rows to be removed
                            rowMerge=possibleArtifactsTime(mergeIndex);
                            rowSplit=possibleArtifactsTime(splitIndex);
                            rowMergeSplit=[rowMerge,rowSplit];
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%% LRO 2019/04/30%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Try having a cell with the information I need
                            % and saving for correction later
                            cellMergeSplitInfo{indexSegment,1}=[segmentMerge,segmentSplit,rowMergeSplit];
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                        else % if there are multiple simultaneous merges and splits
                            % % % % % % % %                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % % % % % % % %                         Modification LRO 2018/08/15-include counting for
                            % % % % % % % %                         the number of merges and splits corrected
                            % % % % % % % %                         in the case of multiple merges and splits, the
                            % % % % % % % %                         number of segments modified will be the difference
                            % % % % % % % %                         of the merges and splits
                            % % % % % % % %                                                     modificationCount=modificationCount+abs(segmentMerge-segmentSplit);
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % pre-alocate the ratios for the calculation of
                            % costMatrix
                            ratioIntensity=zeros(length(segmentMerge),length(segmentSplit));
                            ratioMeanDisp=zeros(length(segmentMerge),length(segmentSplit));
                            
                            
                            % calculate number of merges and splits
                            numberMerges=length(segmentMerge);
                            numberSplits=length(segmentSplit);
                            
                            
                            % calculate all the possible combinations between merges and splits
                            
                            %%%% Modification LRO 2019/04/16, for the cases
                            %%%% where we have more merges than splits it is
                            %%%% not calculating the indices properly, so I am
                            %%%% adding an if condition to deal with that.
                            
                            [mergeIndices,splitIndices] = meshgrid(1:numberMerges,1:numberSplits);
                            
                            
                            % % % %                             %put indices together
                            indicesMS=[mergeIndices(:),splitIndices(:)];
                            
                            for indexMS=1: size(indicesMS,1)
                                
                                % calculate mean intensity and displacements for all
                                % merges and splits
                                
                                %intensity
                                
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % modification LRO 2018/03/16
                                % included the condition to calculate the
                                % intensity only in the interval between the last event before the
                                % merge and the first after the split.
                                
                                timeBefAftEvent  = timeBeforeAfterEvent(seqOfEvents,segmentMerge(indicesMS(indexMS,1)),segmentSplit(indicesMS(indexMS,2)),eventTime);
                                
                                % calculate mean intensity
                                
                                % take into account that track does not necessarily start at frame 1, in all the times
                                
                                % merge segment
                                if issparse(tracksCoordAmpCG)
                                    tracksCoordAmpCGTmp= tracksCoordAmpCG(segmentMerge(indicesMS(indexMS,1)),8*(timeBefAftEvent(1)-firstFrame)+4:8:8*(eventTime-firstFrame));
                                    meanIntMerge=mean(tracksCoordAmpCGTmp(tracksCoordAmpCGTmp~=0));
                                else
                                    meanIntMerge=nanmean(tracksCoordAmpCG(segmentMerge(indicesMS(indexMS,1)),8*(timeBefAftEvent(1)-firstFrame)+4:8:8*(eventTime-firstFrame)),2);
                                end
                                
                                %split segment
                                if issparse(tracksCoordAmpCG)
                                    tracksCoordAmpCGTmp=tracksCoordAmpCG(segmentSplit(indicesMS(indexMS,2)),8*(eventTime-firstFrame+1)+4:8:8*(timeBefAftEvent(2)-firstFrame+1));
                                    meanIntSplit=mean(tracksCoordAmpCGTmp(tracksCoordAmpCGTmp~=0));
                                else
                                    meanIntSplit=nanmean(tracksCoordAmpCG(segmentSplit(indicesMS(indexMS,2)),8*(eventTime-firstFrame+1)+4:8:8*(timeBefAftEvent(2)-firstFrame+1)),2);
                                end
                                
                                
                                
                                % take the x and y values
                                
                                mergeXvals=tracksCoordAmpCG(segmentMerge(indicesMS(indexMS,1)),8*(timeBefAftEvent(1)-firstFrame)+1:8:8*(eventTime-firstFrame));
                                mergeYvals=tracksCoordAmpCG(segmentMerge(indicesMS(indexMS,1)),8*(timeBefAftEvent(1)-firstFrame)+2:8:8*(eventTime-firstFrame));
                                
                                %calculate displacement
                                
                                dispMergeX=diff(mergeXvals,1,2);
                                dispMergeY=diff(mergeYvals,1,2);
                                
                                
                                % calculate mean displacement
                                
                                meanDispMerge=nanmean(sqrt(dispMergeX.^2+dispMergeY.^2),2);
                                
                                
                                % splits
                                
                                % take the x and y values
                                
                                splitXvals=tracksCoordAmpCG(segmentSplit(indicesMS(indexMS,2)), 8*(eventTime-firstFrame+1)+1:8:8*(timeBefAftEvent(2)-firstFrame+1));
                                splitYvals=tracksCoordAmpCG(segmentSplit(indicesMS(indexMS,2)), 8*(eventTime-firstFrame+1)+2:8:8*(timeBefAftEvent(2)-firstFrame+1));
                                %calculate displacement
                                
                                dispSplitX=diff(splitXvals,1,2);
                                dispSplitY=diff(splitYvals,1,2);
                                
                                % calculate mean displacement
                                
                                meanDispSplit=nanmean(sqrt(dispSplitX.^2+dispSplitY.^2),2);
                                
                                % calculate ratio
                                % check which mean have the largest value
                                
                                %intensity
                                
                                intLarger=max(meanIntMerge,meanIntSplit);
                                intSmaller=min(meanIntMerge,meanIntSplit);
                                
                                %displacement
                                
                                dispLarger=max(meanDispMerge,meanDispSplit);
                                dispSmaller=min(meanDispMerge,meanDispSplit);
                                
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % add condition for a null displacement
                                if dispSmaller==0
                                    dispSmaller=1;
                                    
                                end
                                %calculate ratio
                                ratioIntensity(indicesMS(indexMS,1),indicesMS(indexMS,2))=intLarger/intSmaller;
                                ratioMeanDisp(indicesMS(indexMS,1),indicesMS(indexMS,2))=dispLarger/dispSmaller;
                                
                            end % end for loop merge events
                            
                            
                            % for some cases where the split occour and it is
                            % followed by a end we have ratios with value NaN, and
                            % NaN is not a valid input for the LAP function. To
                            % avoid this problem, in the cases where we have
                            % multiple merges and split, we will take the maximum
                            % value of cost matrix and apply for the NaN value
                            
                            
                            
                            if  all(isnan(ratioMeanDisp(:))) % if any ratio displacement is NaN
                                
                                %if all the displacements are NaN we will take only
                                %the value of intensity for the comparison
                                
                                ratioMeanDisp(:)=1;
                                
                            elseif any(isnan(ratioMeanDisp(:)))
                                
                                indiceNaNDisp= isnan(ratioMeanDisp);
                                newValRatioDisp=max(ratioMeanDisp(:));
                                ratioMeanDisp(indiceNaNDisp)=newValRatioDisp;
                                
                            end
                            
                            % calculate the costMatrix
                            
                            costMatrix=ratioIntensity.*ratioMeanDisp;
                            
                            
                            % after all the combinations are calculated and saved in costMatrix, we can
                            % analyse which is the most likely to be a pair of merge and split
                            [linkMerge2Split, ~] = lap(costMatrix, [],[],1,[]);
                            
                            
                            %here I need to combine all the possibilities from merge and splits. If for
                            %example I have multiples merges but only one split, from the output of lap
                            %I need to have the identification of this split is going with this merge.
                            %maybe the best solution could be combine only the number of merges with
                            %the splits, given by the length of merges and splits
                            
                            goodPairsMerge=linkMerge2Split(1:numberMerges);
                            goodPairs=find(goodPairsMerge<=numberSplits);
                            
                            % put the pairs together
                            
                            %KJ 190507: Fix bug here - should stay as
                            %column vectors
                            %                             pairsFinal=[goodPairs',goodPairsMerge(goodPairs)'];
                            pairsFinal=[goodPairs goodPairsMerge(goodPairs)];
                            
                            % go over all the pairs and reform seqOfEvents and tracks
                            
                            for indexPair=1:size(pairsFinal,1)
                                
                                segMergeTemp=segmentMerge(pairsFinal(indexPair,1));% it gives the segments that is merging based in the list of segments that merge
                                segSplitTemp=segmentSplit(pairsFinal(indexPair,2));% it gives the segments that is splitting based in the list of segments that split
                                
                                % define which are the segments to be removed
                                % it  finds the position in the list of the artifacts to be removed
                                rowMergeSplit=[possibleArtifactsTime(mergeIndex(pairsFinal(indexPair,1))),possibleArtifactsTime(splitIndex(pairsFinal(indexPair,2)))];
                                
                                
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%% LRO 2019/04/30%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % Try having a cell with the information I need
                                % and saving for correction later
                                cellMergeSplitInfo{indexSegment,indexPair}=[segMergeTemp,segSplitTemp,rowMergeSplit];
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                
                                %
                                %                                 %% reformCompTracks
                                %                                   seqOfEvents = reformSeqOfEventsFromMerge2split(seqOfEvents,segMergeTemp,segSplitTemp,rowMergeSplit);
                                %                                 % pairwise segments
                                %
                                %                                 [tracksFeatIndxCG,tracksCoordAmpCG] = pairwiseSegTracks(tracksFeatIndxCG,tracksCoordAmpCG,segMergeTemp,segSplitTemp);
                                %
                                %                                 % save tracks in the final format
                                %                                 tracksReform(iTrack).tracksFeatIndxCG = tracksFeatIndxCG;
                                %                                 tracksReform(iTrack).tracksCoordAmpCG = tracksCoordAmpCG;
                                %                                 tracksReform(iTrack).seqOfEvents = seqOfEvents;
                                
                            end %(for indexPair=1:size(pairsFinal,1))
                            
                        end  %(if length(mergeIndex)==1 && length(splitIndex)==1 ... else ...)
                        
                    end %(if ~isempty(mergeIndex) && ~isempty(splitIndex))
                    
                end %(if segmentRepetition(indexSegment)>1)
                
            end %(for indexSegment=1:length(uniqueSegment))
            
            
            %%%%%%%%%%%%%%%% try reforming here%%%%%%%%%%%%%%%5
            %% reformCompTracks
            %%%%%%%%%%%%Modification LRO 2019/04/30%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %take only the elements of the cell that are not empty
            
            cellMergeSplitInfo=  cellMergeSplitInfo(~cellfun(@isempty, cellMergeSplitInfo));
            
            seqOfEvents = reformSeqOfEventsFromMerge2split(seqOfEvents,cellMergeSplitInfo);
            % pairwise segments
            
            [tracksFeatIndxCG,tracksCoordAmpCG] = pairwiseSegTracks(tracksFeatIndxCG,tracksCoordAmpCG,cellMergeSplitInfo);
            
            % save tracks in the final format
            tracksReform(iTrack).tracksFeatIndxCG = tracksFeatIndxCG;
            tracksReform(iTrack).tracksCoordAmpCG = tracksCoordAmpCG;
            tracksReform(iTrack).seqOfEvents = seqOfEvents;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %increase index
            
            timeRepetition=timeRepetition(2:end);
            
        end %(while timeRepetition)
        
    end %(if lengthSeqOfEvents > 4)
    
end %(for iTrack = 1: length(compTracks))

%% regroup tracks

tracksReform = groupSegmentsComptracksAfterReformation(tracksReform);

end

