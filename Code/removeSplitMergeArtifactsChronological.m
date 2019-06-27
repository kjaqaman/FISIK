function tracksReform = removeSplitMergeArtifactsChronological(tracksIn,thresholdTime)
%REMOVESPLITMERGEARTIFACTSCHRONOLOGICAL removes potentially artifactual merges and splits, resulting for instance from detection false positives
%
% SYNOPSIS tracksReform = removeSplitMergeArtifactsChronological(tracksIn,thresholdTime)
%
% INPUT
%          tracksIn       : Tracks in the format of the output of
%                           trackCloseGapsKalmanSparse (or
%                           u-track for GUI version).
%      thresholdTime      : Structure with the following fields
%                            threshSt2M :   minimum allowed start-to-merge
%                            time. Optional. Default=1.
%                            threshSp2M  :   minimum allowed split-to-merge
%                            time. Optional. Default=1.
%                            threshSp2E  :   minimum allowed split-to-end
%                            time. Optional. Default=0.
%
% OUTPUT
%          tracksReform   : Tracks after removing artifacts.
%
% Luciana de Oliveira, August 2017
% Modified Luciana de Oliveira, December 2017, August 2018
% Major bug fix (clean up of loops and if statement ends), Khuloud Jaqaman, April 2019
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modification LRO 2018/08/15- count the number of segments that were
% modified with reformat tracks
modificationCountMerge=0;
modificationCountSplit=0;
numberStart2Merge=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%get threshold information

if isfield(thresholdTime,'threshSt2M')
    threshSt2M = thresholdTime.threshSt2M;
else
    threshSt2M = 1;
end

if isfield(thresholdTime,'threshSp2M')
    threshSp2M = thresholdTime.threshSp2M;
else
    threshSp2M = 1;
end

if isfield(thresholdTime,'threshSp2E')
    threshSp2E = thresholdTime.threshSp2E;
else
    threshSp2E = 0;
end

% load compTracks
compTracks=tracksIn;
%initiate the reform compTrack
tracksReform=tracksIn;
% go over all compTracks
% iTrackNew=length(compTracks);
% add a new field in reformTrack to take account of all the modifications
% done in the compTracks
tracksReform(1).oldTracksInfo=[];

for iTrack = 1: length(compTracks)
    
    seqOfEventsIn = compTracks(iTrack).seqOfEvents;
    
    %calculate length seqOfEvents
    lengthSeqOfEvents= size(seqOfEventsIn,1);
    
    % if the seqOfEvents have more than one segment
    if lengthSeqOfEvents > 2
        
        %copy sequence of events
        seqOfEvents = seqOfEventsIn;
        
        %find indices for which a start of a new track is not due a birth and
        % an end is not due a death
        eventRows = find(~isnan(seqOfEvents(:,4)));
        
        % initiate the vector for the information of segments that were modified
        segInden=zeros(length(eventRows),1); % segIden is the identification of the
        % segment that was already analysed in modified in seqOfEvents
        
        %initiate the index that will take account the segments that are included
        %in segIden
        indexSeg=1;
        
        %go over all the events
        for indexEvent=1:length(eventRows)
            
            % to be sure of analysing only the rows that some event is happening
            iEvent=eventRows(indexEvent);
            
            % determine if the event was already analysed
            segFlag=find(segInden==iEvent, 1);
            
            if isempty(segFlag)
                
                % determine the kind of event
                kindOfEvent=seqOfEvents(iEvent,2);
                
                if kindOfEvent==2
                    
                    %% Merging
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Modification LRO 2018/08/15-need to include a counter
                    % for the number of start to merge
                    numberStart2Merge=numberStart2Merge+1;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    
                    
                    % merge part
                    
                    % for the merging we will have two possible artifacts:
                    % (a) if there is a split right after the merge
                    % (b) if the segment is just borm. This condition will only be tested  if
                    % there is no split after the merge.
                    
                    
                    %find the two merging segments and the time of merging
                    segment1 = seqOfEvents(iEvent,3);
                    segment2 = seqOfEvents(iEvent,4);
                    timeMerge = seqOfEvents(iEvent,1);
                    
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %                 check if the merge is for a segment that just appeared
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    %                 find when each of the segments started
                    iStart1 = seqOfEvents(:,3)==segment1 & ...
                        isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==1;
                    timeStart1 = seqOfEvents(iStart1,1);
                    iStart2 = seqOfEvents(:,3)==segment2 & ...
                        isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==1;
                    timeStart2 = seqOfEvents(iStart2,1);
                    
                    %             get the start-to-merge time for both segments
                    timeStart12Merge = timeMerge - timeStart1;
                    timeStart22Merge = timeMerge - timeStart2;
                    
                    %                 if either time is equal or smaller than the threshold, discard the merge
                    if any([timeStart12Merge timeStart22Merge]<=threshSt2M)
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Modification LRO 2018/08/15-include counting for
                        % the number of merges and splits modified by the
                        % reshape of segments
                        modificationCountMerge=modificationCountMerge+1;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        
                        
                        
                        
                        %                 replace the merge as the end of the segment
                        seqOfEvents(iEvent,4) = NaN;
                        seqOfEvents(iEvent,1) = timeMerge-1;
                        
                        
                        %                     save thes segment in the vector segInden
                        segInden(indexSeg)=iEvent;
                        
                        %                 increase indexSeg
                        indexSeg=indexSeg+1;
                        
                    end
                    
                    
                else
                    
                    %% Splitting
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % check if the split is followed by a merge
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    %find the two splitting segments and the time of splitting
                    segment1 = seqOfEvents(iEvent,3);
                    segment2 = seqOfEvents(iEvent,4);
                    timeSplit = seqOfEvents(iEvent,1);
                    
                    
                    %check whether these two segments merge with each other again
                    iMerge = find(any(seqOfEvents(:,3:4)==segment1,2) & ...
                        any(seqOfEvents(:,3:4)==segment2,2) & seqOfEvents(:,2)==2);
                    
                    %if they merge ...
                    if ~isempty(iMerge)
                        
                        
                        %calculate split to merge time
                        timeMerge = seqOfEvents(iMerge,1);
                        timeSplit2Merge = timeMerge - timeSplit;
                        
                        %if the split to merge time is lower or equal timeThreshold
                        if timeSplit2Merge <=threshSp2M %== 1%timeThreshold
                            
                            
                            %Replace the segment that split as a new birth
                            seqOfEvents(iEvent,4) = NaN;
                            % Replace the merge as the end of the segment
                            seqOfEvents( iMerge,4) = NaN;
                            
                            
                            % save these segments in the vector segInden
                            segInden(indexSeg)=iEvent;
                            segInden(indexSeg+1)=iMerge;
                            %increase indexSeg
                            indexSeg=indexSeg+2;
                        end % if the split is followed by a merge
                        
                    else
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % check if the split is followed by the end of the segment
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        %find when each of the segments ends
                        iEnd1 = seqOfEvents(:,3)==segment1 & ...
                            isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==2;
                        timeEnd1 = seqOfEvents(iEnd1,1);
                        iEnd2 = seqOfEvents(:,3)==segment2 & ...
                            isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==2;
                        timeEnd2 = seqOfEvents(iEnd2,1);
                        
                        %get the split-to-end time for both segments
                        timeSplit2End1 = timeEnd1 - timeSplit;
                        timeSplit2End2 = timeEnd2 - timeSplit;
                        
                        %if either time is equal to 0 frame, discard the split
                        if any([timeSplit2End1 timeSplit2End2]==threshSp2E)
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Modification LRO 2018/08/15-include counting for
                            % the number of merges and splits modified by the
                            % reshape of segments
                            modificationCountSplit=modificationCountSplit+1;
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            
                            %Replace the segment that split as a new birth
                            seqOfEvents(iEvent,4) = NaN;
                            
                            % save thes segment in the vector segInden
                            segInden(indexSeg)=iEvent;
                            
                            %increase indexSeg
                            indexSeg=indexSeg+1;
                            
                        end
                        
                    end %(if ~isempty(iMerge) ... else ...)
                    
                end %(if kindOfEvent==2 ... else ...)
                
            end %(if isempty(segFlag))
            
        end %(for indexEvent=1:length(eventRows))
        
        tracksReform(iTrack).seqOfEvents = seqOfEvents;
        
    end %(if lengthSeqOfEvents > 2)
    
end %(for iTrack = 1: length(compTracks))

%% regroup segments
tracksReform = groupSegmentsComptracksAfterReformation(tracksReform);

%% ~~~ the end ~~~


