function [compTracksAlt,compTracksDef] = convTrackFormatDefault2Alt(compTracks,simultFlag)
%CONVTRACKFORMATDEFAULT2ALT converts compound tracks from default format to alternative format
%
%SYNOPSIS compTracksAlt = convTrackFormatDefault2Alt(compTracks)
%
%INPUT  compTracks   : Compound tracks, in the format of tracksFinal as
%                      output by trackCloseGapsKalman.
%       simultFalg   : 1 to call the function
%                      "removeSimultaneousMergeSplit" at the beginning of
%                      this function, 0 otherwise.
%                      Optional. Default: 0.
%                      HOWEVER, THIS FUNCTION ASSUMES THAT SIMULTANEOUS
%                      MERGES/SPLITS HAVE BEEN TAKEN CARE OF, AND RESULTS
%                      ARE NOT GUARANTEED TO BE CORRECT IF THIS HAS NOT
%                      BEEN DONE.
%
%OUTPUT compTracksAlt: Compound tracks in alternative format, where a
%                      merge/split leads to new track segments, without
%                      continuation of the merging/splitting segments
%       compTracksDef: Compound tracks in default format. Same as input if
%                      simultFlag = 0. If simultFlag = 1, then these are
%                      the tracks after handling simultaneous merges and
%                      splits.
%
%
%
%Khuloud Jaqaman, February 2009, May 2019
% Luciana de Oliveira, April 2019-add conditions for multiple events at the
% same time. MODIFIED HEAVILY by KJ in May 2019.
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

%get number of compound tracks
numTracks = length(compTracks);

if nargin < 2
    simultFlag = 0;
end

%take care of simultaneous merging and splitting events
if simultFlag
    compTracks = removeSimultaneousMergeSplit(compTracks);
end

%% Output

%initialization
compTracksAlt = compTracks;

if nargout >= 2
    compTracksDef = compTracks;
end

%% Conversion

%go over all compound tracks
for iTrack = 1 : numTracks
    
    %get this compound track's information
    seqOfEvents = compTracks(iTrack).seqOfEvents;
    tracksFeatIndx = compTracks(iTrack).tracksFeatIndxCG;
    tracksCoordAmp = compTracks(iTrack).tracksCoordAmpCG;
    
    %KJ 190508: this is remnants of the past - no longer used, but kept to avoid too
    %many code changes
    doubleFreq = 0;
    
    seqOfEventsAlt = seqOfEvents;
    tracksFeatIndxAlt = tracksFeatIndx;
    %082114 - Robel Yirdaw
    %Convert sparse matrix back to full form to avoid bottleneck at
    %lines 169 - 172.
    if (issparse(tracksCoordAmp))
        tracksCoordAmpAlt = full(tracksCoordAmp);
    else
        tracksCoordAmpAlt = tracksCoordAmp;
    end
    
    %shift time in seqOfEventsAlt to make the track start at frame 1 (or 0.5
    %if doubling sampling frequency)
    seqOfEventsAlt(:,1) = seqOfEventsAlt(:,1) - seqOfEvents(1,1) + 1;
    
    %get number of segments and events in default format
    numSegmentsDef = size(tracksFeatIndx,1);
    numSegmentsAlt = numSegmentsDef;
    numEventsDef = size(seqOfEvents,1);
    
    %find all merging and splitting events
    msEvents = find(~isnan(seqOfEvents(:,4)));
    numMSEvents = length(msEvents);
    
    %vector to indicate whether M/S event has been taken care of already
    %needed for the case of simultaneous merges or splits
    event2do = ones(numMSEvents,1);
    
    %reserve memory for matrix storing segment correspondence
    alt2defSegCorr = zeros(numMSEvents,2);
    
    %transpose matrices to speed up calculations
    tracksFeatIndxAlt = tracksFeatIndxAlt';
    tracksCoordAmpAlt = tracksCoordAmpAlt';
    
    %go over merging and splitting events
    for iEventTmp = 1 : numMSEvents
        
        if event2do(iEventTmp)
            
            %get event index
            iEvent = msEvents(iEventTmp);
            
            %get time and type of event, and "original" segment that got merged
            %with or split from
            eventTime = seqOfEventsAlt(iEvent,1);
            eventType = seqOfEventsAlt(iEvent,2);
            originalSegment = seqOfEventsAlt(iEvent,4);
            
            %KJ 190508: get all events of the same type involving the same
            %original segment at the same time
            rowsEvent = find( seqOfEventsAlt(msEvents,1)==eventTime & ...
                seqOfEventsAlt(msEvents,2)==eventType & ...
                seqOfEventsAlt(msEvents,4)==originalSegment );
            
            %indicate that all of these events will be taken care of
            event2do(rowsEvent) = 0;
            
            %KJ 190508: also make sure that there are no events of the other
            %type with the same original segment at the same time - this
            %situation should be taken care of first using the
            %removeSimultaneousMergeSplit code
            rowsOtherEvent = find( seqOfEventsAlt(msEvents,1)==eventTime & ...
                seqOfEventsAlt(msEvents,2)==setdiff([1 2],eventType) & ...
                seqOfEventsAlt(msEvents,4)==originalSegment , 1);
            if ~isempty(rowsOtherEvent)
                error(['convTrackFormatDefault2Alt: Please eliminate ' ...
                    'simultaneous merges and splits first using the function ' ...
                    'removeSimultaneousMergeSplit before converting to ' ...
                    'alternative format!'])
            end
            
            %add 1 to number of segments
            numSegmentsAlt = numSegmentsAlt + 1;
            
            %store segment correspondence
            segmentCorrespond = originalSegment;
            while segmentCorrespond > numSegmentsDef
                segmentCorrespond = alt2defSegCorr(alt2defSegCorr(:,1)==segmentCorrespond,2);
            end
            alt2defSegCorr(rowsEvent,:) = repmat([numSegmentsAlt segmentCorrespond(1)],length(rowsEvent),1);
            
            %separate chunk after merging or splitting from original segment ...
            
            %update feature connectivity matrix
            tracksFeatIndxAlt((doubleFreq+1)*eventTime:end,numSegmentsAlt) = tracksFeatIndxAlt(...
                (doubleFreq+1)*eventTime:end,originalSegment);
            tracksFeatIndxAlt(1:(doubleFreq+1)*eventTime-1,numSegmentsAlt) = 0;
            tracksFeatIndxAlt((doubleFreq+1)*eventTime:end,originalSegment) = 0;
            
            %update matrix of coordinates and amplitudes
            tracksCoordAmpAlt(8*((doubleFreq+1)*eventTime-1)+1:end,numSegmentsAlt) = ...
                tracksCoordAmpAlt(8*((doubleFreq+1)*eventTime-1)+1:end,originalSegment);
            tracksCoordAmpAlt(1:8*((doubleFreq+1)*eventTime-1),numSegmentsAlt) = NaN;
            tracksCoordAmpAlt(8*((doubleFreq+1)*eventTime-1)+1:end,originalSegment) = NaN;
            
            %update sequence of events based on event type
            switch eventType
                
                case 1 %if it's a split
                    
                    %add new event for splitting of original segment into the
                    %new separated segment
                    seqOfEventsAlt(end+1,:) = [eventTime eventType numSegmentsAlt ...
                        originalSegment];
                    
                case 2 %if it's a merge
                    
                    %update index of segment merged into
                    seqOfEventsAlt(msEvents(rowsEvent),4) = numSegmentsAlt;
                    
                    %add new event for merging of original segment into the new
                    %separated segment
                    seqOfEventsAlt(end+1,:) = [eventTime eventType ...
                        originalSegment numSegmentsAlt];
                    
                    %update feature connectivity and coordinate matrices if
                    %sampling frequency is doubled
                    if doubleFreq
                        
                        %get index of merging segment
                        mergingSegment = seqOfEventsAlt(iEvent,3);
                        
                        %update its feature connectivity matrix and matrix of
                        %coordinates and amplitudes
                        tracksFeatIndxAlt((doubleFreq+1)*eventTime-1,mergingSegment) = ...
                            tracksFeatIndxAlt((doubleFreq+1)*eventTime-2,mergingSegment);
                        tracksCoordAmpAlt(8*((doubleFreq+1)*eventTime-1)-7:...
                            8*((doubleFreq+1)*eventTime-1),mergingSegment) = tracksCoordAmpAlt(...
                            8*((doubleFreq+1)*eventTime-1)-15:...
                            8*((doubleFreq+1)*eventTime-1)-8,mergingSegment);
                        
                    end
                    
            end %(switch eventType)
            
            %if the original segment is involved in later events,
            %replace it with the new additional segment
            %KJ 190508: to take care of simultaneous merges/splits, start
            %from the row in seqOfEvents right after the simultaneous
            %events
            indx = find(seqOfEventsAlt(msEvents(rowsEvent(end))+1:numEventsDef,3)==originalSegment);
            seqOfEventsAlt(msEvents(rowsEvent(end))+indx,3) = numSegmentsAlt;
            indx = find(seqOfEventsAlt(msEvents(rowsEvent(end))+1:numEventsDef,4)==originalSegment);
            seqOfEventsAlt(msEvents(rowsEvent(end))+indx,4) = numSegmentsAlt;
            
        end %(if event2do(iEventTmp))
        
    end %(for iEventTmp = 1 : numMSEvents)
    
    %transpose matrices back
    tracksFeatIndxAlt = tracksFeatIndxAlt';
    tracksCoordAmpAlt = tracksCoordAmpAlt';
    
    %sort sequence of events in ascending chronological order
    seqOfEventsAlt = sortrows(seqOfEventsAlt,1);
    
    %get unique event times
    eventTimes = unique(seqOfEventsAlt(:,1));
    
    %go over unique event times
    for iEvent = 1 : length(eventTimes)
        
        %find events happening at this event time
        indxEventsAtEventTime = find(seqOfEventsAlt(:,1)==eventTimes(iEvent));
        numEventsAtEventTime = length(indxEventsAtEventTime);
        
        %if there are more than 2 [1=initiation or termination (OK case);
        %2 = merge or split (OK case); 3 = merge or split AND initiation or
        %termination (possible "problem" case), etc. (all possible
        %"problem" cases)]
        if numEventsAtEventTime > 2
            
            %sort events such that splits/starts come before merges/ends
            seqOfEventsTmp = sortrows(seqOfEventsAlt(indxEventsAtEventTime,:),2);
            
            %extract rows documenting starts, splits, ends and merges
            seqOfEventsTmpStart = seqOfEventsTmp(seqOfEventsTmp(:,2)==1&...
                isnan(seqOfEventsTmp(:,4)),:);
            seqOfEventsTmpSplit = seqOfEventsTmp(seqOfEventsTmp(:,2)==1&...
                ~isnan(seqOfEventsTmp(:,4)),:);
            seqOfEventsTmpEnd = seqOfEventsTmp(seqOfEventsTmp(:,2)==2&...
                isnan(seqOfEventsTmp(:,4)),:);
            seqOfEventsTmpMerge = seqOfEventsTmp(seqOfEventsTmp(:,2)==2&...
                ~isnan(seqOfEventsTmp(:,4)),:);
            
            %make sure that every merge/split is documented on two
            %consecutive rows
            seqOfEventsTmpSplit = sortrows(seqOfEventsTmpSplit,4);
            seqOfEventsTmpMerge = sortrows(seqOfEventsTmpMerge,4);
            
            %put documentation back in sequence of events
            seqOfEventsAlt(indxEventsAtEventTime,:) = [seqOfEventsTmpStart; ...
                seqOfEventsTmpSplit; seqOfEventsTmpMerge; seqOfEventsTmpEnd];
            
        end
        
    end
    
    %add back time offset to get real times
    seqOfEventsAlt(:,1) = seqOfEventsAlt(:,1) + seqOfEvents(1,1) - 1;
    
    %store in output structure
    %082014 - Robel Yirdaw
    %Save tracksCoordAmpCG as a sparse matrix if originally sparse
    compTracksAlt(iTrack).seqOfEvents = seqOfEventsAlt;
    compTracksAlt(iTrack).tracksFeatIndxCG = tracksFeatIndxAlt;
    if (issparse(tracksCoordAmp))
        tracksCoordAmpAlt(isnan(tracksCoordAmpAlt)) = 0;
        compTracksAlt(iTrack).tracksCoordAmpCG = sparse(tracksCoordAmpAlt);
    else
        compTracksAlt(iTrack).tracksCoordAmpCG = tracksCoordAmpAlt;
    end

    %KJ 190508: Not sure this is necessary; keeping it for now, just in case
    %     % because of the multiple interactions with the same segment and same time,
    %     % it is possible to have empty rows at the end of the conversion, so remove
    %     % the empty rows
    %     rowsGood=sum(alt2defSegCorr,2)~=0;
    %     alt2defSegCorr=alt2defSegCorr(rowsGood,:);
    
    compTracksAlt(iTrack).alt2defSegmentCorrespond = alt2defSegCorr;
    
end %(for iTracks = 1 : numTracks)



%% ~~~ the end ~~~
