function tracksReform = groupSegmentsComptracksAfterReformation( tracksIn)
%groupSegmentsComptracksAfterReformation will rearrange compTracks to
%have only the segments that interact
%
% SYNOPSIS  tracksReform = groupSegmentsComptracksAfterReformation( tracksIn)
%
% INPUT
%               tracksIn       :  compTracks after reformation
%
%
% OUTPUT
%          tracksReform     : tracks with the right compTracks, because
%          after reformation some segments can finish not being part of the
%          same compTrack anymore.
%
%  Luciana de Oliveira, February 2018.
% Modified the description of the function April, 2019.
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

% This function is divided in two parts:
% 1) identify which segments are part of the same compTrack
% 2) reform tracksFeatIndxCG, tracksCoordAmpCG and seqOfEvents.

% load compTracks
compTracks=tracksIn;
%initiate the reform compTrack
tracksReform=tracksIn;
% go over all compTracks
iTrackNew=length(compTracks);
% add a new field in reformTrack to take account of all the modifications
% done in the compTracks
tracksReform(1).oldTracksInfo=[];


%% gropup segments

for iTrack = 1: length(compTracks)
    %  for iTrack = 15;
    
    %initiate the flag for segments that were already checked
    
    segCheck=[];
    
    %     fprintf('\nProcessing track=%s ',int2str(iTrack));
    
    % load seqOfEvents and determine the number of
    % segments
    
    seqOfEvents = compTracks(iTrack).seqOfEvents;
    
    %load tracksFeatIndxCG and tracksCoordAmpCG
    tracksFeatIndxCG= compTracks(iTrack).tracksFeatIndxCG;
    tracksCoordAmpCG= compTracks(iTrack).tracksCoordAmpCG;
    
    %calculate number of segments
    numberOfSegments=size(tracksCoordAmpCG,1);
    
    %calculate length seqOfEvents
    
    lengthSeqOfEvents= size(seqOfEvents,1);
    
    % if the seqOfEvents have more than one segment
    
    if lengthSeqOfEvents > 2
        %% regroup segments
        % initiate the cell array that will contain the segments that are
        % part of the same group
        
        
        compTrackSegCell=cell(1,numberOfSegments);
        
        % check the interactions for each segment in this compTrack
        
        for segIndex=1:numberOfSegments
            
            % check if this segment was already checked
            
            if ~ismember(segIndex,segCheck)
                
                %  copy seqOfEvents, because I will remove the rows
                %  that are already checked
                
                seqOfEventsNew=seqOfEvents;
                
                %initiate a new group of segments
                compTrackSeg=[];
                
                % determine in which rows this segment appears
                rowsSegment= seqOfEventsNew(:,3)==segIndex|seqOfEventsNew(:,4)==segIndex;
                
                % determine the number of segments that are interacting
                segmentNumber= unique (seqOfEventsNew(rowsSegment,3:4));
                
                %removing NaNs
                segmentNumber(isnan(segmentNumber))=[];
                
                % save these segments to their group
                
                compTrackSeg=[compTrackSeg,segmentNumber];
                
                
                % remove these rows from seqOfEvents, I need to do
                % that, otherwise I have an infinite while loop
                
                seqOfEventsNew(rowsSegment,:)=[];
                
                clear rowsSegment
                
                % check for the segments that are part of this group which
                % are the other segments that they are interacting
                
                % initiate the flag for the segments interacting with the
                % segment but were already cheched in the previous loop.
                
                iSegChecked=[];
                
                while~isempty(segmentNumber)
                    
                    %get the first seg
                    
                    iSeg = segmentNumber(1);
                    
                    
                    % this analysis need to be done only for the segments
                    % that are not segIndex
                    
                    if iSeg~=segIndex
                        if iSeg~=segCheck
                            % determine in which rows this segment appears
                            rowsSegment= seqOfEventsNew(:,3)==iSeg|seqOfEventsNew(:,4)==iSeg;
                            
                            % determine the number of segments that are interacting
                            segmentNumberNew= unique (seqOfEventsNew(rowsSegment,3:4));
                            %removing NaNs
                            segmentNumberNew(isnan(segmentNumberNew))=[];
                            
                            %update segment number list
                            segExtra=setdiff(segmentNumberNew,segmentNumber);
                            segmentNumber=[segmentNumber;segExtra];
                            
                            % if this segment is still not part of the compTrack
                            if any(~ismember(segmentNumber,compTrackSeg))
                                
                                
                                % save this segments to group this segments in the compTrack
                                compTrackSeg=[compTrackSeg;segmentNumber];
                                
                            end
                        end
                    end
                    
                    %save the segments that where already checked from the
                    %list of segments
                    segCheck=[segCheck;setdiff(iSeg,segCheck )];
                    
                    % some segments can interact multiple times with the same segment. so
                    % update compTracks
                    compTrackSegCell{segIndex}=unique(compTrackSeg);
                    
                    %remove this segment from the list
                    segmentNumber=segmentNumber(2:end);
                    
                    
                end
            end
            
        end
        
        
        % remove the empty cells
        compTrackSegCell=compTrackSegCell(~cellfun('isempty',compTrackSegCell));
        
        %% reform compTrack
        
        
        for indexGroup=1:length(compTrackSegCell)
            
            %load the segments that are part of goup 1
            segmentGroup=compTrackSegCell{indexGroup};
            
            
            % the group that have segment 1 will stay in the original
            % position in compTrack, the other will be add at the end of the
            % compTrack
            
            % check if the group have the first segment
            if any (segmentGroup==1)
                iTrackCurr= iTrack;
            else
                %increase iTrackNew
                iTrackNew=iTrackNew+1;
                
                iTrackCurr= iTrackNew;
                
            end
            
            %% reform tracksFeatIndxCG and tracksCoordAmpCG
            
            % take only the segments that are part of the group segmentGroup
            
            tracksReform( iTrackCurr).tracksFeatIndxCG = tracksFeatIndxCG(segmentGroup,:);
            tracksReform( iTrackCurr).tracksCoordAmpCG =tracksCoordAmpCG(segmentGroup,:);
            
            %% reform seqOfevents
            
            % determine which are the rows that remain for this group
            rowsGood= ismember(seqOfEvents(:,3),segmentGroup)|ismember(seqOfEvents(:,4),segmentGroup);
            tracksReform(iTrackCurr).seqOfEvents=seqOfEvents(rowsGood,:);
            
            % replace number of segments
            if sum(rowsGood)==2
                tracksReform( iTrackCurr).seqOfEvents(:,3)=1;
            else
                seqOfEventsIndex= tracksReform(iTrackCurr).seqOfEvents(:,3:4);
                for i=1:length(segmentGroup)
                    seqOfEventsIndex(seqOfEventsIndex==segmentGroup(i))=i;
                end
                %update seqOfEvents
                tracksReform( iTrackCurr).seqOfEvents(:,3:4)=seqOfEventsIndex;
            end
            
            %% reform tracksFeatIndxCG and tracksCoordAmpCG to give only information for the frames that they appear
            
            % with the reform some tracks can start now at some point
            % different of frame 1, but we need the information of tracks
            % only after they appear and until they disappear
            
            %calculate if there is difference between the time that the
            %original and the new track initiates or finishes
            
            % calculate the initiation times
            
            timeStartOri=compTracks(iTrack).seqOfEvents(1,1);
            timeStartNew=tracksReform( iTrackCurr).seqOfEvents(1,1);
            
            % calculate the end times
            
            timeEndOri=compTracks(iTrack).seqOfEvents(end,1);
            timeEndNew=tracksReform( iTrackCurr).seqOfEvents(end,1);
            
            % calculate the difference between these times
            timeDiffStart= timeStartNew-timeStartOri;
            timeDiffEnd= timeEndOri-timeEndNew;
            
            % reform tracksFeatIndxCG and tracksCoordAmpCG
            
            
            tracksReform( iTrackCurr).tracksFeatIndxCG =tracksReform( iTrackCurr).tracksFeatIndxCG (:,timeDiffStart+1:end-timeDiffEnd);
            tracksReform( iTrackCurr).tracksCoordAmpCG =tracksReform( iTrackCurr).tracksCoordAmpCG(:,8*timeDiffStart+1:end-8*timeDiffEnd);
            
            
            % add a new structure field to the compTracks, that will take
            % account of all the alterations from the original compTrack
            tracksReform(iTrackCurr).oldTracksInfo={iTrack,segmentGroup};
            
        end
    else
        % It was only one track, just copy this info in the field oldTracksInfo
        tracksReform(iTrack).oldTracksInfo={iTrack,1};
    end
end
end



