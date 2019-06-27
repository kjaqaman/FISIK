function [clustHistoryAll,clustHistoryMerged] = ...
    clusterHistoryFromCompTracks_aggregState_motion(tracksAggregStateDef,tracksAggregStateAlt,diffusionAnalysisAlt)
%CLUSTERHISTORYFROMCOMPTRACKS_AGGREGSTATE determines the size and lifetime of all clusters that formed during a simulation
%
%   The function uses the information contained in seqOfEvents and
%   aggregState, two fields from the output of aggregStateFromCompTracks,
%   in the field defaultFormatTracks.
%
%   SYNOPSIS: [clustHistoryAll,clustHistoryMerged] = ...
%    clusterHistoryFromCompTracks_aggregState(tracksAggregStateDef,infoTime)
%
%   INPUT:
%       tracksAggregStateDef:
%                         the structure of track information including
%                         aggregState in default format.
%
%
%       diffusionAnalysis: can be either
%
%   OUTPUT:
%       clustHistoryAll:  a 1D cell with rows = number of tracks in
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
%                        10) Diffusion mode.
%                        11) Diffusion coefficient.
%
%       clustHistoryMerged: Same information as in clustHistoryAll but with
%                         all cells merged into one 2D array, i.e.
%                         individual track information is lost.
%
%   Robel Yirdaw, 09/19/13
%       modified, 02/20/14
%       modified, 04/08/14
%       modified, 11/18/14
%       modified, May/June 2015 (Khuloud Jaqaman)
%       modified, September/October 2016 (Luciana de Oliveira)
%       modified, May 2018 (Tony Vega)
%       modified, October 2018 (Luciana de Oliveira)- if an event start in
%       the first frame, it has a NaN in the event starting the segment. In                                                                                          -
%       the same way if it is ended in the last frame.                                                                                                               -
%       modified, October 2018 (Luciana de Oliveira)- now the event that                                                                                             -
%       start/end and segment is completly determined from seqOfEvents
%       modified, November 2018 (Luciana de Oliveira)-incorporate information
%       of segments from alternative format
% *** NOTE (KJ):
%       As a result of modifications, columns in output have been
%       reshuffled. Thus output is not backward compatible with Robel's
%       functions, unless the latter are modified accordingly.
%
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

%Nothing to do

%% Cluster history collection

%Determine number of compTracks generated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%040814 - now passing defaultFormatTracks directly for memory reasons.
% tracksSize=tracksAggregStateDef(:).tracksFeatIndxCG;
% validTracks=tracksSize>1;
% tracksAggregStateDefTemp=tracksAggregStateDef(validTracks);
[numCompTracks,~] = size(tracksAggregStateDef);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Cluster history from all compTracks will be saved
clustHistoryAll = cell(numCompTracks,1);


%For each compTrack
for compTrackIter=1:numCompTracks
    
    %seqOfEvents for current compTrack
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %040814 - now passing defaultFormatTracks directly for memory reasons.
    seqOfEvents = tracksAggregStateDef(compTrackIter,1).seqOfEvents;

    if nargin>1 % add diffusion information
        %get seqOfEvents in the alternative format
        seqOfEventsAlt = tracksAggregStateAlt(compTrackIter,1).seqOfEvents;
        
        %load alt2defSegCorrenpondence
        alt2defSegmentCorrespond=tracksAggregStateAlt(compTrackIter).alt2defSegmentCorrespond;
    end
    
    if nargin>2 % add diffusion information
        %get diffusion stuff, right now assuming diff Mode, can fix later
        diffMode = diffusionAnalysisAlt(compTrackIter).diffMode;
        % get diffusion coefficient
        diffCoef=diffusionAnalysisAlt(compTrackIter).diffCoef;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %111814
    %For experimental data, the first iteration (frame) in seqOfEvents
    %can be some value other than 1.  This is occurs becuase an object
    %can appear or disappear at any time while in simulation, all of
    %the objects are present at frame 1. This affects how aggregState
    %is accessed since the frame # is supposed to correspond to a
    %column in aggregState. So, need to shift frame numbers in
    %seqOfEvents, if necessary, when accessing aggregState. The very
    %first frame # is the shift amount.
    frameShift = seqOfEvents(1,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %091813
    %aggregState will be used to get cluster sizes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %040814 - now passing defaultFormatTracks directly for memory reasons.
    aggregState = tracksAggregStateDef(compTrackIter,1).aggregState;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Modification Luciana de Oliveira, 10/2016
    
    % Mofication Luciana de Oliveira 12/2017-- add all coluns with
    % information about the cluster, before there was only some.
    
    %If the seqOfEvents has only one track, I need to account this cluster
    %size in the clust history
    
    if size(seqOfEvents,1)==2
        
        %Initialize cluster history
        switch nargin
            case 1 % only default information
                clustHistoryTemp = NaN(1,8);
            case 2 % add segment alternative number
                clustHistoryTemp = NaN(1,9);
            case 3 % add diffusion information
                clustHistoryTemp = NaN(1,11);
        end
        
        % track number
        clustHistoryTemp(1,1)=seqOfEvents(1,3);
        
        %cluster size
        clustHistoryTemp(1,2)= aggregState(seqOfEvents(1,1) - frameShift + 1,...
            seqOfEvents(1,1) - frameShift + 1);
        
        %start time
        clustHistoryTemp(1,3)=seqOfEvents(1,1);
        
        %end time
        clustHistoryTemp(1,4)=seqOfEvents(2,1);
        
        %Cluster lifetime
        clustHistoryTemp(1,5) = ...
            clustHistoryTemp(1,4) - clustHistoryTemp(1,3)+1;
        
        %Type of event starting cluster (1=dissoc., 2=assoc., 0=unknown)
        %Modification LRO 20181030
        clustHistoryTemp(1,6) = ~isnan(seqOfEvents(1,4))*seqOfEvents(1,2);

        %Type of event ending cluster (1=dissoc., 2=assoc., 0=unknown)
        clustHistoryTemp(1,7) = ~isnan(seqOfEvents(2,4))*seqOfEvents(2,2);
        
        if nargin>1 % add segment number alternative format
            
            %sergment number alternative format, in this case is the same as
            %the default format
            clustHistoryTemp(1,9) = clustHistoryTemp(1,1);
            
        end
        
        if nargin>2 % add diffusion information
            
            % diffusion mode
            clustHistoryTemp(1,10) = diffMode;
            % diffusion coefficient
            clustHistoryTemp(1,11) = diffCoef;
            
        end
        
        %Save current clustHistory in collection
        clustHistoryAll{compTrackIter,1} = clustHistoryTemp;
        
    else
        
        %Iteration points of first and last events
        firstEventIter = find(~isnan(seqOfEvents(:,4)),1,'first');
        lastEventIter = find(~isnan(seqOfEvents(:,4)),1,'last');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MODIFICATION LUCIANA DE OLIVEIRA
        % Part of function I add to calculate all cluster sizes for all tracks.
        % Now the function has as output a cluster history with all segments
        % history, I need it for the calculations of bootstraping. In the
        % original function clustHistory have only the history for the clusters
        % that are changing. Now it considers also the segments history for the
        % segments that starts in the first frame and end in the last frame.
        % The function also considers the segments that don't interact with
        % others.
        
        %Determine beginning segments
        trackIndexBegin=find(isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==1);
        
        %Determine changing segments
        trackIndexChanging= find(~isnan(seqOfEvents(:,4)));
        
        %Determine ending segments
        trackIndexEnd=find(isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==2);
        
        %total number of segments in clustHistory
        totSegments = length(trackIndexBegin) + (lastEventIter-firstEventIter+1) + ...
            sum(seqOfEvents(firstEventIter:lastEventIter,2)==1);
        
        %Initialize cluster history
        switch nargin
            case 1 % only default information
                clustHistoryTemp = NaN(totSegments,8);
            case 2 % add segment alternative number
                clustHistoryTemp = NaN(totSegments,9);
            case 3 % add diffusion information
                clustHistoryTemp = NaN(totSegments,11);
        end
        
        clustHistIndxbegin = 1;
        
        % cluster history for the segments that start with a real birth
        
        for trackIndexNew = 1 : length(trackIndexBegin)
            
            % track number
            clustHistoryTemp(clustHistIndxbegin,1) = seqOfEvents(trackIndexBegin(trackIndexNew),3);
            
            %Getting cluster size from entries in aggregState
            clustHistoryTemp(clustHistIndxbegin,2) = ...
                aggregState(seqOfEvents(trackIndexBegin(trackIndexNew),3),...
                seqOfEvents(trackIndexBegin(trackIndexNew),1) - frameShift + 1);
            
            % when it starts
            clustHistoryTemp(clustHistIndxbegin,3) = seqOfEvents(trackIndexBegin(trackIndexNew),1);
            
            % when it ends
            segmentRow = find(seqOfEvents(:,3)==seqOfEvents(trackIndexBegin(trackIndexNew),3) & seqOfEvents(:,2)==2 | seqOfEvents(:,4)==seqOfEvents(trackIndexBegin(trackIndexNew),3));
            clustHistoryTemp(clustHistIndxbegin,4) = seqOfEvents(segmentRow(1),1)-1; %Modification LRO 2018/11/02- the time ends need to be time event-1
            
            %Type of event starting cluster (1=dissoc., 2=assoc., 0=unknown)
            %Modification LRO 20181030
            clustHistoryTemp(clustHistIndxbegin,6) = ~isnan(seqOfEvents(trackIndexBegin(trackIndexNew),4)) * seqOfEvents(trackIndexBegin(trackIndexNew),2);
            
            %Type of event ending cluster (1=dissoc., 2=assoc., 0=unknown)
            %Modification LRO 20181030
            clustHistoryTemp(clustHistIndxbegin,7) = ~isnan(seqOfEvents(segmentRow(1),4)) * seqOfEvents(segmentRow(1),2);
            
            % size of resulting cluster
            clustHistoryTemp(clustHistIndxbegin,8) = ...
                aggregState(seqOfEvents(segmentRow(1),4),...
                seqOfEvents(segmentRow(1),1) - frameShift + 1);
            
            %increase in the clustHistIndx
            clustHistIndxbegin = clustHistIndxbegin + 1;
            
        end  %for track index
        
        % Initiate the cluster index for the clusters that are changing
        clustHistIndx=clustHistIndxbegin;
        
        %Process events in seqOfEvents
        for eventIndx=1:length(trackIndexChanging)
            
            %If event is a dissociation
            if (seqOfEvents(trackIndexChanging(eventIndx),2) == 1)
                
                %The dissociating receptor's (cluster's) segment number
                clustHistoryTemp(clustHistIndx,1) = seqOfEvents(trackIndexChanging(eventIndx),3);
                

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %091713
                %Getting cluster size from entries in aggregState
                
                %Use the row to access seqOfEvents and obtain
                %the segment number in column 3 which corresponds to the
                %newly split segment.  Then access aggregState with this
                %number as the row and the current iteration value (column
                %1 in seqOfEvents of the same row) as the column to
                %get the size of the cluster that split
                %
                %Modified 111814 to accomodate seqOfEvents that start at
                %frames > 1
                
                % cluster size
                clustHistoryTemp(clustHistIndx,2) = ...
                    aggregState(seqOfEvents(trackIndexChanging(eventIndx),3),...
                    seqOfEvents(trackIndexChanging(eventIndx),1) - frameShift + 1);
                
                %Starting iteration point of this cluster                
                clustHistoryTemp(clustHistIndx,3) = seqOfEvents(trackIndexChanging(eventIndx),1);
                
                %Type of event starting cluster (1=dissoc., 2=assoc., 0=unknown)
                %Modification LRO 20181030
                clustHistoryTemp(clustHistIndx,6) = ~isnan(seqOfEvents(trackIndexChanging(eventIndx),4))*seqOfEvents(trackIndexChanging(eventIndx),2);
                
                %Increment saved event count
                clustHistIndx = clustHistIndx + 1;
                
            else
                
                %Event is association
                %Find the cluster that is ending in clustHistory.  If it
                %doesn't exist, then this event is at the start of the
                %simulation and its beginning time is unknown. Only the
                %cluster it is joining will be saved below as a new cluster.
                endingTrackIndx = find(clustHistoryTemp(:,1)==seqOfEvents(trackIndexChanging(eventIndx),3),1,'last');
                
                if (~isempty(endingTrackIndx))
                    
                    %Ending iteration point of this cluster
                    
                    % modification, LRO-11/06/2018- there are cases that
                    % was considering a end point as a merge, now it is
                    % checked
                    if ~isnan(seqOfEvents(trackIndexChanging(eventIndx),4))
                        clustHistoryTemp(endingTrackIndx,4) = seqOfEvents(trackIndexChanging(eventIndx),1)-1;
                    else
                        clustHistoryTemp(endingTrackIndx,4) = seqOfEvents(trackIndexChanging(eventIndx),1);
                    end
                    
                    %Cluster lifetime
                    clustHistoryTemp(endingTrackIndx,5) = ...
                        clustHistoryTemp(endingTrackIndx,4) - clustHistoryTemp(endingTrackIndx,3)+1;
                    
                    %Type of event ending cluster (1=dissoc., 2=assoc., 0=unknown)
                    %Modification LRO 20181030
                    clustHistoryTemp(endingTrackIndx,7) = ~isnan(seqOfEvents(trackIndexChanging(eventIndx),4))*seqOfEvents(trackIndexChanging(eventIndx),2);
                    
                    %021714 - added column #8 for size of resulting cluster
                    %Modified 111814 to accomodate seqOfEvents from
                    %experimental data which can have NaN entries for
                    %column #4 corresponding to segments appearing or
                    %disappearing. The segment in column #4 is unknown if
                    %NaN, so in this case, the final size is unknown.
                    if (~isnan(seqOfEvents(trackIndexChanging(eventIndx),4)))
                        clustHistoryTemp(endingTrackIndx,8) = ...
                            aggregState(seqOfEvents(trackIndexChanging(eventIndx),4),...
                            seqOfEvents(trackIndexChanging(eventIndx),1) - frameShift + 1);
                    end
                    
                end
                
            end %If event 1 or 2
            
            %The cluster from/with which the above is occuring is also
            %changing. Update accordingly.
            %Find the changing segment in clustHistory.  If it doesn't exist, then
            %this event is at the start of the simulation and its beginning
            %time is unkown.  The newly created cluster will be saved below.
            changingTrackIndx = find(clustHistoryTemp(:,1)==seqOfEvents(trackIndexChanging(eventIndx),4),1,'last');
            
            if (~isempty(changingTrackIndx))
                
                %Ending iteration point of this cluster
                clustHistoryTemp(changingTrackIndx,4) = seqOfEvents(trackIndexChanging(eventIndx),1)-1;
                
                %Cluster lifetime
                clustHistoryTemp(changingTrackIndx,5) = ...
                    clustHistoryTemp(changingTrackIndx,4) - clustHistoryTemp(changingTrackIndx,3)+1;

                %Type of event ending cluster (1=dissoc., 2=assoc., 0=unknown)
                %Modification LRO 20181030
                clustHistoryTemp(changingTrackIndx,7) = ~isnan(seqOfEvents(trackIndexChanging(eventIndx),4))*seqOfEvents(trackIndexChanging(eventIndx),2);
                
                %021714 - added column #8 for size of resulting cluster
                %Modified 111814
                clustHistoryTemp(changingTrackIndx,8) = ...
                    aggregState(seqOfEvents(trackIndexChanging(eventIndx),4),...
                    seqOfEvents(trackIndexChanging(eventIndx),1) - frameShift + 1);
                
            end
            
            %A dissociation or association has occured involving the cluster
            %with segment number on column 4 of seqOfEvents.  This starts a new
            %cluster - save its segment number, current size and starting
            %iteration point in clustHistory.
            %Modified 111814 to accomodate seqOfEvents from experimental
            %data which can have NaN entries for column #4 corresponding to
            %segments appearing or disappearing. The segment in column #4
            %is unknown if NaN. So in this case, a new cluster can not be
            %started.
            if (~isnan(seqOfEvents(trackIndexChanging(eventIndx),4)))
                
                %The receptor's (cluster's) segment number
                clustHistoryTemp(clustHistIndx,1) = seqOfEvents(trackIndexChanging(eventIndx),4);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %091713
                %Getting cluster size from entries in aggregState; in this case
                %the row is in column 4 of seqOfEvents since this an update for
                %the changing track
                clustHistoryTemp(clustHistIndx,2) = ...
                    aggregState(seqOfEvents(trackIndexChanging(eventIndx),4),...
                    seqOfEvents(trackIndexChanging(eventIndx),1) - frameShift + 1);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %Starting iteration point of this cluster
                clustHistoryTemp(clustHistIndx,3) = seqOfEvents(trackIndexChanging(eventIndx),1);
                
                %Type of event starting cluster (1=dissoc., 2=assoc., 0=unknown)
                %Modification LRO 20181030
                clustHistoryTemp(clustHistIndx,6) = ~isnan(seqOfEvents(trackIndexChanging(eventIndx),4))*seqOfEvents(trackIndexChanging(eventIndx),2);
                
                %Increment saved event count
                clustHistIndx = clustHistIndx + 1;
                
            end
            
            
        end %eventIndx loop
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Modification Luciana de Oliveira
        % taking the cluster info for the segments that finish in the last frame
        
        %2016/11/01 LRO: I modify the size of clustHistTemp and then now is
        %necessary to remove rows that are not filled with values
        
        %remove rows with all NaNs
        iRowNaN= sum(isnan(clustHistoryTemp),2)==size(clustHistoryTemp,2);
        clustHistoryTemp(iRowNaN,:)=[];
        
        % 2016/11/01 LRO: find the end of the segment if it is in the tracks
        % that are in the end of the frames.
        
        segmentsNotEnd=find(isnan(clustHistoryTemp(:,4)));
        
        for segIndex=1:length(segmentsNotEnd)
            
            segmentEnd=find(clustHistoryTemp(segmentsNotEnd(segIndex),1)==seqOfEvents(trackIndexEnd,3));
            
            if ~isempty(segmentEnd)
                clustHistoryTemp(segmentsNotEnd(segIndex),4)=seqOfEvents(trackIndexEnd(segmentEnd(end)),1);
                
                %Cluster lifetime
                clustHistoryTemp(segIndex,5) = ...
                    clustHistoryTemp(segIndex,4) - clustHistoryTemp(segIndex,3)+1;
                
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%% modification LRO
        %%%% calculate life time for those that still don't have this value
        clustHistoryTemp(:,5)=clustHistoryTemp(:,4)-clustHistoryTemp(:,3)+1;
        
        % % % %
        % % % %         %%%% here is the addition to avoid lifetime=0 coming from old version
        % % % %         %%%% of the function--Luciana de Oliveira--2017/12/26
        % % % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %remove clusters which start and end at exactly the same time point -
        %these would be not detectable with the sub-sampling time step
        clustHistoryTemp = clustHistoryTemp(clustHistoryTemp(:,5)>0,:);
        % % % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Modification LRO 2018/10/30 to accommodate the new terminology for
        % unknown event==0.
        listBeginsAndEnds=isnan(clustHistoryTemp(:,6:7));
        clustHistoryTemp(listBeginsAndEnds(:,1),6)=0;
        clustHistoryTemp(listBeginsAndEnds(:,2),7)=0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % optional information
        
        if nargin>1 % add sergment number alternative format
            
            % initialize alternative segments as the copy of default
            clustHistoryTemp(:,9)=clustHistoryTemp(:,1);
            
            % find the alternative segment number equivalent to the segment
            % in the default format
            for indexSegments=1:size(clustHistoryTemp,1)
                
                % take segment number and time start from the clustHistory
                % (default format)
                segDefNumber=clustHistoryTemp(indexSegments,1);
                startTimeSegDef=clustHistoryTemp(indexSegments,3);
                eventType=clustHistoryTemp(indexSegments,6);
                
                % Find this segment in the alt2default.
                rowSegAlt=alt2defSegmentCorrespond(:,2)==segDefNumber;
                segAltNumber=alt2defSegmentCorrespond(rowSegAlt,1);
                
                % if this segment is not changed from default to
                % alternative, the segAltNumber is empty and it should just
                % maintain the same number as in default
                if ~isempty(segAltNumber)
                    
                    
                    % use time start to find the segment equivalent and
                    % segment number in the alternative format to check the
                    % valeu in the alternative format
                    for indexAlt=1:length(segAltNumber)
                        % need to have the difference for merges and splits
                        %it is a merge
                        
                        if eventType==2
                            
                            rowAltConfirm=find(seqOfEventsAlt(:,1)==startTimeSegDef&seqOfEventsAlt(:,4)==segAltNumber(indexAlt), 1);
                            
                            if ~isempty(rowAltConfirm)
                                
                                clustHistoryTemp(indexSegments,9)=segAltNumber(indexAlt);
                            end
                            
                        elseif eventType==1||eventType==0
                            
                            rowAltConfirm=find(seqOfEventsAlt(:,1)==startTimeSegDef&seqOfEventsAlt(:,3)==segAltNumber(indexAlt), 1);
                            
                            if ~isempty(rowAltConfirm)
                                
                                clustHistoryTemp(indexSegments,9)=segAltNumber(indexAlt);
                            end
                            
                        end
                        
                    end
                    
                end
            end
            
        end %(if nargin>1)
        
        if nargin>2 % add diffusion information
            
            for indexSegments=1:size(clustHistoryTemp,1)
                % diffusion mode
                clustHistoryTemp(indexSegments,10) = diffMode(clustHistoryTemp(indexSegments,9));
                % diffusion coefficient
                clustHistoryTemp(indexSegments,11) = diffCoef(clustHistoryTemp(indexSegments,9));
            end
        end
        
    end %(if size(seqOfEvents,1)==2 ... else ...
    
    
    %Save current clustHistory in collection
    
    clustHistoryAll{compTrackIter,1} = clustHistoryTemp;
    
    clear clustHistoryTemp
    
end

clustHistoryMerged = cat(1,clustHistoryAll{:,1});

end





