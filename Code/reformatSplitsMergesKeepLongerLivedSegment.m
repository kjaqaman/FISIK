function [tracksReformat] = reformatSplitsMergesKeepLongerLivedSegment(tracksIn)
% REFORMATSPLITSMERGESKEEPLONGERLIVEDSEGMENT reformats compound tracks so that shorter-lived segments merge with or split from longer-lived segments
%
% SYNOPSIS [tracksReformat] = reformatSplitsMergesKeepLongerLivedSegment(tracksIn)
%
% INPUT
%          tracksIn       : Tracks in the format of the output of
%                           trackCloseGapsKalmanSparse (or
%                           u-track for GUI version).
% OUTPUT
%          tracksReformat : Tracks after reformating.
%
% Khuloud Jaqaman, April 2019
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

%initiate the reform compTrack
tracksReformat=tracksIn;

% go over all the tracks
for iTrack = 1: length(tracksReformat)
    
    %load track information
    seqOfEvents= tracksReformat(iTrack).seqOfEvents;
    tracksFeatIndxCG = tracksReformat(iTrack).tracksFeatIndxCG;
    tracksCoordAmpCG = tracksReformat(iTrack).tracksCoordAmpCG;
    
    %correction for indexing
    seqOfEventsCorr = seqOfEvents(:,1)-seqOfEvents(1,1)+1;
    
    %find rows reporting merges or splits
    eventRows = find(~isnan(seqOfEvents(:,4)));
    
    %go over all the events    
    for indexEvent = 1 : length(eventRows)
        
        % to be sure of analysing only the rows that some event is happening
        iEvent = eventRows(indexEvent);
        
        % determine the kind of event
        kindOfEvent=seqOfEvents(iEvent,2);
        
        if kindOfEvent == 1 %if the event is a split
            
            %find the two splitting segments and the splitting time
            segment1 = seqOfEvents(iEvent,3);
            segment2 = seqOfEvents(iEvent,4);
            timeEventCorr = seqOfEventsCorr(iEvent,1);
            
            %find the end times of the two segments
            timeEnd1 = seqOfEvents(seqOfEvents(:,3)==segment1 & seqOfEvents(:,2)==2, 1);
            timeEnd2 = seqOfEvents(seqOfEvents(:,3)==segment2 & seqOfEvents(:,2)==2, 1);
            
            %if the splitting segment (segment1) ends after the segment it
            %splits from (segment2), flip the two segments' futures
            if timeEnd1 > timeEnd2
                
                %exchange segment numbers in all events after the split
                [seqOfEventsTemp,dummy] = deal(seqOfEvents(iEvent+1:end,3:4));
                seqOfEventsTemp(dummy==segment2) = segment1;
                seqOfEventsTemp(dummy==segment1) = segment2;
                seqOfEvents(iEvent+1:end,3:4) = seqOfEventsTemp;
                
                %exchange track information between the two segments from
                %the split onwards
                
                %feature indices
                tracksFeatIndxSeg1 = tracksFeatIndxCG(segment1,timeEventCorr:end);
                tracksFeatIndxCG(segment1,timeEventCorr:end) = tracksFeatIndxCG(segment2,timeEventCorr:end);
                tracksFeatIndxCG(segment2,timeEventCorr:end) = tracksFeatIndxSeg1;
                
                %coordinates and amplitudes
                tracksCoordAmpSeg1 = tracksCoordAmpCG(segment1,8*(timeEventCorr-1)+1:end);
                tracksCoordAmpCG(segment1,8*(timeEventCorr-1)+1:end) = tracksCoordAmpCG(segment2,8*(timeEventCorr-1)+1:end);
                tracksCoordAmpCG(segment2,8*(timeEventCorr-1)+1:end) = tracksCoordAmpSeg1;
                
            end
            
        else %if the event is a merge
            
            %find the two merging segments and the merging time
            segment1 = seqOfEvents(iEvent,3);
            segment2 = seqOfEvents(iEvent,4);
            timeEventCorr = seqOfEventsCorr(iEvent,1);
            
            %find the start times of the two segments
            timeStart1 = seqOfEvents(seqOfEvents(:,3)==segment1 & seqOfEvents(:,2)==1, 1);
            timeStart2 = seqOfEvents(seqOfEvents(:,3)==segment2 & seqOfEvents(:,2)==1, 1);
            
            %if the merging segment (segment1) started before the
            %segment it is merging with (segment2), flip the merge
            if timeStart1 < timeStart2
                                
                %exchange segment numbers in all events from the merge
                %onwards
                [seqOfEventsTemp,dummy] = deal(seqOfEvents(iEvent:end,3:4));
                seqOfEventsTemp(dummy==segment2) = segment1;
                seqOfEventsTemp(dummy==segment1) = segment2;
                seqOfEvents(iEvent:end,3:4) = seqOfEventsTemp;
                                
                %exchange track information between the two segments from
                %the merge onwards
                
                %feature indices
                tracksFeatIndxSeg1 = tracksFeatIndxCG(segment1,timeEventCorr:end);
                tracksFeatIndxCG(segment1,timeEventCorr:end) = tracksFeatIndxCG(segment2,timeEventCorr:end);
                tracksFeatIndxCG(segment2,timeEventCorr:end) = tracksFeatIndxSeg1;
                
                %coordinates and amplitudes
                tracksCoordAmpSeg1 = tracksCoordAmpCG(segment1,8*(timeEventCorr-1)+1:end);
                tracksCoordAmpCG(segment1,8*(timeEventCorr-1)+1:end) = tracksCoordAmpCG(segment2,8*(timeEventCorr-1)+1:end);
                tracksCoordAmpCG(segment2,8*(timeEventCorr-1)+1:end) = tracksCoordAmpSeg1;
            
            end
            
        end %(if kindOfEvent == 1 ... else ...)
        
    end %(for indexEvent = 1 : length(eventRows))
    
    % save tracks in the final format
    tracksReformat(iTrack).tracksFeatIndxCG = tracksFeatIndxCG;
    tracksReformat(iTrack).tracksCoordAmpCG = tracksCoordAmpCG;
    tracksReformat(iTrack).seqOfEvents = seqOfEvents;
    
end %(for iTrack = 1: length(tracksReformat))


