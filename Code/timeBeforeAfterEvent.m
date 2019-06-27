function  timeBefAftEvent  = timeBeforeAfterEvent(seqOfEvents,segmentMerge,segmentSplit,timeMerge,segmentContinuous,timeSplit,firstFrame )
% timeBeforeAfterEvent calculates the time before and after the event for
% the pairwsing of segments.
% timeBefAftEvent  = timeBeforeAfterEvent(seqOfEvents,timeMerge,timeSplit )
%
% INPUT
%       seqOfEvents   :        Matrix with number of rows equal to number
%                              of events happening in a compound track and 4
%                              columns:
%                              1st: Frame where event happens;
%                              2nd: 1 = start of track segment, 2 = end of track segment;
%                              3rd: Index of track segment that ends or starts;
%                              4th: NaN = start is a birth and end is a death,
%                              number = start is due to a split, end
%                              is due to a merge, number is the index
%                              of track segment for the merge/split.
%
%       segmentMerge  :        segment that is merging in the merge event to
%                              the contineous segment
%
%       segmentSplit  :        segment that is spliting in the split event
%                              from the contineous segment
%
%       timeMerge  :        segment that is merging in the merge event to
%                              the contineous segment
%
%       timeSplit  :        segment that is spliting in the split event
%                              from the contineous segment
%
%       firstFrame:         first frame that a event happens in a compTrack
%
%
% OUTPUT
%       timeBefAftEvent     :   vector with 4 values:
%                               1) last time that the segment merge has
%                                  been involved in an event before the
%                                  merge event.
%                               2) first time that the segment merge has
%                               been involved in an event after the split
%                               event.
%
%                               3) first time that the segment split has
%                               been involved in an event before the merge
%                               event.
%
%                               4) first time that the segment split has
%                               been involved in an event after the split
%                               event.
%
%
%
% Luciana de Oliveira, March 2018
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
lastFrame= seqOfEvents(end,1);
% for the simultaneous merge to split there is only one time
% of event, in that case
if nargin<7
    timeSplit=timeMerge;
end

% determine which are the events in this seqOfEvents

notNanEvents=find(~isnan(seqOfEvents(1:end,4)));


% check if there is an event before the merge for both merge and split

% take only the events before the merge
possibleBefMerge= seqOfEvents(notNanEvents,1)<timeMerge;
possibleAftSplit= find(seqOfEvents(notNanEvents,1)>timeSplit);
% take only the seqOfEvents part that has no NaNs and is before and
% after the merge 2 split event
seqOfEventsBefore = seqOfEvents(notNanEvents(possibleBefMerge),:);
seqOfEventsAfter  = seqOfEvents(notNanEvents(possibleAftSplit),:);
if nargin<7
    
    
    %%%%%%%%%%%%% Modification LRO 2019/04/17%%%%%%%%%%
    % it was not taking the proper interval
    eventBefMerge1=find(seqOfEventsBefore(:,3)==segmentMerge);
    eventBefMerge2=find(seqOfEventsBefore(:,4)==segmentMerge);
    eventBefMerge=[eventBefMerge1;eventBefMerge2];
    eventAftSplit1=find(seqOfEventsAfter(:,3)==segmentSplit);
     eventAftSplit2=find(seqOfEventsAfter(:,4)==segmentSplit);
    eventAftSplit=[eventAftSplit1;eventAftSplit2];
    
    
    if ~isempty(eventBefMerge)
        timeBefMerge=seqOfEventsBefore(max(eventBefMerge),1);
    else
        timeBefMerge=1;
    end
    
    if ~isempty(eventAftSplit)
        timeAftSplit=seqOfEvents(notNanEvents(possibleAftSplit(min(eventAftSplit))),1);
    else
        timeAftSplit=lastFrame;
    end
    
    % output time
    
    timeBefAftEvent=[timeBefMerge,timeAftSplit];
elseif nargin==7
    % check if there is an event before the merge for both merge and split
    [eventBefMerge,~]=find(seqOfEventsBefore(:,3:4)==segmentMerge);
    [eventAftSplit,~]=find(seqOfEventsAfter(:,3:4)==segmentSplit);
    [eventBefCont,~]=find(seqOfEventsBefore(:,3:4)==segmentContinuous);
    [eventAftCont,~]=find(seqOfEventsAfter(:,3:4)==segmentContinuous);
    
    if eventBefMerge
        timeBefMerge=seqOfEventsBefore(eventBefMerge(end),1);
    else
        timeBefMerge=firstFrame;
    end
    
    if eventAftSplit
        timeAftSplit=seqOfEvents(notNanEvents(possibleAftSplit(eventAftSplit(1))),1);
    else
        timeAftSplit=lastFrame;
    end
    
    if eventBefCont
        timeBefCont=seqOfEventsBefore(eventBefCont(end),1);
    else
        timeBefCont=firstFrame;
    end
    
    if eventAftCont
        timeAftCont=seqOfEventsAfter(eventAftCont(1),1);
    else
        timeAftCont=lastFrame;
    end
    
    % output time
    timeBefAftEvent=[timeBefMerge,timeAftSplit,timeBefCont,timeAftCont];
    
end

end

