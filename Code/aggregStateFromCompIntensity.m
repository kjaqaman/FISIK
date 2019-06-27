function [compTracksOut] = aggregStateFromCompIntensity(compTracks)
%AGGREGSTATEFROMCOMPINTENSITY recovers particle oligomeric states directly from their intensity.
%
%SYNOPSIS compTracks = aggregStateFromCompIntensity(compTracks)
%
%INPUT  compTracks   : Compound tracks, in the format of tracksFinal as
%                      output by receptorAggregationSimpleOptimized (it is
%                      a field of the structure receptorInfoLabeled).
%
%OUTPUT compTracks   : -Structure with the 2 fields: "defaultFormatTracks" and
%                       "alternativeFormatTracks".
%                      -Both contain the fields "tracksFeatIndxCG",
%                       "tracksCoordAmpCG", "seqOfEvents" and "aggregState".
%                      -"alternativeFormatTracks" also contains the field
%                       "alt2defSegmentCorrespond".
%                      -"defaultFormatTracks" is the format of the output of
%                       "trackCloseGapsKalmanSparse" (or u-track for the GUI version).
%                      -"alternativeFormatTracks" is the format generated by
%                       the function "convFormatDefault2Alt", where tracks
%                       do not continue through merges and splits, but a
%                       merge consists of 2 track segments merging to form a
%                       3rd segment, and a split consists of 1 track segment
%                       splitting into 2 different segments.
%                      -The field "aggregState" has the same dimensions as
%                       "tracksFeatIndxCG", and indicates the estimated
%                       number of units (e.g. receptors) within each
%                       detected particle/feature.
%
%
%
%IMPORTANT: This code assumes that unit intensity has mean 1 and std 0
%
%Luciana de Oliveira, April 2019, based on aggregStateFromCompTracksMIQP
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

%convert tracks from default format (output of trackCloseGapsKalman) to
%alternative format where there is no track continuation through a merge or
%a split
compTracksDef = compTracks;
fprintf('\nIn aggregState...entering convTrackFormatDefault2Alt...');
convDef2AltTic = tic;
compTracks = convTrackFormatDefault2Alt(compTracksDef);
convDef2AltTime = toc(convDef2AltTic);
fprintf('\nDone with convTrackFormatDefault2Alt. Elapsed time is %g seconds.\n',convDef2AltTime);



%% Calculation

%go over all compound tracks
for iTrack = 1 : numTracks
    %         for iTrack = 11798
    
    fprintf('\nIn aggregState...at iTrack %d',iTrack);
    
    %get this compound track's information
    seqOfEvents = compTracks(iTrack).seqOfEvents;
    tracksAmp = compTracks(iTrack).tracksCoordAmpCG(:,4:8:end);
    doubleFreq = mod(seqOfEvents(1,1)*2,2)==1;
    aggregStateMat  = full(tracksAmp);
    
    %% store aggregation state matrix in compound tracks structure with default format
    
    %get number of segments in default format
    numSegmentsDef = size(compTracksDef(iTrack).tracksFeatIndxCG,1);
    
    %copy aggregation state matrix
    aggregStateDef = aggregStateMat(1:numSegmentsDef,:);
    
    %copy out segment correspondence array
    segmentCorrespond = compTracks(iTrack).alt2defSegmentCorrespond;
    
    %go over additional segments in alternative format and store their
    %aggregation state in the original segment location
    for iSegment = 1 : size(segmentCorrespond,1)
        segmentNew = segmentCorrespond(iSegment,1);
        segmentOld = segmentCorrespond(iSegment,2);
        aggregStateDef(segmentOld,:) = max([aggregStateDef(segmentOld,:);...
            aggregStateMat(segmentNew,:)]);
    end
    
    %Moved from above with conversion to sparse added
    if (issparse(compTracks(iTrack).tracksCoordAmpCG))
        aggregStateMat(isnan(aggregStateMat)) = 0;
        compTracks(iTrack).aggregState = sparse(aggregStateMat);
        %Also converting to sparse the original line below
        aggregState = aggregStateDef(:,1+doubleFreq:(1+doubleFreq):end);
        aggregState(isnan(aggregState)) = 0;
        compTracksDef(iTrack).aggregState = sparse(aggregState);
    else
        compTracks(iTrack).aggregState = aggregStateMat;
        compTracksDef(iTrack).aggregState = aggregStateDef(:,1+doubleFreq:(1+doubleFreq):end);
    end
    
end %(for iTrack = 1 : numTracks)

%store results in output structure
compTracksOut = struct('defaultFormatTracks',compTracksDef,'alternativeFormatTracks',...
    compTracks);

%% ~~~ the end ~~~

end



