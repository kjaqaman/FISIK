function compTracks = convReceptClust2CompTracksPhotobleaching(clust2recept,recept2clust,...
    receptorTraj,receptorIntensity)
%CONVRECEPTCLUST2COMPTRACK generates compound tracks out of receptor clustering information over time
%
%SYNOPSIS compTracks = convReceptClust2CompTracks(clust2recept,recept2clust,receptorTraj)
%
%INPUT  clust2recept: Equivalent to clust2receptAssign as output by
%                     receptorAggregationSimple.
%       recept2clust: Equivalent to recept2clustAssign as output by
%                     receptorAggregationSimple.
%       receptorTraj: Equivalent to receptorTraj as output by
%                     receptorAggregationSimple.
%       receptorIntensity: (Number of receptors) - by - (number of
%                          time points) array of receptor intensities.
%                          Optional. Default: All 1.
%
%OUTPUT compTracks  : Compound tracks, in the format of tracksFinal as
%                     output by trackCloseGapsKalman.
%
%Khuloud Jaqaman, January 2009
%
%Modified 08/20/14, Robel Yirdaw
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

%get maximum number of clusters, maximum cluster size and number of frames
[numClustMax,sizeClustMax,numFrames] = size(clust2recept);

%get number of receptors system dimensionality
[numReceptors,dimension,numFrames] = size(receptorTraj);

%assign receptor intensities if not input
if nargin < 4 || isempty(receptorIntensity)
    receptorIntensity = ones(numReceptors,numFrames);
end

%% First find interacting clusters and construct a global cluster table
%each global cluster will make one compound track

%find number of clusters in 1st frame
clust2receptFrame1 = clust2recept(:,:,1);
numClustFrame1 = length(find(clust2receptFrame1(:,1)~=0));

%initialize array of global clusters
clust2receptGlobal = zeros(numClustFrame1,sizeClustMax);

%initialize array indicating which clusters to look at
lookAtCluster = ones(numClustFrame1,1);

%go over all clusters in first frame
for iClust = 1 : numClustFrame1
    
    if lookAtCluster(iClust) == 1

        %find receptors in current cluster in first frame
        receptClustWholeMovieNew = clust2receptFrame1(iClust,:);
        receptClustWholeMovieNew = receptClustWholeMovieNew(receptClustWholeMovieNew~=0)';
        listGrew = 1;

        while listGrew

            %copy list of receptors
            receptClustWholeMovie = receptClustWholeMovieNew;

            %go over all frames that the receptor is present and find which other receptors ever cluster with
            %these receptors
                        
            % modification LRO 2019/05/15 to accomodate receptors that are
            % not ending in the last frame
           
            receptorsAdditional = []; %fix by KJ; take out of the for-loop
            for indexRecept = 1 : length(receptClustWholeMovie)
                
                % find the last time that the receptor is present
                lastFrame=find(recept2clust(receptClustWholeMovie(indexRecept),:)~=0,1,'last');
                
                %go over all times that receptor exists and add any
                %interacting partners
                for iFrame = 1 : lastFrame
                    tmp = clust2recept(recept2clust(receptClustWholeMovie(indexRecept),iFrame),:,iFrame);
                    tmp = tmp(:);
                    receptorsAdditional = [receptorsAdditional; tmp(tmp~=0)];
                end
                
            end
            
            %add these receptors to list of receptors
            receptClustWholeMovieNew = unique([receptClustWholeMovie; receptorsAdditional]);

            %check whether new list of receptors is longer than old list
            listGrew = length(receptClustWholeMovieNew) > length(receptClustWholeMovie);

        end

        %find which clusters these receptors belong to in the first frame
        clustIDFrame1 = unique(recept2clust(receptClustWholeMovieNew,1));

        %tell code not to look at these clusters again
        lookAtCluster(clustIDFrame1) = 0;

        %store global cluster information in list
        numReceptInClust = length(receptClustWholeMovieNew);
        clust2receptGlobal(iClust,1:numReceptInClust) = receptClustWholeMovieNew;
        
    end

end

%remove all repetition and all zeros
clust2receptGlobal = unique(clust2receptGlobal,'rows');
clust2receptGlobal = clust2receptGlobal(sum(clust2receptGlobal,2)~=0,:);

%get number of global clusters = number of compound tracks
numClustGlobal = size(clust2receptGlobal,1);

%% Then convert the global clusters to the format of the output of trackCloseGapsKalman

%initialize the compound tracks (output of trackCloseGapsKalman)
compTracks = repmat(struct('tracksFeatIndxCG',[],'tracksCoordAmpCG',[],...
    'seqOfEvents',[]),numClustGlobal,1);

%go over all global receptor clusters
for iClust = 1 : numClustGlobal

    %get participating receptors in current global cluster
    clusterMembers = clust2receptGlobal(iClust,:);
    clusterMembers = clusterMembers(clusterMembers~=0);
    
    %go over all frames with the receptor present (Modification LRO 2019/05/16)
    % find the last time that each receptor is present, take the largest time
    % when there are multiple receptors in the same cluster
    %define lastFrame
    %KJ modifications - previous stuff seemed unncessarily complicated
    lastFrame = NaN(length(clusterMembers),1);
    for indexMembers = 1 : length(clusterMembers)
        lastFrame(indexMembers) = find(recept2clust(clusterMembers(indexMembers),:)~=0,1,'last');
    end 
    lastFrame = max(lastFrame);
    
    %fill in the local-cluster numbers in "tracksFeatIndxCG"
    %local-cluster numbers are the equivalent of feature indices (from
    %detection)
    clustersFrame1 = unique(recept2clust(clusterMembers,1));
    numClustInit = length(clustersFrame1);
    tracksFeatIndxCG = [clustersFrame1 zeros(numClustInit,lastFrame-1)]; %fix KJ
    
    %initialize sequence of events and final number of clusters
    seqOfEvents = [ones(numClustInit,2) (1:numClustInit)' NaN(numClustInit,1)];
    numClustFinal = numClustInit;
    
    for iFrame = 1 : lastFrame-1
    
        %find the local-clusters in current frame
        clustersCurrentFrame = tracksFeatIndxCG(:,iFrame);
        indxNonZero = find(clustersCurrentFrame~=0); %the non-empty rows
        
        %find the local-clusters in next frame
        clustersNextFrame = unique(recept2clust(clusterMembers,iFrame+1));
        
        % remove zero receptors
        clustersNextFrame=clustersNextFrame(clustersNextFrame~=0); 
        
        %get the cluster sizes
        clustSizeCurrentFrame = [];
        clustSizeCurrentFrame(indxNonZero,:) = clust2recept( ...
            clustersCurrentFrame(indxNonZero),:,iFrame);
        [tmpDim1,tmpDim2] = size(clustSizeCurrentFrame);
        clustSizeCurrentFrame = [clustSizeCurrentFrame; zeros(size(...
            clustersCurrentFrame,1)-tmpDim1,tmpDim2)];
        clustSizeCurrentFrame = sum(clustSizeCurrentFrame~=0,2);
        clustSizeNextFrame = clust2recept(clustersNextFrame,:,iFrame+1); %%% here need to be changed LRO 
        clustSizeNextFrame = sum(clustSizeNextFrame~=0,2);
        
        %determine which clusters have only one receptor and which have
        %more
        clustersBigCurrent = find(clustSizeCurrentFrame>1);
        clustersOneCurrent = find(clustSizeCurrentFrame==1);
        clustersBigNext = find(clustSizeNextFrame>1);
        
        %for clusters with only one receptor in current frame, map what
        %cluster they go to in next frame
        clusterMapCurrent2Next = zeros(size(clustersCurrentFrame));
        receptTmp = clust2recept(clustersCurrentFrame(clustersOneCurrent),1,iFrame);
        clusterMapCurrent2Next(clustersOneCurrent,1) = recept2clust(...
            receptTmp,iFrame+1);
        
        %for clusters with more than one receptor in current frame, map
        %what clusters they go to and store splitting information
        for iClustBig = clustersBigCurrent'
            
            %find where receptors of current cluster go
            receptTmp = clust2recept(clustersCurrentFrame(iClustBig),:,iFrame);
            receptTmp = receptTmp(receptTmp~=0);
            clustMapTmp = recept2clust(receptTmp,iFrame+1);
            clustMapTmp = unique(clustMapTmp);
            clustMapTmp= clustMapTmp(clustMapTmp~=0);
            
            if length(clustMapTmp) == 1 %if receptors stay in the same cluster
                
                %store that cluster number
                clusterMapCurrent2Next(iClustBig,1) = clustMapTmp(1);
                
            else %if they split
                
                %store the first cluster number as the continuation
                clusterMapCurrent2Next(iClustBig,1) = clustMapTmp(1);
                clustMapTmp = clustMapTmp(2:end);
                
                %store the rest of the clusters as splits
                numSplits = length(clustMapTmp);
                numClustFinal = numClustInit + numSplits;
                clusterMapCurrent2Next(numClustInit+1:numClustFinal,1) = clustMapTmp;
                
                %store splitting information in seqOfEvents
                seqOfEventsTmp = [(iFrame+1)*ones(numSplits,1) ones(numSplits,1) ...
                    (numClustInit+1:numClustFinal)' iClustBig*ones(numSplits,1)];
                seqOfEvents = [seqOfEvents; seqOfEventsTmp];
                numClustInit = numClustFinal;
                
            end %(if length(clustMapTmp) == 1 ... else ...)
            
        end %(for iClustBig = clustersBigCurrent')

        %for clusters with more than one receptor in next frame, map
        %what clusters they come from and store merging information
        for iClustBig = clustersBigNext'

            %find where receptors of cluster come from
            receptTmp = clust2recept(clustersNextFrame(iClustBig),:,iFrame+1);
            receptTmp = receptTmp(receptTmp~=0);
            clustMapTmp = recept2clust(receptTmp,iFrame);
            clustMapTmp = unique(clustMapTmp);

            %if receptors come from more than one cluster, this indicates a
            %merge
            %in this case, remove all mappings except for the first mapped
            %cluster, and designate the rest as merges
            if length(clustMapTmp) > 1

                %find all mentions of iClustBig
                clustBigIndx = find(clusterMapCurrent2Next==clustersNextFrame(iClustBig));
                
                %remove all mentions but the first
                clusterMapCurrent2Next(clustBigIndx(2:end),1) = 0;
                
                %store merging information in seqOfEvents
                numMerges = length(clustBigIndx) - 1;
                seqOfEventsTmp = [(iFrame+1)*ones(numMerges,1) 2*ones(numMerges,1) ...
                    clustBigIndx(2:end) clustBigIndx(1)*ones(numMerges,1)];
                seqOfEvents = [seqOfEvents; seqOfEventsTmp];

            end %(if length(clustMapTmp) == 1)

        end %(for iClustBig = clustersBigCurrent')
        
        %store cluster ("feature") connectivity information
        tracksFeatIndxCG(1:numClustFinal,iFrame+1) = clusterMapCurrent2Next;

    end %(for iFrame = 1 : lastFrames - 1)
    
    %sort sequence of events in increasing order of time
    [~,indxSort] = sort(seqOfEvents(:,1));
    seqOfEvents = seqOfEvents(indxSort,:);

    %go over unique event times
    eventTimes = unique(seqOfEvents(:,1));
    for iEvent = 1 : length(eventTimes)

        %find events happening at this event time
        indxEventsAtEventTime = find(seqOfEvents(:,1)==eventTimes(iEvent));
        numEventsAtEventTime = length(indxEventsAtEventTime);

        %if there is more than 1, sort events such that splits and starts
        %come before merges and ends
        if numEventsAtEventTime > 1
            seqOfEventsTmp = sortrows(seqOfEvents(indxEventsAtEventTime,:),2);
            seqOfEvents(indxEventsAtEventTime,:) = seqOfEventsTmp;
        end
    end
    
    %check whether any track segments are empty
    indxEmpty = find(max(tracksFeatIndxCG,[],2)==0);

    %if there are empty track segments, that means they split from one
    %track segment and merge with another at the same time point
    %in this case, switch the merging and splitting segments to avoid empty
    %segments
    if ~isempty(indxEmpty)
        for iSegment = indxEmpty'

            %determine the split event row and time, and the segment that got
            %split from
            splitRow = find(seqOfEvents(:,2)==1&seqOfEvents(:,3)==iSegment);
            splitTime = seqOfEvents(splitRow,1);
            splitSegment = seqOfEvents(splitRow,4);

            %determine the merge event row and time, adn the segment that got
            %mered with
            mergeRow = find(seqOfEvents(:,2)==2&seqOfEvents(:,3)==iSegment);
            mergeTime = seqOfEvents(mergeRow,1);
            mergeSegment = seqOfEvents(mergeRow,4);
            
            %if the mergeRow is not immediately after the splitRow,
            %push the mergeRow up to just below the split row
            if mergeRow ~= (splitRow + 1)
                seqOfEvents(splitRow+1:mergeRow,:) = seqOfEvents(...
                    [mergeRow splitRow+1:mergeRow-1],:);
                mergeRow = splitRow + 1;
            end
            
            
            %switch segment numbers in sequence of events at times of
            %splitting and merging
            seqOfEvents(splitRow,3:4) = seqOfEvents(splitRow,4:-1:3);
            seqOfEvents(mergeRow,3:4) = seqOfEvents(mergeRow,4:-1:3);
           
            %replace segment numbers before splitting
            seqOfEventsBS = seqOfEvents(1:splitRow-1,:);
            seqOfEventsBS(seqOfEventsBS(:,3)==splitSegment,3) = iSegment;
            seqOfEventsBS(seqOfEventsBS(:,4)==splitSegment,4) = iSegment;
            seqOfEvents(1:splitRow-1,:) = seqOfEventsBS;
            
            %replace segment numbers after merging
            seqOfEventsAM = seqOfEvents(mergeRow+1:end,:);
            seqOfEventsAM(seqOfEventsAM(:,3)==mergeSegment,3) = iSegment;
            seqOfEventsAM(seqOfEventsAM(:,4)==mergeSegment,4) = iSegment;
            seqOfEvents(mergeRow+1:end,:) = seqOfEventsAM;

            %move feature indices to their new segment
            tracksFeatIndxCG(iSegment,1:splitTime-1) = tracksFeatIndxCG(...
                splitSegment,1:splitTime-1);
            tracksFeatIndxCG(splitSegment,1:splitTime-1) = 0;
            tracksFeatIndxCG(iSegment,mergeTime:end) = tracksFeatIndxCG(...
                mergeSegment,mergeTime:end);
            tracksFeatIndxCG(mergeSegment,mergeTime:end) = 0;

        end
    end
    
    %%%%%%%%%% LRO 2019/05/14%%%%%%%%%%%
    %Need to find the segments that end with a real end (not by merge)
    
    % first identify the segments that are in this compTrack
    segComp=unique(seqOfEvents(:,3:4));
    segComp=segComp(~isnan(segComp));
    
    % take only those that are not ended by a merge
    for indexSeg = 1 : length(segComp)
        
        %find if the segment end has been recorded (i.e. it's a merge)
        segFlag=find(seqOfEvents(:,2)==2 & seqOfEvents(:,3)==segComp(indexSeg), 1);
        
        if isempty(segFlag) % means that this segment end has not been recorded
            
            %find in which frame it ends adn store information
            frameEnd=find(tracksFeatIndxCG(segComp(indexSeg),:)~=0,1,'last');
            seqOfEventsTmp = [frameEnd 2 segComp(indexSeg) NaN];
            seqOfEvents = [seqOfEvents; seqOfEventsTmp];
            
        end
    end
    
    %%%% modification 2019/05/14 not sure if I still need this, but for now I
    %%%% will leave it here
    [~,indxSort] = sort(seqOfEvents(:,1));
    seqOfEvents = seqOfEvents(indxSort,:);
    
    %collect coordinates and amplitudes in tracksCoordAmpCG
    tracksCoordAmpCG = NaN(size(tracksFeatIndxCG,1),8*lastFrame); %fix KJ
    for iFrame = 1 : lastFrame %fix KJ
        clusterIndx = tracksFeatIndxCG(:,iFrame);
        clusterIndxGood = find(clusterIndx~=0);
        receptorIndx = clust2recept(clusterIndx(clusterIndxGood),:,iFrame);
        coordMat = receptorTraj(receptorIndx(:,1),:,iFrame);
        intensityMat = zeros(size(receptorIndx,1),1);
        for iReceptor = 1 : size(receptorIndx,1)
            indxTmp = receptorIndx(iReceptor,:);
            indxTmp = indxTmp(indxTmp~=0);
            intensityMat(iReceptor) = sum(receptorIntensity(indxTmp,iFrame));
        end
        tracksCoordAmpCG(clusterIndxGood,(iFrame-1)*8+1:(iFrame-1)*8+dimension) = coordMat;
        tracksCoordAmpCG(clusterIndxGood,(iFrame-1)*8+dimension+1) = 0;
        tracksCoordAmpCG(clusterIndxGood,(iFrame-1)*8+4) = intensityMat;
        tracksCoordAmpCG(clusterIndxGood,(iFrame-1)*8+5:(iFrame-1)*8+8) = 0;
    end
    
    %store this compound track's information
    %082014 - Robel Yirdaw
    %Saving tracksFeatIndxCG and tracksCoordAmpCG as sparse matrices
    compTracks(iClust).tracksFeatIndxCG = sparse(tracksFeatIndxCG);
    tracksCoordAmpCG(isnan(tracksCoordAmpCG)) = 0;
    compTracks(iClust).tracksCoordAmpCG = sparse(tracksCoordAmpCG);
    compTracks(iClust).seqOfEvents = seqOfEvents;

end %(for iClust = 1 : numClustGlobal)

%% ~~~ the end ~~~

