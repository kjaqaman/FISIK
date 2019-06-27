function [tracksFeatIndxCG,tracksCoordAmpCG] = pairwiseSegTracks(tracksFeatIndxCG,tracksCoordAmpCG,cellMergeSplitInfo)
% This function pairwise the segments from simultaneos merges and splits
% artifacs
% SYNOPSIS
%
% tracksReform = pairwiseSegTracks( tracksFeatIndxCG,tracksCoordAmpCG, segmentMerge,segmentSplit)
%
% INPUT
%        tracksFeatIndxCG : tracks intensity coming from the specific track
%        number
%        tracksCoordAmpCG : tracks coordinates coming from the specific track
%        number
%        segmentMerge     : segment that is merging in the merge event to
%        the contineous segment
%        segmentSplit     : segment that is spliting in the split event
%        from the contineous segment
% OUTPUT
%          tracksReform     : tracks after pairwesing the merges and the
%          splits, creating longer segments
%
%
% Luciana de Oliveira, January 2018 
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
              %%%%%%%%%%%%%%%% Modification LRO 2019/04/30%%%%%%%%%%%%%
              % to acomodate multiple events in the same time, I add the
              % condition of reforming only after checking all the events
              %call values from cell Merge Split Info
for indexInfo=1:size(cellMergeSplitInfo,1)
    for indexInfo2=1:size(cellMergeSplitInfo,2)
        
              
segmentMerge=cellMergeSplitInfo{indexInfo,indexInfo2}(1);
segmentSplit=cellMergeSplitInfo{indexInfo,indexInfo2}(2);
                % take the rows for the segment that merge and that split
                tracksFeatIndxCGMerge=tracksFeatIndxCG(segmentMerge,:);
                tracksFeatIndxCGSplit=tracksFeatIndxCG(segmentSplit,:);
                
                tracksCoordAmpCGMerge=tracksCoordAmpCG(segmentMerge,:);
                tracksCoordAmpCGSplit=tracksCoordAmpCG(segmentSplit,:);
                
                % put together the two piece of tracks, that is, take the non nan
                % info from the segment that split and replace it in the segment that
                % merge.
                                              
                tracksFeatIndxCGMerge=max(tracksFeatIndxCGMerge,tracksFeatIndxCGSplit);
                tracksCoordAmpCGMerge=max(tracksCoordAmpCGMerge,tracksCoordAmpCGSplit);
                
                % replace the merge segment for the new combined one and remove the
                % split segment
                tracksFeatIndxCG(segmentMerge,:)=tracksFeatIndxCGMerge;
                tracksCoordAmpCG(segmentMerge,:)= tracksCoordAmpCGMerge;
                
                % replace by NaNs the segment split
                tracksFeatIndxCG(segmentSplit,:)=NaN;
                tracksCoordAmpCG(segmentSplit,:)= NaN;
                
    end
end
% Remove rows with only Nans
        tracksFeatIndxCG(~any(~isnan(tracksFeatIndxCG), 2),:)=[];
        tracksCoordAmpCG(~any(~isnan( tracksCoordAmpCG), 2),:)=[];
end

