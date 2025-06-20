function [receptorInfoLabeledFixedID] = genReceptorInfoLabeledFixedID(receptorInfoAll,labelRatio,...
    intensityQuantum,compTrackFlag, fracNumFluor, labelOverlap)
%GENRECEPTORINFOLABELEDFIXEDID sub-samples simulated receptor trajectories and outputs information for labeled subset, including compound tracks
%
%SYNPOSIS receptorInfoLabeledFixedID = genReceptorInfoLabeledFixedID(receptorInfoAll,...
%    labelRatio,intensityQuantum,compTrackFlag,fracNumFluor)
%
%INPUT  receptorInfoAll : Output of receptorAggregationSimpleOptimized, or
%                         receptorAggregationVEGFR2ModeDistCorr.
%       labelRatio      : Vector of labeling ratios. Each entry corresponds
%                         to a different channel.
%       intensityQuantum: Row vector with 2 values, namely the mean and std of
%                         individual fluorophore intensity.
%       compTrackFlag    : compTrackFlag = 1 indicates
%                         user wants results with compound tracks,
%                         compTrackInd = 0 indicates user does not want
%                         compund tracks. 
%                         Optional. Default: 1.
%       fracNumFluor    : 3D array gives the fraction of labeled
%                         receptors that are tagged with fluorophores equal
%                         to the row value (iterate through the columns).
%                         Optional. Default: 1.
%       labelOverlap    : 1 to allow labeling overlap between channels, 0
%                         otherwise. ONLY needed when labeling more than
%                         one channel, i.e. length(labelRatio) >= 2.
%                         Optional. Default: 0.
%
%OUTPUT receptorInfoLabeledFixedID: Structure array with number of elements =
%                         number of different labeling ratios. For each
%                         labeling ratio, fields are similar to
%                         receptorInfoAll, but for the labeled
%                         receptors only. It has two additional fields:
%               .compTracks : The receptor trajectories and
%                             interactions over time, as converted by
%                             convReceptClust2CompTracks. Same format
%                             as output of trackCloseGapsKalmanSparse
%                             (or u-track for GUI version of tracker).
%               .labelRatio : The labeling ratio used to subsample.
%
%Code started by Robel Yirdaw in 2014, initially as a copy-paste from
%receptorAggregationSimple_new.
%
%Modified May 2015, Khuloud Jaqaman
%Modified May 2020 by Zachariah Malik for matrix fracNumFluor
%Modified August 2020 by Zachariah Malik to generate
%receptorInfoLabeledNoCT
%%%% NOTE: Final few issues; receptorInfoLabeledNoCT with incomplete
%%%% labeling is not perfect, if you plug it into
%%%% clusterHistoryFromReceptorInfo you will observe invisible
%%%% interactions.
% % Modified, Dec 2022 by Fnu Bilal, for fixed receptor ID
%Mofified Oct 2023, Jesus Vega-Lugo to allow the labeling of multiple
%channels without overlapping labels i.e. one recepto getting labeled on both channels
%
%NOTE: NO LONGER SURE WHAT fracNumFluor DOES OR WHAT IT MEANS.
%
% Copyright (C) 2025, Jaqaman Lab - UTSouthwestern
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


%% Input/Output

if nargin < 4 || isempty(compTrackFlag)
    compTrackFlag = 1;
end

if nargin < 5 || isempty(fracNumFluor)
    fracNumFluor = 1;
end

if nargin < 5 || isempty(labelOverlap)
    labelOverlap = 0;
end
% Compute these values later; to be used if using parfor
%maxNumFluor = length(fracNumFluor);
%fluorCumSum = cumsum(fracNumFluor);

%Modification LRO 20190603-include fixed random seeds
% LRO 2019/07/18 just for testing of profiling, use only one random seed
% rng(1,'twister')

timeStep = receptorInfoAll.simParam.timeStep;
receptorTraj = receptorInfoAll.receptorTraj;
recept2clustAssign = receptorInfoAll.recept2clustAssign;
clust2receptAssign = receptorInfoAll.clust2receptAssign;

%20200720 Obtain diffHist
if isfield(receptorInfoAll,'diffHist')
    diffHist = receptorInfoAll.diffHist;
else
    diffHist = [];
end

[numReceptors,~,numIterSim] = size(receptorTraj);

%09/05/14 (ryirdaw) - modified to allow a vector labelRatio
%numLabelRatio = length(labelRatio);
%receptorInfoLabeled(numLabelRatio,1) = struct('receptorTraj',[],...
%     'recept2clustAssign',[],...
%     'clust2receptAssign',[],...
%     'compTracks',[],...
%     'labelRatio',[])

%20200716 Include diffHist cell array (ZMALIK)
% numLabelRatio = length(labelRatio);
% numFracFluor = length(fracNumFluor(1,1,:));
% receptorInfoLabeledWithCT(numLabelRatio*numFracFluor,1) =...
%     struct('receptorTraj',[],...
%     'recept2clustAssign',[],...
%     'clust2receptAssign',[],...
%     'compTracks',[],...
%     'labelRatio',[],...
%     'fracRecept2NumFluor',[],...
%     'indxLabeled',[]);
% %20200813 New output
% receptorInfoLabeledNoCT(numLabelRatio*numFracFluor,1) =...
%     struct('receptorTraj',[],...
%     'recept2clustAssign',[],...
%     'clust2receptAssign',[],...
%     'labelRatio',[],...
%     'fracRecept2NumFluor',[],...
%     'diffHist',[]);

%20230217 receptor label with or without compound tracks (Bilal)
numLabelRatio = length(labelRatio);
numFracFluor = length(fracNumFluor(1,1,:));
receptorInfoLabeledFixedID(numLabelRatio*numFracFluor,1) =...
    struct('receptorTraj',[],'receptorIntensity',[],...
    'recept2clustAssign',[],'clust2receptAssign',[],...
    'clustPosition',[],'clustIntensity',[],...
    'compTracks',[],...
    'labelRatio',[],'indxLabeled',[],...
    'fracRecept2NumFluor',[],'diffHist',[]);

%% Labeling and sub-sampling
% Initiatlize receptorInfoLabeled index

%%%%%%%%%%%%

%Nov 2023: Jesus Vega-Lugo added to avoid labeling overlap when labeling 2 or more
%channels. It also guaratees getting the expected number of labeled
%receptors based on the labelRatio.

%assign labels randomly
labelFlag = rand(numReceptors,numLabelRatio) <= labelRatio;

indxLabeled = cell(1,numLabelRatio);
indxNotLabeled = cell(1,numLabelRatio);

for iCh = 1:numLabelRatio
    
    indxLabeled{iCh} = find(labelFlag(:,iCh)==1);
    indxNotLabeled{iCh} = find(labelFlag(:,iCh)==0);
    
    expRecepLabeled = round(numReceptors*labelRatio(iCh));
    numRecepLabeled = numel(indxLabeled{iCh});
    
    %make sure you get the expected number of labeled receptors per channel
    if numRecepLabeled > expRecepLabeled

        numRecepToElim = numRecepLabeled - expRecepLabeled;
        idxToElim = randperm(numRecepLabeled,numRecepToElim);
        indxNotLabeled{iCh} = sort(vertcat(indxNotLabeled{iCh},indxLabeled{iCh}(idxToElim)));
        indxLabeled{iCh}(idxToElim) = [];
        
    elseif numRecepLabeled < expRecepLabeled
        
        numRecepToAdd = expRecepLabeled - numRecepLabeled;
        idxToAdd = randperm(numel(indxNotLabeled{iCh}),numRecepToAdd);
       
        indxLabeled{iCh} = sort(vertcat(indxLabeled{iCh},indxNotLabeled{iCh}(idxToAdd)));
         indxNotLabeled{iCh}(idxToAdd) = [];
         
    end
  
end

%make sure there is no label overlap 
if numLabelRatio > 1 && ~labelOverlap
    for iCh = 1:numLabelRatio - 1 
       
        [~,~,idxLabelOverlap] = intersect(indxLabeled{iCh},indxLabeled{iCh+1});

        indxNotLabeled{iCh+1} = sort(vertcat(indxNotLabeled{iCh+1},indxLabeled{iCh+1}(idxLabelOverlap)));

        indxLabeled{iCh+1}(idxLabelOverlap) = [];
        newLabel = randsample(setdiff(indxNotLabeled{iCh},indxLabeled{iCh+1}),numel(idxLabelOverlap),false);
        indxLabeled{iCh+1} = sort(vertcat(indxLabeled{iCh+1},newLabel));

        [~,idxToElim,~] = intersect(indxNotLabeled{iCh+1},indxLabeled{iCh+1});
        indxNotLabeled{iCh+1}(idxToElim) = [];
       
    end
end

%%%%%%%%%%%%%%%%%%

rILindx = 1;
for fNFindx = 1:numFracFluor
    
    for iCh = 1:numLabelRatio
        
        maxNumFluor = length(fracNumFluor(fracNumFluor(:,1,fNFindx) > 0));
        fluorCumSum = cumsum(fracNumFluor(:,1,fNFindx));
        
        %if labeling ratio is less than one ...
        if labelRatio(iCh) < 1
            
            %Get receptorTraj for labeled receptors
            receptorTrajLabeled = receptorTraj(indxLabeled{iCh},:,:);
            
            %%%%%%%%ZMALIK 2020/March/10
            %%%%%%%%Now more than one fluorophore may tag a receptor
            %assign the intensities of the labeled receptors and initialize
            %variables
            %%%%%%%%ZMALIK 2020/May/07
            %%%%%%%%Adjusted for matrix fracNumFluor
            receptorNumFluorFlag = rand(length(indxLabeled{iCh}),numIterSim);
            receptorIntensityLabeled = randn(length(indxLabeled{iCh}),numIterSim)*...
                intensityQuantum(2);
            for iNumFluor = 1:maxNumFluor
                if iNumFluor == 1
                    receptorNumFluorFlag(receptorNumFluorFlag < ...
                        fluorCumSum(iNumFluor) & receptorNumFluorFlag >= 0) ...
                        = iNumFluor;
                else
                    receptorNumFluorFlag(receptorNumFluorFlag < ...
                        fluorCumSum(iNumFluor) & receptorNumFluorFlag >= ...
                        fluorCumSum(iNumFluor - 1)) = iNumFluor;
                end
                receptorIntensityLabeled(receptorNumFluorFlag == iNumFluor) = ...
                    receptorIntensityLabeled(receptorNumFluorFlag == iNumFluor) ...
                    + iNumFluor*intensityQuantum(1);
            end
            
            receptorIntensityLabeled(receptorIntensityLabeled < eps) = eps;
            
            % Modification 2022/12/28, Bilal, Store only the labeled receptors and update the
            % cluster2receptor and receptor2cluster accordingly
            
            %modify the cluster-to-receptor assignments to include only labeled
            %receptors
            clust2receptAssignLabeled = clust2receptAssign;
            
            % identify which are the receptors that are maintained after
            % labeling
            labelRecep = ismember(clust2receptAssignLabeled,  indxLabeled{iCh}');
            
            % replace by zero the receptors that are not labeled, Modified by Bilal
            
            clust2receptAssignLabeled=clust2receptAssignLabeled.*labelRecep;
            recept2clustAssignLabeled = recept2clustAssign;
            
            clust2receptAssignNotLabeled = clust2receptAssignLabeled(indxNotLabeled{iCh},:,:);
            
            for iIter = 1 : numIterSim
                clustSizeNotLabeled = sum(clust2receptAssignNotLabeled(:,:,iIter),2);
                
                indx = find(clustSizeNotLabeled~=0);
                
                clustToBack = clust2receptAssignNotLabeled(indx,:,iIter);  % move them together
                
                [sizeClustToBack, clusterSize]= size(clustToBack);
                
                for iCluster =1: sizeClustToBack
                    
                    %get receptors belonging to this cluster
                    clusterMembers = find(clustToBack(iCluster, 1:clusterSize)~=0);
                    
                    clust2receptAssignLabeled(clustToBack(iCluster, clusterMembers(1)),:,iIter) = [ ...
                        clustToBack(iCluster, clusterMembers) zeros(1, clusterSize-length(clusterMembers))];
                    
                    recept2clustAssignLabeled((clustToBack(iCluster, clusterMembers)),iIter) = clustToBack(iCluster, clusterMembers(1));
                    
                end
            end
            
            % Modofication, 2022/12/30, Bilal, update the labeled receptors and remove
            % the unlabeled once
            clust2receptAssignLabeled(indxNotLabeled{iCh},:) = 0;
            recept2clustAssignLabeled = recept2clustAssignLabeled(indxLabeled{iCh},:);
            
            % Modification 2023/02/09, Bilal, use only the labeled
            % receptors and remove the not labeled ones
            clust2receptAssignLabeled = clust2receptAssignLabeled(indxLabeled{iCh},:,:);
         
            for iLabeled = 1 : length(indxLabeled{iCh})
                clust2receptAssignLabeled(clust2receptAssignLabeled==indxLabeled{iCh}(iLabeled)) = iLabeled;
                recept2clustAssignLabeled(recept2clustAssignLabeled==indxLabeled{iCh}(iLabeled))= iLabeled;
            end
                   
            % eliminate the unlabeled receptors
            
            %make sure that zeros come after receptor indices (in both of
            %them)
            clust2receptAssignLabeled(clust2receptAssignLabeled==0)=NaN;
            clust2receptAssignLabeled = sort(clust2receptAssignLabeled,2);
            clust2receptAssignLabeled(isnan(clust2receptAssignLabeled))=0;
            
            %remove empty rows and columns from clust2receptAssign,
            %Bilal
            cluster2receptor = max(clust2receptAssignLabeled,[],3);
            columnSum = sum(cluster2receptor);
            clust2receptAssignLabeled = clust2receptAssignLabeled(:,columnSum~=0,:);
            
            %KJ 231211: Document cluster positions and intensity per frame
            clustPosition = receptorTrajLabeled;
            for iIter = 1 : numIterSim
                clustPosition(clust2receptAssignLabeled(:,1,iIter)==0,:,iIter) = NaN;
            end
            clustIntensity = zeros(size(receptorIntensityLabeled));
            for iIter = 1 : numIterSim
                tmp = clust2receptAssignLabeled(:,:,iIter);
                indxSomething = find(tmp~=0);
                intensitySomething = receptorIntensityLabeled(tmp(indxSomething),iIter);
                tmp(indxSomething) = intensitySomething;
                clustIntensity(:,iIter) = sum(tmp,2);
            end
            clustIntensity(clustIntensity==0) = NaN;
            
            %prepare output
            receptorInfoLabeledFixedID(rILindx).receptorTraj = receptorTrajLabeled;
            receptorInfoLabeledFixedID(rILindx).receptorIntensity = receptorIntensityLabeled;
            
            receptorInfoLabeledFixedID(rILindx).recept2clustAssign = recept2clustAssignLabeled;
            receptorInfoLabeledFixedID(rILindx).clust2receptAssign = clust2receptAssignLabeled;
            
            receptorInfoLabeledFixedID(rILindx).clustPosition = clustPosition;
            receptorInfoLabeledFixedID(rILindx).clustIntensity = clustIntensity;
            
            receptorInfoLabeledFixedID(rILindx).labelRatio = labelRatio(iCh);
            receptorInfoLabeledFixedID(rILindx).fracRecept2NumFluor = fracNumFluor(:,fNFindx);
            receptorInfoLabeledFixedID(rILindx).indxLabeled = indxLabeled{iCh};
            
            receptorInfoLabeledFixedID(rILindx).diffHist = []; %FIX THIS!!!!
            
            receptorInfoLabeledFixedID(rILindx).compTracks = [];
            
            %put labeled receptor trajectories and clusters into the format of the
            %output of trackCloseGapsKalman
            if compTrackFlag == 1
                compTracksLabeled = convReceptClust2CompTracksFixedID(clust2receptAssignLabeled,...
                    recept2clustAssignLabeled,receptorTrajLabeled,receptorIntensityLabeled);
                receptorInfoLabeledFixedID(rILindx).compTracks = compTracksLabeled;
            end
            
        else %if all receptors are labeled
            
            receptorNumFluorFlag = rand(numReceptors,numIterSim);
            receptorIntensity = randn(numReceptors,numIterSim)*...
                intensityQuantum(2);
            for iNumFluor = 1:maxNumFluor
                if iNumFluor == 1
                    receptorNumFluorFlag(receptorNumFluorFlag < ...
                        fluorCumSum(iNumFluor) & receptorNumFluorFlag >= 0) ...
                        = iNumFluor;
                else
                    receptorNumFluorFlag(receptorNumFluorFlag < ...
                        fluorCumSum(iNumFluor) & receptorNumFluorFlag >= ...
                        fluorCumSum(iNumFluor - 1)) = iNumFluor;
                end
                receptorIntensity(receptorNumFluorFlag == iNumFluor) = ...
                    receptorIntensity(receptorNumFluorFlag == iNumFluor) ...
                    + iNumFluor*intensityQuantum(1);
            end
            
            receptorIntensity(receptorIntensity < eps) = eps;
                        
            %KJ 231211: Document cluster positions and intensity per frame
            clustPosition = receptorTraj;
            for iIter = 1 : numIterSim
                clustPosition(clust2receptAssign(:,1,iIter)==0,:,iIter) = NaN;
            end
            clustIntensity = zeros(size(receptorIntensity));
            for iIter = 1 : numIterSim
                tmp = clust2receptAssign(:,:,iIter);
                indxSomething = find(tmp~=0);
                intensitySomething = receptorIntensity(tmp(indxSomething),iIter);
                tmp(indxSomething) = intensitySomething;
                clustIntensity(:,iIter) = sum(tmp,2);
            end
            clustIntensity(clustIntensity==0) = NaN;
            
            %prepare output
            receptorInfoLabeledFixedID(rILindx).receptorTraj = receptorInfoAll.receptorTraj;
            receptorInfoLabeledFixedID(rILindx).receptorIntensity = receptorIntensity;
            
            receptorInfoLabeledFixedID(rILindx).recept2clustAssign = receptorInfoAll.recept2clustAssign;
            receptorInfoLabeledFixedID(rILindx).clust2receptAssign = receptorInfoAll.clust2receptAssign;
            
            receptorInfoLabeledFixedID(rILindx).clustPosition = clustPosition;
            receptorInfoLabeledFixedID(rILindx).clustIntensity = clustIntensity;
            
            receptorInfoLabeledFixedID(rILindx).labelRatio = labelRatio(iCh);
            receptorInfoLabeledFixedID(rILindx).fracRecept2NumFluor = fracNumFluor(:,fNFindx);
            receptorInfoLabeledFixedID(rILindx).indxLabeled = (1:numReceptors)';
            
            receptorInfoLabeledFixedID(rILindx).diffHist = diffHist;
            
            receptorInfoLabeledFixedID(rILindx).compTracks = [];
            
            %put labeled receptor trajectories and clusters into the format of the
            %output of trackCloseGapsKalman
            if compTrackFlag == 1
                compTracksLabeled = convReceptClust2CompTracksFixedID(clust2receptAssign,...
                    recept2clustAssign,receptorTraj,receptorIntensity);
                receptorInfoLabeledFixedID(rILindx).compTracks = compTracksLabeled;
            end
            
        end
        
        % Increase the index
        rILindx = rILindx + 1;
        
    end % for each labelRatio element
    
end % for each column of fracNumFluor

end