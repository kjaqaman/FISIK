function receptorInfoLabeled = genReceptorInfoLabeledPhotobleaching(receptorInfoAll,...
    labelRatio,intensityQuantum,fracSurvive)
%GENRECEPTORINFOLABELEDPHOTOBLEACHING sub-samples simulated receptor trajectories, with exponential decay over time, and outputs information for labeled subset, including compound tracks
%
%SYNPOSIS receptorInfoLabeled = genReceptorInfoLabeledPhotobleaching(receptorInfoAll,...
%   labelRatio,intensityQuantum,fracSurvive)
%
%INPUT  receptorInfoAll : Output of receptorAggregationSimpleOptimized.
%       labelRatio      : Labeling ratio.
%       intensityQuantum: Row vector with 2 values, namely the mean and std of
%                         individual fluorophore intensity.
%       fracSurvive     : Fraction of surviving labels between consecutive
%                         time points. It is related to to the
%                         photobleaching rate, k, by the equation:
%                                   fracSurvive = exp(-kdt)
%                         where dt is the simulation time step.
%
%OUTPUT receptorInfoLabeled: Structure array with number of elements =
%                         number of different labeling ratios. For each
%                         labeling ratio, fields are similar to
%                         receptorInfoAll, but for the labaled
%                         receptors only. It has two additional fields:
%               .compTracks : The receptor trajectories and
%                             interactions over time, as converted by
%                             convReceptClust2CompTracks. Same format
%                             as output of trackCloseGapsKalmanSparse
%                             (or u-track for GUI version of tracker).
%               .labelRatio : The labeling ratio used to subsample.
%
%Luciana de Oliveira, April 2019, based on genReceptorInfoLabeled
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

%% NOTE
%Photobleaching currently implemented only for labelRatio < 1


%% Input/Output

receptorTraj = receptorInfoAll.receptorTraj;
recept2clustAssign = receptorInfoAll.recept2clustAssign;
clust2receptAssign = receptorInfoAll.clust2receptAssign;

[numReceptors,~,numIterSim] = size(receptorTraj);

%09/05/14 (ryirdaw) - modified to allow a vector labelRatio
numLabelRatio = length(labelRatio);
receptorInfoLabeled(numLabelRatio,1) = struct('receptorTraj',[],...
    'recept2clustAssign',[],...
    'clust2receptAssign',[],...
    'compTracks',[],...
    'labelRatio',[]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Labeling and sub-sampling

for lRindx=1:numLabelRatio
    % the first frame need to have all the receptors labeled from a
    % particular labeled fraction, the decay will start from the second
    % frame and will be just in the labeled receptors.
    
    %if labeling ratio is less than one ...
    if labelRatio(lRindx) < 1
        
        %Modification LRo 2019/05/08: labeled fraction will be dependent on time
        %calculate distribution for the first frame, that is taken from the
        %full distribution
        labelFlag = rand(numReceptors,1) <= labelRatio(lRindx);
        indxLabeled = find(labelFlag==1);
        indxNotLabeled = find(labelFlag==0);
        
        %extract the trajectories of the labeled receptors
        receptorTrajLabeled = receptorTraj(indxLabeled,:,:);
        
        %assign the intensities of the labeled receptors
        receptorIntensityLabeled = intensityQuantum(1) + ...
            randn(length(indxLabeled),numIterSim) * intensityQuantum(2);
        receptorIntensityLabeled(receptorIntensityLabeled < eps) = eps;
        
        %extract the cluster assignments of the labeled receptors
        recept2clustAssignLabeled = recept2clustAssign(indxLabeled,:);
        
        %modify the cluster-to-receptor assignments to include only labeled
        %receptors
        clust2receptAssignLabeled = clust2receptAssign;
        %convert receptors that are not labeled to zero
        for iNotLabeled = indxNotLabeled'
            clust2receptAssignLabeled(clust2receptAssignLabeled==iNotLabeled) = 0;
        end
        %update the indices of the labeled receptors
        for iLabeled = 1 : length(indxLabeled)
            clust2receptAssignLabeled(clust2receptAssignLabeled==indxLabeled(iLabeled)) = iLabeled;
        end
        %make sure that zeros come after receptor indices
        clust2receptAssignLabeled = sort(clust2receptAssignLabeled,2,'descend');
        
        %remove empty clusters (which include unlabeled receptors)
        %modify receptor-to-cluster assignments accordingly
        for iIter = 1 : numIterSim
            clustSize = sum(clust2receptAssignLabeled(:,:,iIter)~=0,2);
            indxFull = find(clustSize~=0);
            indxEmpty = find(clustSize==0);
            clust2receptAssignLabeled(:,:,iIter) = clust2receptAssignLabeled(...
                [indxFull;indxEmpty],:,iIter);
            for iFull = 1 : length(indxFull)
                recept2clustAssignLabeled(recept2clustAssignLabeled(:,iIter)...
                    ==indxFull(iFull),iIter) = iFull;
            end
        end
        
        %remove empty rows and columns from clust2receptAssign
        cluster2receptor = max(clust2receptAssignLabeled,[],3);
        columnSum = sum(cluster2receptor);
        clust2receptAssignLabeled = clust2receptAssignLabeled(:,columnSum~=0,:);
        rowSum = sum(cluster2receptor,2);
        clust2receptAssignLabeled = clust2receptAssignLabeled(rowSum~=0,:,:);
        
        %% Now calculate for all remaining frames
        
        numLabeled=length(indxLabeled );
        
        %start a new flag for those surviving photobleaching
        labelFlag=ones(numLabeled,1);
        
        for iterIndex=2:numIterSim-1
            
            %choose the surviving labeled receptors after photobleaching
            labelFlag = (rand(numLabeled,1) <= fracSurvive(lRindx)) .* labelFlag;
            indxNotLabeled = find(labelFlag==0);
            
            %replace by NaN the receptors that are not labeled from this
            %frame onwards
            receptorTrajLabeled(indxNotLabeled,:,iterIndex:end) = NaN; %%%%% LRO modification 05/07/2019 include the frame here
            
            %replace the intensity of receptors that are not labeled by
            %zero
            receptorIntensityLabeled(indxNotLabeled,iterIndex:end) = 0;
            
            %remove the cluster assignments of photobleached receptors
            recept2clustAssignLabeled(indxNotLabeled,iterIndex:end) = 0;  %%%%% LRO modification 05/07/2019 include the frame here
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%% Modification LRO modification 05/07/2019 include the frame here
            
            %take the clust2receptAssign from the current frame until the
            %end
            clust2receptAssignLabeledTmp=clust2receptAssignLabeled(:,:,iterIndex:end);
            
            %convert receptors that are not labeled to zero
            for iNotLabeled = indxNotLabeled'
                
                %replace by zero the receptors that are not labeled from
                %the current frame until the end
                clust2receptAssignLabeledTmp(clust2receptAssignLabeledTmp==iNotLabeled)=0;
                
            end
            
            % replace this information in the original clust2receptAssign
            clust2receptAssignLabeled(:,:,iterIndex:end) = clust2receptAssignLabeledTmp;
            
        end %for each time step
        
        %make sure that zeros come after receptor indices
        clust2receptAssignLabeled = sort(clust2receptAssignLabeled,2,'descend');
        
        convTracksStartTime = tic;
        
        %put labeled receptor trajectories and clusters into the format of the
        %output of trackCloseGapsKalman
        compTracksLabeled = convReceptClust2CompTracksPhotobleaching(clust2receptAssignLabeled,...
            recept2clustAssignLabeled,receptorTrajLabeled,receptorIntensityLabeled);
        
        convTracksETime = toc(convTracksStartTime);
        fprintf('\n============\nTime for convReceptClust2CompTracks is %g seconds. \n============\n',convTracksETime);
        
        
        receptorInfoLabeled(lRindx).receptorTraj = receptorTrajLabeled;
        receptorInfoLabeled(lRindx).recept2clustAssign = recept2clustAssignLabeled;
        receptorInfoLabeled(lRindx).clust2receptAssign = clust2receptAssignLabeled;
        receptorInfoLabeled(lRindx).compTracks = compTracksLabeled;
        receptorInfoLabeled(lRindx).labelRatio = labelRatio(lRindx);
        
    else %if all receptors are labeled
        
        receptorInfoLabeled(lRindx).receptorTraj = receptorInfoAll.receptorTraj;
        receptorInfoLabeled(lRindx).recept2clustAssign = receptorInfoAll.recept2clustAssign;
        receptorInfoLabeled(lRindx).clust2receptAssign = receptorInfoAll.clust2receptAssign;
        receptorInfoLabeled(lRindx).labelRatio = labelRatio(lRindx);
        
        %assign receptor intensities
        receptorIntensity = intensityQuantum(1) + ...
            randn(numReceptors,numIterSim) * intensityQuantum(2);
        receptorIntensity(receptorIntensity < eps) = eps;
        
        convTracksStartTime = tic;
        
        %put labeled receptor trajectories and clusters into the format of the
        %output of trackCloseGapsKalman
        compTracksLabeled = convReceptClust2CompTracks(clust2receptAssign,...
            recept2clustAssign,receptorTraj,receptorIntensity);
        
        convTracksETime = toc(convTracksStartTime);
        fprintf('\n============\nTime for convReceptClust2CompTracks is %g seconds. \n============\n',convTracksETime);
        
        receptorInfoLabeled(lRindx).compTracks = compTracksLabeled;
        
    end %(if labelRatio < 1 ... else ...)
    
end % for each labelRatio element

end