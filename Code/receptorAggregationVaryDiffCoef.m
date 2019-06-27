function [receptorInfoAll,receptorInfoLabeled,timeIterArray,errFlag] ...
    = receptorAggregationVaryDiffCoef(modelParam,simParam)
%RECEPTORAGGREGATIONVARYDIFFCOEF simulates the free diffusion and interactions of molecules in 2D, where diffusion coefficient varies between molecules
%
%SYNOPSIS [receptorInfoAll,receptorInfoLabeled,timeIterArray,errFlag] ...
%    = receptorAggregationVaryDiffCoef(modelParam,simParam)
%
%INPUT  modelParam: Structure with the fields:
%           diffCoef        : Mean diffusion coefficient (microns^2/s)
%           stdDiff         : Standard deviation of diffusion
%                                 coefficient, as fraction of the mean (unitless).
%           receptorDensity : Receptor density (# molecules/microns^probDim).
%           aggregationProb : Probability of association if a receptor
%                             bumps into another receptor or receptor
%                             complex.
%           aggregationDist : Distance between 2 receptors to consider them
%                             bumping into each other, thus potentially
%                             associating (microns).
%           dissociationRate: Rate of dissociation of a receptor from a
%                             receptor complex (/s).
%           labelRatio      : Receptor labeling ratio.
%           intensityQuantum: Row vector with 2 values, namely the mean and std of
%                             individual fluorophore intensity.
%           initPositions   : Receptor initial positions. If supplied, they
%                             will be used. If not, random positions will
%                             be chosen.
%       simParam: Structure with the fields:
%           probDim         : System dimensionality (1, 2 or 3). Default: 2.
%           observeSideLen  : Observation side length (microns). Either one
%                             value, used for all probDims, or a row
%                             vector with a value for each probDim.
%                             Default: 1 in all probDims.
%           timeStep        : Simulation time step (s).
%                             Default: 0.01/max(diffCoef,dissociationRate).
%           simTime         : Total time of simulation (s).
%                             Default: 100 * timeStep.
%           initTime        : Initialization time, for system to reach
%                             steady state.
%                             Default: 100 * timeStep.
%           randNumGenSeeds : Random number generator seed.
%                             Default: 100.
%                 Whole structure optional. Individual fields also optional.
%
%OUTPUT receptorInfoAll     : Structure with fields:
%           .receptorTraj        : (Number of receptors) - by - (probDim)
%                                  - by - (number of iterations) array storing
%                                  receptor positions.
%           .recept2clustAssign  : (Number of receptors) - by - (number of
%                                  iterations) array storing the cluster
%                                  assignment of each receptor.
%           .clust2receptAssign  : (Maximum number of clusters) - by - (maximum
%                                  cluster size) - by - (number of iterations)
%                                  array storing the receptor members of each
%                                  cluster.
%            .modelParam    : The imput model parameters
%            .simParam      : The imput simulation parameters
%       receptorInfoLabeled : The same as receptorInfoAll, but for labaled
%                             receptors only. It also has the additional
%                             field:
%           .compTracks          : The receptor trajectories and
%                                  interactions over time, as converted by
%                                  convReceptClust2CompTracks. Same format
%                                  as output of trackCloseGapsKalmanSparse
%                                  (or u-track for GUI version of tracker).
%       timeIterArray       : Vector storing time corresponding to each
%                             iteration.
%       errFlag             : 0 if function executes normally, 1 otherwise.
%
%Luciana de Oliveira, April 2019, based on receptorAggregationSimpleOptimized
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

%% NOTE:
%
%While the function was originally designed for 1D, 2D and 3D, most of
%the recent updates were with only 2D in mind. The function and simulation
%concepts have only been tested extensively in 2D. Therefore, please test
%extensively before running any 1D or 3D simulations.

%% Output
tic
errFlag = 0;
receptorInfoAll = [];
%09/05/14 (ryirdaw)
%need to block the following otherwise conversion from struct to double
%error at at the end
%receptorInfoLabeled = struct;
timeIterArray = [];
%assocStats = [];
%collProbStats = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modification LRO 2017/11/30
% Initiate the counts for the possible associations
% % % % % % % % % % % % % % % % % % % % totalPossibleAssPath=0;
totalPossAssSure=0;
% % % % % % % % % % % % % % % % % % % % % % % % % % totalPossAssPathCirc=0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Input

%check if correct number of arguments were used when function was called
if nargin < 1
    disp('--receptorAggregationSimple: Too few input arguments');
    errFlag  = 1;
    return
end

%extract model parameters from modelParam

%diffusion coefficient
if isfield(modelParam,'diffCoef')
    diffCoef = modelParam.diffCoef;
    
else
    disp('--receptorAggregationSimple: Please supply diffusion coefficient');
    errFlag = 1;
end


%add a standard deviation for diffusion coeficient
if isfield(modelParam,'stdDiff')
    stdDiff = modelParam.stdDiff;
    
else
    disp('--receptorAggregationSimple: Please supply std diffusion coefficient');
    errFlag = 1;
end

%receptor density
if isfield(modelParam,'receptorDensity')
    receptorDensity = modelParam.receptorDensity;
else
    disp('--receptorAggregationSimple: Please supply receptor density');
    errFlag = 1;
end

%aggregation probability
if isfield(modelParam,'aggregationProb')
    aggregationProb = modelParam.aggregationProb;
else
    disp('--receptorAggregationSimple: Please supply aggregation probability');
    errFlag = 1;
end

%aggregation distance
if isfield(modelParam,'aggregationDist')
    aggregationDist = modelParam.aggregationDist;
else
    disp('--receptorAggregationSimple: Please supply aggregation distance');
    errFlag = 1;
end

%dissociation rate
if isfield(modelParam,'dissociationRate')
    dissociationRate = modelParam.dissociationRate;
else
    disp('--receptorAggregationSimple: Please supply dissociation rate');
    errFlag = 1;
end

%labeling ratio
if isfield(modelParam,'labelRatio')
    labelRatio = modelParam.labelRatio;
else
    disp('--receptorAggregationSimple: Please supply labeling ratio');
    errFlag = 1;
end

%receptor initial positions
if isfield(modelParam,'initPositions')
    initPositions = modelParam.initPositions;
else
    initPositions = [];
end

%intensity quantum
if isfield(modelParam,'intensityQuantum')
    intensityQuantum = modelParam.intensityQuantum;
else
    disp('--receptorAggregationSimple: Please supply intensity quantum');
    errFlag = 1;
end

%exit if there are missing model parameters
if errFlag == 1
    disp('--receptorAggregationSimple: Please supply missing variables!');
    return
end

%check model parameters ...

%some must be positive
%09/05/14 - modified to accomodate a vector labelRatio
if any([receptorDensity aggregationDist (labelRatio(:)') intensityQuantum(1)] <= 0)
    disp('--receptorAggregationSimple: Receptor density, aggregation distance, labeling ratio and intensity quantum should be positive');
    errFlag = 1;
    return
end


%and some must be non-negative
%03/25/14 - modified to accomodate a vector aggregationProb
if any([diffCoef (aggregationProb') dissociationRate intensityQuantum(2)] < 0)
    disp('--receptorAggregationSimple: Diffusion coefficient, aggregation probability and dissociation rate should be non-negative');
    errFlag = 1;
    return
end

%extract simulation parameters from simParam

%if simParam wasn't supplied at all
if nargin < 2 || isempty(simParam)
    
    probDim = 2;
    observeSideLen = ones(1,probDim);
    timeStep = 0.01 / max(diffCoef,dissociationRate);
    simTime = 100 * timeStep;
    initTime = 100 * timeStep;
    randNumGenSeeds = [100 100];
    
else
    
    %system probDimality
    if isfield(simParam,'probDim')
        probDim = simParam.probDim;
    else
        probDim = 2;
    end
    
    %observation side length
    if isfield(simParam,'observeSideLen')
        observeSideLen = simParam.observeSideLen;
        if length(observeSideLen) == 1
            observeSideLen = observeSideLen * ones(1,probDim);
        end
    else
        observeSideLen = ones(1,probDim);
    end
    
    %time step
    if isfield(simParam,'timeStep')
        timeStep = simParam.timeStep;
    else
        timeStep = 0.01 / max(diffCoef,dissociationRate);
    end
    
    %simulation time
    if isfield(simParam,'simTime')
        simTime = simParam.simTime;
    else
        simTime = 100 * timeStep;
    end
    
    %initialization time
    if isfield(simParam,'initTime')
        initTime = simParam.initTime;
    else
        initTime = 100 * timeStep;
    end
    
    %random number generator seeds
    if isfield(simParam,'randNumGenSeeds')
        randNumGenSeeds = simParam.randNumGenSeeds;
    else
        randNumGenSeeds = [100 100];
    end
    
end

%determine number of iterations to perform
totalTime = simTime + initTime;
numIterSim = ceil( simTime / timeStep ) + 1;
numIterInit = ceil( initTime / timeStep );
numIterations = numIterSim + numIterInit;

%store time correponding to each iteration
timeIterArray = (0 : numIterSim - 1)' * timeStep;


%calculate dissociation probability
dissociationProb = dissociationRate * timeStep;

%initialize random number generators
% rand('twister',randNumGenSeeds(1));
% randn('state',randNumGenSeeds(2));
rng(randNumGenSeeds(1),'twister')

%% receptor initial positions and clustering

%calculate observation region size
obsRegionSize = prod(observeSideLen);

if isempty(initPositions)
    
    %calculate number of receptors
    numReceptors = round(obsRegionSize * receptorDensity);
    
    %initialize receptor positions
    initPositions = rand(numReceptors,probDim) .* repmat(observeSideLen,numReceptors,1);
    
else
    
    %get number of receptors
    numReceptors = size(initPositions,1);
    
end

%%%%%%%%%LRO 2019/04/30%%%%%%%%%%%%%%%
%%%Modificaiton: the diffusion coefficient is a distribution and each
%%%receptor have a different value

diffCoef=max(normrnd(diffCoef, diffCoef*stdDiff,[numReceptors,1]),eps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%assuming that the displacement taken in a timeStep is normally distributed
%with mean zero, get standard deviation of step distribution from diffusion
%coefficient
% modification: now the stepStd is an individual value for each
% receptor
stepStd = sqrt( 2 * diffCoef * timeStep );

%Modification LRO 2019/04/30 - since it has an individual value for
%each receptor, for the aggregation distance uses the average value
averageStepStd=mean(stepStd);

%adjust aggregationDist to account for the finite simulation time step and
%the expected receptor displacement in that time step
aggregationDistCorr = max(aggregationDist,sqrt(probDim)*averageStepStd*2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Starting with all receptors as monomers.  Set cluster2receptor and
%receptor2cluster as 1D vector of 1:numReceptors.  Then clusterSize will be
%numReceptors long with all values of 1. Note all are column vectors.
receptor2cluster = (1:numReceptors)';
cluster2receptor = (1:numReceptors)';
clusterSize = ones(numReceptors,1);

[numClusters,maxClustSize] = size(cluster2receptor);

%% Main simulation body

%reserve memory for output vectors
receptorTraj = zeros(numReceptors,probDim,numIterations);
recept2clustAssign = zeros(numReceptors,numIterations);
clust2receptAssign = zeros(numReceptors,maxClustSize,numIterations);

%store initial information
receptorTraj(:,:,1) = initPositions;
recept2clustAssign(:,1) = receptor2cluster;
clust2receptAssign(1:numClusters,1:maxClustSize,1) = cluster2receptor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Diagnostic quantities (ryirdaw)
assocStats = struct('numSureAssoc',NaN(numIterations,1),...
    'numPotColl',NaN(numIterations,1),'numColl',NaN(numIterations,1),...
    'numPotColl_Assoc',NaN(numIterations,1),...
    'sureAssocCountBySize',NaN(numReceptors,numIterations),...
    'numCollProbPairs',NaN(numIterations,1),...
    'numPathAssoc',NaN(numIterations,1),'numPathCirAssoc',NaN(numIterations,1),'numberPossibleAssociations',NaN(numIterations,1),'numberPossibleAssociationsPath',NaN(numIterations,1));


collProbStatStruct = struct('collisionProb',NaN,'pwDist',NaN,...
    'primaryNodeRadius',NaN,'partnerNodeRadii',NaN);
collProbStats(numIterations,1) = collProbStatStruct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

progressText(0,'Simulation');

%iterate in time
for iIter = 2 : numIterations
    %fprintf('\niIter = %d\n',iIter);
    
    %% Dissociation
    
    %allow receptors in clusters to dissociate in current time point
    [cluster2receptor,newReceptor2Cluster,clusterSize] = receptorDissociationAlg(...
        cluster2receptor,receptor2cluster,clusterSize,dissociationProb);
    %end
    
    %% Association flags based on outcome of dissociation
    
    %initialize
    aggregationProbVec = ones(numReceptors,1);
    
    %09/05/13 (ryirdaw)
    %The call to receptorDissociationAlg above could have dissociated
    %receptors from existing clusters.  These free receptors should not be
    %allowed to associate in receptorAggregationAlg below. Set the
    %aggregation probability for these receptors to zero.
    %Note: receptor2cluster vector from receptorDissociationAlg is being
    %recieved now (above, newReceptor2Cluster).
    
    %KJ, 2018/11/07: The above description is erroneous. Both dissociating
    %receptors and the clusters they come from should not associate in this
    %step. Check actual implementation.
    prevReceptor2cluster=receptor2cluster;
    
    if (max(newReceptor2Cluster) > max(receptor2cluster))
        %A dissociation has occured. To confirm and identify the
        %receptors involved, determine the cluster size of each receptor
        %now and compare with size values from previous iterations. The
        %best way to do this is to tally the number of other receptors each
        %receptor is associated with.
        %Store the new and previous sizes on two columns for each receptor.
        
        % Dev 2018-Nov-29: Changed cluster-receptor mappings to use
        % accumarray.
        
        accumulatedNewClusters=accumarray(newReceptor2Cluster,1,[length(newReceptor2Cluster),1]);
        newClusterSizes=accumulatedNewClusters(newReceptor2Cluster);
        
        accumulatedPrevClusters=accumarray(prevReceptor2cluster,1,[length(prevReceptor2cluster),1]);
        prevClusterSizes=accumulatedPrevClusters(prevReceptor2cluster);
        
        
        
        sizeNewPrev_=[newClusterSizes prevClusterSizes];
        sizeNewPrev=sizeNewPrev_;
        
        %For those receptors who have dissociated set the
        %aggregationProbVec to 0.  NOTE: if the other receptors remain
        %clustered, their aggregationProbVec must stay as 1.
        aggregationProbVec( (sizeNewPrev(:,1) - sizeNewPrev(:,2) < 0) ) = 0;
        
        %Reassign receptor2cluster to the new set reflecting the
        %dissociation.  NOTE: the position vector also shows the
        %dissociation that has occured.
        receptor2cluster = newReceptor2Cluster;
        
    end
    
    %% New receptor/cluster positions
    
    %get indices of clusters with more than one receptor
    clustersBig = find(clusterSize>1);
    
    %get receptor positions at previous time point
    positionsOld = receptorTraj(:,:,iIter-1);
    
    %generate receptor displacements
    receptorDisp = stepStd.*randn(numReceptors,probDim);
    
    %assign receptors in a cluster the displacement of the receptor with
    %the smallest index
    for iCluster = clustersBig'
        
        %get receptors belonging to this cluster
        clusterMembers = cluster2receptor(iCluster,1:clusterSize(iCluster));
        
        %%%%%%%%%%%%Modification LRO-2019/04/30%%%%%%%%%%%%%%%%%%%
        % incluede here also the colletive value of diffusion coeficient
        % and update the std for those receptors
        
        
        stepStdCluster=mean(stepStd(clusterMembers));
        
        clustDisp=stepStdCluster.*randn(1,probDim);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %assign all receptors in this cluster the displacement of that
        %receptor
        receptorDisp(clusterMembers,:) = repmat(clustDisp,...
            clusterSize(iCluster),1);
        
    end
    
    %calculate the new receptor positions at the current time point
    positionsNew = positionsOld + receptorDisp;
    
    %make sure that receptors stay inside the region of interest
    correctionBoundaryLow = min(positionsNew,0);
    positionsNew = positionsNew - 2 * correctionBoundaryLow;
    correctionBoundaryUp = max(positionsNew - repmat(observeSideLen,numReceptors,1),0);
    positionsNew = positionsNew - 2 * correctionBoundaryUp;
    
    %% Association
    
    try
        
        numClustPre = length(cluster2receptor(:,1));
        
        [cluster2receptor,receptor2cluster,clusterSize,receptPositions,...
            aggregationProbVec,sureAssocCount, numberPossibleAssociations] = ...
            receptorAggregationAlg_maxWeightedMatching_sureCollOptimized(positionsNew,...
            aggregationDistCorr,aggregationProbVec,aggregationProb,receptor2cluster,...
            cluster2receptor,clusterSize);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Modification LRO 2017/11/30 to calculate the number of
        % possible associations now I add a new output for the function
        totalPossAssSure=totalPossAssSure+numberPossibleAssociations;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        assocStats.numberPossibleAssociations(iIter,1)=numberPossibleAssociations;
        
        numClustPost = length(cluster2receptor(:,1));
        assocStats.numSureAssoc(iIter,1) = numClustPre - numClustPost;
        
        assocStats.sureAssocCountBySize(1:length(sureAssocCount),iIter) = sureAssocCount;
        
        clear numClustPre numClustPost sureAssocCount
        
    catch newAssocFunExcep
        fprintf('\nError at iIter = %d\n',iIter);
        disp(newAssocFunExcep.message);
        pause;
        return;
    end
    
    [numClusters,maxClustSize] = size(cluster2receptor);
    
    %store new receptor information
    receptorTraj(:,:,iIter) = positionsNew;
    recept2clustAssign(:,iIter) = receptor2cluster;
    clust2receptAssign(1:numClusters,1:maxClustSize,iIter) = cluster2receptor;
    
    progressText((iIter-1)/(numIterations-1),'Simulation');
    
end %(for iIter = 2 : numIterations)

%% Post-processing

%remove the initialization period from simulation
receptorTraj = receptorTraj(:,:,numIterInit+1:end);
recept2clustAssign = recept2clustAssign(:,numIterInit+1:end);
clust2receptAssign = clust2receptAssign(:,:,numIterInit+1:end);

%remove empty rows and columns from clust2receptAssign
cluster2receptor = max(clust2receptAssign,[],3);
columnSum = sum(cluster2receptor);
clust2receptAssign = clust2receptAssign(:,columnSum~=0,:);
rowSum = sum(cluster2receptor,2);
clust2receptAssign = clust2receptAssign(rowSum~=0,:,:);

% %put receptor trajectories and clusters into the format of the output of
% %trackCloseGapsKalman
% compTracks = convReceptClust2CompTracks(clust2receptAssign,...
%     recept2clustAssign,receptorTraj);

%put information in receptorInfoAll
receptorInfoAll = struct('receptorTraj',receptorTraj,'recept2clustAssign',...
    recept2clustAssign,'clust2receptAssign',clust2receptAssign,'simParam',simParam,'modelParam',modelParam);

elapsedTime = toc;
%% Receptor labeling and sub-sampling
tic
%KJ (150528): call function to label and sub-sample
receptorInfoLabeled = genReceptorInfoLabeled(receptorInfoAll,...
    labelRatio,intensityQuantum);
elapsedTime = toc;
end


%% ~~~ the end ~~~


%% subfunction

function [cluster2receptor,receptor2cluster,clusterSize] = receptorDissociationAlg(...
    cluster2receptor,receptor2cluster,clusterSize,dissociationProb)

%get number of receptors and number of clusters
numReceptors = length(receptor2cluster);
numClusters = length(clusterSize);
maxClusterSize = max(clusterSize);

%copy some input variables for modification and output
cluster2receptorTmp = [cluster2receptor; zeros(numReceptors-numClusters,maxClusterSize)];
receptor2clusterTmp = receptor2cluster;
clusterSizeTmp = [clusterSize; zeros(numReceptors-numClusters,1)];

%find clusters with more than one receptor
clustersBig = find(clusterSize > 1);

%06/27/13 (1 of 2)
%{
dissociationProb used on clusters instead of receptors (above).
Probability can be determined in two ways:  a single value for all current
clusters or each cluster has its own value. For a dissociating cluster, a
receptor will be randomly picked and removed from the cluster (below).
Dissociation happens one receptor at a time.
%}
%Each cluster has its own value
clusterDissociateFlag = rand(numClusters,1) < dissociationProb;

%go over these clusters
numClustersInit = numClusters; %this is for storing receptors that dissociate
for iCluster = clustersBig'
    
    %get receptors belonging to this cluster
    clusterMembers = cluster2receptor(iCluster,1:clusterSize(iCluster));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %06/27/13 (2 of 2)
    %For a cluster that is dissociating, determine which receptor to
    %dissociate by randomly picking one of the receptors.
    dissociateFlag = zeros(numel(clusterMembers),1);
    if (clusterDissociateFlag(iCluster))
        %Current cluster is dissociating. Pick a random integer between 1
        %and number of receptors in cluster.
        recept2Dissociate = randi(numel(clusterMembers),1);
        %Set the flag for picked receptor to 1.
        dissociateFlag(recept2Dissociate) = 1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %if all receptors want to dissociate, assign the first one a flag
    %of 0 (just for algorithmic reasons)
    if all(dissociateFlag)
        dissociateFlag(1) = 0;
    end
    
    %remove receptors that want to dissocate from cluster
    receptorsReject = clusterMembers(dissociateFlag==1);
    
    % Dev 2018 Nov 29: Modified approach to use ismember instead of
    % setdiff.
    if ( ~isempty(clusterMembers) && ~isempty(receptorsReject))
        clusterMembersTrim_ = sort(clusterMembers(~ismember(clusterMembers,receptorsReject)));
    else
        clusterMembersTrim_ = sort(clusterMembers);
    end
    numMembers_ = length(clusterMembersTrim_);
    
    numMembers=numMembers_;
    clusterMembersTrim=clusterMembersTrim_;
    
    
    %append to the end of the cluster vector those receptors that
    %dissociated
    numClustersFinal = numClustersInit + length(receptorsReject);
    cluster2receptorTmp(numClustersInit+1:numClustersFinal,1) = receptorsReject;
    clusterSizeTmp(numClustersInit+1:numClustersFinal) = 1;
    
    %keep the receptors that did not dissociate in their proper cluster
    cluster2receptorTmp(iCluster,:) = 0;
    cluster2receptorTmp(iCluster,1:numMembers) = clusterMembersTrim';
    clusterSizeTmp(iCluster) = numMembers;
    
    %update vector storing for every receptor its cluster number
    receptor2clusterTmp(receptorsReject) = (numClustersInit+1:numClustersFinal)';
    numClustersInit = numClustersFinal;
    
end

%remove empty row and columns and get final number of clusters
columnSum = sum(cluster2receptorTmp);
cluster2receptorTmp = cluster2receptorTmp(:,columnSum~=0);
rowSum = sum(cluster2receptorTmp,2);
cluster2receptorTmp = cluster2receptorTmp(rowSum~=0,:);
clusterSizeTmp = clusterSizeTmp(rowSum~=0);

%copy temporary variable into output variables
cluster2receptor = cluster2receptorTmp;
receptor2cluster = receptor2clusterTmp;
clusterSize = clusterSizeTmp;


end
