function [receptorInfoAll,receptorInfoLabeled,timeIterArray,errFlag] ...
    = receptorAggregationVEGFR2ModeDistCorrFixedID(modelParam,simParam)
%RECEPTORAGGREGATIONVEGFR2MODEDISTCORR simulates the diffusion and interactions of molecules in 2D, with molecules in distinct mobility modes and inter-dependence between interactions and mobility modes
%
%SYNOPSIS [receptorInfoAll,receptorInfoLabeled,timeIterArray,errFlag] ...
%    = receptorAggregationVEGFR2ModeDistCorrFixedID(modelParam,simParam)
%
%INPUT  modelParam: Structure with fields:
%           mobilityModeCell: (Number of modes)x(max cluster size + 1) cell
%                             array. Each row corresponds to a unique
%                             mobility mode. The first column gives
%                             the fraction of receptors living in that
%                             mode. Each column after that gives the
%                             diffusion coefficient distribution for that
%                             mode, per cluster size (oligomeric state).
%                             Diffusion coefficient in (micron^2/s).
%           receptorDensity : Receptor density (# molecules/microns^probDim).
%           aggregationProb : (Maximum cluster size + 1)x(number of
%                             mobility modes)x(number of mobility modes)
%                             array. Probability of association if a
%                             monomer bumps into another monomer or
%                             cluster. Dependent on cluster size and
%                             mobility mode. The first index refers to the
%                             final size of the newly formed cluster. The
%                             second index refers to the mobility mode of
%                             the associating cluster (size >= 1). The
%                             third index refers to the mobility mode of
%                             the associating monomer (size = 1).
%           aggregationDist : Distance between 2 receptors/clusters to
%                             consider them bumping into each other, thus
%                             poised for associating (microns).
%           dissociationRate: (maximum cluster size)x(number of mobility
%                             modes) array. Rate of dissociation of a
%                             monomer from a cluster (/s). For example,
%                             dissociationRate(2,1) is the rate of dissociation
%                             from a size 2 cluster in mobility mode 1.
%           labelRatio      : Receptor labeling ratio (fraction).
%           intensityQuantum: Row vector with 2 values, namely the mean and
%                             std of individual fluorophore intensity.
%           initPositions   : Receptor initial positions. If supplied, they
%                             will be used. If not, simulation will start
%                             with random positions.
%           mobModeProbAfterMerge : 4D array that gives the probability of
%                                   selection of a mobility mode based on:
%                                   The new cluster size and the mobility
%                                   modes of the two mergers. First index
%                                   is new cluster size. Second index is
%                                   the mode of the newly formed cluster.
%                                   The third index is the mode of the
%                                   merging cluster (size >= 1). The fourth
%                                   index is the mode of the merging
%                                   monomer (size = 1).
%           mobModeProbAfterSplit : 4D array that gives the probability of
%                                   selection of a mobility mode based on:
%                                   The old cluster size and its mobility.
%                                   First index is the size of the cluster
%                                   before splitting. Second index is the
%                                   mobility mode of the cluster before
%                                   splitting. The third index is the mode
%                                   of the split cluster (size >= 1).
%                                   Fourth index is the mode of the split
%                                   monomer (size = 1).
%           diffCoefScaleAfterMerge:4D array that scales the diffusion
%                                   coefficient of a cluster that was
%                                   formed as a merge. First index
%                                   is new cluster size. Second index is
%                                   the mode of the newly formed cluster.
%                                   The third index is the mode of the
%                                   merging cluster (size >= 1). The fourth
%                                   index is the mode of the merging
%                                   monomer (size = 1). Optional, if not
%                                   input, then the scale is set to one,
%                                   that is, diffusion coefficients remain
%                                   unchanged.
%
%           KJ 240229: I am not sure how relevant the below note from Zach is.
%           Probably safer to just follow the straightforward convention
%           that 1 is Mode 1 (restricted), 2 is Mode 2 (mobile), etc.
%           KJ 240304: Whichever mobility convention is followed, all
%           variables dealing with mobility must follow the same ordering
%           convention. This includes definition of mobility modes.
%           NOTE: The analysis code (as of 20200716) orders mobility mode
%           from most mobile to least mobile. It is highly recommended to
%           follow this convention when building modelParam.
%
%       simParam: Structure with fields:
%           probDim         : System dimensionality (1, 2 or 3).
%                             Default: 2.
%           observeSideLen  : Observation side length (microns). Either one
%                             value, used for all probDims, or a row
%                             vector with a value for each dimension.
%                             Default: 1 in all dimensions.
%           timeStep        : Simulation time step (s).
%                             Default: 0.01/max(diffCoef,dissociationRate).
%           simTime         : Total time of simulation (s).
%                             Default: 100 * timeStep.
%           initTime        : Initialization time, for system to reach
%                             steady state.
%                             Default: 100 * timeStep.
%           randNumGenSeeds : Random number generator seed.
%                             Default: 100.
%           compTrackFlag   : Flag (0/1) indicating whether compound tracks
%                             are constructed.
%                             Default: 1.
%           Whole structure optional. Individual fields also optional.
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
%            .diffHist      : 2x(number of time points) cell array showing
%                             the diffusion coefficients (row 1) and
%                             diffusion modes (row 2) of the simulated
%                             receptors at each time point. 
%       receptorInfoLabeled : The same as receptorInfoAll, but for labaled
%                             receptors only. If input compTrackFlag = 1,
%                             it then has the additional field:
%           .compTracks          : The receptor trajectories and
%                                  interactions over time, as converted by
%                                  convReceptClust2CompTracks. Same format
%                                  as output of trackCloseGapsKalmanSparse
%                                  (or u-track for GUI version of tracker).
%                             FOR NOW, THE FIELD diffHist IS LEFT EMPTY.
%       timeIterArray       : Vector storing time corresponding to each
%                             iteration.
%       errFlag             : 0 if function executes normally, 1 otherwise.
%
% 2020 April
% Zachariah Malik, based on receptorAggregationVEGFR2ConstClustLab
% 2020 January (ZMalik): Now can specify diffusion coefficient
% distributions and mobility modes for individual particles.
% 2020 February (ZMalik): Input PDF for picking probability mode after
% splitting instead of CDF.
% 2020 April (ZMalik) : Achieve consistent cluster labelling. Effective
% molecular radius is now specified per individual molecule and used for
% aggregation.
% 2020 July (ZMalik): Arrange mobilityModeCell on a per cluster basis as
% well. Scale diffusion coefficients after merge events.
% 2023 Jan, Fnu Bilal: run code for fixedID receptors.
% 2024 Feb Khuloud Jaqaman: Cleaned up code and made consistent with
% receptorAggregationSimpleOptimizedFixedID.
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


%% Output

errFlag = 0;
receptorInfoAll = [];
timeIterArray = [];

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
    disp('--receptorAggregationVEGFR2LessSimple: Too few input arguments');
    errFlag  = 1;
    return
end

%extract model parameters from modelParam

% Mobility Mode Cell Array
if isfield(modelParam,'mobilityModeCell')
    mobilityModeCell = modelParam.mobilityModeCell;
else
    disp('--receptorAggregationVEGFR2LessSimple: Please supply mobility mode Cell Array');
    errFlag = 1;
end

%receptor density
if isfield(modelParam,'receptorDensity')
    receptorDensity = modelParam.receptorDensity;
else
    disp('--receptorAggregationVEGFR2LessSimple: Please supply receptor density');
    errFlag = 1;
end

%aggregation probability
if isfield(modelParam,'aggregationProb')
    aggregationProb = modelParam.aggregationProb;
else
    disp('--receptorAggregationVEGFR2LessSimple: Please supply aggregation probability');
    errFlag = 1;
end

%aggregation distance
if isfield(modelParam,'aggregationDist')
    aggregationDist = modelParam.aggregationDist;
else
    disp('--receptorAggregationVEGFR2LessSimple: Please supply aggregation distance');
    errFlag = 1;
end

%dissociation rate
if isfield(modelParam,'dissociationRate')
    dissociationRate = modelParam.dissociationRate;
else
    disp('--receptorAggregationVEGFR2LessSimple: Please supply dissociation rate');
    errFlag = 1;
end

%Probability of picking a mobility mode after merging
if isfield(modelParam,'mobModeProbAfterMerge')
    mobModeProbAfterMerge = modelParam.mobModeProbAfterMerge;
else
    disp('--receptorAggregationVEGFR2LessSimple: Please supply the probability of picking a mobility mode after merging');
    errFlag = 1;
end

%Probability of picking a mobility mode after splitting
if isfield(modelParam,'mobModeProbAfterSplit')
    mobModeProbAfterSplit = modelParam.mobModeProbAfterSplit;
else
    disp('--receptorAggregationVEGFR2LessSimple: Please supply the probabilities of picking a mobility mode after splitting');
    errFlag = 1;
end

%Diffusion Coefficient scale after merging
if isfield(modelParam,'diffCoefScaleAfterMerge')
    diffCoefScaleAfterMerge = modelParam.diffCoefScaleAfterMerge;
else
    diffCoefScaleAfterMerge = ones(length(dissociationRate(:,1)),length(dissociationRate(1,:)),...
        length(dissociationRate(1,:)),length(dissociationRate(1,:)));
end

%labeling ratio
if isfield(modelParam,'labelRatio')
    labelRatio = modelParam.labelRatio;
else
    disp('--receptorAggregationVEGFR2LessSimple: Please supply labeling ratio');
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
    disp('--receptorAggregationVEGFR2LessSimple: Please supply intensity quantum');
    errFlag = 1;
end

%exit if there are missing model parameters
if errFlag == 1
    disp('--receptorAggregationVEGFR2LessSimple: Please supply missing variables!');
    return
end

%check model parameters ...

%some must be positive
%09/05/14 - modified to accomodate a vector labelRatio
if any([receptorDensity aggregationDist (labelRatio(:)') intensityQuantum(1)] <= 0)
    disp('--receptorAggregationVEGFR2LessSimple: Receptor density, aggregation distance, labeling ratio and intensity quantum should be positive');
    errFlag = 1;
    return
end

%and some must be non-negative
%2020/Feb/06 - modified for 3D array aggregationProb and matrix
%dissociationRate
%03/25/14 - modified to accomodate a vector aggregationProb
if any(intensityQuantum(2) < 0)
    disp('--receptorAggregationVEGFR2LessSimple: aggregation probability should be non-negative');
    errFlag = 1;
    return
end

if any(dissociationRate< 0, 'all' )  % Modification, Bilal, 2023/03/13, the condition must be inside
    disp('--receptorAggregationVEGFR2LessSimple: dissociation rates may not be negative')
    errFlag = 1;
    return
end

if any(aggregationProb< 0, 'all')  % Modification, Bilal, 2023/03/13, the condition must be inside
    disp('--receptorAggregationVEGFR2LessSimple: aggregation probabilities may not be negative')
    errFlag = 1;
    return
end

if any([modelParam.mobModeProbAfterMerge,modelParam.mobModeProbAfterSplit] < 0,'all') % Modification, Bilal, 2023/03/13, the condition must be inside
    disp('--receptorAggregationVEGFR2LessSimple: mobility mode probabilities may not be negative')
    errFlag = 1;
    return
end

%extract simulation parameters from simParam

%if simParam wasn't supplied at all
if nargin < 2 || isempty(simParam)
    
    probDim = 2;
    observeSideLen = ones(1,probDim);
    timeStep = 0.01 / max(diffCoefDist,dissociationRate);
    simTime = 100 * timeStep;
    initTime = 100 * timeStep;
    randNumGenSeeds = [100 100];
    compTrackFlag =1; 
    
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
        %Modified: 2020/Jan/29, timeStep determined from
        %dissociationRate alone, since diffCoefDist is not an array
        %anymore. Also, now dissociationRate is a matrix.
        timeStep = 0.01 / max(dissociationRate,[],'all');
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
    
    %compound track flag
    if isfield(simParam,'compTrackFlag')
        compTrackFlag = simParam.compTrackFlag;
    else
        compTrackFlag = 1;
    end
    
end

%determine number of iterations to perform
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

%% receptor initial positions, clustering, and mode determination

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


%%%%%%%%%%%%%%%ZMALIK 2020/Jan/28%%%%%%%%%%%%%%%
% Let's determine the mobility modes of our receptors and the
% mobility merge/split partitions

% Compute number of modes
numModes = length(mobilityModeCell(:,1));

% For future use in the aggregation algorithm
mobProbMergePartition = cumsum(mobModeProbAfterMerge,2);

%%% EDIT: 2020/Feb/27 Now we input the probabilities of picking mobility
%%% modes after a split. Here we will produce the partition.
%%% ZMalik
% For future use in the dissociation algorithm
for oldMobMode = 1:numModes
    mobModeProbAfterSplitFromOldMobMode = mobModeProbAfterSplit(:,...
        oldMobMode,:,:);
    mobModeProbAfterSplitFromOldMobModeSqueeze = squeeze(...
        mobModeProbAfterSplitFromOldMobMode);
    mobModeProbAfterSplitFromOldMobModeHorzCat = ...
        mobModeProbAfterSplitFromOldMobModeSqueeze; %Initialize
    for iConc = 2:numModes
        mobModeProbAfterSplitFromOldMobModeHorzCat = horzcat(...
            mobModeProbAfterSplitFromOldMobModeHorzCat(:,:,1),...
            mobModeProbAfterSplitFromOldMobModeSqueeze(:,:,iConc));
    end
    mobModeProbAfterSplitFromOldMobModePartition = cumsum(...
        mobModeProbAfterSplitFromOldMobModeHorzCat,2);
    pairPerm = 1; %Initialize
    for receptNewMob = 1:numModes
        for clustNewMob = 1:numModes
            mobProbSplitPartition(:,oldMobMode,...
                clustNewMob,receptNewMob) = ...
                mobModeProbAfterSplitFromOldMobModePartition(:,pairPerm);
            pairPerm = pairPerm + 1;
        end
    end
end


% Partition the line segment from 0 to 1 in a cumsum.
% That is, modePartition(i) - modePartition(i-1) is
% the probability of being in mode i.
modeFrac = cell2mat(mobilityModeCell(:,1));
modePartition = cumsum(modeFrac);

% Now we assign each particle an initial mobility mode and diffusion
% coefficient
initMobMode = rand(numReceptors,1);
initDiffCoef = zeros(numReceptors,1);
for i = 2 : numModes
    receptIndOfModei = initMobMode <= modePartition(i) & initMobMode > ...
        modePartition(i-1);
    initMobMode(receptIndOfModei) = i;
    initDiffCoef(receptIndOfModei) = random(mobilityModeCell{i,2},...
        nnz(receptIndOfModei),1);
end
% Treat the first mobility mode as a special case
receptIndOfModei = initMobMode <= modePartition(1);
initMobMode(receptIndOfModei) = 1;
initDiffCoef(receptIndOfModei) = random(mobilityModeCell{1,2},...
    nnz(receptIndOfModei),1);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute stepStd from diffCoef.
% Each cluster is assigned a stepStd value.
%initStepStd = sqrt( 2 * initDiffCoef * timeStep );

%Modification LRO 2019/04/30 - since it has an individual value for
%each receptor, for the aggregation distance uses the average value
%averageinitStepStd=mean(initStepStd);

%adjust aggregationDist to account for the finite simulation time step and
%the expected receptor displacement in that time step
%aggregationinitDistCorr = max(aggregationDist,sqrt(probDim)*averageinitStepStd*2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain the initial column vectors relating to diffusion for use.
diffCoef = initDiffCoef;
mobMode = initMobMode;

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

%Reserve memory for output cell array. The first row gives a column
%vector of diffusion coefficients at each time interval. The second row
%gives the mobility mode assignment at that time.
diffHist = cell(2,numIterations);

%store initial information
receptorTraj(:,:,1) = initPositions;
recept2clustAssign(:,1) = receptor2cluster;
clust2receptAssign(1:numClusters,1:maxClustSize,1) = cluster2receptor;
diffHist{1,1} = diffCoef;
diffHist{2,1} = mobMode;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Diagnostic quantities (ryirdaw)
assocStats = struct('numSureAssoc',NaN(numIterations,1),...
    'numPotColl',NaN(numIterations,1),'numColl',NaN(numIterations,1),...
    'numPotColl_Assoc',NaN(numIterations,1),...
    'sureAssocCountBySize',NaN(numReceptors,numIterations),...
    'numCollProbPairs',NaN(numIterations,1),...
    'numPathAssoc',NaN(numIterations,1),'numPathCirAssoc',NaN(numIterations,1),'numberPossibleAssociations',NaN(numIterations,1),'numberPossibleAssociationsPath',NaN(numIterations,1));


% collProbStatStruct = struct('collisionProb',NaN,'pwDist',NaN,...
%     'primaryNodeRadius',NaN,'partnerNodeRadii',NaN);
% collProbStats(numIterations,1) = collProbStatStruct;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

progressText(0,'Simulation');

%iterate in time
for iIter = 2 : numIterations
    %fprintf('\niIter = %d\n',iIter);
    
    %% Dissociation
    % 2020/Feb/05 ZMalik, made a try/catch block for dissociation in
    % the same fashion as the association try/catch block
    %        try
    %allow receptors in clusters to dissociate in current time point
    [cluster2receptor,receptor2cluster,newClusterSize,diffCoef,mobMode] = ...
        receptorDissociationAlg(cluster2receptor,receptor2cluster,...
        clusterSize,dissociationProb,diffCoef,mobMode,mobilityModeCell,...
        mobProbSplitPartition);
    
    %        catch newDissFunExcep
    %             fprintf('\nError at iIter = %d\n',iIter);
    %             disp(newDissFunExcep.message);
    %             pause;
    %             return;
    %        end
    %% Association flags based on outcome of dissociation
    
    %initialize
    %aggregationProbVec = ones(numReceptors,1);
    
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
    %prevReceptor2cluster=receptor2cluster;
    
    prevClusterSize = clusterSize;
    
    % We cannot allow association if a cluster is empty or if we have a
    % cluster/receptor pair that was involved in dissociation
    aggregationProbVec = prevClusterSize == newClusterSize;
    aggregationProbVec(newClusterSize == 0) = 0;
    
    % Reassign newClusterSize back to clusterSize
    clusterSize = newClusterSize;
    %       if (max(newReceptor2Cluster) > max(receptor2cluster))
    %A dissociation has occured. To confirm and identify the
    %receptors involved, determine the cluster size of each receptor
    %now and compare with size values from previous iterations. The
    %best way to do this is to tally the number of other receptors each
    %receptor is associated with.
    %Store the new and previous sizes on two columns for each receptor.
    
    % Dev 2018-Nov-29: Changed cluster-receptor mappings to use
    % accumarray.
    
    %           accumulatedNewClusters=accumarray(newReceptor2Cluster,1,[length(newReceptor2Cluster),1]);
    %           newClusterSizes=accumulatedNewClusters(newReceptor2Cluster);
    
    %           accumulatedPrevClusters=accumarray(prevReceptor2cluster,1,[length(prevReceptor2cluster),1]);
    %           prevClusterSizes=accumulatedPrevClusters(prevReceptor2cluster);
    
    
    
    %           sizeNewPrev_=[newClusterSizes prevClusterSizes];
    %           sizeNewPrev=sizeNewPrev_;
    
    %For those receptors who have dissociated set the
    %aggregationProbVec to 0.  NOTE: if the other receptors remain
    %clustered, their aggregationProbVec must stay as 1.
    %           aggregationProbVec( (sizeNewPrev(:,1) - sizeNewPrev(:,2) < 0) ) = 0;
    
    %Reassign receptor2cluster to the new set reflecting the
    %dissociation.  NOTE: the position vector also shows the
    %dissociation that has occured.
    %           receptor2cluster = newReceptor2Cluster;
    
    %       end
    
    %% New receptor/cluster positions
    %%%%%%%%%%%%%%%%%%Modification ZMalik 2020/Feb/03%%%%%%%%%%%%%%%%%%%%%%
    % Now stpStd corresponds to each cluster. Then we will assign the receptors
    % their proper displacements.
    %get receptor positions at previous time point
    positionsOld = receptorTraj(:,:,iIter-1);
    
    %Initialize receptor displacement vector:
    receptorDisp = zeros(numReceptors,probDim);
    
    %Compute stpStd
    stpStd = sqrt( 2 * diffCoef * timeStep );
    
    %Generate cluster displacement
    
    for iCluster = 1:numClusters
        % Get receptors belonging to this cluster
        clusterMembers = cluster2receptor(iCluster,1:clusterSize(iCluster));
        stepStdCluster = stpStd(iCluster); % Get step standard for the cluster
        clustDisp = stepStdCluster.*randn(1,probDim);
        % Assign all receptors in this cluster the same displacement
        receptorDisp(clusterMembers,:) = repmat(clustDisp,...
            clusterSize(iCluster),1);
    end
    
    
    effectiveMolRadius = max(aggregationDist/2,sqrt(2)*stpStd);
    % ABOVE IS WHAT WE WANT
    % aggregationDistCorr = max(aggregationDist,sqrt(2)*stpStd);
    %get indices of clusters with more than one receptor
    %clustersBig = find(clusterSize>1);
    %generate receptor displacements
    %receptorDisp = stepStd.*randn(numReceptors,probDim);
    
    %assign receptors in a cluster the displacement of the receptor with
    %the smallest index
    %for iCluster = clustersBig'
    
    %get receptors belonging to this cluster
    %    clusterMembers = cluster2receptor(iCluster,1:clusterSize(iCluster));
    
    %%%%%%%%%%%%Modification LRO-2019/04/30%%%%%%%%%%%%%%%%%%%
    % include here also the collective value of diffusion coeficient
    % and update the std for those receptors
    
    %stepStdCluster=mean(stepStd(clusterMembers));
    
    %clustDisp=stepStdCluster.*randn(1,probDim);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %assign all receptors in this cluster the displacement of that
    %receptor
    %receptorDisp(clusterMembers,:) = repmat(clustDisp,...
    %    clusterSize(iCluster),1);
    
    %end
    
    %calculate the new receptor positions at the current time point
    positionsNew = positionsOld + receptorDisp;
    
    %make sure that receptors stay inside the region of interest
    correctionBoundaryLow = min(positionsNew,0);
    positionsNew = positionsNew - 2 * correctionBoundaryLow;
    correctionBoundaryUp = max(positionsNew - repmat(observeSideLen,numReceptors,1),0);
    positionsNew = positionsNew - 2 * correctionBoundaryUp;
    
    %% Association
    
    %        try
    
    numClustPre = length(cluster2receptor(:,1));
    
    [cluster2receptor,receptor2cluster,clusterSize,receptPositions,...
        aggregationProbVec,sureAssocCount, numberPossibleAssociations,...
        mobMode,diffCoef,effectiveMolRadius] = receptorAggregationAlg_VEGFR2ModeDistCorrFixedID(...
        positionsNew,effectiveMolRadius,aggregationDist,aggregationProbVec,aggregationProb,...
        receptor2cluster,cluster2receptor,clusterSize,mobMode,...
        diffCoef,mobilityModeCell,mobProbMergePartition,diffCoefScaleAfterMerge);
    
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
    
    %        catch newAssocFunExcep
    %            fprintf('\nError at iIter = %d\n',iIter);
    %            disp(newAssocFunExcep.message);
    %            pause;
    %            return;
    %        end
    
    [numClusters,maxClustSize] = size(cluster2receptor);
    
    %store new receptor information
    receptorTraj(:,:,iIter) = receptPositions; %200128 ZM Fixed bug based on CA fix from receptorAggregationSimple
    recept2clustAssign(:,iIter) = receptor2cluster;
    clust2receptAssign(1:numClusters,1:maxClustSize,iIter) = cluster2receptor;
    
    %Store new diffusion and mobility mode information
    diffHist{1,iIter} = diffCoef;
    diffHist{2,iIter} = mobMode;
    
    progressText((iIter-1)/(numIterations-1),'Simulation');
    
end %(for iIter = 2 : numIterations)

%% Post-processing

%remove the initialization period from simulation
receptorTraj = receptorTraj(:,:,numIterInit+1:end);
recept2clustAssign = recept2clustAssign(:,numIterInit+1:end);
clust2receptAssign = clust2receptAssign(:,:,numIterInit+1:end);
diffHist = diffHist(:,numIterInit+1:end);

%Now we keep empty rows and columns from clust2receptAssign
%remove empty rows and columns from clust2receptAssign
%cluster2receptor = max(clust2receptAssign,[],3);
%columnSum = sum(cluster2receptor);
%clust2receptAssign = clust2receptAssign(:,columnSum~=0,:);
%rowSum = sum(cluster2receptor,2);
%clust2receptAssign = clust2receptAssign(rowSum~=0,:,:);

% %put receptor trajectories and clusters into the format of the output of
% %trackCloseGapsKalman
% compTracks = convReceptClust2CompTracks(clust2receptAssign,...
%     recept2clustAssign,receptorTraj);

%put information in receptorInfoAll
receptorInfoAll = struct('diffHist',{diffHist},...
    'receptorTraj',receptorTraj,'recept2clustAssign',...
    recept2clustAssign,'clust2receptAssign',clust2receptAssign,...
    'simParam',simParam,'modelParam',modelParam);

%% Receptor labeling and sub-sampling

%KJ (150528): call function to label and sub-sample
%20200813 ZM: Option of which type of receptorInfoLabeled we want
%     [receptorInfoLabeledWithCT,receptorInfoLabeledNoCT] ...
%         = genReceptorInfoLabeled(receptorInfoAll,labelRatio,...
%         intensityQuantum,fracNumFluor,compTrackInd);


% Modification, Bilal, 20230320, call the fixed ID function

receptorInfoLabeled = genReceptorInfoLabeledFixedID(receptorInfoAll,...
    labelRatio,intensityQuantum,compTrackFlag);

end


%% ~~~ the end ~~~


%% subfunction

function [cluster2receptor,receptor2cluster,clusterSize,diffCoef,mobMode] = ...
    receptorDissociationAlg(cluster2receptorPrev,receptor2clusterPrev,...
    clusterSizePrev,dissociationProb,diffCoefPrev,mobModePrev,...
    mobilityModeCell,mobProbSplitPartition)
%RECEPTORDISSOCIATIONALG Dissociation Algorithm for molecules
%
% SYNOPSIS [cluster2receptor,receptor2cluster,clusterSize,diffCoef,mobMode] = ...
%    receptorDissociationAlg(cluster2receptorPrev,receptor2clusterPrev,...
%    clusterSizePrev,dissociationProb,diffCoefPrev,mobModePrev,...
%    mobilityModeCell,mobProbSplitPartition)
%
%INPUTS
%       cluster2receptorPrev: 2D array of clusters along the rows and
%                             member receptors along the columns.
%       receptor2clusterPrev: 1D array of current cluster labels for each
%                             receptor (rows = receptors).
%       clusterSizePrev     : 1D array of cluster sizes.
%       dissociationProb    : (maximum cluster size)x(number of mobility
%                             modes) array of dissociation probabilities.
%       diffCoefPrev        : 1D array of diffusion coefficients.
%       mobModePrev         : 1D array of mode assignment per cluster
%       mobilityModeCell    : (number of modes)x(2) cell array, each row
%                             corresponds to a unique mobility mode. The
%                             second column gives the fraction of particles
%                             living in that mode.
%       mobProbSplitPartition: Cumulative sum of the array
%                              mobModeProbAfterSplit along the second
%                              dimension.
%OUTPUTS
%       cluster2receptor    : 2D array of clusters along the rows and
%                             member receptors along the columns.
%       receptor2cluster    : 1D array of cluster labels for each receptor.
%       clusterSize         : 1D array of cluster sizes.
%       diffCoef            : 1D array of diffusion coefficients.
%       mobMode             : 1D array of mode assignment per cluster.


%copy some input variables for modification and output
cluster2receptorTmp = cluster2receptorPrev; %% Now consistent in size
receptor2clusterTmp = receptor2clusterPrev;
clusterSizeTmp = clusterSizePrev;
diffCoefTmp = diffCoefPrev;
mobModeTmp = mobModePrev;

%find clusters with more than one receptor
%clustersBig = find(clusterSizePrev > 1); (Not needed anymore, specify
%in input that monomers do not dissociate)

%06/27/13 (1 of 2)
%{
dissociationProb used on clusters instead of receptors (above).
Probability can be determined in two ways:  a single value for all current
clusters or each cluster has its own value. For a dissociating cluster, a
receptor will be randomly picked and removed from the cluster (below).
THIS IS NO LONGER THE CASE. NOW RECEPTOR WITH HIGHEST ID WILL BE REOMOVED
(SEE BELOW).
Dissociation happens one receptor at a time.
%}

%%%%%%%%%%%%%%%%%Zachariah Malik 2020/Jan/22%%%%%%%%%%%%%%%%%%%%%%
% Mod: Now dissociation probability is dependent on oligomer size and
% mobility mode

%Maximum cluster size
maxClusterSize = max(clusterSizePrev);

%Extract number of clusters from input
numClusters = length(clusterSizePrev);

%Number of mobility modes
numModes = length(mobilityModeCell(:,1));

%Initialize flag
clusterDissociateFlag = rand(numClusters, 1);

% Empty clusters cannot dissociate
% Also clusters of size 1 cannot dissociate (KJ 240605)
clusterDissociateFlag(clusterSizePrev <= 1) = 0;

% Flag which clusters (of size > 1; KJ 240605) will dissociate
for i = 1:numModes %Mobility mode i
    for j = 2:maxClusterSize %Cluster of size j %KJ 240605: Fixed bug here. Previously, this loop started at j=1, which does not make sense, as that allows clusters of size 1 to dissociate
        %clusterIndOfSizej = clusterSizePrev == j; %Obtain indices of clusters with size j
        %clusterIndOfModei = mobModePrev == i; %Indices of clusters in mode i, at first time step
        %clusterIndOfSizejAndModei = clusterIndOfSizej.*clusterIndOfModei; %Clusters in both
        %clusterIndOfSizejAndModei = logical(clusterIndOfSizejAndModei);
        %Combined the above four expressions into one compound logical
        %expression.
        clusterIndOfSizejAndModei = clusterSizePrev == j & mobModePrev == i;
        clusterDissociateFlag(clusterIndOfSizejAndModei) = ...
            clusterDissociateFlag(clusterIndOfSizejAndModei) < ...
            dissociationProb(j,i);%Whether given cluster dissociates
    end
end

clusterDissociateFlag = logical(clusterDissociateFlag);

% Now for the diffusion coefficients
% First we will update the diffusion coefficient of every cluster that
% dissociated.

%%%%ZMALIK: Updated for loop 2020/March/31
% Assign each of the nonzero numbers in diffCoefUpdateFlag their
% corresponding mode (and keep track of their mode)
for iClustUpdate = 1:numClusters
    if clusterDissociateFlag(iClustUpdate) == 1
        
        iUpdatePairMobMode = rand;
        iUpdatePairPrevMob = mobModePrev(iClustUpdate);
        iUpdatePairOldSize = clusterSizePrev(iClustUpdate);
        for iNewReceptNewMob = 1:numModes  %Loop through possible receptor modes
            for iClustUpdateNewMob = 1:numModes %Loop through the possible cluster modes
                if iUpdatePairMobMode < mobProbSplitPartition(...
                        iUpdatePairOldSize,iUpdatePairPrevMob,...
                        iClustUpdateNewMob,iNewReceptNewMob)
                    break
                end
            end
            if iUpdatePairMobMode < mobProbSplitPartition(...
                    iUpdatePairOldSize,iUpdatePairPrevMob,...
                    iClustUpdateNewMob,iNewReceptNewMob)
                break
            end
        end
        iClustUpdateMobMode = iClustUpdateNewMob;
        mobModeTmp(iClustUpdate) = iClustUpdateMobMode;
        %             diffCoefTmp(iClustUpdate) = random(mobilityModeCell{...
        %                 iClustUpdateMobMode,iUpdatePairOldSize});
        diffCoefTmp(iClustUpdate) = random(mobilityModeCell{...
            iClustUpdateMobMode,2},1); % Modification, Bilal, 20230321, mobilityModeCell should give only one value
        %Add one to the cluster size to get the right column
        
        % Randomly pick which receptor in this cluster leaves (not the
        % first receptor ie, the receptor with the same label as the
        % cluster)
        clusterMembers = cluster2receptorPrev(iClustUpdate,...
            2:clusterSizePrev(iClustUpdate));
        %         recept2Dissociate = clusterMembers(randi(numel(clusterMembers),1)); %KJ 240606: commented out original procedure, which allowed any receptor (other than the one with smallest ID) to dissociate. This may however cause problem when subsampling, because the dissociating receptor might be the receptor with smallest ID among those labeled. 
        recept2Dissociate = clusterMembers(end); %KJ 240606: to avoid problems when labeling a fraction and subsampling, the receptor with the largest ID is always chosen to leave the cluster. This extends the strategy of not allowing the receptor with the smallest ID to leave the cluster.
        
        clusterSizeTmp(recept2Dissociate) = 1;
        cluster2receptorTmp(recept2Dissociate) = recept2Dissociate;
        receptor2clusterTmp(recept2Dissociate) = recept2Dissociate;
        
        clusterSizeTmp(iClustUpdate) = clusterSizePrev(iClustUpdate)-1;
        cluster2receptorTmp(iClustUpdate,1:length(clusterMembers))...
            = [iClustUpdate, clusterMembers(clusterMembers ~= recept2Dissociate)];
        cluster2receptorTmp(iClustUpdate,...
            length(clusterMembers)+1:maxClusterSize) = 0;
        
        iNewReceptMobMode = iNewReceptNewMob;
        mobModeTmp(recept2Dissociate) = iNewReceptMobMode;
        diffCoefTmp(recept2Dissociate) = random(mobilityModeCell{...
            iNewReceptMobMode,2});
        
    end
end

%copy temporary variables into output variables
cluster2receptor = cluster2receptorTmp;
receptor2cluster = receptor2clusterTmp;
clusterSize = clusterSizeTmp;
diffCoef = diffCoefTmp;
mobMode = mobModeTmp;

end