function [cluster2receptor,receptor2cluster,clusterSize,receptPositions,...
    aggregationProbVec,sureAssocCount,assocFlag] = receptorAggregationAlg_maxWeightedMatching_sureCollOptimized(...
    receptPositionsPrev,aggregationDist,aggregationProbVec,aggregationProb,...
    receptor2clusterPrev,cluster2receptorPrev,clusterSizePrev)
%RECEPTAGGREGATIONALG_MAXWEIGHTEDMATCHING_SURECOLL performs receptor-to-receptor
%   and receptor-to-cluster associations for pairs whose pairwise distance
%   falls below the default value (aggregationDist). These are considered
%   sure collisions. This function replaces receptorAggregationAlg in the 
%   simulation.  maxWeightedMatching is used to perform associations, 
%   limiting all receptors/clusters to associate with only one other 
%   receptor. Cluster-to-cluster associations are not allowed.
%
%   INPUT:
%       receptPositonsPrev:     2D array of current receptor positions
%       aggregationDist:        the distance threshold for association
%       aggregationProbVec:     1D array of receptor aggregation flags
%       aggregationProb:        1D array giving domain from which 
%                               aggregation probabilities will be drawn,
%                               for each cluster size (formation). Element
%                               #1 is assumed to be a default value, used
%                               when aggregationProb is a scalar or when
%                               its size < a given cluster size.
%       receptor2clusterPrev:   1D array of current cluster labels for each
%                               receptor (rows = receptors)
%       cluster2receptorPrev:   2D array of clusters along the rows and
%                               member receptors along the columns
%       clusterSizePrev:        1D array of cluster sizes 
%
%   OUTPUT:
%       cluster2receptor:   2D array of clusters along the rows and
%                           member receptors along the columns
%       receptor2cluster:   1D array of cluster labels for each receptor
%       clusterSize:        1D array of cluster sizes
%       receptPositions:    2D array of recpetor positions 
%       aggregationProbVec: the updated array of receptor association flags
%       sureAssocCount:     number of associations that occured
%       LRO included numberPossibleAssociations 
%       numberPossibleAssociations: count for the number of possible associations
%       in each time step.
%
%   Robel Yirdaw, 10/02/13
%       Modified, 01/28/14
%       Modified, 03/25/14
%       Modified, 07/03/14 (renamed sureColl)
%       Modified, 08/11/14 (sure assoc counts returned)
%       Modified, 2017/11/30 Luciana de Oliveira
%       Modified, 2018 Nov 29 by Devin O'Kelly
%           -Changed calculation of particle interactions to use
%           rangesearch instead of squareform/pdist.
%           -Changed some variable names to reflect their functionality
%           more accurately.
%           -Modified memory deallocation to spend less time thrashing.
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

% assocFlag=0; % try putting it here.


    %Total number of receptors % TOTAL NUMBER OF POSSIBLE CLUSTER SIZES.
    numReceptors = length(aggregationProb)-1;   %%%%%%%%%%%%%%%%%%% temporally need this -1
    
    %Determine association flag for all receptors
    %01/28/14: moved below for nodes instead of receptors
    %aggregFlag = rand(numReceptors,1) < aggregationProb;
    
    %The new vectors initialized to previous
    receptor2cluster = receptor2clusterPrev; % Index list associating a receptor to a cluster.
    cluster2receptor = cluster2receptorPrev; % Inverse of above.
    receptPositions = receptPositionsPrev;   % Last location of receptors.
    clusterSize = clusterSizePrev;           % Number of items in each cluster (redundant?)

    %081114
    sureAssocCount = NaN;  %???
    
    %Current number of nodes
    numNodesInit = sum(cluster2receptorPrev(:,1) ~= 0); %Probably not super efficient, but probably doesn't matter.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %01/28/14
    %Determine association flag for all nodes. A row in cluster2receptor is
    %a node. Using the aggregationProb value assigned to the first member
    %of the cluster (row) when generating the flag for the node. This is
    %necessary because a cluster or receptor may have aggregatinProb value
    %of 0 indicating that it shouldn't undergo association because it has
    %undergone dissociation in the current iteration. Also, setting
    %aggregationProb to 1 for all members of a cluster in the main
    %simulation has been disabled. Thus, the possible values found in
    %aggregationProb are 0 (for a receptor or a set of receptors forming a
    %cluster) or the input value; previously it also included the
    %value of 1 for all receptors of a cluster as well.
    %
    %nodeFlag = rand(numNodesInit,1) <...
    %    aggregationProb(cluster2receptorPrev(1:numNodesInit,1),1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   03/25/14
    %   The function now inputs associationProb and associationProbVec. In
    %   the earlier version, associationProb was actually
    %   associationProbVec. Here, associationProb is a vector and gives the 
    %   probability of forming a cluster of a given size. On the other
    %   hand, associationProbVec is an association flag for each receptor.
    %
    %   Thus, we first set up nodeFlag using associationProbVec. Those with
    %   value of 0, will not undergo association, and this would be due to
    %   having participated in a dissociation event in the current
    %   iteration.  However, this quantity can be used for other purposes.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nodeFlag = aggregationProbVec(cluster2receptorPrev(1:numNodesInit,1)); % 0 if dissociated recently. 
    
    assocFlag=0; %%%%%% testing here, maybe was one case that was out before. LRO 2019/04/02
    
    
    %If there are receptors that can associate, then continue.
    if any(nodeFlag)
        %Current nodes
        nodesInit = cell(numNodesInit,1);
        %Node sizes
        nodeSizes = NaN(numNodesInit,1);        
        %Position vector for current nodes        
        nodePositionsInit = NaN(numNodesInit,2);
        
        %Association flag for nodes
        %01/28/14: moved above
        %nodeFlag = NaN(numNodesInit,1);
        
        %Current largest cluster size
        largestClusterInit = 0;
        
        %Build node information % Dev: Unnecessarily iterative?
        for nodeIter=1:numNodesInit
            %A row in cluster2receptor is a node
            nodeMembers = cluster2receptorPrev(nodeIter,cluster2receptorPrev(nodeIter,:) ~= 0); % SPLIT
            %Assign members
            nodesInit{nodeIter} = nodeMembers; 
            %Determine size of node
            nodeSizes(nodeIter) = length(nodeMembers);            
            %Get positions
            nodePositionsInit(nodeIter,1:2) = receptPositionsPrev(nodeMembers(1,1),1:2);
            
            %Get flag - assuming all members of a cluster have the same
            %flag (1)
            %01/28/14: moved above for nodes instead of receptors
            %nodeFlag(nodeIter,1) = aggregFlag(nodeMembers(1,1),1);
            
            %Save largest cluster size
            if (length(nodeMembers) > largestClusterInit)
                largestClusterInit = length(nodeMembers);
            end
            
            nodeMembers(:)=[]; 
        end %for all initial nodes
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Get id of those that can associate (flag = 1)        
        associatableNodeIndices = find(nodeFlag == 1); % Only monomers can associate.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Get pairwise distances for those receptors with flag 1
%         pwdists=pdist(nodePositionsInit(associatableNodeIndices,:));
%         aggregNodePWDist = squareform(pwdists);
%        
        [closeNeighbors,closeNeighborsDists] = rangesearch( nodePositionsInit(associatableNodeIndices,:), ...
                                                            nodePositionsInit(associatableNodeIndices,:), ...
                                                            aggregationDist);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %In the pairwise distance matrix, find the nodes where distance
        %criteria is satisfied. The nodes are the columns in the matrix and
        %what is stored in aggreNodeIndx_aggregDist is the column index.
        %These are not (necessarily) node ids.
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % modification LRO 2017/11/30- as we need to know how
        % much times two receptors are in the minimum distance
        % for association, I add this step to quantify the number of times
        % two receptors are the min distance
        % aggregNodeIndx_aggregDist is output.
        
        
        numNeighbors=zeros(numel(closeNeighbors),1);
        for k=1:numel(closeNeighbors)
            numNeighbors(k)=numel(closeNeighbors{k})-1;
        end
        
        nodesWithCloseNeighbors = find(numNeighbors); 
        
        % total number of receptors that are in the possible distance for
        % associating
% % % % % %         numberPossibleAssociations=length(nodesWithCloseNeighbors);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
       
        
        
        if  isempty (nodesWithCloseNeighbors)  %% LRO 04022019 flag for those that are not associating 
                        assocFlag=0;
        end
        
        
        if (~isempty(nodesWithCloseNeighbors))
            %Nodes associating
            aggregNodeID = associatableNodeIndices(nodesWithCloseNeighbors);
            %Count of nodes associating
            numNodesAggreg = length(aggregNodeID);
            
            %The number of pairwise associations
            numNodePairsWithinDist = sum(numNeighbors)./2;
            
            %Collect edges and weights. % BAD NAMES
            paramEdges = NaN(numNodePairsWithinDist,2);
            paramWeights = NaN(numNodePairsWithinDist,1);
            edgeCount = 0;
            
            %Iterate through each % INDIRECT
            for aggregNodeIter=1:numNodesAggreg                
                %Get a receptor from list of aggregating nodes. This will
                %be used to access the columns of the pair-wise distance
                %matrix below.
                currNodeID = aggregNodeID(aggregNodeIter);
                
                currentlyAggregatingNodeIndex=nodesWithCloseNeighbors(aggregNodeIter);
                %Get a list of potential partners
                % Find all neighbors within the aggregation distance, don't
                % allow for self-selection. Only move forward in the list
                % (i.e., if nodes 1 and 12 are close together, only test
                % for 1 to interact with 12, and don't try again when it's
                % 12's turn)
                
                [closeNeighborsIndices, indSort] = sort(closeNeighbors{currentlyAggregatingNodeIndex}(2:end).');
                lessFlag=closeNeighborsIndices<currentlyAggregatingNodeIndex;
                closeNeighborsDist= closeNeighborsDists{currentlyAggregatingNodeIndex}(indSort+1);
                closeNeighborsIndices(lessFlag)=[];
                closeNeighborsDist(lessFlag)=[];
                
                
                %The node ids corresponding to the index values above -
                %note this must be taken from aggregNodeID_aggregFlag, not
                %aggregNodeID.      
                partnerNodeID = associatableNodeIndices(closeNeighborsIndices);               
                
                %Cluster associations not allowed
                if ( (nodeSizes(currNodeID) > 1) &&...  % If the current node is a cluster, 
                        any(nodeSizes(partnerNodeID) > 1) ) % and any of the node's partners are a cluster
                    
                    %Get flags of nodes in partnerNodeID with sizes > 1
                    cluster2clusterAssoc = (nodeSizes(partnerNodeID) > 1); 
                    %Remove values for those nodes
                    
                    partnerNodeID(cluster2clusterAssoc) = []; % Don't allow fusion with those other nodes.
                    closeNeighborsIndices(cluster2clusterAssoc) = []; % Ignore from here out. 
                    closeNeighborsDist(cluster2clusterAssoc) = []; % Ignore from here out. 
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %03/25/14
                %Determine the sizes of clusters that can potentially be
                %formed and use these to draw probabilities from the given
                %range, in aggregationProb. If aggregationProb is a scalar,
                %or there are more sizes than values in aggregationProb,
                %the vector will be expanded for the maximum needed size 
                %and the default domain assigned for these sizes.
                %
                if (~isempty(partnerNodeID))
                    %Determine the cluster sizes that can be formed
                    finalSizes = nodeSizes(currNodeID) + nodeSizes(partnerNodeID); % Shouldn't all of the partner nodes be size 1?
                    %If all sizes are not accounted for in aggregationProb,
                    %expand it and assign default value which is element 1.
                    if (length(aggregationProb) < max(finalSizes)) % REDUNDANT. MOVE OUTSIDE.
                        aggregationProb(length(aggregationProb)+1:max(finalSizes),1) =...
                            aggregationProb(1,1);
                    end
                    %Draw a random value for each cluster to be formed,
                    %using its size, and get association flag.
                    isAssociating = (rand(length(partnerNodeID),1) < aggregationProb(finalSizes) );
                    
                    
                    
                     if  isempty ( isAssociating)  %% LRO 04022019 flag for those that are not associating 
                        assocFlag=0;
                      end
                                        
                    %Remove from the list partner nodes not associating
                    
                    partnerNodeID(~isAssociating) = [];
                    closeNeighborsIndices(~isAssociating) = [];
                    closeNeighborsDist(~isAssociating) = [];
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
                    
                %Number of possible partners for the current node
                numPartners = length(partnerNodeID);
                %Determine edges and weights of current node with all
                %potential partners
                if (numPartners > 0)
                    %Determine edges. Note that partnerNodeID can have more
                    %than one value.
                    tempEdges(1:numPartners,1) = currNodeID;
                    tempEdges(1:numPartners,2) =  partnerNodeID;
                    tempEdges((numPartners+1):end,:)=[];
                    
                    %Determine weights
                    tempWeights = 1./closeNeighborsDist;
                    
                    %Accumulate
                    paramEdges(edgeCount+1:(edgeCount+numPartners),1:2) = tempEdges;
                    paramWeights(edgeCount+1:(edgeCount+numPartners)) = tempWeights;
                    
                    %Increment total count of edges (pairs of nodes)
                    edgeCount = edgeCount + numPartners;
                     
                end

            end %for each associating node

            %Next, depending on the number of edges determined, perform
            %selection - if more than one edge found, use
            %maxWeightedMatching to pick edges.  
            selectedEdges = [];
            if (edgeCount == 1)
                %Trim empty rows in edges only
                paramEdges(isnan(paramEdges(:,1)),:) = [];             
                selectedEdges = paramEdges;
            elseif (edgeCount > 1)
                %Trim blank rows in edges and weights
                paramEdges(isnan(paramEdges(:,1)),:) = [];
                paramWeights(isnan(paramWeights(:,1)),:) = [];
                %The node input parameter is the largest node number
                paramNumNodes = max(paramEdges(:));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Peform weighted matching to pick a single pairwise
                %interaction for every receptor                
                selectionFlags = maxWeightedMatching(paramNumNodes,paramEdges,paramWeights);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Per returned flags above, pick the edges
                selectedEdges = paramEdges(selectionFlags,:);
            end
               
            %If edges selected, then associations have occured. Construct
            %the new positions, cluster to receptor and receptor to cluster
            %relationships.
            if (~isempty(selectedEdges))               
                %The new largest cluster size
                largestClusterNew = largestClusterInit;
                %Copy initial nodes to be modified
                nodesNew = nodesInit;
                
                %081114
                %Initialize sure association count array. Rows are new node
                %sizes and the largest possible is the current largest + 1.
                sureAssocCount = zeros(max(nodeSizes) + 1,1);
                sureAssocCount(1) = NaN;
                
                %Build new nodes list
                for newClusterIter=1:length(selectedEdges(:,1))
                    %Get node pair from list
                    tempNode = selectedEdges(newClusterIter,1:2);
                    %Insert new node
                    nodesNew{tempNode(1)} = [nodesNew{tempNode(1)} nodesNew{tempNode(2)}];
                    %The merged node will be set to NaN
                    nodesNew{tempNode(2)}(:) = NaN;
                    
                    %Update positions - receptors in nodes that associated
                    %will get the average position of the two nodes
                    %Average of node positions
                    %newPos = mean(nodePositionsInit(tempNode(1:2),:));
                    %receptPositions(nodesNew{tempNode(1)},1) = newPos(1);
                    %receptPositions(nodesNew{tempNode(1)},2) = newPos(2);
                    
                    %Average of receptor positions
                    %newPos = mean(receptPositionsPrev([nodesInit{tempNode(1:2)}],:));
                    %receptPositions(nodesNew{tempNode(1)},1) = newPos(1);
                    %receptPositions(nodesNew{tempNode(1)},2) = newPos(2);
                    %Or equivalently, but with one less indexing than above
                    newPosX = sum(nodePositionsInit(tempNode(1:2),1).*nodeSizes(tempNode(1:2)))/sum(nodeSizes(tempNode(1:2)));
                    newPosY = sum(nodePositionsInit(tempNode(1:2),2).*nodeSizes(tempNode(1:2)))/sum(nodeSizes(tempNode(1:2)));
                    receptPositions(nodesNew{tempNode(1)},1) = newPosX;
                    receptPositions(nodesNew{tempNode(1)},2) = newPosY;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %070314
                    %Update aggregationProbVec since it is now returned and
                    %will be used by the second association function, for
                    %nodes whose pairwise distance falls between
                    %aggregationDist and aggregationDistCorr, i.e.,
                    %aggregationDist <= pairwise distance <
                    %aggreagtionDistCorr. Setting the flag to 0 here will
                    %prevent those that associated here from associating in
                    %the second function.
                    aggregationProbVec(nodesNew{tempNode(1)}) = 0;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    %The new largest cluster size
                    if (length(nodesNew{tempNode(1)}) > largestClusterNew)
                        largestClusterNew = length(nodesNew{tempNode(1)});
                    end
                    
                    %081114
                    %Accumulate sure association counts. %%%%%%%%%%%%%%%%
                    %probably here I need the flag
                    currNewNodeSize = length(nodesNew{tempNode(1)});  %%%%%%%%%%%%%%%%%%%%%%%%%%% i will be probably be here I need the flag.
                    assocFlag=1;   %%% LRO 04022019 try using this flag for association
                    sureAssocCount(currNewNodeSize) =...
                        sureAssocCount(currNewNodeSize) + 1;
                    
%                     clear tempNode currNewNodeSize
                end %for each new cluster
                
                %Build new receptor2cluster, cluster2receptor and
                %clusterSize variables.  Initialize cluster2receptor with
                %zeros not NaN - unused columns of cluster rows are
                %supposed to have zero values, not NaN.
                %Note: combining the while loop below with the for above
                %is not a good choice specifically because of
                %receptor2cluster. The construction of cluster2receptor can
                %be handled in the for above easily but not
                %receptor2cluster.  A change to a row in c2r means the
                %cluster label for those receptors has to be changed.
                %However, the label for all receptors in clusters
                %subsequent to the changed one will also have to be
                %updated - the cluster labels are in ascending order.
                receptor2cluster = NaN(numReceptors,1);
                cluster2receptor = zeros(numNodesInit,largestClusterNew);
                clusterSize = NaN(numNodesInit,1);
                numNodesNew = 1;
                initNodesIter = 1;
                while (initNodesIter <= numNodesInit)
                    %NaN rows in the new set of nodes will be skipped
                    itNode=nodesNew{initNodesIter};
                    if (~isnan(itNode(:)))
                        clusterSize(numNodesNew) = length(itNode);
                        cluster2receptor(numNodesNew,1:length(itNode)) = ...
                            itNode;
                        receptor2cluster(itNode,1) = numNodesNew;
                        numNodesNew = numNodesNew + 1;
                    end
                    
                    initNodesIter = initNodesIter + 1;                    
                end %while
                
                %Trim empty array elements                
                cluster2receptor(numNodesNew:end,:) = [];
                clusterSize(numNodesNew:end,:) = [];
                
            end %if associations occured
                
        end %if any distance within aggregationDist exists in distance mtx
        
    end %if any node has association flag = 1
        
end %function
