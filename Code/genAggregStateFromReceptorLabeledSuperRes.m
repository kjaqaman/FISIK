function detectionAggregState = genAggregStateFromReceptorLabeledSuperRes(receptorInfoAll,observedFraction)
%GENAGGREGSTATEFROMRECEPTORLABELEDSUPERRES calculates oligomeric state distribution from a snapshot time point
%
%   SYNPOSIS: detectionAggregStateLabelReceptors = genAggregStateFromReceptorLabeledSuperRes(receptorInfoAll,observedFraction)
%
%   INPUT:
%
%       receptorInfoAll  :  Output of receptorAggregationSimpleOptimized
%
%       observedFraction :  Fraction of observed receptors
%
%
%   OUTPUT:
%
%       detectionAggregState:  Vector with dimention 1 x maximum oligomeric
%                              state, storing the number of particles with
%                              oligomeric state 1, 2, 3, etc.
%
%  Luciana de Oliveira, May 2017.
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

% load the receptorInfoAll in the last frame
receptorPerCluster = receptorInfoAll.clust2receptAssign(:,:,end);

% calculate the number of receptors per cluster size
maxClusterSize=size(receptorPerCluster,2);



%% labeling the receptors

% reserve space for the detectionAggregState
detectionAggregState=zeros(1,maxClusterSize);

%label receptors

labelFlag = rand(size(receptorPerCluster)) <= observedFraction;
receptorPerClusterLabeled = labelFlag .* receptorPerCluster;

% remove the rows with only zeros
rowSum = sum(receptorPerClusterLabeled,2);
receptorPerClusterLabeled = receptorPerClusterLabeled(rowSum~=0,:);


%% calculate aggregate state


% convert receptorPerClusterLabeled in a logical matrix
receptorPerClusterLabeled= logical(receptorPerClusterLabeled);

% calculate the density of receptors in each cluster size
densityReceptors=sum(receptorPerClusterLabeled,2);

% calculating the number of receptors in each cluster size
for i=1:maxClusterSize
    detectionAggregState(i) = sum(densityReceptors==i);
end




