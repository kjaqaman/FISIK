function clusterStatsStatic = clusterDensityStatic(detectionAggregState,infoSpace)
%CLUSTERDENSITYSTATIC calculates cluster/oligomer densities from static snapshot data
%
%   SYNPOSIS: clusterStatsStatic = clusterDensityStatic(detectionAggregState,infoSpace)
%
%   INPUT:   
%
%       detectionAggregState:  Output of function genAggregStateFromReceptorLabeledSuperRes
%
%       infoSpace           :  Structure with fields           
%                  
%                    .probDim         :  Problem dimensionality.
%
%                    .areaSideLen    :  Simulation/image side length values,
%                             which can be a single value or a value per
%                             side. In units of interest (e.g. um).
%
%                   
%   OUTPUT:
% 
% 
%           clusterStatsStatic: a struct with the following fields (outputs as in
%           function clusterNumbersFromCompTracks)
% 
%           a) clusterCount (max cluster in in the last frame): number of 
%              clusters per size.
%           b) clusterFrac (max cluster in the last frame): fraction
%              of clusters per size.
%           c) clusterDensity (max cluster in the last frame): density
%              of clusters per size.
%           d) receptorCount : number of receptors
%              in the last frame.
%           e) receptorDensity : density of receptors
%              in the last frame.
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

%get intensity and space information
probDim = infoSpace.probDim;
areaSideLen = infoSpace.areaSideLen;

clusterCount=detectionAggregState;

% to have the same outputs as in the function clusterNumbersFromCompTracks

%Determine area to be used in the calculation of density per cluster
numSideLenVals = length(areaSideLen);
if (numSideLenVals == 1)
    simArea = areaSideLen ^ probDim;
else
    simArea = prod(areaSideLen);
end

maxClustSize=length(clusterCount);

%calculate cluster fractions 
clusterFrac = clusterCount / sum(clusterCount);

%calculate cluster density
clusterDensity = clusterCount / simArea;

%calculate receptor density per cluster size
receptorCountClust=zeros(length(clusterCount),1);
receptorDensityClust=zeros(length(clusterCount),1);

%number of receptors and clusters with size 1 is the same
receptorCountClust(1)=clusterCount(1);
receptorDensityClust(1)=clusterDensity(1);

for clusterIndex=2:maxClustSize
 receptorCountClust(clusterIndex) = clusterCount(clusterIndex).* clusterIndex;
 receptorDensityClust (clusterIndex)= receptorCountClust (clusterIndex)/ simArea;
end

% calculate the total number of receptors in the last frame

receptorCount=sum(receptorCountClust);
receptorDensity=sum(receptorDensityClust);


%Save values for return
clusterStatsStatic.clusterCount = clusterCount';
clusterStatsStatic.clusterFrac = clusterFrac';
clusterStatsStatic.clusterDensity = clusterDensity';
clusterStatsStatic.receptorCount = receptorCount;
clusterStatsStatic.receptorDensity = receptorDensity;

end


