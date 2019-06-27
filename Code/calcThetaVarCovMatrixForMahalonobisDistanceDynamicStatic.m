function [thetaTarget,thetaProbe,vTarget,vProbe]=calcThetaVarCovMatrixForMahalonobisDistanceDynamicStatic(paramMatrixTarget,paramMatrixProbe, systemState)

% This fuction reform the probe and target paramMatrix to have the same size 
% for the calculation of Mahalanobis distance and pValue and also calculate
% the values of theta and v for target and probe
%
%  INPUT:   
%      paramMatrixTarget    : target intermediate statistics array
% 
%      paramMatrixProbe     : probe intermediate statistics array
%     
%      systemState : flag for dynamic or static data
%                             1 - dynamic
%                             2 - static
% 
%  OUTPUT:
%       
%     thetaTarget: mean value of the paramMatrix target intermediate statistics
%     thetaProbe: mean value of the paramMatrix probe intermediate statistics
%     vTarget: variance covariance target paramMatrix
%     vProbe: variance covariance probe paramMatrix
%
% Luciana de Oliveira, december 2016
% Modification LRO, February 2017: The modifications are done to have a
% genetal calculation for static or dynamic data. At this time I included a
% new imput and divided the function in two blocks.
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


% In this version when target and probe param matrix have different sizes 
%it takes off the extra rows of on/off rate and fill zeros in the density



%% INPUT

%%%%% Modification LRO, February 2017

%Test if it is dynamic or static data

if systemState==1
%% dynamic data

%calculates the maximum cluster size

%%
%calculates the maximum cluster size
maxClusterSizeTarget=size(paramMatrixTarget,1)/3;
maxClusterSizeProbe=size(paramMatrixProbe,1)/3;



%target
paramMatrixTargetOffRate = paramMatrixTarget(maxClusterSizeTarget+1:2*maxClusterSizeTarget,:);
paramMatrixTargetDensity = paramMatrixTarget(2*maxClusterSizeTarget+1:end,:);



%probe
paramMatrixProbeOffRate = paramMatrixProbe(maxClusterSizeProbe+1:2*maxClusterSizeProbe,:);
paramMatrixProbeDensity = paramMatrixProbe(2*maxClusterSizeProbe+1:end,:);



%% removing rowns with only NaNs

%cutOffNaN is the maximum number off NaNs permited in a row to this row enter
%in the calculation of Mahalanobis distance

% calculate how many observations we have
numOfMoviesTarget= size(paramMatrixTargetOffRate,2);
numOfMoviesProbe= size(paramMatrixProbeOffRate,2);
% minimum number of data points
cutOffNaN=5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% modification 2018/12/18
% compare target and probe
if maxClusterSizeTarget~=maxClusterSizeProbe
    
    clusterDiff=abs(maxClusterSizeTarget-maxClusterSizeProbe);
    
    if maxClusterSizeTarget>maxClusterSizeProbe
        
        %% off rate
        
        paramMatrixProbeOffRate = [paramMatrixProbeOffRate;nan(clusterDiff,numOfMoviesProbe)];
        
        
        %% density
        
         paramMatrixProbeDensity = [paramMatrixProbeDensity;zeros(clusterDiff,numOfMoviesProbe)];
       
    else
       %% off rate
        paramMatrixTargetOffRate = [paramMatrixTargetOffRate;nan(clusterDiff,numOfMoviesTarget)];
        
        
            
        %% density
        paramMatrixTargetDensity=[paramMatrixTargetDensity; zeros(clusterDiff,numOfMoviesTarget)];
    
       
    end
end

% put the paramMatrix together

% % paramMatrixTarget=[paramMatrixTargetOnRate;paramMatrixTargetOffRate;paramMatrixTargetDensity];
% % paramMatrixProbe=[paramMatrixProbeOnRate;paramMatrixProbeOffRate;paramMatrixProbeDensity];


paramMatrixTarget=[paramMatrixTargetOffRate;paramMatrixTargetDensity];
paramMatrixProbe=[paramMatrixProbeOffRate;paramMatrixProbeDensity];

%% removing NaNs

%% OnRate

% if there is a NaN and a number or a NaN and a NaN in the probe and
% target, this row should be removed.

%target and probe
sumNotNaNT = sum(~isnan(paramMatrixTarget),2);
sumNotNaNP = sum(~isnan(paramMatrixProbe),2);

iRowNaN = sumNotNaNT>=cutOffNaN & sumNotNaNP>=cutOffNaN;


% update paramMatrix Target and probe
paramMatrixTarget=paramMatrixTarget(iRowNaN,:);
paramMatrixProbe=paramMatrixProbe(iRowNaN,:);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% static data

elseif systemState==0    
    
maxClusterSizeTarget=size(paramMatrixTarget,1);
 maxClusterSizeProbe=size(paramMatrixProbe,1);

 if maxClusterSizeTarget ~= maxClusterSizeProbe
     
      %calculates the difference between the sizes
     
     clusterDiffDensity=abs(maxClusterSizeTarget-maxClusterSizeProbe);
    
    if (maxClusterSizeTarget > maxClusterSizeProbe)
    paramMatrixProbe = [paramMatrixProbe;zeros( clusterDiffDensity,size(paramMatrixProbe,2))];
    elseif (maxClusterSizeTarget < maxClusterSizeProbe)
    paramMatrixTarget=[paramMatrixTarget; zeros(clusterDiffDensity,size(paramMatrixTarget,2))];
    end
 end
end 
    
%% calculate theta and cov matrix

%theta: target and probe intermediate statistics data
thetaTarget=nanmean(paramMatrixTarget,2);
thetaProbe=nanmean(paramMatrixProbe,2);

% v: target and probe variance-covariance matrix
%target
vTarget = nancov(paramMatrixTarget');
vProbe =nancov(paramMatrixProbe');

end