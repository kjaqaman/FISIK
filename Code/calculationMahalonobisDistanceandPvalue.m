function [s,p]=calculationMahalonobisDistanceandPvalue(paramMatrixTarget,paramMatrixProbe,flagUseAsIs)

% This fuction calculates the Mahalanobis distance value using the target
%and probe intermediate statistics given.
%
%
%  INPUT:   
%     paramMatrixTarget    : target intermediate statistics array
% 
%     paramMatrixProbe     : probe intermediate statistics array
%
%     flagUseAsIs          : 1 to use parameter matrices as input, 0 to
%                            reformat them first. 
%                            Optional. Default: 1.
%
%
%   OUTPUT:
%       
%      s                    : calculated Mahalanobis distance
%
%      p                    : p value
%
%
% Luciana de Oliveira, July 2016
% Khuloud Jaqaman, February 2025: added flagUseAsIs option.
%
% NOTE: Currently function works only for flagUseAsIs = 1.
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reform of the parameters and calculation of theta and v

flagUseAsIs = 1;

if flagUseAsIs
    
    thetaTarget = nanmean(paramMatrixTarget,2);
    thetaProbe = nanmean(paramMatrixProbe,2);
    vTarget = nancov(paramMatrixTarget');
    vProbe = nancov(paramMatrixProbe');

else

% [thetaTarget,thetaProbe,vTarget,vProbe] = ...
%     calcThetaVaramCovMatrixForMahalonobisDistance(paramMatrixTarget,paramMatrixProbe);

end

%% calculation of Mahalanobis distance

% s: Mahalanobis distance

s = transpose(thetaProbe-thetaTarget)*(pinv(vProbe+vTarget))*(thetaProbe-thetaTarget);

%p Value

% the degrees of freedom are calculated as the number of intermediate statistis
dof= size(thetaTarget,1);
p = 1 - chi2cdf(s,dof);

end