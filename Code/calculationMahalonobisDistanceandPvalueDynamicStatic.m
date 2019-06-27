function [s,p]=calculationMahalonobisDistanceandPvalueDynamicStatic(paramMatrixTarget,paramMatrixProbe,systemState)

% This fuction calculates the Mahalanobis distance value using the target
%and probe intermediate statistics given.
%
%
%  INPUT:
%     paramMatrixTarget: 2D matrix with each column = intermediate
%                        statistics vector of 1 simulation/movie.
%                        Intermediate statistics vector depends on data
%                        type.
%                   For dynamic data, it is the vector of off rates
%                        per cluster/oligomer size (from 1 to maximum size)
%                        followed by densities per cluster/oligomer size
%                        (from 1 to maximum size). This is the
%                        concatenation of the outputs rateOffPerClust and
%                        densityPerClust from clusterOnOffRatesAndDensity,
%                        after filling missing values, as described below.
%                   For static data, it is the vector of densities per
%                        cluster/oligomer size (from 1 to maximum size).
%                        This is the field "clusterDensity" in the output
%                        clusterStatsStatic from clusterDensityStatic,
%                        after filling missing values, as described next.
%                   Filling missing values: If the maximum cluster size is 
%                        not the same between the different
%                        simulations/movies to be concatenated,
%                        NaN should be used to fill in the missing off rate
%                        values and 0 should be used to fill in the missing
%                        density values. Note that both off rate and
%                        density must go up to the same maximum cluster
%                        size. For example, if maximum cluster size is 5, the
%                        intermediate statistics vector is of length 10 in
%                        the case of dynamic data and 5 in the case of
%                        static data.
%     paramMatrixProbe:  Same as paramMatrixTarget, but for probe data.
%     systemState     :  Flag for dynamic or static data
%                             1 - dynamic
%                             2 - static
%
%   OUTPUT:
%
%     s               : Mahalanobis distance.
%     p               : p-value.
%
%
% Luciana de Oliveira, July 2016
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reform of the parameters and calculation of theta and v

[thetaTarget,thetaProbe,vTarget,vProbe]=calcThetaVarCovMatrixForMahalonobisDistanceDynamicStatic(paramMatrixTarget,paramMatrixProbe,systemState);

%% calculation of Mahalanobis distance

% s: Mahalanobis distance

s = transpose(thetaProbe-thetaTarget)*(pinv(vProbe+vTarget))*(thetaProbe-thetaTarget);

%p Value

% the degrees of freedom are calculated as the number of intermediate
% statistis
dof= size(thetaTarget,1);

p = 1 - chi2cdf(s,dof);

end