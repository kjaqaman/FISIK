function [objFunc,residuals] = multiStepFunction(param,x,y)
%MULTISTEPFUNCTION calculates difference between multi-step function and input data
%
%SYNOPSIS [objFunc,residuals] = multiStepFunction(param,x,y)
%
%INPUT  param        : Column vector with multi-step function parameters.
%                      For n steps, first n entries with step locations,
%                      then n+1 entries with y-value between steps.
%       x            : Independent variable of time series.
%       y            : Dependent variable of time series.
%
%OUTPUT objFunc      : Objective function for minimization.
%       residuals    : Difference between model and data. objFunc =
%                      sum(residuals.^2).
%
%Khuloud Jaqaman, July 2014
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

objFunc = [];
residuals = [];

%% Input

%check whether correct number of input arguments was used
if nargin < 2
    disp('--multiStepFunction: Incorrect number of input arguments!');
    return
end

if nargin < 3 || isempty(y)
    y = [];
end

%% Calculation

%get number of steps
numSteps = (length(param)-1)/2;

%get step locations
stepPos = param(1:numSteps);

%get function values between steps
betweenStepsVal = param(numSteps+1:end);

%go over x and calculate function
multiStepVal = betweenStepsVal(1)*ones(size(x));
for iStep = 1 : numSteps
    multiStepVal(x>stepPos(iStep)) = betweenStepsVal(iStep+1);
end

%take difference with data and calculate objective function
if isempty(y)
    residuals = [];
    objFunc = [];
else
    residuals = multiStepVal - y;
    objFunc = residuals' * residuals;
end


%% ~~~ the end ~~~

