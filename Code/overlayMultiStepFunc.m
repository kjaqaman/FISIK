function overlayMultiStepFunc(x,y,stepX,valY,figTitle,goodStep)
%OVERLAYMULTISTEPFUNC overlays a multi-step function on time series
%
%SYNOPSIS overlayMultiStepFunc(x,y,stepX,valY)
%
%INPUT  x            : Independent variable of time series.
%       y            : Dependent variable of time series.
%       stepX        : x-values at which there is a step.
%       valY         : y-values between steps. If there are n steps, there
%                      will be n+1 y-values. First value is before first
%                      step, etc., and last value is after last step.
%       figTitle     : Figure title. Optional.
%       goodStep     : 1 to indicate real step, 0 to indicate gradual
%                      increase/decrease. 
%                      Optional. Default: All ones.
%
%OUTPUT the plot
%
%Khuloud Jaqaman, August 2014
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


if nargin < 5 || isempty(figTitle)
    figure
else
    figure('Name',figTitle,'NumberTitle','off')
end
hold on

if nargin < 6 || isempty(goodStep)
    goodStep = ones(size(stepX));
end

plot(x,y,'marker','.')

numSteps = length(stepX);
if numSteps == 0
    plot([x(1) x(end)],valY*[1 1],'r');
else
    plot([x(1) stepX(1)],valY(1)*[1 1],'r');
    if goodStep(1)
        plot(stepX(1)*[1 1]',valY(1:2),'r');
    else
        plot(stepX(1)*[1 1]',valY(1:2),'k');
    end
    for iStep = 1 : numSteps-1
        plot([stepX(iStep) stepX(iStep+1)],valY(iStep+1)*[1 1],'r');
        if goodStep(iStep)
            plot(stepX(iStep+1)*[1 1]',valY(iStep+1:iStep+2),'r');
        else
            plot(stepX(iStep+1)*[1 1]',valY(iStep+1:iStep+2),'k');
        end
    end
    plot([stepX(end) x(end)],valY(end)*[1 1],'r');
end


%% ~~~ the end ~~~

