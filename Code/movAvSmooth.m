function [newTraj,trend,errFlag] = movAvSmooth(traj,window,boundary,doPlot)
%movAvSmooth detects trends in a series by calculating a 2-sided moving average
%
%SYNPOSIS [newTraj,trend,errFlag] = movAvSmooth(traj,window,boundaryOption)
%
%INPUT  traj    : Observations of time series to be detrended. A structure 
%                 array with field:
%           .observations: 2D array of measurements and their uncertainties.
%                          Missing points should be indicated with NaN.
%                          If there is only one series, it can also be
%                          input directly as a 2D array.
%       window  : Number of points to be used on each side for averaging.
%                 Optional. Default: 2.
%       boundary: Input specifying the way to pad trajetory end points.
%                 1: Use the end point values.
%                 2: Use NaN.
%                 Optional. Default: 1.
%       doPlot  : 1 to make plots of trajectories and trends, 0 otherwise.
%                 Optional. Default: 0.
%
%OUTPUT newTraj: De-trended series.
%       trend  : Estimated trend.
%       Both newTraj and trend will be in the same format as traj.
%       errFlag: 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, October 2004
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

newTraj = [];
trend = [];
errFlag = 0;

%% Input

%check if correct number of arguments were used when function was called
if nargin < 1
    disp('--movAvSmooth: Incorrect number of input arguments.');
    errFlag  = 1;
    return
end

%check input trajectory
if isstruct(traj)
    dataVec = 0;
    if ~isfield(traj,'observations')
        disp('--movAvSmooth: Input traj needs field ''observations''.')
        errFlag = 1;
    end
else
    dataVec = 1;
    tmp = traj;
    clear traj
    traj.observations = tmp;
    clear tmp
end
numTraj = length(traj);

%get trajectory lengths and make sure all are column vectors
trajLength = zeros(numTraj,1);
for iTraj = 1 : numTraj
    [nRow,nCol] = size(traj(iTraj).observations);
    if nCol > 2
        disp('--movAvSmooth: Trajectories should be column vectors.');
        errFlag = 1;
    end
    trajLength(iTraj) = nRow;
end

%check window size
if nargin < 2 || isempty(window)
    window = 2;
else
    if window <= 0 %check that it has accetable value
        disp('--movAvSmooth: Number of points used for averaging on each side should be >= 1.');
        errFlag = 1;
    end
end

%check boundary padding option
if nargin < 3 || isempty(boundary)
    boundary = 1;
else
    if boundary ~=1 && boundary ~= 2
        disp('--movAvSmooth: Input argument "boundary" can take only values 1 and 2.');
        errFlag = 1;
    end
end

%check whether to make plots
if nargin < 4 || isempty(doPlot)
    doPlot = 0;
end
        
%exit if there are problems in input data
if errFlag
    disp('--movAvSmooth: Please fix input data.');
    return
end

%% Trend calculation and subtraction

newTraj = traj;
trend = traj;

%go over all trajetories
for iTraj = 1 : numTraj
    
    %get current trajectory
    trajCurrent = traj(iTraj).observations;
    newTrajCurrent = trajCurrent;
    
    %append trajectory at both ends to allow trend calculation for points
    %less than "window" points away from either end
    switch boundary
        case 1
            trajCurrent1 = [trajCurrent(1,1)*ones(window,1); trajCurrent(:,1); trajCurrent(end,1)*ones(window,1)];
        case 2
            trajCurrent1 = [NaN*ones(window,1); trajCurrent(:,1); NaN*ones(window,1)];
    end
    
    %compute moving average
    trendCurrent = NaN*ones(trajLength(iTraj),1);
    for i = 1 : trajLength(iTraj)
        trendCurrent(i) = nanmean(trajCurrent1(i:i+2*window));
    end
    
    %subtract trend from trajectory
    trajCurrent1 = trajCurrent1(window+1:end-window);
    newTrajCurrent(:,1) = trajCurrent1 - trendCurrent;
    
    %store in output structures
    newTraj(iTraj).observations = newTrajCurrent;
    trend(iTraj).observations = trendCurrent;
    
end

%% Plotting

if doPlot
    for iTraj = 1 : numTraj
        
        if isfield(traj,'time')
            timeVec = traj(iTraj).time(:,1);
        elseif isfield(traj,'timePoint')
            timeVec = traj(iTraj).timePoint;
        else
            timeVec = (1 : trajLength(iTraj))';
        end
        
        figure
        
        %plot original trajectory and trend
        subplot(1,2,1)
        plot(timeVec,traj(iTraj).observations(:,1))
        if size(traj(iTraj).observations,2) == 2
            myErrorbar(timeVec,traj(iTraj).observations(:,1),traj(iTraj).observations(:,2))
        end
        hold on
        plot(timeVec,trend(iTraj).observations,'r')
        legend('Original','Trend')
        xlabel('Time')
        ylabel('Values')
        
        %plot de-trended trajectory
        subplot(1,2,2)
        plot(timeVec,newTraj(iTraj).observations(:,1))
        if size(newTraj(iTraj).observations,2) == 2
            myErrorbar(timeVec,newTraj(iTraj).observations(:,1),newTraj(iTraj).observations(:,2))
        end
        legend('Detrended')
        xlabel('Time')
        ylabel('Values')
        
    end
end

%% Final touches

%put output as vectors if input was a vector
if dataVec
    newTraj = newTrajCurrent;
    trend = trendCurrent;
end

%% %%% ~~ the end ~~ %%%%%

