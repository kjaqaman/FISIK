function [stepSizeDown,stepSizeUp,fracSteps,numStepsT,segmentInt,flagDisappear,...
    stepStatusDown,stepStatusUp] = getPhotoBleachStepSize(tracksFinal,...
    numStepsRange,alpha,doPlot,checkFit,minLength)
%GETPHOTOBLEACHSTEPSIZE determines photobleaching step size by fitting multi-step functions to single-molecule tracks
%
%SYNOPSIS [stepSizeDown,stepSizeUp,fracSteps,numStepsT,segmentInt,flagDisappear,...
%    stepStatusDown,stepStatusUp] = getPhotoBleachStepSize(tracksFinal,...
%    numStepsRange,alpha,doPlot,checkFit,minLength)
%
%INPUT  tracksFinal  : Tracker output. Structure array or matrix.
%       numStepsRange: Row vector indicating minimum and maximum number
%                      steps to attempt to fit.
%                      Optional. Default: [0 5].
%       alpha        : Alpha-value for F-test to determine number of steps.
%                      Optional. Default: 0.05.
%       doPlot       : 1 to plot step size histogram, 0 otherwise.
%                      Optional. Default: 0.
%       checkFit     : 1 to visually check fit for all tracks, 2 to
%                      visually check fit only for tracks where a step
%                      change has been detection, and 0 otherwise.
%                      Optional. Default: 0.
%       minLength    : Minimum length of track to use for analysis.
%                      Optional. Default: 25.
%
%OUTPUT stepSizeDown : Distribution of obtained step sizes going down.
%       stepSizeUp   : Distribution of obtained step sizes going up.
%       fracSteps    : (maximum number of steps)-by-(4) matrix. 
%                      Column 1: Number of steps in a track.
%                      Column 2: Fraction of tracks with that number of
%                      steps, regardless of whether down or up.
%                      Column 3: Fraction of tracks with that number of
%                      steps going down.
%                      Column 4: Fraction of tracks with that number of
%                      steps going up.
%       numStepsT    : (number of tracks)-by-3 array indicating number of
%                      steps found in each track. Columns show total, down,
%                      up.
%       segmentInt   : (number of tracks)-by-(maximum number of steps)
%                      array storing mean track intensity between steps,
%                      starting from last segment (after last step)
%                      and going back to first segment (before first step).
%       flagDisappear: Vector of size = number of tracks, with value 1 if
%                      track ends before end of movie, and 0 otherwise.
%       stepStatusDown:Same size as stepSizeDown, indicating step status (1
%                      = good = real step down; 0 = bad = gradual
%                      decrease).
%       stepStatusUp : Same size as stepSizeUp, indicating step status (1
%                      = good = real step up; 0 = bad = gradual increase).
%
%Khuloud Jaqaman, August 2014
%2024 April: Modified heavily by KJ
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
stepSizeDown = [];
stepSizeUp = [];
fracSteps = [];

%% Input

%check whether correct number of input arguments was used
if nargin < 1
    disp('--getPhotoBleachStepSize: Incorrect number of input arguments!');
    return
end

if nargin < 2 || isempty(numStepsRange)
    numStepsRange = [0 5];
end

if nargin < 3 || isempty(alpha)
    alpha = 0.05;
end

if nargin < 4 || isempty(doPlot)
    doPlot = 0;
end

if nargin < 5 || isempty(checkFit)
    checkFit = 0;
end

if nargin < 6 || isempty(minLength)
    minLength = 25;
end

%% Calculation

%retain only tracks without any merging or splitting
if isstruct(tracksFinal)
    numSegments = getNumSegmentsPerTrack(tracksFinal);
    tracksFinal = tracksFinal(numSegments == 1);
end

%get amplitude matrix from tracksFinal
if isstruct(tracksFinal)
    tracksMat = convStruct2MatNoMS(tracksFinal);
else
    tracksMat = tracksFinal;
end
ampMat = tracksMat(:,4:8:end)';
[numFrames,numTracks] = size(ampMat);

%get track lifetime and end time
trackSEL = getTrackSEL(tracksMat);

%retain only tracks longer than minLength for analysis
% trackLength = zeros(numTracks,1);
% for iTrack = 1 : numTracks
%     trackLength(iTrack) = length(find(~isnan(ampMat(:,iTrack))));
% end
indxGood = find(trackSEL(:,3)>= minLength);
keepTrack = zeros(numTracks,1);
keepTrack(indxGood) = 1;

%indicate which tracks end before end of movie
flagDisappear = trackSEL(:,2) < numFrames;

numSteps = zeros(numTracks,3);
stepSize = NaN(numTracks,numStepsRange(2));
stepStatus = stepSize;
segmentInt = NaN(numTracks,numStepsRange(2)+1);

for iTrack = indxGood'
    
    %fit amplitude time series
    [stepX,valY,goodStep] = fitDataWithMultiStepFunc((1:numFrames)',ampMat(:,iTrack),numStepsRange,alpha,0);
    
    %get step sizes and number of steps
    %also save step status (good or bad)
    numStepsTmp = length(stepX);
    if numStepsTmp > 0
        stepSizeTmp = -diff(valY)';
        stepSize(iTrack,1:numStepsTmp) = stepSizeTmp;
        numSteps(iTrack,:) = [numStepsTmp sum(stepSizeTmp>0) sum(stepSizeTmp<0)];
        stepStatus(iTrack,1:numStepsTmp) = goodStep;
    end
    
    %save intensity of segments between photobleaching steps.
    %start with last segment before disappearance and go back in time
    segmentInt(iTrack,1:numStepsTmp+1) = valY(end:-1:1);
    
    %check fit if requested
    switch checkFit
        
        case 1 %check all tracks, whether with or without intensity step changes
            
            x = (1:numFrames)';
            y = ampMat(:,iTrack);
            indxFirst = find(~isnan(y),1,'first');
            indxLast = find(~isnan(y),1,'last');
            x = x(indxFirst:indxLast);
            y = y(indxFirst:indxLast);
            overlayMultiStepFunc(x,y,stepX,valY,['Track ' num2str(iTrack)],goodStep)
            userEntry = input(['Step Detection for Track ' num2str(iTrack) ' OK? y/n '],'s');
            close
            if strcmp(userEntry,'n')
                keepTrack(iTrack) = 0;
            end
            
        case 2 %check only tracks where a step change in intensity has been detected
            
            if numStepsTmp > 0
                
                x = (1:numFrames)';
                y = ampMat(:,iTrack);
                indxFirst = find(~isnan(y),1,'first');
                indxLast = find(~isnan(y),1,'last');
                x = x(indxFirst:indxLast);
                y = y(indxFirst:indxLast);
                overlayMultiStepFunc(x,y,stepX,valY,['Track ' num2str(iTrack)],goodStep)
                userEntry = input(['Step Detection for Track ' num2str(iTrack) ' OK? y/n '],'s');
                close
                if strcmp(userEntry,'n')
                    keepTrack(iTrack) = 0;
                end
                
            end
    end
    
end

%copy numSteps for output
numStepsT = numSteps;

%keep info for good tracks only
indxGood = find(keepTrack);
numSteps = numSteps(indxGood,:);
stepSize = stepSize(indxGood,:);
stepStatus = stepStatus(indxGood,:);

%re-format for output
stepSizeDown = stepSize(stepSize > 0);
stepStatusDown = stepStatus(stepSize > 0);
stepSizeUp = -stepSize(stepSize < 0);
stepStatusUp = stepStatus(stepSize < 0);
fracSteps = [(histcounts(numSteps(:,1),numStepsRange(1)-0.5:numStepsRange(2)+0.5))' ...
    (histcounts(numSteps(:,2),numStepsRange(1)-0.5:numStepsRange(2)+0.5))' ...
    (histcounts(numSteps(:,3),numStepsRange(1)-0.5:numStepsRange(2)+0.5))'];
fracSteps = fracSteps / size(numSteps,1);
fracSteps = [(numStepsRange(1):numStepsRange(2))' fracSteps];

%% Plotting

if doPlot
    
    figure
    subplot(1,2,1)
    histogram(stepSizeDown)
    title('Histogram of detected steps DOWN')
    subplot(1,2,2)
    histogram(stepSizeUp)
    title('Histogram of detected steps UP')

    figure
    subplot(1,3,1)
    plot(fracSteps(:,1),fracSteps(:,2),'o-')
    title('Fraction of tracks with specified number of steps TOTAL')
    subplot(1,3,2)
    plot(fracSteps(:,1),fracSteps(:,3),'o-')
    title('Fraction of tracks with specified number of steps DOWN')
    subplot(1,3,3)
    plot(fracSteps(:,1),fracSteps(:,4),'o-')
    title('Fraction of tracks with specified number of steps UP')

end


%% ~~~ the end ~~~

