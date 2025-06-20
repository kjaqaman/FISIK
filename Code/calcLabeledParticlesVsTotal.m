function [probPartLabel,probPartOligoStateLabel] = calcLabeledParticlesVsTotal(probMolLabel,probPartOligoStateFull)
%CALCLABELEDPARTICLESVSTOTAL maps labeled fraction from molecules to particles and maps oligomeric state distribution from total population to labeled subset
%
%SYNPOSIS [probPartLabel,probPartOligoStateLabel] = calcLabeledParticlesVsTotal(probMolLabel,probPartOligoStateFull)
%
%INPUT 
%   probMolLabel          : Probability of labeling a molecule (i.e. 
%                           fraction of molecules that get labeled).
%   probPartOligoStateFull: Vector of probabilities for a particle (i.e.
%                           molecular complex) to be a monomer, dimer, etc.
%                           in the full population of molecules/particles.
%                           Vector length reflects maximum oligomeric
%                           state. Sum of probabilities must equal 1.
%
%OUTPUT
%   probPartLabel         : Probability of a particle (i.e. molecular
%                           complex) getting labeled.
%   probPartOligoStateLabel:Same as probPartOligoStateFull but only for the
%                           labeled subset of particles.
%
%Khuloud Jaqaman, Oct. 2024
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


%% Input

%get maximum oligomeric state
maxOligo = length(probPartOligoStateFull);

%% Calculations

%contribution to each oligomeric state in labeled subset from each oligomeric state in full population
termJI = zeros(maxOligo);
for iOligo = 1 : maxOligo %labeled
    for jOligo = iOligo : maxOligo %full
        termJI(jOligo,iOligo) = binopdf(iOligo,jOligo,probMolLabel)*probPartOligoStateFull(jOligo);
    end
end

%total for each oligomeric state in labeled subset
termIOligo = sum(termJI);

%fraction of labeled particles
probPartLabel = sum(termIOligo);

%fraction of labeled particles with each apparent oligomeric state
probPartOligoStateLabel = termIOligo / probPartLabel;

end