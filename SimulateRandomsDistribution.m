function [randomsSinogram] = SimulateRandomsDistribution(transEmMap,RandMthd,projectorType,pixelInfo)
%SIMULATERANDOMSDISTRIBUTION [randomsSinogram] = SimulateRandomsDistribution(transEmMap,RandMthd,projectorType,pixelInfo)
%   This function simulates a distribution of random background radiation. For now, this is a simple
% 	constant-valued sinogram, or an estimate straight from data. This could be simulated properly (?).

warning('SimulateRandomsDistribution is currently a stub');
% ! Should these be scaled according to breathing pattern, or is that done later?

% Figure out how many random sinograms are needed
if iscell(transEmMap)
	nToSimulate = numel(transEmMap);
else
	nToSimulate = 1;
end

% Simulate according to the method suggested
switch RandMthd
case 'none'
	randomsSinogram = cell(1,nToSimulate);
	for it = 1:nToSimulate
		randomsSinogram{it} = zeros(pixelInfo.sino);
	end
case 'simple'
	constSimpleRands = 10;
	randomsSinogram = cell(1,nToSimulate);
	for it = 1:nToSimulate
		randomsSinogram{it} = constSimpleRands*(1./nToSimulate)*ones(pixelInfo.sino);
	end
otherwise
	error('Unrecognised randoms simulation method requested - has it been implemented yet?');
end

end

