function [scatterSinograms] = SimulateScatterDistribution(trueSino,attFactorSino,randSino,ScatterMthd,projectorType,pixelInfo)
%SIMULATESCATTERDISTRIBUTION [scatterSinograms] = SimulateScatterDistribution(trueSino,attFactorSino,randSino,ScatterMthd,projectorType,pixelInfo)
%   This function simulates a scatter distribution, using one of several methods:
%	(a) Simple: Using a simple symmetric monoexponential filter on the true distribution
% 	(b) e7: The scatter estimate generated by e7tools is used in place of a direct sim
% 	(c) Single Scatter Simulation (SSS) algorithm - analytical generation of a scatter 
%	    distribution of photons scattered once prior to detection [Watson 2000]

% TO-DO
% 1. If it can be shown that scatter is ~ independent of position, can streamline the below
% 2. Implement calls to different subfunctions (simple, SSS, montecarlo)
% 3. Revise exactly which datasets are required for each method
% 4. Is transEmMap really needed? Surely trueSino is better; it avoids unnecessary reprojection
% 5. Is this parallelisable?

warning('SimulateScatterDistribution is currently a stub');
% ! Should these be scaled according to breathing pattern, or is that done later?

% Figure out how many scatter sinograms are needed 
if iscell(transEmMap)
	nToSimulate = numel(transEmMap);
else
	nToSimulate = 1;
end

% Simulate according to the method suggested
switch ScatterMthd
case 'none'
	scatterSinograms = cell(1,nToSimulate);
	for it = 1:nToSimulate
		scatterSinograms{it} = zeros(pixelInfo.sino);
	end
% Simple convolution-based scatter simulation (Bailey 1993)
case 'simple'
	scatterSinograms = cell(1,nToSimulate);
	for it = 1:nToSimulate
		fprintf(' - Simulating scatter sinograms (simple): Dynamic sinogram %d/%d',it,nToSimulate);
		scatterSinograms{it} = ConvScatterEstimation(trueSino{it},attFactorSino{it},randSino{it},'Bailey93');
	end
% Single scatter simulation (Watson 2000):
case 'SSS'
	error('SSS not currently implemented');
	%scatterSinograms = cell(1,nToSimulate);
	%for it = 1:nToSimulate
	%	fprintf('== Simulating scatter sinograms (SSS): Dynamic sinogram %d/%d',it,nToSimulate);
	%	scatterSinograms{it} = SingleScatterSimulation(trueSino{it},attFactorSino{it},randSino{it},'Watson00');
	%end
% Monte Carlo simulation (very slow)
case 'MonteCarlo'
	error('Monte Carlo scatter simulation not currently implemented');
	% scatterSinograms = cell(1,nToSimulate);
	% for it = 1:nToSimulate
	% 	fprintf('== Simulating scatter sinograms (Monte Carlo): Dynamic sinogram %d/%d',it,nToSimulate);
	% 	scatterSinograms{it} = MonteCarloScatterEstimation(trueSino{it},attFactorSino{it},randSino{it},'Zaidi01');
	% end
% Error handling:
otherwise
	error('Unrecognised scatter simulation method requested - has it been implemented yet?');
end
w
if iscell(transEmMap)
	for it = 1:numel(transEmMap)
		scatterSinograms{it} = zeros(pixelInfo.sino);
	end
else
	scatterSinograms{1} = zeros(pixelInfo.sino);
end

end

