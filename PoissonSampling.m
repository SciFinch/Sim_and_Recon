function [noisySinograms] = PoissonSampling(meanDistbnSino,totCounts,dynamicWeights,rescaleFlag)
%POISSONSAMPLING [sim.noisyPrompts] = PoissonSampling(sim.prompts,totCounts,dynamicWeights,rescaleFlag);
%   The purpose of this function is to sample a mean distribution in a
%   sinogram using a poisson number generator. This has the capacity to be
%   very slow, so a placeholder might involve MATLAB's native PoissRnd
%   function (for now), although I suspect this isn't exactly what's
%   needed.
% 	rescaleFlag allows user-specification on whether resulting sinograms
%	will retain their original total values or the new scaled values from
% 	the noise generation process [the former is more convenient but might
%	not be Poisson (check this!)]
%
%   The user specifies the total number of counts across ALL sinograms
%   contained in the cells of meanDistbnSino, which will be split according
%   to the fractions within dynamicWeights (this is because the number of
%   counts across dynamic sinograms will generally vary due to motion,
%   radioactive decay, washout, etc.).

% Compatibility checks:
if iscell(meanDistbnSino)
	nToMakeNoisy = numel(meanDistbnSino);
elseif size(nToMakeNoisy,4) == 1
	nToMakeNoisy = 1;
	meanDistbnSino{1} = meanDistbnSino;
else
	error('Expected meanDistbnSino in cell format');
end

if numel(dynamicWeights) ~= nToMakeNoisy
    error('the number of weights provided must equal the number of sinograms');
end

if abs(sum(dynamicWeights) - 1) > 1E-4
    error('dynamicWeights must sum to unity');
end

% Perform noise simulation
warning('PoissonSampling is currently MATLAB-based')
warning('Need to come back to PoissonSampling and ensure consistency with emission scalefactors');

% Find total counts in mean-value sinograms
currTotCounts = 0;
for it = 1:nToMakeNoisy
	currTotCounts = currTotCounts + sum(meanDistbnSino{it}(:));
end

noisySinograms = cell(1,nToMakeNoisy);
for g = 1:nToMakeNoisy
	% Current sinogram [normalised so that sum(dynamicSino) = 1]:
	toMakeNoisy = meanDistbnSino{it}/currTotCounts;

	% Now scale to desired number of counts [this is only separate to the above for readability]
	% Note that this depends not just on total counts, but the fraction of time spent in the position
	toMakeNoisy = toMakeNoisy * (totCounts * dynamicWeights(g));

	% Run noise generation algorithm
	currNoisySino = poissrnd(toMakeNoisy);

	% Optionally rescale noisy sinogram so that the reconstructed values remain the same
	if rescaleFlag == 1
		% In this case, sinograms will be sparsely non-zero, but non-zero entries will be very
		% large and non-integer - this is basically a 'scaled Poisson' noise realisation, whatever that is.
	    noisySinograms{g} = currNoisySino * currTotCounts/(totCounts * dynamicWeights(g));
	else
		% This will have Poisson noise in a more realistic sense - non-zero entries will be
		% small integers, but reconstructed values might not match those in the FDG map(s).
	    noisySinograms{g} = currNoisySino;
	end
end

end

