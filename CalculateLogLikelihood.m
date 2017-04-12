function [logLikelihoodVal] = CalculateLogLikelihood(promptSinograms,modelSinograms)
% CALCULATELOGLIKELIHOOD [logLikelihoodVal] = CalculateLogLikelihood(promptSinogram,modelSinogram)
% 	This is a simple function which calculates the Poisson Log-likelihood used in PET
%	reconstruction.

% Consistency checks:
if iscell(promptSinograms) && iscell(modelSinograms)
	nDynamics = numel(promptSinograms);
elseif size(promptSinograms,4) == 1 && size(modelSinograms,4) == 1
	nDynamics = 1;
	promptSinograms{1} = promptSinograms;
	modelSinograms{1} = modelSinograms;
else
	error('Input sinograms are inconsistent');
end

% Calculate the bin values of the log-likelihood for each gate and collect
likelihoodSinogram = zeros(size(promptSinograms{1}));
for it = 1:nDynamics
	currLikelihoodSino = promptSinograms{it}.*log(modelSinograms{it}) - modelSinograms{it};
	currLikelihoodSino(isinf(currLikelihoodSino)|isnan(currLikelihoodSino)) = 0;
	likelihoodSinogram = likelihoodSinogram + currLikelihoodSino;
end
% sum over all intensities
logLikelihoodVal = sum(likelihoodSinogram(:));
warning('Likelihood model for multiple frames has not been verified');

end