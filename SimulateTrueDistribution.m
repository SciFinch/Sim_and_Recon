function [trueSinograms] = SimulateTrueDistribution(transEmMaps,projectorType,pixelInfo)
%Usage: [trueSinograms] = SimulateTrueDistribution(transEmMaps,projectorType,pixelInfo);
% This function finds an estimate for the true emission distribution, as observed by the scanner.
% This is effectively a wrapper function for FwdProject, but handles any additional processing.

% Notify User
fprintf(' - Simulating true sinogram(s)\n')

% Find number of images to be projected
if iscell(transEmMaps)
	nProjs = numel(transEmMaps);
else
	nProjs = 1;
	thisImgCell{1} = transEmMaps;
	transEmMaps = thisImgCell;
	clear thisImgCell;
end

% Forward project emission maps:
trueSinograms = FwdProject(transEmMaps,projectorType,pixelInfo);

end

