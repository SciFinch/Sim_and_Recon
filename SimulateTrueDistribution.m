function [trueSinograms] = SimulateTrueDistribution(transEmMaps,projectorType,pixelInfo)
%SIMULATETRUEDISTRIBUTION [trueSinograms] = SimulateTrueDistribution(transEmMaps,projectorType,pixelInfo);
% This function finds an estimate for the true emission distribution, as observed by the scanner.
% This is effectively a wrapper function for FwdProject, but handles any additional processing.

% Notify User
fprintf('== Simulating true sinogram(s)\n')

% Find number of images to be projected
if iscell{transEmMaps}
	nProjs = numel(transEmMaps);
else
	nProjs = 1;
end

% Forward project emission maps:
trueSinograms = cell(1,nProjs);
for it = 1:nProjs
	trueSinograms{it} = FwdProject(transEmMaps{it},projectorType,pixelInfo);
end

end

