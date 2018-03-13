function [sensitivityImage] = CalculateSensitivity(attAndNormFactors,projectorType,interpType,pixelInfo)
%CALCULATESENSITIVITY [sensitivityImage] = CalculateSensitivity(attAndNormFactors,dHF,dAP,dRL,projectorType,interpType,pixelInfo)
%   This function calculates the combined sensitivity for all time bins in a dynamic dataset.
%	This includes motion correction by default (although zeroed motion fields can be used to
%	use this function without correction).
%	NOTE that the motion fields are back/inverse transformations: 
%	motion position -> reference position
%	This function uses the extended-FOV functionality automatically (thus output images
%	will be in padded-image size).

% Note: depending on the algorithm, attAndNormFactors may or may not exist in
% multiple positions (since attFactors can depend on motion position).

% Compatibility checks:
if iscell(attAndNormFactors)
	nFactorImages = numel(attAndNormFactors);
else
	nFactorImages = 1;
	attAndNormFactors{1} = attAndNormFactors;
end
if nFactorImages > 1
    error('Expected only one input sinogram. Did you mean to use CalculateSensitivityWithMotion instead?')
end

%% Backproject the sensitivity factor sinograms(s)
% NOTE: The BckProject function automatically pads the resulting image(s)
sensitivityImage = BckProject(attAndNormFactors,projectorType,interpType,pixelInfo);


end

