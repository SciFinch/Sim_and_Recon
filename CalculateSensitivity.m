function [sensitivityImage] = CalculateSensitivity(attAndNormFactors,dHF,dAP,dRL,projectorType,interpType,pixelInfo)
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
if iscell(dHF)
	nTransformations = numel(dHF);
else
	nTransformations = 1;
	dHF{1} = dHF;
	dAP{1} = dAP;
	dRL{1} = dRL;
end

if iscell(attAndNormFactors)
	nFactorImages = numel(attAndNormFactors);
else
	nFactorImages = 1;
	attAndNormFactors{1} = attAndNormFactors;
end

%% Backproject the sensitivity factor sinograms(s)
% NOTE: The BckProject function automatically pads the resulting image(s)
if (nFactorImages == 1) && (nTransformations > 1)
	% Transform same image many times (only 1 backprojection required)
	temp = BckProject(attAndNormFactors,projectorType,interpType,pxinfo);
	for it = 1:nTransformations
		normImages{it} = temp;
	end
elseif (nFactorImages == nTransformations) && (nTransformations > 1)
	% Need to backproject each sens factor sinogram prior to transformation
	normImages = BckProject(attAndNormFactors,projectorType,interpType,pxinfo);
else
	error('Inconsistent number of transformations relative to number of sensitivity sinograms');
end

%% Transform and generate sensitivity image
% (Note that TransformImage can operate with a cell of multiple images)
sensitivityContributions = TransformImage(normImages,pxinfo,dHF,dAP,dRL,interpType,0);

sensitivityImage = zeros(pixelInfo.pxSizePadded);
for it = 1:nTransformations
	sensitivityImage = sensitivityImage + sensitivityContributions{it};
end


end

