function [emWeightSinos,muWeightSinos] = CalcWeightedSinos(transEmImages,transMuImages,motionModel,respSigs,pixelInfo,projectorType)
%CALCWEIGHTEDSINOS [emWeightSinos,muWeightSinos] = CalcWeightedSinos(transEmMaps,transMuMaps,motionModel,respSigs,pixelInfo,projectorType)
%
% This function finds the derivatives of the transformed emission estimates and attenuation
% coefficient images with respect to the parameter of the motion model, prior to transforming
% these into sinogram space. These resulting 'weight sinograms' are then used to 
% in the ME-MCIR respiratory signal estimation update.
%

nDynamics = numel(respSigs);

%% Find first derivatives of the motion model
[dHFds,dAPds,dRLds] = EstimateFieldsFromModelDeriv(respSigs,motionModel,pixelInfo);
warning('Gradient directions have not been checked')

%% Find emission image gradients wrt motion model parameter:
for it = 1:nDynamics
	% Find emission image spatial gradients:
	[em_APgrad,em_RLgrad,em_HFgrad] = gradient(transEmImages{it},pixelInfo.pxdims(1),pixelInfo.pxdims(2),pixelInfo.pxdims(3));

	% Chain rule:
	weightIms{it} = (dRLds{it}.*em_APgrad + dAPds{it}.*em_RLgrad + dHFds{it}.*em_HFgrad);
	weightIms{it}(isnan(weightIms{it})) = 0;
end

% Find sinograms from weight images:
emWeightSinos = FwdProject(weightIms,projectorType,pixelInfo);
clear weightIms em_*;

%% Find attenuation coefficient image gradients wrt motion model parameter:
for it = 1:nDynamics
	% Find attenuation coefficient image spatial gradients
	[mu_APgrad,mu_RLgrad,mu_HFgrad] = gradient(transMuImages{it},pixelInfo.pxdims(1),pixelInfo.pxdims(2),pixelInfo.pxdims(3));

	% Chain rule:
	weightMus{it} = (dRLds{it}.*mu_APgrad + dAPds{it}.*mu_RLgrad + dHFds{it}.*mu_HFgrad);
	weightMus{it}(isnan(weightMus{it})) = 0;
end

% Find sinograms from weight mu-maps:
muWeightSinos = FwdProject(weightMus,projectorType,pixelInfo);
clear weightMus mu_*;

end

