function [dHF,dAP,dRL] = EstimateFieldsFromModel(respSig,motionModel,regStr,pixelInfo)
%EstimateFieldsFromModel [dHF,dAP,dRL] = EstimateFieldsFromModel(respSig,model,pixelInfo)
%   Given the inputs for a motion model, this function produces the
%   transformation fields for each direction. This relies upon the model
%   being in struct format, and each direction being labelled with respect
%   to the subject's axes
%
%	NOTE: Fields will be automatically padded with zeros to match extFOV format

% TO DO
% 1. Make multiple models available
% 2. Build a model class and infrastructure to save/load in a specific way

if size(respSig,1) > size(respSig,2)
	error('Respiratory signal values are expected in [1 x N] format')
end
nTransformations = numel(respSig); %check this is the right direction

% Calculate the order of the motion model polynomials
switch motionModel.type
    case 'none'
        error('no-model case not yet implemented')
    case 'linear'
        polyOrder = 1;
    case 'quadratic'
        polyOrder = 2;
    case 'cubic'
        polyOrder = 3;
    otherwise
        error('unrecognised model type');
end

% streamlining 
switch regStr
	case 'fwd'
		RLcoeffs = motionModel.coeffs.RL.fwd;
		APcoeffs = motionModel.coeffs.AP.fwd;
		HFcoeffs = motionModel.coeffs.HF.fwd;
	case 'bck'
		RLcoeffs = motionModel.coeffs.RL.bck;
		APcoeffs = motionModel.coeffs.AP.bck;
		HFcoeffs = motionModel.coeffs.HF.bck;
	case 'inv'
		RLcoeffs = motionModel.coeffs.RL.inv;
		APcoeffs = motionModel.coeffs.AP.inv;
		HFcoeffs = motionModel.coeffs.HF.inv;
	otherwise
		error('Unrecognised registration string. Must be one of {"fwd","bck","inv"}');
end

% Perform estimation for each voxel, in each direction
% (note, padding is performed *after* field estimation)
vandermondeMat = bsxfun(@power,respSig',0:polyOrder);

rawRL = vandermondeMat*RLcoeffs;
rawAP = vandermondeMat*APcoeffs;
rawHF = vandermondeMat*HFcoeffs;

for t = 1:nTransformations 
    dRL{t} = reshape( rawRL(t,:), pixelInfo.pxSize );
    dAP{t} = reshape( rawAP(t,:), pixelInfo.pxSize );
    dHF{t} = reshape( rawHF(t,:), pixelInfo.pxSize );
end
clear vandermondeMat rawRL rawAP rawHF;

% Pad motion fields to match extFOV
dRL = ToggleImagePadding(dRL,pixelInfo);
dHF = ToggleImagePadding(dHF,pixelInfo);
dAP = ToggleImagePadding(dAP,pixelInfo);

end

