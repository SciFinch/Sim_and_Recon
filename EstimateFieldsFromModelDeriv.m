function [dHF,dAP,dRL] = EstimateFieldsFromModelDeriv(respSig,model,pixelInfo)
%EstimateFieldsFromModelDeriv [dHF,dAP,dRL] = EstimateFieldsFromModelDeriv(respSig,model,pixelInfo)
%   Given the inputs for a motion model, this function produces the DERIVATIVE
%   transformation fields for each direction. This relies upon the model
%   being in struct format, and each direction being labelled with respect
%	NOTE: Fields will be automatically padded with zeros to match extFOV format

% TO DO
% 1. Incorporate this into EstimateFieldsFromModelDeriv somehow

%warning: this only works for pairs of  fwd/inv models, and currently only
%works for the basic 2nd-order polynomial model used in the ME-MCIR work

nTransformations = numel(respSig); %check this is the right direction
warning('Directions need checking in EstimateFieldsFromModelDeriv');

% streamlining for parallelisation
RLcoeffs = model.coeffs.RL.fwd;
APcoeffs = model.coeffs.AP.fwd;
HFcoeffs = model.coeffs.HF.fwd;

% Perform derivative of 2nd-order polynomial estimation for each voxel, 
% in each direction
vandermondeMat = bsxfun(@power,respSig',-1:1);
vandermondeMat(:,3) = 2*vandermondeMat(:,3);
vandermondeMat = vandermondeMat(:,2:end);

%TODO: Generalise this to any function, and tie in with
%EstimateFieldsFromModel

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

