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

% streamlining for parallelisation
RLcoeffs = model.RL;
APcoeffs = model.AP;
HFcoeffs = model.HF;

% Perform 2nd-order polynomial estimation for each voxel, in each direction
parfor t = 1:nTransformations 
    currRespSig = respSig(t);
    vandeMat = repmat([0 1 2*currRespSig]',[1 prod(pixelInfo.pxSize)]);
    
    dRL{t} = reshape( sum(vandeMat.*RLcoeffs), pixelInfo.pxSize );
    dAP{t} = reshape( sum(vandeMat.*APcoeffs), pixelInfo.pxSize );
    dHF{t} = reshape( sum(vandeMat.*HFcoeffs), pixelInfo.pxSize );
end

% Pad motion fields to match extFOV
dRL = ToggleImagePadding(dRL,pxinfo);
dHF = ToggleImagePadding(dHF,pxinfo);
dAP = ToggleImagePadding(dAP,pxinfo);


end

