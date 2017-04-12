function [dHF,dAP,dRL] = EstimateFieldsFromModel(respSig,model,pixelInfo)
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

nTransformations = numel(respSig); %check this is the right direction

% streamlining for parallelisation
RLcoeffs = model.RL;
APcoeffs = model.AP;
HFcoeffs = model.HF;

% Perform 2nd-order polynomial estimation for each voxel, in each direction
parfor t = 1:nTransformations 
    currRespSig = respSig(t);
    vandeMat = repmat([1 currRespSig currRespSig^2]',[1 prod(pxsize)]);
    
    dRL{t} = reshape( sum(vandeMat.*RLcoeffs), pxinfo );
    dAP{t} = reshape( sum(vandeMat.*APcoeffs), pxinfo );
    dHF{t} = reshape( sum(vandeMat.*HFcoeffs), pxinfo );
end

% Pad motion fields to match extFOV
dRL = ToggleImagePadding(dRL,pxinfo);
dHF = ToggleImagePadding(dHF,pxinfo);
dAP = ToggleImagePadding(dAP,pxinfo);


end

