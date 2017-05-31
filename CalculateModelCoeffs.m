function [coeffs_RL, coeffs_AP, coeffs_HF] = CalculateModelCoeffs(dRL,dAP,dHF,respSigVals,modelTypeStr,pixelInfo)
% Usage: [coeffs_RL, coeffs_AP, coeffs_HF] = CalculateModelCoeffs(dRL,dAP,dHF,respSigVals,modelTypeStr,pixelInfo)
%
%	Given a co-labelled set of respiratory signal values and the corresponding 
%	motion fields, this function will find a set of coefficients of the 
%	specified type. This works for fwd, bck, or inv fields.
%	This uses the linear Vandermonde method for finding polynomial
%	coefficients.
%	Note: Coefficients are saved with the polynomial labels:
%			P(x) = a_0 x^0 + a_1 x^1 + ... + a_M x^M
%		  i.e., coeffs_RL(:,1) are the zeroeth coeffs,
%				coeffs_RL(:,2) are the first coeffs, etc.
%	Warning: Currently, only scalar surrogates are accepted for 
%			 model formation. These are expected in [1 x N]
%			 format.
%

% error checking
if size(respSigVals.forModel,1) > size(respSigVals.forModel,2)
	error('Respiratory signal values are expected in [1 x N] format')
end

if iscell(dHF)
	nMotionFields = numel(dHF);
elseif size(dHF,4) == 1
	warning('Expected motion fields in cell format');
	nMotionFields = 1;
	dRL{1} = dRL;
	dAP{1} = dAP;
	dHF{1} = dHF; 
else
	error('Expected motion fields in cell format');
end

% Check type of model to use
switch modelTypeStr
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

% Regression according to the Vandermonde method
vandermondeMat = bsxfun(@power,respSigVals.forModel',0:polyOrder);
B = pinv(vandermondeMat)*vandermondeMat;
C = B\pinv(vandermondeMat); %<-this is used to find coeffs of P(x), given y
clear vandermondeMat B;


% AP displacements
% - Sampled y values:
currMFs = zeros(prod(pixelInfo.pxSize),nMotionFields);
for tt = 1:nMotionFields
    currMFs(:,tt) = dAP{tt}(:);
end
% - Find least-squares coefficients, given sampled x values
coeffs_AP = C*currMFs';
clear currMFs;

% RL displacements
% - Sampled y values:
currMFs = zeros(prod(pixelInfo.pxSize),nMotionFields);
for tt = 1:nMotionFields
    currMFs(:,tt) = dRL{tt}(:);
end
% - Find least-squares coefficients, given sampled x values
coeffs_RL = C*currMFs';
clear currMFs;

% HF displacements
% - Sampled y values:
currMFs = zeros(prod(pixelInfo.pxSize),nMotionFields);
for tt = 1:nMotionFields
    currMFs(:,tt) = dHF{tt}(:);
end
% - Find least-squares coefficients, given sampled x values
coeffs_HF = C*currMFs';
clear currMFs;

end