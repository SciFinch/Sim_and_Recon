function [motionModel] = FormMotionModel(datasetStr,modelTypeStr,pixelInfo,saveFlag)
% Usage: [motionModel] = FormMotionModel(datasetStr,modelTypeStr,pixelInfo,saveFlag)
%
%	This function forms a motion model, given a set of motion fields.
%	These fields are loaded automatically and used in conjunction
%	with sampled respiratory signal values (scalar) to find a set
% 	of coefficients which parameterise the motion of each voxel in
%	the field of view.
%


fprintf('Generating Motion Model for dataset: %s\n',datasetStr);

motionModel.type = modelTypeStr;
[respSigVals] = ReadRespsigs(datasetStr);
motionModel.minExpectedInput = min(respSigVals.forModel);
motionModel.maxExpectedInput = max(respSigVals.forModel);

%% Calculate forward model
% (these should always exist...)
fprintf('- Loading motion fields for forward model\n');
[dHF,dAP,dRL] = ReadMotionFields(datasetStr,'fwd',pixelInfo);
fprintf('- Calculating coefficients for forward model\n');
[motionModel.coeffs.RL.fwd,motionModel.coeffs.AP.fwd, motionModel.coeffs.HF.fwd] = ...
    CalculateModelCoeffs(dRL,dAP,dHF,respSigVals,modelTypeStr,pixelInfo);
clear dHF dAP dRL;

%% Inverse and Backward models
% Note: Depending on the dataset, either or both of these may exist.
% It's best to check, and handle cases where only one exists.
fprintf('- Loading motion fields for backward model\n');
[dHF,dAP,dRL] = ReadMotionFields(datasetStr,'bck',pixelInfo);
if isempty(dHF{1})
    warning('- No backwards MFs found. Skipping backward model');
    % Stub (so that the fields exist)
    motionModel.coeffs.RL.bck = [];
    motionModel.coeffs.AP.bck = [];
    motionModel.coeffs.HF.bck = [];
else
    fprintf('- Calculating coefficients for backward model\n');
    [motionModel.coeffs.RL.bck,motionModel.coeffs.AP.bck, motionModel.coeffs.HF.bck] = ...
        CalculateModelCoeffs(dRL,dAP,dHF,respSigVals,modelTypeStr,pixelInfo);
    clear dHF dAP dRL;
end

fprintf('- Loading motion fields for inverse model\n');
[dHF,dAP,dRL] = ReadMotionFields(datasetStr,'inv',pixelInfo);
if isempty(dHF{1})
    warning('- No inverse MFs. Skipping inverse model');
    % Stub (so that the fields exist)
    motionModel.coeffs.RL.inv = [];
    motionModel.coeffs.AP.inv = [];
    motionModel.coeffs.HF.inv = [];
else
    fprintf('- Calculating coefficients for inverse model\n');
    [motionModel.coeffs.RL.inv,motionModel.coeffs.AP.inv, motionModel.coeffs.HF.inv] = ...
        CalculateModelCoeffs(dRL,dAP,dHF,respSigVals,modelTypeStr,pixelInfo);
    clear dHF dAP dRL;
end

%% Saving (optional)
if saveFlag == 1
    fprintf('- Saving model\n');
    datasetInfoStruct = GetDatasetInfo(datasetStr);
    if exist(datasetInfoStruct.motionModels.all,'file')
        warning('- Overwriting existing motion model...')
    end
    save(datasetInfoStruct.motionModels.all,'motionModel','-v7.3');       
end