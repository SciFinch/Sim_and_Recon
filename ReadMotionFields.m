function [dHF,dAP,dRL] = ReadMotionFields(datasetStr,regStr,pixelInfo)
%READMOTIONFIELDS [dHF,dAP,dRL] = ReadMotionFields(datasetStr,pixelInfo)
%   Given a dataset identifier, datasetStr, this function will load up the
%   volumes containing the motion fields for a given dataset. This is done
%   in a function so that the output can be standardised - the output will
%   be 3 cell arrays containing the motion transformations.
%
%   regStr specifies whether the motion transformations should be fwd
%   (i.e., ref -> position), or bck/inv (i.e, position -> ref). In each
%   case, the label for the cell array is the position label.
%   [Note: In many cases, bck != inv transformations. Generally, for PET
%   reconstruction, it is best to use inverse fields, otherwise
%   inconsistencies between fwd and bck projection can lead to artefacts.]
%
%   NOTE: Before including a new dataset in this function, make sure it is
%   properly and rigorously tested to identify the correct orientation of
%   the vectors relative to the body, as labelled.

datasetInfoStruct = GetDatasetInfo(datasetStr);

nDynamics = datasetInfoStruct.nDynamics;

% Find the relevant filenames to load
switch regStr
    case 'fwd'
        filenames.RL = datasetInfoStruct.motionFields.RL.fwd;
        filenames.HF = datasetInfoStruct.motionFields.HF.fwd;
        filenames.AP = datasetInfoStruct.motionFields.AP.fwd;
        orientationFlag = datasetInfoStruct.flags.orientationChecked.MF.fwd;
    case 'bck'
        filenames.RL = datasetInfoStruct.motionFields.RL.bck;
        filenames.HF = datasetInfoStruct.motionFields.HF.bck;
        filenames.AP = datasetInfoStruct.motionFields.AP.bck;
        orientationFlag = datasetInfoStruct.flags.orientationChecked.MF.bck;        
    case 'inv'
        filenames.RL = datasetInfoStruct.motionFields.RL.inv;
        filenames.HF = datasetInfoStruct.motionFields.HF.inv;
        filenames.AP = datasetInfoStruct.motionFields.AP.inv;
        orientationFlag = datasetInfoStruct.flags.orientationChecked.MF.inv;
    otherwise
        error('Unable to load MFs: Unrecognised registration string provided');
end

if orientationFlag == false
    warning('The motion fields being loaded have not had their orientation checked!')
end

% Load motion fields
dRL = cell(1,nDynamics);
dHF = cell(1,nDynamics);
dAP = cell(1,nDynamics);

% Check that you're getting something compatible with the current project
[dHF{1},imageSize,imagePixelDims] = AutoLoadImage(filenames.HF{1});
if sum(abs(pixelInfo.imageSize - imageSize)) > 1E4
    warning('Motion fields not same size as current image volumes');
end
if sum(abs(pixelInfo.imagePixelDims - imagePixelDims)) > 1E4
    warning('Motion fields not same pixel scale as current image volumes');
end

% Finally, load all fields
for it = 1:nDynamics
    [dRL{it},~,~] = AutoLoadImage(filenames.RL{it});
    [dHF{it},~,~] = AutoLoadImage(filenames.HF{it});
    [dAP{it},~,~] = AutoLoadImage(filenames.AP{it});
end

% Check that the scaling of the motion fields is sensible
if median(dHF{it}) < norm(imagePixelDims)
    warning('Motion fields seem scaled incorrectly, since median displacement < pixel size');
end

end

