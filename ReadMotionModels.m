function [mmHF,mmAP,mmRL] = ReadMotionModels(datasetStr,regStr,pixelInfo)
%READMOTIONMODELS [mmHF,mmAP,mmRL] = ReadMotionModels(datasetStr,regStr,pixelInfo)
%   As with the ReadMotionFields function, this function returns motion
%   transformations in the format of a motion model. The main purpose of
%   this function is to tidy away all the commenting/explanation required
%   for individual datasets, and to standardise the output for use in
%   simulation/reconstruction code.

% TO-DOs:
% 1. Handle outputs of each motion model mat file
% 2. Recognise which type of motion model data needs to be loaded
% 3. Incorporate outputs into mmRL, mmAP, mmRL
% 4. Check and standardise orientation of each dataset
% 5. Consider standardising this whole process 

datasetInfoStruct = GetDatasetInfo(datasetStr);
warning('ReadMotionModels is currently a stub');

switch datasetStr
case 'CK1'
    load(datasetInfoStruct.motionModels.all)
otherwise
    error('Dataset not yet included in ReadMotionModels');
end
 
mmHF = []; %zeros(pixelInfo.pxSizePadded);
mmAP = []; %zeros(pixelInfo.pxSizePadded);
mmRL = []; %zeros(pixelInfo.pxSizePadded);
    
end

