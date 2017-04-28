function [motionModel] = ReadMotionModels(datasetStr)
%READMOTIONMODELS [motionModel] = ReadMotionModels(datasetStr)
%   As with the ReadMotionFields function, this function returns motion
%   transformations in the format of a motion model. The main purpose of
%   this function is to tidy away all the commenting/explanation required
%   for individual datasets, and to standardise the output for use in
%   simulation/reconstruction code.

% TO-DOs:
% 1. Handle outputs of each motion model mat file
% 2. Recognise which type of motion model data needs to be loaded
% 3. Check and standardise orientation of each dataset

datasetInfoStruct = GetDatasetInfo(datasetStr);
warning('ReadMotionModels is currently a stub');

load(datasetInfoStruct.motionModels.all)

    
end
%'motionModel'