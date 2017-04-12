function [respSigs] = ReadRespsigs(datasetStr)
%READRESPSIGS [respSigs] = ReadRespsigs(datasetStr)
%   A cleanup and standardisation function for outputting the respiratory
%   signals used in motion modelling and correction with various datasets.
%   This is stratified into 3 types: all values, values used for training a
%   motion model for that dataset, and values to be used for testing a
%   motion model for that dataset.

% TO-DO:
% 1. Load respSigs, and assign values specifically for each dataset (using switch)
% 2. For respsigAll, check nDynamics == numel(respsigAll)


warning('ReadRespsigs is currently [mostly] a stub');

datasetInfoStruct = GetDatasetInfo(datasetStr);

%nDynamics = datasetInfoStruct.nDynamics;
load(datasetInfoStruct.respSig.fn);



% Assign relevant respiratory values
respSigs.all = [0];
respSigs.training = [0];
respSigs.testing = [0];

end

