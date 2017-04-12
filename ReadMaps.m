function [emMap,muMap] = ReadMaps(datasetStr,pixelInfo)
%READMAPS [emMap,muMap] = ReadMaps(datasetStr,pixelInfo);
%   Given a dataset identifier, datasetStr, this function will load up the
%   image volumes containing the distributions of tracer and attenuation
%   coefficient for emMap and muMap respectively. These are then designed
%   to be used for PET simulation. This way, the outputs for different
%   datasets can be standardised with verbose notes here, rather than
%   cluttering the master script they are used for.

% Possible to-dos:
% 1. include padding for specific datasets (could use ToggleImagePadding)
% 2. Better zeroing
% 3. Ensuring orientations are standardised
% 4. Checking attmaps are in per millimetre

datasetInfoStruct = GetDatasetInfo(datasetStr);

nDynamics = datasetInfoStruct.nDynamics;

% Find the relevant filenames to load
emission_fn = datasetInfoStruct.maps.emission.fn;
mumap_human_fn = datasetInfoStruct.maps.attenuation.human;
mumap_hardware_fn = datasetInfoStruct.maps.attenuation.hardware;

if isempty(mumap_hardware_fn)
	flags.muHardware = 0;
else
	flags.muHardware = 1;
end

% Load motion fields
[emMap,~,~] = AutoLoadImage(emission_fn);
[muMap_human,~,~] = AutoLoadImage(mumap_human_fn);

if flags.muHardware
	[muMap_hardware,~,~] = AutoLoadImage(mumap_hardware_fn);
else
	muMap_hardware = 0*muMap_human;
end

muMap = muMap_human + muMap_hardware;

% Some maps need neatening prior to use in simulation
switch datasetStr
case 'CK1'
	%baseMumap_human(:,:,1:2) = 0;
	%baseMumap_human(:,:,126:127) = 0;
	muMap(:,1:53,:) = 0;
	muMap(:,285:end,:) = 0;
	muMap(1:125,:,:) = 0;
	muMap(225:end,:,:) = 0;
	emMap(:,1:53,:) = 0;
	emMap(:,285:end,:) = 0;
	emMap(1:125,:,:) = 0;
	emMap(225:end,:,:) = 0;
end

% Assumng orientation [AP RL HF]
if ~(sum(sum(emMap(1,:,:))) + sum(sum(emMap(end,:,:))) + sum(sum(emMap(:,1,:))) + sum(sum(emMap(:,end,:))))
	warning('Edges of transaxial plane are non-zero');
end

end

