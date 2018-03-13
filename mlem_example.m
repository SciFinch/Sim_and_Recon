%% Main execution file for project
addpath('C:\Users\db12\Documents\Sim_and_Recon\branches\quicktest');
% Useful for testing:
%imagesc(squeeze(imageCell{1}(:,110,1:127))'); axis image; colormap;

%% Initialisation
% general
datainfo = 'simple';
interpType = 'linear';
projectorType = 1;
flags.fullSimulation = true; % if false, either CKX dataset, or load sims

% dimensional data
% NOTE: extFOV is for extended field of view - keeps track of all counts during proc'ing of each gate. The +0.25*sizeZ is hard-coded 
% (rather than soft-coding with regular .FOV in each function) to maintain consistency with any changes made later.
pxinfo.pxSize = [344 344 127]; % size of projected image
pxinfo.padSize = [0 0 floor(0.25*pxinfo.pxSize(3))]; % how much to pad image in each direction
pxinfo.pxSizePadded = pxinfo.pxSize + 2*pxinfo.padSize; % size of padded image
pxinfo.pxdims = [2.0445 2.0445 2.03125]; % dimensions of projected image voxels in mm
pxinfo.sino = [344 252 837]; % size of (3D, span-11) sinogram

% simulation specifications
nGates = 1;
scatterMthd = 'conv';
randMthd = 'simple';
totCounts = 75E6;
dynamicWeights = ones(1,nGates)/nGates;

% reconstruction specifications
nIterations = 20;
wait_its = 0;
regionID = 'all';

%% Open Files for Simulation

% Read in base images for simulation
% NOTE TO SELF: - make sure you standardise orientation wrt projtr; arraydim.
%               - play around with data formats to reduce resource usage
[emMap,muMap] = ReadMaps(datainfo,pxinfo);
if strcmp(datainfo,'simple')
    muMap = muMap/10;
    warning('Edit the simple dataset muMaps to be in per mm');
end

%% Simulation
% The following functions will simulate PET data using the dynamic MR volumes and 
% segmented UTE images:

% Generate true sinograms
sim.trues = SimulateTrueDistribution(emMap,projectorType,pxinfo);

% Generate attenuation factor sinograms
sim.AFs = SimulateAttFactors(muMap,projectorType,pxinfo);
warning('ensure that cm conversion is included in all atten procs')

% Generate normalisation factor sinograms
sim.norm = SimulateNormFactors(projectorType,pxinfo,datainfo);
for it = 1:nGates
    temp{it} = sim.norm{1};
end
sim.norm = temp; clear temp;
warning('Need a better way to handle sim.norm for multiple gates')
warning('Need to revist S & R sims - c.f.: dynPET sims')
% Generate randoms sinograms
sim.randoms = SimulateRandomsDistribution(emMap,randMthd,projectorType,pxinfo);

% Generate scatter sinograms
sim.scatters = SimulateScatterDistribution(sim.trues,sim.AFs,sim.randoms,scatterMthd,projectorType,pxinfo);

% Calculate the distribution of prompt counts, on average
for g = 1:nGates
    sim.prompts{g} = sim.norm{g}.*sim.AFs{g}.*sim.trues{g} + sim.scatters{g} + sim.randoms{g};
end

% Simulate the Poisson sampling process
rescaleFlag = 1;
[sim.noisyPrompts] = PoissonSampling(sim.prompts,totCounts,dynamicWeights,rescaleFlag);

% Option for loading/saving rather than generating every time
warning('save/load simulation stub here')

%% Reconstruction

tic
% Initialise for Reconstruction:
dataSinograms = sim.noisyPrompts;
est.AFs = sim.AFs;
est.scatters = sim.scatters;
est.randoms = sim.randoms; 

% - Find normalisation image [can relate to SimulateNormFactors above]
est.norm = sim.norm;

% - initial estimates [including soft-restart from specified point]
% ...
emEst = ones(pxinfo.pxSizePadded);

% - calculate new sensitivity image
SensitivityImage = CalculateSensitivity(est.AFs,projectorType,interpType,pxinfo);
    
% MLEM reconstruction
for itNum = 1:nIterations
    % Project transformed maps
    est.trues = SimulateTrueDistribution(emEst,projectorType,pxinfo);
     
    % Calculate ratio sinograms and generate update images
    ratioSinos = cell(1,nGates);
    for g = 1:nGates
        est.prompts = est.norm{g}.*est.AFs{g}.*est.trues{g} + est.scatters{g} + est.randoms{g};
        ratioSinos{g} = dataSinograms{g}./est.prompts;
    end
    ratioImages = BckProject(ratioSinos,projectorType,interpType,pxinfo);
   
    % Sum over updateImages and apply sensitivity correction
    updateImage = zeros(pxinfo.pxSizePadded);
    for g = 1:nGates, updateImage = updateImage + ratioImages{g}./SensitivityImage{g}; end
    
    % update image variable
    emEst = emEst.*updateImage;
    emEst(isnan(emEst)|isinf(emEst)) = 0;
    imagesc(squeeze(emEst(:,180,:))'); title(itNum); axis image; drawnow;
end
toc
% imagesc(squeeze(est.trues{1}(:,180,:))'); axis image; drawnow;
% imagesc(squeeze(updateImage(:,180,:))'); axis image; drawnow;
% imagesc(squeeze(ratioImages{1}(:,180,:))'); axis image; drawnow;
% imagesc(squeeze(emEst(:,180,:))'); axis image; drawnow;
% imagesc(squeeze(ratioSinos{1}(:,180,:))'); axis image; drawnow;
