% Attempt to structure this code properly

%% Initialisation
% general
datainfo = 'stub';
interpType = 'linear';
projectorType = 3;
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

% reconstruction specifications
nIterations = 1;
wait_its = 0;
regionID = 'all';

%% Open Files for Simulation

% Read in base images for simulation
% NOTE TO SELF: - make sure you standardise orientation wrt projtr; arraydim.
%               - play around with data formats to reduce resource usage
[emMap,muMap] = ReadMaps(datainfo,pxinfo);

% Read in the motion fields used to simulate motion positions
[dHF,dAP,dRL] = ReadMotionFields(datainfo,'fwd',pxinfo);
[dHFinv,dAPinv,dRLinv] = ReadMotionFields(datainfo,'inv',pxinfo);

%% Simulation
% The following functions will simulate PET data using the dynamic MR volumes and 
% segmented UTE images:

% Generate motion transformed distribution maps
[transEmMap,transMuMap] = TransformMaps(emMap,muMap,pxinfo,dHF,dAP,dRL,interpType,0);

% Generate true sinograms
sim.trues = SimulateTrueDistribution(transEmMap,projectorType,pxinfo);

% Generate attenuation factor sinograms
sim.AFs = SimulateAttFactors(transEmMap,projectorType,pxinfo);
warning('ensure that cm conversion is included in all atten procs')

% Generate normalisation factor sinograms
sim.norm = SimulateNormFactors(projectorType,pxinfo,datainfo);

% Generate scatter sinograms
sim.scatters = SimulateScatterDistribution(transEmMap,scatterMthd,projectorType,pxinfo);

% Generate randoms sinograms
sim.randoms = SimulateRandomsDistribution(transEmMap,randMthd,projectorType,pxinfo);

% Calculate the distribution of prompt counts, on average
for g = 1:nGates
    sim.prompts{g} = sim.norm{g}.*sim.AFs{g}.*sim.trues{g} + sim.scatters{g} + sim.randoms{g};
end

% Simulate the Poisson sampling process
[sim.noisyPrompts] = PoissonSampling(sim.prompts,totCounts,dynamicWeights);

% Option for loading/saving rather than generating every time
warning('save/load simulation stub here')

%% Open Files for Reconstruction

% Load motion models
[phi.HF,phi.AP,phi.RL] = ReadMotionModels(datainfo,'fwd',pixelInfo);
[phi.HFinv,phi.APinv,phi.RLinv] = ReadMotionModels(datainfo,'inv',pixelInfo);

% Load respiratory signals
[respSigs] = ReadRespsigs(datainfo);

%% Reconstruction

% Initialise for Reconstruction:
dataSinograms = sim.noisyPrompts;

% - Generate ROI image
roiImage = GenerateROI(datainfo,regionID,pixelInfo);

% - Find normalisation image [can relate to SimulateNormFactors above]
% ...

% - initial estimates [including soft-restart from specified point]
% ...
emEst = zeros(pxinfo.pxSizePadded);
respSigEst = zeros(1,nGates) + min(respSigs.all) + 0.5*range(respSigs.all);

% CHECK ALL OF THE FOLLOWING AGAINST BEST-WORKING MEMCIR
% Don't forget to include extended field of view stuff
for itNum = 1:nIterations
    
    % IMAGE UPDATE:
    
    % Estimate motion fields
    [dHF,dAP,dRL] = EstimateFieldsFromModel(respSigEst,phi,pxinfo);
    [dHFinv,dAPinv,dRLinv] = EstimateFieldsFromModel(respSigEst,phi_inv,pxinfo);
    
    % Transform emission and attenuation map estimates
    [transEmMap,transMuMap] = TransformMaps(emMap,muMap,pxinfo,dHF,dAP,dRL,interpType,0);
    
    % Project transformed maps
    est.AFs = SimulateAttFactors(transMuMap,projectorType,pxinfo);
    est.trues = SimulateTrueDistribution(transEmMap,projectorType,pxinfo);
    
    % Calculate new sensitivity image
    SensitivityImage = CalculateSensitivity(est.AFs,dHFinv,dAPinv,dRLinv,projectorType,pxinfo);
    
    % Calculate ratio sinograms and generate update images
    ratioSinos = cell(1,nGates);
    for g = 1:nGates
        est.prompts{g} = est.norm{g}.*est.AFs{g}.*est.trues{g} + est.scatters{g} + est.randoms{g};
        ratioSinos{g} = dataSinograms{g}./prompts{g};
    end
    ratioImages = BckProject(ratioSinos,projectorType,pxinfo);
    updateImages = TransformImage(ratioImages,pxinfo,dHFinv,dAPinv,dRLinv,interpType,0);
    
    % Sum over updateImages and apply sensitivity correction
    temp = zeros(pxinfo.pxSizePadded);
    for g = 1:nGates, temp = temp + updateImages{g}; end
    updateImage = temp./SensitivityImage; clear temp;
    
    % update image variable
    emEst = emEst.*updateImg;
    
    % RESPSIG UPDATE

    transEmMap = TransformImage(emMap,pxinfo,dHF,dAP,dRL,interpType,0);
    est.trues = SimulateTrueDistribution(transEmMap,projectorType);
    ratioSinos = cell(1,nGates);
    for g = 1:nGates
        est.prompts{g} = est.norm{g}.*est.AFs{g}.*est.trues{g} + est.scatters{g} + est.randoms{g};
        ratioSinos{g} = dataSinograms{g}./prompts{g};
    end
    
    % Find the motion-model-derivative--weighted sinograms
    [emWeightSinos,muWeightSinos] = CalcWeightedSinos(transEmMap,transMuMap,phi,respSigEst,pxinfo);
    
    % Combined weights for update sinograms [this will be data-model--dependent]
    for g = 1:ngates
        totWeightSino{g} = est.AFs{g}.*emWeightSinos{g} - est.prompts{g}.*muWeightSinos{g};
    end

    % Calculate the respsig additive-update sinograms
    RespSigUpdate = zeros(1,nGates);
    for g = 1:nGates
        RespSigUpdateSinos{g} = (ratioSinos{g}-1).*totWeightSino{g};
        RespSigUpdateSinos{g}(isnan(RespSigUpdateSinos{g})|isinf(RespSigUpdateSinos{g})) = 0;

        RespSigUpdate(g) = sum(RespSigUpdateSinos{g}(:).*roiSino(:));
    end

    if k <= wait_its
        stepsize = 0;
    else
        stepsize = 0;
    end
    
    % update respiratory signal estimates
    respSigEst = respSigEst + stepsize*RespSigUpdate;    
    
end

