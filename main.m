%% Main execution file for project
addpath('C:\Users\db12\Sim_and_Recon\branches\TestProjectors');
% Useful for testing:
%imagesc(squeeze(imageCell{1}(:,110,1:127))'); axis image; colormap;

%% Initialisation
% general
datainfo = 'simple';
interpType = 'linear';
projectorType = 2;
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
nGates = 4;
scatterMthd = 'conv';
randMthd = 'none';
totCounts = 75E6;
dynamicWeights = ones(1,nGates)/nGates;

% reconstruction specifications
%nIterations = 1;
%wait_its = 0;
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

% Read in the motion fields used to simulate motion positions
[dHF,dAP,dRL] = ReadMotionFields(datainfo,'fwd',pxinfo);
[dHFinv,dAPinv,dRLinv] = ReadMotionFields(datainfo,'inv',pxinfo);

%% Simulation
% The following functions will simulate PET data using the dynamic MR volumes and 
% segmented UTE images:

% Generate motion transformed distribution maps
[transEmMap,transMuMap] = TransformMaps(emMap,muMap,pxinfo,dHF,dAP,dRL,interpType,0);
%clear emMap muMap dHF dAP dRL dHFinv dAPinv dRLinv;

% Generate true sinograms
sim.trues = SimulateTrueDistribution(transEmMap,projectorType,pxinfo);

% Generate attenuation factor sinograms
sim.AFs = SimulateAttFactors(transMuMap,projectorType,pxinfo);
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
sim.randoms = SimulateRandomsDistribution(transEmMap,randMthd,projectorType,pxinfo);

% Generate scatter sinograms
sim.scatters = SimulateScatterDistribution(sim.trues,sim.AFs,sim.randoms,scatterMthd,projectorType,pxinfo);

% Calculate the distribution of prompt counts, on average
for g = 1:nGates
    sim.prompts{g} = sim.norm{g}.*sim.AFs{g}.*sim.trues{g} + sim.scatters{g} + sim.randoms{g};
end

% Simulate the Poisson sampling process
[sim.noisyPrompts] = PoissonSampling(sim.prompts,totCounts,dynamicWeights,0);

% Option for loading/saving rather than generating every time
warning('save/load simulation stub here')

%% Open Files for Reconstruction

% Load motion models
motionModel = ReadMotionModels(datainfo);

% Load respiratory signals
respSigs = ReadRespsigs(datainfo);

%% Reconstruction

% Initialise for Reconstruction:
dataSinograms = sim.noisyPrompts;

% - Generate ROI image
roiImage = GenerateROI(datainfo,regionID,pxinfo); %TODO
roiSino = FwdProject(roiImage,projectorType,pxinfo);

% - Find normalisation image [can relate to SimulateNormFactors above]
% est.norm = SimulateNormFactors(projectorType,pxinfo,datainfo);
% for it = 1:nGates
%     temp{it} = est.norm{1};
% end
% est.norm = temp; clear temp;
est.norm = sim.norm;

%
est.randoms = sim.randoms;
est.scatters = sim.scatters;

% - initial estimates [including soft-restart from specified point]
% ...
muEst = ToggleImagePadding(muMap,pxinfo);

figure(1);
nIterations = 70;
wait_its = 6;

startIt = 1;
if startIt == 1
    % estimates:
    emEst = ones(pxinfo.pxSizePadded);
    respSigEst = zeros(1,nGates) + min(respSigs.all) + 0.5*range(respSigs.all);
    % recording:
    emEstRec = zeros([pxinfo.pxSizePadded nGates]);
    respSigRec = nan*zeros(nGates,nIterations);
else
    % set estimates so that recon can begin from startIt:
    emEst  = emEstRec(:,:,:,startIt-1);
    respSigEst = respSigRec(:,startIt-1)';
end


% CHECK ALL OF THE FOLLOWING AGAINST BEST-WORKING MEMCIR
% Don't forget to include extended field of view stuff
for itNum = startIt:nIterations
    
    % IMAGE UPDATE:
    
    % Estimate motion fields
    [dHF,dAP,dRL] = EstimateFieldsFromModel(respSigEst,motionModel,'fwd',pxinfo);
    [dHFinv,dAPinv,dRLinv] = EstimateFieldsFromModel(respSigEst,motionModel,'inv',pxinfo);
    
    % Transform emission and attenuation map estimates
    transEmEst = TransformImage(emEst,pxinfo,dHF,dAP,dRL,interpType,0);
    transMuEst = TransformImage(muEst,pxinfo,dHF,dAP,dRL,interpType,0);
    
%         figure(4);
%         for t = 1:nGates
%             imagesc(squeeze(est.trues{t}(:,110,1:127))'); title(t); axis image;
%           %  imagesc(squeeze(transEmEst{t}(:,177,:))'); title(t); axis image;
%             drawnow;
%             pause(0.5);
%         end
%     
    % Project transformed maps
    est.AFs = SimulateAttFactors(transMuEst,projectorType,pxinfo);
    est.trues = SimulateTrueDistribution(transEmEst,projectorType,pxinfo);
    
    % Calculate new sensitivity image
    sensitivityImage = CalculateSensitivity(est.AFs,dHFinv,dAPinv,dRLinv,projectorType,interpType,pxinfo);
    sensitivityImage(isinf(sensitivityImage)|isnan(sensitivityImage)) = 0;

    % Calculate ratio sinograms and generate update images
    ratioSinos = cell(1,nGates);
    for g = 1:nGates
        est.prompts{g} = est.norm{g}.*est.AFs{g}.*est.trues{g} + est.scatters{g} + est.randoms{g};
        ratioSinos{g} = dataSinograms{g}./(est.prompts{g} + eps);
        ratioSinos{g}(isnan(ratioSinos{g})|isinf(ratioSinos{g})) = 0;
    end
    ratioImages = BckProject(ratioSinos,projectorType,interpType,pxinfo);
    updateImages = TransformImage(ratioImages,pxinfo,dHFinv,dAPinv,dRLinv,interpType,0);
    
    % Sum over updateImages and apply sensitivity correction
    temp = zeros(pxinfo.pxSizePadded);
    for g = 1:nGates, temp = temp + updateImages{g}; end
    updateImage = temp./sensitivityImage; clear temp;
    
    % update image variable
    emEst = emEst.*updateImage;
    emEst(isinf(emEst)|isnan(emEst)) = 0;
    
    emEstRec(:,:,:,itNum) = emEst;
    
    subplot(121);
    imagesc(squeeze(emEst(:,175,:))'); title(itNum);
    axis image;
    drawnow;
    
    % RESPSIG UPDATE
    if itNum > wait_its
        
        transEmEst = TransformImage(emEst,pxinfo,dHF,dAP,dRL,interpType,0);
        est.trues = SimulateTrueDistribution(transEmEst,projectorType,pxinfo);
        ratioSinos = cell(1,nGates);
        for g = 1:nGates
            est.prompts{g} = est.norm{g}.*est.AFs{g}.*est.trues{g} + est.scatters{g} + est.randoms{g};
            ratioSinos{g} = dataSinograms{g}./est.prompts{g};
        end
        
        % Find the motion-model-derivative--weighted sinograms
        [emWeightSinos,muWeightSinos] = CalcWeightedSinos(transEmEst,transMuEst,motionModel,respSigEst,pxinfo,projectorType);
        
        % Combined weights for update sinograms [this will be data-model--dependent]
        for g = 1:nGates
            totWeightSino{g} = est.AFs{g}.*emWeightSinos{g} - est.prompts{g}.*muWeightSinos{g};
        end
        
        % Calculate the respsig additive-update sinograms
        RespSigUpdate = zeros(1,nGates);
        for g = 1:nGates
            RespSigUpdateSinos{g} = (ratioSinos{g}-1).*totWeightSino{g};
            RespSigUpdateSinos{g}(isnan(RespSigUpdateSinos{g})|isinf(RespSigUpdateSinos{g})) = 0;
            
            RespSigUpdate(g) = sum(RespSigUpdateSinos{g}(:).*roiSino{1}(:));
        end
        
    end
    if itNum <= wait_its
        stepsize = 0;
        RespSigUpdate = zeros(size(respSigEst));
    else
        if itNum == wait_its+1
            stepsize = 1/sum(abs(RespSigUpdate));
        end
    end
    
    % update respiratory signal estimates
    respSigEst = respSigEst + stepsize*RespSigUpdate;
    
    respSigRec(:,itNum) = respSigEst;
    
    subplot(122);
    plot(respSigRec'); xlim([1 nIterations]); 
    ylim([min(respSigs.all) - 0.2*range(respSigs.all), min(respSigs.all) + 1.2*range(respSigs.all)]); 
    drawnow;
end

% %
% 
% figure;
% for itNum = 1:nIterations
%     subplot(121);
%     imagesc(squeeze(emEstRec(:,175,:,itNum))'); title(itNum);
%     axis image;
%     drawnow;
%     subplot(122);
%     plot(respSigRec'); xlim([1 nIterations]);
%     ylim([min(respSigs.all) - 0.2*range(respSigs.all), min(respSigs.all) + 1.2*range(respSigs.all)]);
%     drawnow;
%     pause(0.1);
% end