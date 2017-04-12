% Attempt to simulate data from a real data reconstruction
% Considerations:
% - how will the noise & low res of the real data be taken into account?
% - can you be sure that the real data is in a stationary position?
% - how will the attenuation effects be taken into account?

%% initialise PET tools
% CONFIGURE PATHS
% - Check what OS I am running on:
if(strcmp(computer(), 'GLNXA64'))
    os = 'linux';
    pathBar = '/';
    sepEnvironment = ':';
elseif(strcmp(computer(), 'PCWIN') || strcmp(computer(), 'PCWIN64'))
    os = 'windows';
    pathBar = '\';
    sepEnvironment = ';';
else
    disp('OS not compatible');
    return;
end

% APIRL PATH
apirlPath = 'C:\Users\db12\apirl_withBinaries_win64\trunk\';
apirlBinariesPath = 'C:\Users\db12\apirl\tags\APIRL1.3.1_win64_vCPU\';
addpath(genpath([apirlPath pathBar 'matlab']));
setenv('PATH', [getenv('PATH') sepEnvironment apirlBinariesPath pathBar 'build' pathBar 'bin']);
setenv('LD_LIBRARY_PATH', [getenv('LD_LIBRARY_PATH') sepEnvironment apirlBinariesPath pathBar 'build' pathBar 'bin']);

PET.scanner = 'mMR';
PET.method =  'otf_siddon_cpu'; %options: {'otf_siddon_cpu','otf_siddon_gpu'}
PET.PSF.type = 'none';
PET.radialBinTrim = 0;
PET.Geom = '';
PET.nSubsets = 1; %by default 21
PET.method_for_randoms = 'from_ML_singles_matlab'; %'from_e7_binary_interfile';
PET.method_for_scatter = 'from_e7_binary_interfile';

PET.sinogram_size.span = 11; % Any span, 0 for multislice 2d, -1 for 2d.
%PET.sinogram_size.span = 1; % Any span, 0 for multislice 2d, -1 for 2d.
PET = classGpet(PET);

% other initialisation stuff

pxsize.FOV = [344 344 127];
pxdim.FOV = [2.0445 2.0445 2.03125];
pxsize.sino = [344 252 837];
roiVec = -7:7;
%% Read in base images 
%load('C:\...')
%baseImg = giplread('C:\Users\db12\Repository_WIN\Collaborative\Christoph Data 310516\Dan Recons\includingMotion\MCIR_mfs.gipl');
baseImg = giplread('C:\Users\db12\Repository_WIN\RealDataAsSim\FDGmap.gipl');
baseImg = baseImg/1000;

% Mu-Maps for attenuation
% note: umaps in /cm
%! Warning: real data attmaps have values < 0 (wtf??), so there's a
%question whether to include these here or elsewhere
%[baseMumap_hardware, refMuMap] = interfileReadSiemensImage('C:\Users\db12\Repository_WIN\Collaborative\Christoph Data 310516\PET_ACQ_51_20160201134327_umap_hardware_00.v.hdr');
%baseMumap_human = interfileReadSiemensImage('C:\Users\db12\Repository_WIN\Collaborative\Christoph Data 310516\PET_ACQ_51_20160201134327_umap_human_00_RPE.v.hdr');
[baseMumap_hardware, refMuMap] = interfileReadSiemensImage('C:\Users\db12\Repository_WIN\Collaborative\Christoph Data 310516\PET_ACQ_51_20160201134327_umap_hardware_00.v.hdr');
baseMumap_human = interfileReadSiemensImage('C:\Users\db12\Repository_WIN\Collaborative\Christoph Data 310516\PET_ACQ_51_20160201134327_umap_human_00_RPE.v.hdr');

% warning('Testing Mu Maps, CURRENTLY ZERO')
% baseMumap_hardware = baseMumap_hardware*0;
% baseMumap_human = baseMumap_human*0;


% Need to mask so that only body visible (no edge-of-field effects)
baseImg(:,:,1:2) = 0;
baseImg(:,:,126:127) = 0;
baseImg(:,1:53,:) = 0;
baseImg(:,285:end,:) = 0;
baseImg(1:125,:,:) = 0;
baseImg(225:end,:,:) = 0;
baseImg(baseImg<eps) = 0;
imagesc(squeeze(baseImg(170,:,:))',[0 0.1]); axis image;

%baseMumap_human(:,:,1:2) = 0;
%baseMumap_human(:,:,126:127) = 0;
baseMumap_human(:,1:53,:) = 0;
baseMumap_human(:,285:end,:) = 0;
baseMumap_human(1:125,:,:) = 0;
baseMumap_human(225:end,:,:) = 0;

baseMumap_human(baseMumap_human<0)=0;

%% Motion data
% load motion information
%load('C:\Users\db12\Repository_WIN\Collaborative\Christoph Data 310516\MotionModel and RespSigs\RespSigs.mat');
%respSigObserved = load('C:\Users\db12\Repository_WIN\Collaborative\Christoph Data 310516\MotionModel and RespSigs\RespSigs.mat');
respSigObserved = load('C:\Users\db12\Repository_WIN\Collaborative\Christoph Data 310516\NavSig.mat');
respSigObserved = respSigObserved.rsig;

load('C:\Users\db12\Repository_WIN\Collaborative\Christoph Data 310516\MotionModel and RespSigs\MotionModels.mat');

pn = 'C:\Users\db12\Repository_WIN\Collaborative\Christoph Data 310516\Motion Fields\';
fn = 'CKdat_gMR_g';
gRef = num2str(1);
pxdim.FOV = PET.image_size.voxelSize_mm;
pxsize.FOV = PET.image_size.matrixSize;
xvec = (1:pxsize.FOV(1))*pxdim.FOV(1);
yvec = (1:pxsize.FOV(2))*pxdim.FOV(2);
zvec = (1:pxsize.FOV(3))*pxdim.FOV(3);
[x,y,z] = ndgrid(xvec,yvec,zvec);


gRef = num2str(1);
for g = 1:8
    gCurr = num2str(g);
    fnCurr = [fn gRef '_to_g' gCurr '_m'];
    % forward transformations
    temp = giplread([pn fnCurr 'x.gipl']);
    temp = padarray(padarray(ipermute(temp,[3 2 1]),[31 31 0],'pre'),[28 28 0],'post');
    xf(:,:,:,g) = temp;
    temp = giplread([pn fnCurr 'y.gipl']);
    temp = padarray(padarray(ipermute(temp,[3 2 1]),[31 31 0],'pre'),[28 28 0],'post');
    yf(:,:,:,g) = temp;
    temp = giplread([pn fnCurr 'z.gipl']);
    temp = padarray(padarray(ipermute(temp,[3 2 1]),[31 31 0],'pre'),[28 28 0],'post');
    zf(:,:,:,g) = temp;
    % inverse (~back) transformations
    temp =  giplread([pn fnCurr 'x_inv_x1000.gipl'])/1000;
    temp = padarray(padarray(ipermute(temp,[3 2 1]),[31 31 0],'pre'),[28 28 0],'post');
    xb(:,:,:,g) = temp;
    temp =  giplread([pn fnCurr 'y_inv_x1000.gipl'])/1000;
    temp = padarray(padarray(ipermute(temp,[3 2 1]),[31 31 0],'pre'),[28 28 0],'post');
    yb(:,:,:,g) = temp;
    temp =  giplread([pn fnCurr 'z_inv_x1000.gipl'])/1000;
    temp = padarray(padarray(ipermute(temp,[3 2 1]),[31 31 0],'pre'),[28 28 0],'post');
    zb(:,:,:,g) = temp;
end


%% Simulation
% note that you can look directly at the e7 recon file dump for some of
% this stuff if you don't want to try and simulate it from scratch.

% sim details
nGates = 8;

% if not scatters/rands:
nCounts = 5E+7;
%nCounts = 10E+10;
% if scatters/rands
% nCounts = 7.5E7;
[c,~] = histcounts(respSigObserved,nGates);
relSignal = c/sum(c); clear c; %proportion of total counts in each gate

fdgMap = zeros([size(baseImg),nGates]);
attMap = zeros([size(baseImg),nGates]);
for g = 1:nGates
    fdgMap(:,:,:,nGates+1-g) = interpn(x,y,z,baseImg,...
        x+zf(:,:,:,g),y+xf(:,:,:,g),z+yf(:,:,:,g),'linear',0);

    % might need some mass preservation here
    attMap(:,:,:,nGates+1-g) = interpn(x,y,z,baseMumap_human,...
        x+zf(:,:,:,g),y+xf(:,:,:,g),z+yf(:,:,:,g),'linear',0);
end
warning('MFs for Sim currently not in default orientation or order')
%attMap(attMap<1E-5) = 0;
%warning('DataSim: Values <1E-5 in the mu-map are being zeroed')
% for g = nGates:-1:1
%     imagesc(squeeze(fdgMap(:,210,:,g))',[0 80]); axis image;
%     drawnow;
%     pause(0.2);
% end

% 
% for g = 1:nGates
%     RespSig = respSigObserved(g);
%     VandeMat = repmat([1 RespSig RespSig^2]',[1 prod(pxsize.FOV)]);
%     
%     fwdMF_x(:,:,:,g) = reshape( sum(VandeMat.*phi.x), pxsize.FOV );
%     fwdMF_y(:,:,:,g) = reshape( sum(VandeMat.*phi.y), pxsize.FOV );
%     fwdMF_z(:,:,:,g) = reshape( sum(VandeMat.*phi.z), pxsize.FOV );
% end
% fdgMap_phi = zeros([size(baseImg),nGates]);
% for g = 1:nGates
%     fdgMap_phi(:,:,:,g) = interpn(x,y,z,baseImg,...
%         x+fwdMF_y(:,:,:,g),y+fwdMF_x(:,:,:,g),z+fwdMF_z(:,:,:,g),'linear',0);
% end
% for g = 1:nGates
%     imagesc(squeeze(fdgMap_phi(:,210,:,g))',[0 80]); axis image;
%     drawnow;
%     pause(0.2);
% end
% for g = 1:nGates
%     subplot(2,1,1); imagesc(squeeze(fdgMap(:,210,:,nGates + 1 - g))',[0 80]); axis image;
%     subplot(2,1,2); imagesc(squeeze(fdgMap_phi(:,210,:,g))',[0 80]); axis image;
%     drawnow;
%     pause(0.2);
% end
% for g = 1:nGates
%     subplot(2,1,1); imagesc(squeeze(fdgMap(:,210,:,g))',[0 80]); axis image;
%     subplot(2,1,2); imagesc(squeeze(fdgMap_phi(:,210,:,g))',[0 80]); axis image;
%     drawnow;
%     pause(0.2);
% end
% phi.x = fwd_model.x; %RL
% phi.y = fwd_model.y; %AP
% phi.z = fwd_model.z; %HF

%%

% apply resolution filter [basic for now]
% more realistic: depends on position, known values, in sinogram space[?]
% kernSize = 5;
% blurFwhm = 1; %in px, not mm
% for g = 1:nGates
%     fdgMap(:,:,:,g) = smooth3(fdgMap(:,:,:,g),'gaussian',kernSize,blurFwhm);
%     attMap(:,:,:,g) = smooth3(attMap(:,:,:,g),'gaussian',kernSize,blurFwhm);
% end
warning('blur disabled')

% generate true sinograms
simTrues = zeros([344 252 837 nGates]);
for g = 1:nGates
    proj = PET.P(fdgMap(:,:,:,g));
    simTrues(:,:,:,g) = proj;
end
clear proj;

% generate attenuation sinograms
simACFs = ones(size(simTrues));
for g = 1:nGates
    % PET.ACF expects the attenuation maps in "per cm" format
    proj = PET.ACF(squeeze(attMap(:,:,:,g))+baseMumap_hardware,refMuMap);
    simACFs(:,:,:,g) = proj;
end
clear proj;
% imagesc(squeeze(simACFs(:,10,1:127,1))'); axis image;
simAFs = 1./simACFs;
simAFs(isinf(simAFs)) = 0;

% generate scatter sinograms
simScatts = zeros(size(simTrues));

% generate random sinograms
simRands = zeros(size(simTrues));

% include normalisation
%ncfs = PET.NCF('C:\Users\db12\Repository_WIN\Collaborative\Christoph Data 310516\norm.n'); % time-invariant.
ncfs = ones(pxsize.sino);%PET.NCF('C:\Users\db12\Repository_WIN\Collaborative\Christoph Data 310516\norm.n'); % time-invariant.
NFs = 1./ncfs;
NFs(isinf(NFs)|isnan(NFs)) = 0;

% generate simulated prompts:
simPrompts = repmat(NFs,[1 1 1 nGates]).*simTrues.*simAFs + simScatts + simRands;
%simPrompts(simPrompts<eps)=0;
%imagesc(squeeze(simPrompts(:,10,1:127))'); axis image;

% add noise to prompt sinograms
for g = 1:nGates
    proj = simPrompts(:,:,:,g);
    currCounts = sum(proj(:));
    scalefactor = nCounts.*relSignal(g);
    proj = (currCounts/scalefactor)*poissrnd(proj*scalefactor/currCounts);
    simPrompts(:,:,:,g) = proj;
end 
clear proj
%clear fdgMap attMap simTrues simAFs simScatts simRands;

% IF storing prompt sinograms:
for g = 1:nGates
    currSino = simPrompts(:,:,:,g);
    comprSimPrompts{g} = sparse(double(currSino(:)));
end
save('C:\Users\db12\Repository_WIN\RealDataAsSim\SimData_wAtt.mat','comprSimPrompts','-v7.3');

%% Recon

projNum = 3; %IMPORTANT

nIterations = 40;
wait_its = 6;


load('C:\Users\db12\Repository_WIN\Collaborative\Christoph Data 310516\MotionModel and RespSigs\RespSigs.mat');

respSigObserved = imgDerivSig;
%respSigObserved = navDerivSig;

% phi.x = fwd_model.x; %RL
% phi.y = fwd_model.y; %AP
% phi.z = fwd_model.z; %HF
% phi_inv.x = bck_model.x;
% phi_inv.y = bck_model.y;
% phi_inv.z = bck_model.z;
warning('Model orientation currently changed from default')
phi.x = fwd_model.y; %RL
phi.y = fwd_model.x; %AP
phi.z = fwd_model.z; %HF
phi_inv.x = bck_model.y;
phi_inv.y = bck_model.x;
phi_inv.z = bck_model.z;


% roi image: [note: might need to define for each dataset uniquely]
roiImage = ones(pxsize.FOV);
roiImage(:,:,1:25) = 0;
roiImage(:,:,95:end) = 0;
roiImage(1:145,:,:) = 0;
roiImage(210:end,:,:) = 0;
roiImage(:,1:110,:) = 0;
roiImage(:,165:end,:) = 0;
% roiSino = PET.P(roiImage);

% Works...ish
% roiImage(:,:,1:3) = 0;
% roiImage(:,:,125:127) = 0;
% roiImage(:,1:75,:) = 0;
% roiImage(:,265:end,:) = 0;
% roiImage(1:125,:,:) = 0;
% roiImage(225:end,:,:) = 0;
%imagesc(squeeze(activityEst(170,:,:))'.*squeeze(roiImage(170,:,:))'); axis image;

% roiImage(:,:,1:3) = 0;
% roiImage(:,:,125:127) = 0;
% roiImage(:,1:120,:) = 0;
% roiImage(:,230:end,:) = 0;
% roiImage(1:125,:,:) = 0;
% roiImage(225:end,:,:) = 0;
%imagesc(squeeze(activityEst(170,:,:))'.*squeeze(roiImage(170,:,:))'); axis image;
%imagesc(squeeze(activityEst(:,170,:))'.*squeeze(roiImage(:,170,:))'); axis image;

 roiSino = PET.P(roiImage);


fovMask = double(PET.ones());
[mx,my,mz] = gradient(fovMask);
edgeMask = round(sqrt(mx.^2 + my.^2 + mz.^2));
%edgeMask(:,:,1) = 1;
%edgeMask(:,:,end) = 1;
 normImage = PET.PT(ones(pxsize.sino));
 
% Starting Estimate for Image
startIt = 1;
if startIt == 1
    activityEst = ones(pxsize.FOV);
    respsigEst = repmat(min(respSigObserved)+0.5*range(respSigObserved),[1 nGates]);
    
    respsigRec = zeros([nGates nIterations])*nan;
    activityRec = zeros([pxsize.FOV nIterations]);
    
    motionEstimationControlParam = 0;
    
else
    activityEst = activityRec(:,:,:,startIt-1);
    respsigEst = respsigRec(:,startIt-1)';
end

%activityRec = cell(1,3);
promptSino = simPrompts;
sliceNum_AP = 170;

%figure;
RespSigUpdate = zeros([1,nGates]);
for k = startIt:nIterations
    tic
    disp(['Iteration #' num2str(k) ' beginning']);
    % Estimate current motion
    if k > wait_its || k == 1
        if k > 1
            oldFwdMF_x = fwdMF_x;
            oldFwdMF_y = fwdMF_y;
            oldFwdMF_z = fwdMF_z;
        end
        
        parfor g = 1:nGates
            RespSig = respsigEst(g);
            VandeMat = repmat([1 RespSig RespSig^2]',[1 prod(pxsize.FOV)]);
            
            fwdMF_x(:,:,:,g) = reshape( sum(VandeMat.*phi.x), pxsize.FOV );
            fwdMF_y(:,:,:,g) = reshape( sum(VandeMat.*phi.y), pxsize.FOV );
            fwdMF_z(:,:,:,g) = reshape( sum(VandeMat.*phi.z), pxsize.FOV );
            
            bckMF_x(:,:,:,g) = reshape( sum(VandeMat.*phi_inv.x), pxsize.FOV );
            bckMF_y(:,:,:,g) = reshape( sum(VandeMat.*phi_inv.y), pxsize.FOV );
            bckMF_z(:,:,:,g) = reshape( sum(VandeMat.*phi_inv.z), pxsize.FOV );
        end
        %  clear RespSig VandeMat g;
        
        % Find current attenuation position estimates
        kernSize_recon = 5;
        blurFwhm_recon = 1; %in px, not mm
        attMap = zeros([pxsize.FOV,nGates]);
        attSino = ones([pxsize.sino,nGates]);
        for g = 1:nGates
            temp = interpn(x,y,z,baseMumap_human,...
                x+fwdMF_x(:,:,:,g),y+fwdMF_y(:,:,:,g),z+fwdMF_z(:,:,:,g),'linear');
            temp = temp + baseMumap_hardware;
          %  temp = smooth3(temp,'gaussian',kernSize_recon,blurFwhm_recon);
            temp(temp<1E-5) = 0;
            attMap(:,:,:,g) = temp;
            
            proj = PET.ACF(squeeze(attMap(:,:,:,g)),refMuMap);
            attSino(:,:,:,g) = NFs.*(1./proj);
        end
        attSino(isinf(attSino))=0;
        subplot(5,4,1); imagesc(squeeze(attSino(:,10,1:127,nGates))'); axis image; title('AF Sino'); drawnow;
                
        % Find correct sensitivity, including motion
        % S = sum_g M_g' A' X_g' N 1  == sum_g M_g' A' X_g N 1
        SensitivityImg = zeros(pxsize.FOV);
        for g = 1:nGates
            % Can remove the following if you are convinced it is
            % cancelling properly when forming update image:
 %           normImage = PET.PT(attSino(:,:,:,g));  % THIS DIDN'T TAKE 1/10 FACTOR INTO ACCOUNT!
            
            normImage = PET.Sensitivity(attSino(:,:,:,g));
                       
            % Calculate motion contribution to sensitivity:
            motionNormImage = interpn(x,y,z,normImage,...
                x+bckMF_x(:,:,:,g),y+bckMF_y(:,:,:,g),z+bckMF_z(:,:,:,g),'linear');
            SensitivityImg = SensitivityImg + motionNormImage;
        end
        SensitivityImg = SensitivityImg.*fovMask;
        subplot(5,4,2); imagesc(squeeze(SensitivityImg(sliceNum_AP,:,:))'); axis image; title('Sensitivity'); drawnow;
     end
    %  clear UnitSino bckprojUnitSino motionUnitSino;
    
    % Estimate current images
    updateImage = zeros(pxsize.FOV);
    for g = 1:nGates
        currGateData = promptSino(:,:,:,g);
        subplot(5,4,3); imagesc(squeeze(currGateData(:,10,1:127))'); axis image; title('Curr. Data'); drawnow;
        
        fwdImage = interpn(x,y,z,activityEst,...
            x+fwdMF_x(:,:,:,g),y+fwdMF_y(:,:,:,g),z+fwdMF_z(:,:,:,g),'linear');
        %  fwdImage(fwdImage<0)=0;
        fwdImage(isnan(fwdImage)|isinf(fwdImage))=0;
        fwdImage(:,:,1:3) = 0; fwdImage(:,:,(end-2):end) = 0;
        fwdImage = fwdImage.*fovMask;
        subplot(5,4,4); imagesc(squeeze(fwdImage(sliceNum_AP,:,:))'); axis image; title('Fwd Image (imgup)'); drawnow;

        %  meanSino = meanSino.*attEstimate(:,:,:,g) + scatterSino(:,:,:,g) + randSino;
        meanSino = PET.P(fwdImage);
        meanSino = meanSino.*attSino(:,:,:,g);
        %  meanSino(meanSino<eps & meanSino>0) = eps;
        subplot(5,4,5); imagesc(squeeze(meanSino(:,10,1:127))'); axis image; title('Data Model (imgup)'); drawnow;
        
        ratioSino = attSino(:,:,:,g).*single(currGateData./(meanSino));
        ratioSino(isinf(ratioSino)) = 0;
        ratioSino(isnan(ratioSino)) = 0;
        %    ratioSino(ratioSino<0)=0;
        %     clear currGateData frame fwdImage meanSino;
        subplot(5,4,6); imagesc(squeeze(ratioSino(:,10,1:127))'); axis image; title('Ratio Sino (imgup)'); drawnow;
        
        % Can remove the attSino here if you're sure it's being cancelled
        % out in the sensitivity image (remove that one too if so):
        bckRatioSino = PET.PT(ratioSino);
        subplot(5,4,7); imagesc(squeeze(fovMask(sliceNum_AP,:,:).*bckRatioSino(sliceNum_AP,:,:))'); axis image; title('Ratio Image'); drawnow;
        
        % this is to prevent weird out of FOV intensities being dragged
        % into the accepted FOV:
        bckRatioSino = bckRatioSino.*fovMask; 
        
        % Move current contribution to update image into ref position
        updateImage = updateImage + interpn(x,y,z,bckRatioSino,x+bckMF_x(:,:,:,g),y+bckMF_y(:,:,:,g),z+bckMF_z(:,:,:,g),'linear');
        subplot(5,4,8); imagesc(squeeze(fovMask(sliceNum_AP,:,:).*updateImage(sliceNum_AP,:,:)./SensitivityImg(sliceNum_AP,:,:))'); axis image; title(['Update (' num2str(g) ' gates)']); drawnow;
        %     clear bckRatioSino;
    end
    activityEst = activityEst.*(updateImage./(SensitivityImg));
    
    % Is the following really necessary? Consider <0 = 0
    activityEst(isnan(activityEst)|isinf(activityEst))=0;
%    activityEst = activityEst.*fovMask;
 %   subplot(5,4,9); imagesc(squeeze(activityEst(170,:,1:127))'); axis image;
    %  activityEst(:,:,1:2) = 0; activityEst(:,:,(end-2):end) = 0;
    activityEst(:,:,1) = 0; activityEst(:,:,end) = 0;
    % activityEst(activityEst>200) = 200;
    %  clear updateImage;
    subplot(5,4,9); imagesc(squeeze(activityEst(sliceNum_AP,:,:))',[0 quantile(activityEst(:),0.999)]); axis image; title(['Image Estimate (' num2str(k) ')']);
    
    % Estimate respiratory position
    if k > wait_its
        for g = 1:nGates
            currGateData = promptSino(:,:,:,g);
            currMuMap = attMap(:,:,:,g);
 %           currMuMap = currMuMap - baseMumap_hardware;
            currMuMap = currMuMap/10; %needs to me in per mm
            
            fwdImage = interpn(x,y,z,activityEst,...
                x+fwdMF_x(:,:,:,g),y+fwdMF_y(:,:,:,g),z+fwdMF_z(:,:,:,g),'linear');
            %   fwdImage(fwdImage<0)=0;
            fwdImage(isnan(fwdImage)|isinf(fwdImage))=0;
%            fwdImage(:,:,1:3) = 0; fwdImage(:,:,(end-2):end) = 0;
   %         fwdImage = fwdImage.*roiImage;
%            fwdImage(fwdImage<=1E-5)=0;
%            fwdImage(fwdImage>10)=10;
            
            subplot(5,4,10); imagesc(squeeze(fwdImage(sliceNum_AP,:,:))'); axis image; title('Fwd Img (respup)');  drawnow;
            
            currMuMap(isnan(currMuMap)|isinf(currMuMap))=0;
 %           currMuMap(:,:,1:3) = 0; currMuMap(:,:,(end-2):end) = 0;
 %           currMuMap = currMuMap.*roiImage;
            
            %  meanSino = meanSino.*attEstimate(:,:,:,g) + scatterSino(:,:,:,g) + randSino;
            meanSino = PET.P(fwdImage);
            meanSino = meanSino.*attSino(:,:,:,g);

            ratioSino = single(currGateData./(meanSino));
            ratioSino(isnan(ratioSino)|isinf(ratioSino))=0;
        %    clear currGateData;
            subplot(5,4,11); imagesc(squeeze(ratioSino(:,10,1:127))'); axis image; title('Ratio Sino (respup)');  drawnow;
            
            [gradFwdImageX,gradFwdImageY,gradFwdImageZ] = gradient(fwdImage,pxdim.FOV(1),pxdim.FOV(2),pxdim.FOV(3));
            subplot(5,4,12); imagesc(squeeze(gradFwdImageZ(sliceNum_AP,:,:))'); axis image; title('Z-gradient (emission)');  drawnow;
            
            [gradFwdMuX,gradFwdMuY,gradFwdMuZ] = gradient(currMuMap,pxdim.FOV(1),pxdim.FOV(2),pxdim.FOV(3));
            subplot(5,4,13); imagesc(squeeze(gradFwdMuZ(sliceNum_AP,:,:))'); axis image; title('Z-gradient (attenuation)');  drawnow;

            RespSig = respsigEst(g);
            VandeMatrix = repmat([0 1 2*RespSig]',[ 1 prod(pxsize.FOV) ]);
            derivMF_x = reshape(sum(VandeMatrix.*phi.x),pxsize.FOV);
            derivMF_y = reshape(sum(VandeMatrix.*phi.y),pxsize.FOV);
            derivMF_z = reshape(sum(VandeMatrix.*phi.z),pxsize.FOV);
   %         clear VandeMatrix fwdImage;
            
            gradDotImage = gradFwdImageX.*derivMF_x ...
                + gradFwdImageY.*derivMF_y ...
                + gradFwdImageZ.*derivMF_z;
            gradDotImage(isnan(gradDotImage)|isinf(gradDotImage))=0;
   %         clear gradFwdImage*;
            gradDotImage = gradDotImage.*roiImage;   
            subplot(5,4,14); imagesc(squeeze(gradDotImage(sliceNum_AP,:,:))'); axis image; title('Weighted Grad (emiss)'); drawnow;
            
            gradDotMu = gradFwdMuX.*derivMF_x ...
                + gradFwdMuY.*derivMF_y ...
                + gradFwdMuZ.*derivMF_z;
            gradDotMu(isnan(gradDotMu)|isinf(gradDotMu))=0;
            gradDotMu = gradDotMu.*roiImage;
   %         clear gradFwdImage*;
            subplot(5,4,15); imagesc(squeeze(gradDotMu(sliceNum_AP,:,:))'); axis image; title('Weighted Grad (atten)'); drawnow;
            
     %       gradDotImage = gradDotImage.*roiImage;
            fwdGradDotImage = PET.P(gradDotImage);
            subplot(5,4,16); imagesc(squeeze(fwdGradDotImage(:,10,1:127))'); axis image; title('Weight Sino (emiss)'); drawnow;
            
     %       gradDotMu = gradDotMu.*roiImage;
            fwdGradDotMu = PET.P(gradDotMu);
            subplot(5,4,17); imagesc(squeeze(fwdGradDotMu(:,10,1:127))'); axis image; title('Weight Sino (atten)'); drawnow;

            combWeights = attSino(:,:,:,g).*fwdGradDotImage - meanSino.*fwdGradDotMu;
            subplot(5,4,18); imagesc(squeeze(combWeights(:,10,1:127))'); axis image; title('Combined weights'); drawnow;
            
            RespSigUpdateSino = (ratioSino-1).*(combWeights);
            RespSigUpdateSino(isnan(RespSigUpdateSino)|isinf(RespSigUpdateSino)) = 0;
            subplot(5,4,19); imagesc(squeeze(RespSigUpdateSino(:,10,1:127))'); axis image; title('RespSig update'); drawnow;
%            RespSigUpdate(g) =  sum( roiSino(:).*RespSigUpdateSino(:) );
            RespSigUpdate(g) =  sum( RespSigUpdateSino(:) );
            %  RespSigUpdate(g) =  sum( sum ( sum( RespSigUpdateSino(~isnan(roiSino)))));
  %          clear ratioSino fwdGradDotImage frame;
        end
        
        %Check the maximum change displacements to try and identify when
        %convergence is below the tolerance
        if k == (wait_its + 1)
            motionEstimationControlParam = repmat(0.25/mean(abs(RespSigUpdate)),[1 nGates]);
        end
        
        if k > 1 && k > wait_its
            absDispDiff = sqrt((oldFwdMF_x - fwdMF_x).^2 +...
                (oldFwdMF_y - fwdMF_y).^2 +  (oldFwdMF_z - fwdMF_z).^2);
            
            maxDispDiff = zeros(1,nGates);
            for g = 1:nGates
                maxDispDiff(g) = max(max(max(absDispDiff(:,:,:,g))));
                if maxDispDiff(g) < 0.5*pxdim.FOV(3)
                    motionEstimationControlParam(g) = motionEstimationControlParam(g)*1;%motionEstimationControlParam(g)*0.25;
                end
            end
        end
        
    end
      
    % Updatmotion estimateep:
    
    %stepSize = (0.93)^(k-wait_its+1)*motionEstimationControlParam; % @@@
    %stepSize = 3.2*10^(-5);%motionEstimationControlParam; % @@@
    
    stepSize = motionEstimationControlParam; % @@@
    respsigEst = respsigEst + stepSize.*RespSigUpdate;
    
    respsigEst(respsigEst>=(max(respSigObserved)+0.2*range(respSigObserved))) = (max(respSigObserved)+0.2*range(respSigObserved));
    respsigEst(respsigEst<=(min(respSigObserved)-0.2*range(respSigObserved))) = (min(respSigObserved)-0.2*range(respSigObserved));
    clear stepSize
    subplot(5,4,20); plot(1:k,respsigRec(:,1:k)); title('RespSig update'); drawnow; ylim([min(respSigObserved)-0.2*range(respSigObserved), max(respSigObserved)+0.2*range(respSigObserved)]); xlim([1 nIterations])    ;

    % Store step:
    respsigRec(:,k) = respsigEst(:);
    activityRec(:,:,:,k) = squeeze(activityEst);      
%     % Display step    :
%     sliceNum = 170    ;
%     toDisplay = squeeze(activityEst(sliceNum,:,:))'    ;
%     subplot(1,2,1); plot(1:k,respsigRec(:,1:k)); ylim([min(respSigObserved)-0.2*range(respSigObserved), max(respSigObserved)+0.2*range(respSigObserved)])    ;
%     xlabel('Iteration Number, k'); ylabel('RespSig Estimate (mm)')    ;
%     legend    ;
%     subplot(1,2,2); imagesc(toDisplay,[0 40])    ;
%     title(['Iteration ' num2str(k) '/' num2str(nIterations)])    ;
%     axis image    ;
%     drawnow    ;
%     toc
    
    % if ~mod(k,5)
    saveas(gcf,['C:\Users\db12\Dropbox\ResultWatch\' 'SIM_emandatt_recon_realdata_' num2str(nGates) 'G_k' num2str(k) '.png']);
    disp(['Iteration ' num2str(k) ' complete']);
    %  end
end% % 
% %imagesc(squeeze(test(:,:,70))',[0 40]);% 
% figure;
% for it = 1:40
%     sliceNum = 170;
%     toDisplay = squeeze(activityRec(sliceNum,:,:,it))';
%     imagesc(toDisplay,[0 50]);
%     title(['Iteration ' num2str(it) '/' num2str(nIterations)]);
%     axis image;
%     drawnow;
%     pause(0.5);
% end% % 
% for sliceNum = 120:250
%     %   sliceNum = 170;
%     toDisplay = squeeze(fwdImage(sliceNum,:,:))';
%     imagesc(toDisplay);
%     title(['slice ' num2str(sliceNum)]);
%     axis image;
%     drawnow;
%     pause(0.1);
% end
% NORMALISATION CURRENTLY OFF

%anf = acfs .* ncfs;
%anf(anf~=0) = 1./anf(anf~=0);
sensImage = PET.Sensitivity(mean(repmat(NFs,[1 1 1 nGates]).*simAFs,4));
%additive = 0*(randoms + scatter).*ncfs.*acfs;
additive = zeros(pxsize.sino);
recon = PET.ones();
%recon = PET.OPOSEM(sinogram,additive,sensImage,recon, 1);
recon = PET.OPOSEM(mean(simPrompts,4),additive,sensImage,recon,20);

figure; imagesc(squeeze(recon(170,:,:))'); axis image;