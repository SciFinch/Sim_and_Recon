%% This code is used to generate the simple data used to test this project
addpath('C:\Users\db12\Dropbox\Work_Docs\MATLAB Scripts');
addpath('C:\Users\db12\Sim_and_Recon\branches\TestProjectors');
simpleData_basePath = 'C:\Users\db12\Dropbox\Work_Docs\SimpleData\';
datainfo = GetDatasetInfo('simple');

%% User definitions:
phantom_type = 'sphere';
%phantom_type = 'nema';
nGates = 4;
flags.visualise = true;

% Respiratory signal values: See MM definition part
warning('Correct orientation of AP and RL not checked yet')

%% Construct phantom
% dimensional data
pxinfo.pxSize = [344 344 127]; % size of projected image
pxinfo.padSize = [0 0 floor(0.25*pxinfo.pxSize(3))]; % how much to pad image in each direction
pxinfo.pxSizePadded = pxinfo.pxSize + 2*pxinfo.padSize; % size of padded image
pxinfo.pxdims = [2.0445 2.0445 2.03125]; % dimensions of projected image voxels in mm
pxinfo.sino = [344 252 837]; % size of (3D, span-11) sinogram

% Form a numerical phantom:
switch phantom_type
    case 'cylinder'
%         imvol = phantom('Modified Shepp-Logan',nx,[1 0.8 0.9 0 0 0]);
%         
%         R = 0.15;
%         r = 0.5;
%         for ii = 1:8
%             x = r*cos(respSigs(ii));
%             y = r*sin(respSigs(ii));
%             imvol = imvol + ...
%                 phantom('Modified Shepp-Logan',nx,[2 R-0.015*ii R-0.015*ii x y 0]);
%         end
%         imvol = repmat(imvol,[1 1 nz/2]);
%         imvol = padarray(imvol,[0 0 nz/4],0,'both');
%         
    case 'sphere'
        % define spheres (in millimetres), centre of image == [0 0 0]
        % - radii (mm):
        sphRad = [72 36 12 12 12 12 12 12 12 12 12 12];
        
        % - relative intensities:
        % (these will be included into the image additively)
        sphAddIntens = [1 2 2 2 2 2 2 2 2 2 2 2];
        
        % - fractional positions of each sphere
        fracs = [[0 0 0];
            [0 0 0]
            [.25 .25 .25]
            [.25 .25 -.25]
            [.25 -.25 .25]
            [-.25 .25 .25]
            [-.25 -.25 .25]
            [.25 -.25 -.25]
            [-.25 .25 -.25]
            [-.25 -.25 -.25]
            [0 0 .75]
            [0 0 -.75]];
        
        if numel(sphRad) ~= numel(sphAddIntens) || numel(sphRad) ~= size(fracs,1)
            error('Inconsistent number of radii vs intensities');
        end
        
        % calculate positions in mm:
        sphPos = zeros(numel(sphRad),3);
        for sph = 1:numel(sphRad)
            sphPos(sph,:) = fracs(sph,:).*(pxinfo.pxdims.*pxinfo.pxSize/2);
        end
        
        imvol = zeros(pxinfo.pxSize);
        for sph = 1:numel(sphRad)
            
            % Centre of current sub-sphere:
            sphCent = sphPos(sph,:);
            
            % Co-ordinate space:
            % (defined by body orientation)
            ap_vec = linspace(-(pxinfo.pxdims(1)*pxinfo.pxSize(1))/2,(pxinfo.pxdims(1)*pxinfo.pxSize(1))/2,pxinfo.pxSize(1));
            rl_vec = linspace(-(pxinfo.pxdims(2)*pxinfo.pxSize(2))/2,(pxinfo.pxdims(2)*pxinfo.pxSize(2))/2,pxinfo.pxSize(2));
            hf_vec = linspace(-(pxinfo.pxdims(3)*pxinfo.pxSize(3))/2,(pxinfo.pxdims(3)*pxinfo.pxSize(3))/2,pxinfo.pxSize(3));
            
            [ap,rl,hf] = ndgrid(ap_vec - sphCent(1),rl_vec - sphCent(2),hf_vec - sphCent(3));
            clear ap_vec rl_vec hf_vec;
            
            % Radial ordinate of current sphere:
            %    rr = sqrt( (ap - sphCent(1)).^2 + (rl - sphCent(2)).^2 + (hf - sphCent(3)).^2 );
            rr = sqrt( (ap).^2 + (rl).^2 + (hf).^2 );
            clear ap rl hf;
            
            % Draw current sphere (could toggle opacity/transparency here):
            imvol(rr <= sphRad(sph)) = imvol(rr <= sphRad(sph)) + sphAddIntens(sph);
            clear rr;
        end
    otherwise
        error('phantom type undefined');
end

% Attenuation (Mu) Map (in per-centimetres):
mumap = imvol;
mumap(imvol==1) = 0.11; % plastic
mumap(imvol==2) = 0.11; % plastic
mumap(imvol==3) = 0.097; % water

% Emission (Em) Map (in arbitrary units):
emmap = imvol;
emmap(imvol==1) = 1; % relatively cold
emmap(imvol==3) = 2; % warm
emmap(imvol==2) = 5; % hot

% Save numerical phantoms:
giplwrite(datainfo.maps.emission.fn,emmap*100,'short_noscale',pxinfo.pxdims);
giplwrite(datainfo.maps.attenuation.human,mumap*10000,'short_noscale',pxinfo.pxdims);
giplwrite(datainfo.maps.attenuation.hardware,0,'short_noscale',pxinfo.pxdims);

if flags.visualise
    figure(1);
    set(1,'OuterPosition',[762.6,34.6,780.8,830.4]);
    subplot(211);
    ap_vec = linspace(-(pxinfo.pxdims(1)*pxinfo.pxSize(1))/2,(pxinfo.pxdims(1)*pxinfo.pxSize(1))/2,pxinfo.pxSize(1));
    rl_vec = linspace(-(pxinfo.pxdims(2)*pxinfo.pxSize(2))/2,(pxinfo.pxdims(2)*pxinfo.pxSize(2))/2,pxinfo.pxSize(2));
    hf_vec = linspace(-(pxinfo.pxdims(3)*pxinfo.pxSize(3))/2,(pxinfo.pxdims(3)*pxinfo.pxSize(3))/2,pxinfo.pxSize(3));
    [ap,rl,hf] = ndgrid(ap_vec,rl_vec,hf_vec);
    [faces,verts,colors] = isosurface(ap,rl,hf,imvol,0,ap);
    patch('Vertices', verts, 'Faces', faces, ...
        'FaceVertexCData', colors, ...
        'FaceColor','interp', ...
        'edgecolor', 'interp');
    title('Outline of 3D phantom');
    xlabel('AP position');
    ylabel('RL position');
    zlabel('HF position');
    view(30,15);
    axis vis3d image;
    colormap copper;
    subplot(212);
    imagesc(squeeze(imvol(:,end/2,:))'); axis image;
    title('Phantom Cross-Section (coronal slice)');
    clear ap_vec rl_vec hf_vec ap rl hf;
end

%% Motion definition
% In this case, the respiratory signal will be the translation of the
% centre of the phantom. This could be extended to a 'navigator', and could
% also be extended to non-rigid models.

nTimeSteps = datainfo.nDynamics; % currently = 4 (hardcoded)
model_type = 'linear';
% model_type = 'quadratic';

% Generate respiratory signals:
respSig_all = linspace(-12,12,nTimeSteps); % CoM in mm
save(datainfo.respSig.fn,'respSig_all');

% Motion model creation
% Generate transformations to simulate for phantom data (not for recon)
% note: currently only a rigid up/down translation
motionModel = FormMotionModel('simple',model_type,pxinfo,true);

%% Transform phantom 
% n.b.: derives transformations from motion model directly, which makes
% this unsuitable for proper respSig estimation experiments since it is 
% a circular argument

[dHF,dAP,dRL] = EstimateFieldsFromModel(respSig_all,motionModel,'fwd',pxinfo);

PaddedEmMap = ToggleImagePadding(emmap,pxinfo);
PaddedMuMap = ToggleImagePadding(mumap,pxinfo);
[transEmMaps,transMuMaps] = ...
    TransformMaps(PaddedEmMap,PaddedMuMap,pxinfo,dHF,dAP,dRL,'linear',0);

if flags.visualise
    imsize = pxinfo.pxSizePadded;
    ap_vec = linspace(-(pxinfo.pxdims(1)*imsize(1))/2,(pxinfo.pxdims(1)*imsize(1))/2,imsize(1));
    % rl_vec = linspace(-(pxinfo.pxdims(2)*pxinfo.imsize(2))/2,(pxinfo.pxdims(2)*imsize(2))/2,imsize(2));
    hf_vec = linspace(-(pxinfo.pxdims(3)*imsize(3))/2,(pxinfo.pxdims(3)*imsize(3))/2,imsize(3));
    figure(3);
    set(3,'OuterPosition',[762.6,34.6,780.8,830.4]);
    subplot(212);
    toShow = zeros(imsize);
    for tt = 1:nTimeSteps, toShow = toShow + transEmMap{tt}; end
    toShow = toShow/nTimeSteps;
    imagesc(ap_vec,hf_vec,squeeze(toShow(:,end/2,:))',[0 max(imvol(:))]);
    title('Mean distribution');
    colormap copper;
    axis image;
    clear toShow;
    pause(0.5);
    subplot(211);
    for tt = 1:nTimeSteps
        imagesc(ap_vec,hf_vec,squeeze(transEmMap{tt}(:,end/2,:))',[0 max(imvol(:))]);
        title(sprintf('Transformed position %d/%d',tt,nTimeSteps));
        xlabel('AP direction');
        ylabel('HF direction');
        colormap copper;
        axis image;
        drawnow;
        pause(0.5);
    end
    clear ap_vec rl_vec hf_vec;
end

%% Project phantoms
projectorType = 3;
trueSinograms = SimulateTrueDistribution(transEmMaps,projectorType,pxinfo);

for rr = 1:2:252
    imagesc(squeeze(trueSinograms{1}(:,rr,1:127))',[0 0.8*max(trueSinograms{1}(:))]);  axis image;
    title(179/(252-rr)); drawnow; pause(0.1);
end
