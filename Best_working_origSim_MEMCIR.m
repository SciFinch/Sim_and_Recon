%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       BEST WORKING VERSION
% This is the best working version of ME-MCIR using simple simulations and
% multi-2D projectors.
% Included:
%   - Attenuation
%
%
% The following is a list of versions saved in this file previously:
%   - Best_working_basicSim_MEMCIR.m (used as base code)
%   [C:\Users\db12\Documents\MATLAB\testingMCIR_extdFOV_3D_MEMCIR_AC__hiI_supression.m on 22/02/2017 (DB)]
%
%
%NOTE: There is a trivial sign error in this version of the algorithm. It
%does not affect the algorithm's ability to function, but should probably
%be ironed out in future versions. See original file for details. This only
%affects th_est vs th_obs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
%interpType = 'nearest';
interpType = 'linear';

%InitialiseApirl;

%% Simulation
phantom_type = 'sphere';
nx = 64;
nz = 32;
enx = nx + round(0.5*nx);
enz = nz + round(0.5*nz);
nx_s = round((enx-nx)/2) + 1;
nz_s = round((enz-nz)/2) + 1;
nx_f = enx - round((enx-nx)/2);
nz_f = enz - round((enz-nz)/2);

% Respiratory signal values:
th = linspace((2*pi)/8,2*pi,8);

% Form a numerical phantom:
if strcmp(phantom_type,'cylinder')
    imvol = phantom('Modified Shepp-Logan',nx,[1 0.8 0.9 0 0 0]);

    R = 0.15;
    r = 0.5;
    for ii = 1:8
        x = r*cos(th(ii));
        y = r*sin(th(ii));
        imvol = imvol + ...
            phantom('Modified Shepp-Logan',nx,[2 R-0.015*ii R-0.015*ii x y 0]);
    end
    imvol = repmat(imvol,[1 1 nz/2]);
    imvol = padarray(imvol,[0 0 nz/4],0,'both');
elseif strcmp(phantom_type,'sphere')
    imvol = zeros(nx,nx,nz)+0.0;
    
    % define image features
    sphRad = [nx/6,nx/16,nx/16,nx/16,nx/16,nx/16,nx/16,nx/16,nx/16];
    sphPos = [[nx nx nz]'/2,[nx nx nz]'/2,[nx/3 nx/3 nz/3]',[2*nx/3 nx/3 nz/3]',...
        [nx/3 2*nx/3 nz/3]',[nx/3 nx/3 2*nz/3]',[2*nx/3 2*nx/3 nz/3]',...
        [2*nx/3 nx/3 2*nz/3]',[nx/3 2*nx/3 2*nz/3]',[nx/3 2*nx/3 2*nz/3]'];
    
    sphAddIntens = [1 2 2 2 2 2 2 2 2];
    
    for sph = 1:numel(sphRad)
        
    sphCent = sphPos(:,sph);
    for xx = 1:nx
        for yy = 1:nx 
            for zz = 1:nz
                if (xx - sphCent(1))^2 + (yy - sphCent(2))^2 + (zz - sphCent(3))^2 <= sphRad(sph)^2
                    imvol(xx,yy,zz) = imvol(xx,yy,zz) + sphAddIntens(sph);
                end
            end
        end
    end
    end
else
    error('phantom type undefined');
end

imvol = padarray(imvol,round([(enx-nx)/2 (enx-nx)/2 (enz-nz)/2]),0,'both');
imagesc(squeeze(imvol(:,round((enx-1)/2),:))'); axis image;

%% Motion definition
T = 6;
A = 6;%nz/4-1; %2 works
for t = 1:T
    th_obs(t) = pi*(t-T/2)/T;
    dz(t) = A*sin(th_obs(t));
end
dx = zeros(1,T);
dy = zeros(1,T);

attvol = imvol;
attvol(imvol==1) = 0.096; %per cm
attvol(imvol==2) = 0.15; %per cm
attvol(imvol==3) = 0.05; %per cm
attvol = attvol/10; % per cm to per mm
imagesc(squeeze(attvol(:,round((enx-1)/2),:))'); axis image;


[x,y,z] = ndgrid(1:enx,1:enx,1:enz);
for t = 1:T
    tImvol(:,:,:,t) = interpn(x,y,z,imvol,x+dx(t),y+dy(t),z+dz(t),'linear',0);
    tAttvol(:,:,:,t) = interpn(x,y,z,attvol,x+dx(t),y+dy(t),z+dz(t),'linear',0);
end

% figure(1);
% subplot(1,3,1); imagesc(squeeze(mean(tImvol(:,round((nx-1)/2),:,:),4))'); axis image; title('Motion Affected Ground Truth')
% subplot(1,3,2); imagesc(squeeze(imvol(:,round((nx-1)/2),:))'); axis image; title('Ground Truth')
% subplot(1,3,3); imagesc(squeeze(tImvol(:,round((nx-1)/2),:,1))'); axis image; title('Motion Affected Ground Truth')

for t = 1:T
    for zz = 1:nz
        [temp,xp] = radon(tImvol(nx_s:nx_f,nx_s:nx_f,nz_s+zz-1,t));
        trueS(:,:,zz,t) = temp;
        
        [temp,xp] = radon(tAttvol(nx_s:nx_f,nx_s:nx_f,nz_s+zz-1,t));
        attS(:,:,zz,t) = exp(-temp);
    end
end
clear temp;
data = trueS.*attS;

% Test motion correction:
% for zz=1:nz
%     imagesc(squeeze(data(:,:,zz,T))'); title(zz); drawnow; pause(0.1);
% end

%nCounts = 10E10;
%nCounts = nz * 50E6/500;
nCounts = 50E6/500;
currCounts = sum(data(:));
data = data*nCounts/currCounts;
data = poissrnd(data);
data = data*currCounts/nCounts;

%% Reconstruction

%f = f/M'A'X'1 M'A'X' m/XAMf
nIterations = 200;
beta = 0.05*0.000001;
imEst = ones(enx,enx,enz);
th_est = pi*ones(T,1);
thRec = [];%NaN*zeros(T,nIterations);
endplane = 0;
repflag = 0;
epsilon = 0;%1E-4;
for k = 1:nIterations
    
    % PREPARATION
    % Find current motion states
    for t = 1:T
        dx(t) = 0;
        dy(t) = 0;
        dz(t) = A*sin(th_est(t));
    end
    
    % Estimate attenuation effects
    for t = 1:T
        for zz = 1:nz
            tAttEst(:,:,:,t) = interpn(x,y,z,attvol,x+dx(t),y+dy(t),z+dz(t),interpType,0);
            [temp,xp] = radon(tAttEst(nx_s:nx_f,nx_s:nx_f,nz_s+zz-1,t));
            attEstS(:,:,zz,t) = exp(-temp);
        end
    end
    clear temp;
    
    % Find current sensitivity image
    sensImg = zeros(size(imEst));  
    for t = 1:T
        for zz = 1:nz
            temp(:,:,zz) = iradon(attEstS(:,:,zz,t),0:179,interpType,'none');
        end
        temp = temp(2:(end-1),2:(end-1),:);
        normImg = padarray(temp,[(enx-nx)/2 (enx-nx)/2 (enz-nz)/2],0,'both');
                      
%        sensImg = sensImg + interpn(x+dx(t),y+dy(t),z+dz(t),normImg,x,y,z,interpType,eps);
        sensImg = sensImg + interpn(x,y,z,normImg,x+dx(t),y+dy(t),z+dz(t),interpType,eps);
        clear temp;
    end
    subplot(3,3,1); imagesc(squeeze(sensImg(:,round((enx-1)/2),:))'); title('Sensitivity'); axis image; drawnow;
    subplot(3,3,2); imagesc(squeeze(sensImg(nx_s:nx_f,round((enx-1)/2),nz_s:nz_f))'); title('Sensitivity (crop)'); axis image; drawnow;
    
    % ESTIMATE IMAGE
    subplot(3,3,4); imagesc(squeeze(imEst(:,round((enx-1)/2),:))'); axis image; title(['Estimate ' num2str(k-1)]); drawnow;
    subplot(3,3,5); imagesc(squeeze(imEst(nx_s:nx_f,round((enx-1)/2),nz_s:nz_f))'); axis image; title(['Estimate (crop)']); drawnow;
    updateIm = zeros(enx,enx,enz);
    for t = 1:T
%        transImg = interpn(x,y,z,imEst,x+dx(t),y+dy(t),z+dz(t),interpType,0);
        transImg = interpn(x+dx(t),y+dy(t),z+dz(t),imEst,x,y,z,interpType,0);
        
        for zz = 1:nz
            [temp,xp] = radon(transImg(nx_s:nx_f,nx_s:nx_f,nz_s+zz-1));
            imEstS(:,:,zz) = temp;
        end
        clear temp;
        datamodel = imEstS.*attEstS(:,:,:,t);
        
        ratioSino = attEstS(:,:,:,t).*squeeze(data(:,:,:,t))./(datamodel + epsilon);
        % get a high intensity here on zplane 27 during g=1, iteration 16, nn interp 
        for zz = 1:nz
            temp(:,:,zz) = iradon(ratioSino(:,:,zz),0:179,interpType,'none');
        end
        temp = temp(2:(end-1),2:(end-1),:);
        temp = padarray(temp,[(enx-nx)/2 (enx-nx)/2 (enz-nz)/2],0,'both');
%        temp = interpn(x+dx(t),y+dy(t),z+dz(t),temp,x,y,z,interpType,0);
        temp = interpn(x,y,z,temp,x+dx(t),y+dy(t),z+dz(t),interpType,0);
              
        updateIm = updateIm + temp;
        clear temp;
    end
   % updateIm = updateIm;
    subplot(3,3,3); imagesc(squeeze(mean(datamodel(:,round((nx-1)/2),:),4))'); title('Mean Data Model'); drawnow;
    subplot(3,3,7); imagesc(squeeze(mean(updateIm(:,round((enx-1)/2),:),4))'); axis image; title(['Update ' num2str(k)]); drawnow;
    subplot(3,3,8); imagesc(squeeze(updateIm(nx_s:nx_f,round((enx-1)/2),nz_s:nz_f))'); axis image; title(['Update ' num2str(k) ' (crop)']); drawnow;
       
    imEst = (imEst./(sensImg + epsilon)).*updateIm;
    if endplane ~= 0
        if repflag == 0
            imEst(:,:,1:endplane) = 0; imEst(:,:,(end-(endplane-1)):end)=0;
        else
            imEst(:,:,1:endplane) = repmat(imEst(:,:,endplane+1),[1 endplane]);
            imEst(:,:,(end-(endplane-1)):end)= repmat(imEst(:,:,end-endplane),[1 endplane]);
        end
    end
    imEst(isnan(imEst))=0;
    
    % ESTIMATE MOTION PARAMETERS
    if k > 6
        for t = 1:T
            th = th_est(t);
            
%            transImg = interpn(x,y,z,imEst,x+dx(t),y+dy(t),z+dz(t),interpType,0);
            transImg = interpn(x+dx(t),y+dy(t),z+dz(t),imEst,x,y,z,interpType,0);
            toproj = transImg;
            toproj(isnan(toproj))=0;
            for zz = 1:nz
                [temp,xp] = radon(toproj(nx_s:nx_f,nx_s:nx_f,nz_s+zz-1));
                imEstS(:,:,zz) = temp;
            end
            clear temp;
            datamodel = imEstS.*attEstS(:,:,:,t);
            
            ratioSino = squeeze(data(:,:,:,t))./(datamodel + epsilon);
            ratioSino(isnan(ratioSino))=0;
            
            % Motion model derivative:
            dXdth = 0;
            dYdth = 0;
            dZdth = A*cos(th);            
            
            % project model-weighted derivative of emission image
            [im_xgrad,im_ygrad,im_zgrad] = gradient(transImg);
            weightIm = (dXdth*im_xgrad+dYdth*im_ygrad+dZdth*im_zgrad);
            weightIm(isnan(weightIm))=0;
            for zz = 1:nz
                [temp,xp] = radon(weightIm(nx_s:nx_f,nx_s:nx_f,nz_s+zz-1));
                weightEmSino(:,:,zz) = temp;
            end
            clear temp;
            
            % project model-weighted derivative of mu-map
            [mu_xgrad,mu_ygrad,mu_zgrad] = gradient(tAttEst(:,:,:,t));
            weightMu = (dXdth*mu_xgrad+dYdth*mu_ygrad+dZdth*mu_zgrad);
            weightMu(isnan(weightMu))=0;
            for zz = 1:nz
                [temp,xp] = radon(weightMu(nx_s:nx_f,nx_s:nx_f,nz_s+zz-1));
                weightMuSino(:,:,zz) = temp;
            end
            clear temp;
            
            totWeightS = attEstS(:,:,:,t).*weightEmSino - ...
                datamodel.*weightMuSino;

            thUpdateSino = (1-ratioSino).*totWeightS;
            thUpdateSino(isnan(thUpdateSino)|isinf(thUpdateSino)) = 0;
          
            subplot(3,3,6); imagesc(squeeze(totWeightS(:,round((enx-1)/2),:))'); axis image; title('Weight Sino'); drawnow;
            
            thUpdate(t) = sum(thUpdateSino(:));
        end
        th_est = th_est + beta*thUpdate';
    end
    thRec = [thRec th_est];

    subplot(3,3,9); plot(thRec'); axis square; title('Update Contribn'); drawnow;
end