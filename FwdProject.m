function [sinograms] = FwdProject(imgToProject,projectorType,pxinfo)
% FWDPROJECT [sinograms] = FwdProject(imgToProject,projectorType,pxinfo)
%	This function will take volumetric images and calculate a 3D sinogram for each
% 	using a projector function, specified by projectorType.
%	NOTE: if images are provided with an extended FOV, they will be de-padded automatically.


% Ensure input is cell format, otherwise make compatible:
if iscell(imgToProject)
	nToProject = numel(imgToProject);
else
	nToProject = 1;
	thisImgCell{1} = imgToProject;
    imgToProject = thisImgCell;
    clear thisImgCell;
end

% If padded, crop to make compatible with projector sizes:
currSize = size(imgToProject{1});
if abs(sum(currSize(1:3) - pxinfo.pxSizePadded)) < 1E-4
    % Padded image, therefore crop:
    toProject_unpadded = ToggleImagePadding(imgToProject,pxinfo);
else
    % Not padded, therefore continue:
    toProject_unpadded = imgToProject;
end
clear currSize;

% Project the input images, one at a time
% (cannot parallelise at this point since APIRL and CECR are already threaded)
for it = 1:nToProject
	% Identify which projector method will be used, and apply it
	switch projectorType
	% Matlab's basic projectors:
	case 1
%		if abs(size(toProject_unpadded{it},3) - pxinfo.sino(3)) > 1E-4
		if abs(size(toProject_unpadded{it},3) - 127) > 1E-4
			error('Matlab projector requires number of transaxial image slices == those in the sinogram');
        end
        
        theta = linspace(0,179,252);
		parproject = zeros([pxinfo.sino(1:2) 127]);
		parfor zz = 1:pxinfo.pxSize(3)
	    	[temp,xp] = radon(toProject_unpadded{it}(:,:,zz),theta);
            [~,lowR] = max(xp==-344/2);
            [~,highR] = max(xp==+344/2);
	    	parproject(:,:,zz) = temp(lowR+1:highR,:);
        end
        sinograms{it} = parproject;
	    clear temp parproject;
        
	% CECR projectors:
	case 2
        if it == 1
            InitialiseCECR;
        end
		%sino=cecrfwd(image,scannerModelFile,scannerModelNumber,totalThreadNumber,totalSubsetNumber,subsetIndex,verbose_level);
        sinograms{it} = cecrfwd(single(toProject_unpadded{it}),input_file,5002,32,1,0,0);
	
    % APIRL projectors:
	case 3
        if it == 1
            InitialiseApirl;
        end
		sinograms{it} = PET.P(toProject_unpadded{it});

	% Handle projectorType not_in {1,2,3}
	otherwise
		error('Unrecognised projectorType requested')
	end

end


end