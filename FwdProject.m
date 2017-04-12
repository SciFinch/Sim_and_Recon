function [sinogram] = FwdProject(imgToProject,projectorType,pxinfo)
% FWDPROJECT [sinogramm] = FwdProject(imgToProject,projectorType,pxinfo)
%	This function will take volumetric images and calculate a 3D sinogram for each
% 	using a projector function, specified by projectorType.
%	NOTE: if images are provided with an extended FOV, they will be de-padded automatically.


% Ensure input is cell format, otherwise make compatible:
if iscell(imgToProject)
	nToProject = numel(imgToProject);
else
	nToProject = 1;
	imgToProject{1} = imgToProject;
end


% Project the input images, one at a time
% (cannot parallelise at this point since APIRL and CECR are already threaded)
for it = 1:nToProject

	% If padded, crop to make compatible with projector sizes:
	currSize = size(imgToProject{it});
	if abs(sum(currSize(1:3) - pxinfo.pxSizePadded)) < 1E-4
		% Padded image, therefore crop:
		toProject_unpadded = ToggleImagePadding(imgToProject{it},pxinfo);
	else
		% Not padded, therefore continue:
		toProject_unpadded = imgToProject{it};
	end
	clear currSize;


	% Identify which projector method will be used, and apply it
	switch projectorType
	% Matlab's basic projectors
	case 1
		if abs(size(toProject_unpadded,3) - pxinfo.sino(3)) > 1E-4
			error('Matlab projector requires number of transaxial image slices == those in the sinogram');
		end

		sinogram = zeros(pxinfo.sino);
		parfor zz = 1:pxinfo.pxSize(3)
	    	[temp,~] = radon(toProject_unpadded(:,:,zz));
	    	sinogram{it}(:,:,zz) = temp;
		end
	    clear temp;

	% CECR projectors
	case 2
		InitialiseCECR;

		%sino=cecrfwd(image,scannerModelFile,scannerModelNumber,totalThreadNumber,totalSubsetNumber,subsetIndex,verbose_level);
        sinogram{it} = cecrfwd(single(toProject_unpadded),input_file,5002,32,1,0,0);
	% APIRL projectors
	case 3
		InitialiseApirl;

		sinogram{it} = PET.P(toProject_unpadded);

	% Handle projectorType not_in {1,2,3}
	otherwise
		error('Unrecognised projectorType requested')
	end

end


end