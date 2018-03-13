function [bckProjImage] = BckProject(sinoToProject,projectorType,interpType,pxinfo,subsetNum)
% BCKPROJECT [bckProjImage] = BckProject(sinoToProject,projectorType,interpType,pxinfo,subsetNum)
%	This function will take a sinogram and calculate a 3D image
% 	using a back-projector function, specified by projectorType.
%	Remember that this is NOT the inverse of the FwdProject function.

nSubsets = 21;
nAnglesPerSubset = 12;
if nargin < 5
    subsetNum = -1;
end

if iscell(sinoToProject)
	nToProject = numel(sinoToProject);
else
	nToProject = size(sinoToProject,4);
    for it = 1:nToProject
		temp{it} = sinoToProject(:,:,:,4);
    end
    sinoToProject = temp;
    clear temp;
end

% For each sinogram to be projected:
for it = 1:nToProject

	% Identify which projector method will be used, and apply it:
	switch projectorType
	% Matlab's basic projectors
	case 1
        theta = linspace(0,179,pxinfo.sino(2));
        
        if subsetNum > 0
            subsetIdx = GetSubsetIndices(subsetNum); %shuffled
            angles_to_project = theta(subsetIdx);
            
            if size(sinoToProject{it},2) ~= nAnglesPerSubset
                error('Expected a subset sinogram');
            end
        else
            angles_to_project = theta;
        end
        
		temp = zeros(pxinfo.pxSize);
		parfor zz = 1:pxinfo.pxSize(3)
			temp(:,:,zz) = iradon(sinoToProject{it}(:,:,zz),angles_to_project,interpType,'none',pxinfo.pxSize(1));
		end
		bckProjImage{it} = temp;
		clear temp;

	% CECR projectors
	case 2
        if subsetNum > 0
            error('subset initialisation not defined for CECR projectors');
        end
		InitialiseCECR;
        
        % Note: the projectors are insensitive to a global scale factor in
        % the images, so a factor of 1000 is used to avoid a bug where
        % input values smaller than 1E-5 produce an output of 0

		%img=cecrbck(sino,scannerModelFile,scannerModelNumber,totalThreadNumber,totalSubsetNumber,subsetIndex,verbose_level);
		bckProjImage{it} = cecrbck(1000*single(sinoToProject{it}),input_file,5002,32,1,0,0);
        bckProjImage{it} = bckProjImage{it}/1000;
        
	% APIRL projectors
	case 3
        if subsetNum > 0
            error('subset initialisation not defined for APIRL projectors');
        end        
		InitialiseApirl;

		bckProjImage{it} = PET.PT(sinoToProject{it});

	% Handle projectorType not_in {1,2,3}
	otherwise
		error('Unrecognised projectorType requested')
	end
end

% Pad images to padded FOV size:
bckProjImage = ToggleImagePadding(bckProjImage,pxinfo);

end