function [bckProjImage] = BckProject(sinoToProject,projectorType,interpType,pxinfo)
% BCKPROJECT [bckProjImage] = BckProject(sinoToProject,projectorType,interpType,pxinfo)
%	This function will take a sinogram and calculate a 3D image
% 	using a back-projector function, specified by projectorType.
%	Remember that this is NOT the inverse of the FwdProject function.


if iscell(sinoToProject)
	nToProject = numel(sinoToProject);
else
	nToProject = 1;
	sinoToProject{1} = sinoToProject;
end

% For each sinogram to be projected:
for it = 1:nToProject

	% Identify which projector method will be used, and apply it:
	switch projectorType
	% Matlab's basic projectors
	case 1
		temp = zeros(pxinfo.pxSize);
		for zz = 1:nz
			temp(:,:,zz) = iradon(sinoToProject(:,:,zz),0:179,interpType,'none');
		end
		bckProjImage{it} = temp(2:(end-1),2:(end-1),:);
		clear temp;

	% CECR projectors
	case 2
		InitialiseCECR;

		%img=cecrbck(sino,scannerModelFile,scannerModelNumber,totalThreadNumber,totalSubsetNumber,subsetIndex,verbose_level);
		bckProjImage{it} = cecrbck(single(sinoToProject),input_file,5002,32,1,0,0);

	% APIRL projectors
	case 3
		InitialiseApirl;

		bckProjImage{it} = PET.PT(sinoToProject);

	% Handle projectorType not_in {1,2,3}
	otherwise
		error('Unrecognised projectorType requested')
	end

	bckProjImage{it} = ToggleImagePadding(bckProjImage{it},pxinfo);
end

warning('BckProject is currently stubbed: does not pad images')


end