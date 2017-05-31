function padFlag = ispadded(image,pixelInfo)
% Usage: padFlag = ispadded(image,pixelInfo)
% Helper function - determines whether an image is padded or not, 
% returns 1 or 0. If input is a set of images stored in a cell,
% this function returns the a single value for all images.
% If images are a mix of padded and unpadded images, this function
% returns -1.

if ~iscell(image)
	% Simple case: simply check whether image is padded
	currSize = size(image);
	padFlag = (abs(sum(currSize(1:3) - pixelInfo.pxSizePadded)) < 1E-4);
else
	% If is a cell, check all images
	nImages = numel(image);
	pad_perImage = zeros(1,nImages);
	for it = 1:nImages
		currSize = size(image{it});
		pad_perImage(it) = (abs(sum(currSize(1:3) - pixelInfo.pxSizePadded)) < 1E-4);
	end

	% If not all images agree, return -1:
	if mean(pad_perImage) == 1.0
		padFlag = 1;
	elseif mean(pad_perImage) == 0.0
		padFlag = 0;
	else
		padFlag = -1;
        warning('Inconsistent padding in images');
	end	
end