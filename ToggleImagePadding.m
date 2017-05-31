function [resizedImage] = ToggleImagePadding(imgToResize,pxinfo)
% TOGGLEIMAGEPADDING [resizedImage] = ToggleImagePadding(imgToResize,pxinfo)
%	This function is designed to extend or reduce the amount of zero padding around
% 	an image in a way that is compatible with the 'extendedFOV' aspect of the 
%	ME-MCIR algorithm. This has been converted into a function to ensure that it
% 	is standardised throughout the algorithm. Generally, images are only resized along
%	the third array dimension (corresponding to the transaxial direction) since this
%	is the only direction along which the human body leaves the FOV. This function
%	does support the other two directions too, and also supports a set of images
% 	to be resized, in cell format.
% 	Note that the pad/crop functionality is determined automatically by comparing
% 	the input image volume(s) with the sizes stored in pxinfo.

% Dev notes:
% 1) this currently only supports zero-padding as a default. Would be useful
% 	 to extend this to padding with anything from single numbers to whole
%	 image planes.
% 2) the function assumes that the padding is symmetric, and that when cropping,
%    end planes are simply removed. This could be made more sophisiticated
%	 by having a user-specified image centre, which could be extended in the future
%    into support for scanner bed positions

% Check whether image is cell array, and how many resizes are required:
if iscell(imgToResize)
	nToResize = numel(imgToResize);
else
	nToResize = 1;
	thisImageCell{1} = imgToResize;
    imgToResize = thisImageCell;
    clear thisImageCell;
end

% Determine whether the image is to be padded or cropped
flags.isPadded = CheckPadding(imgToResize,pxinfo);

% Apply cropping/padding to each input image:
for it = 1:nToResize
    switch flags.isPadded
        % If currently padded, crop:
        case true
            % Calculate original pixel subscripts in terms of padded subscripts:
            x_start = pxinfo.padSize(1) + 1;
            x_end = pxinfo.pxSizePadded (1) - pxinfo.padSize(1);
            y_start = pxinfo.padSize(2) + 1;
            y_end = pxinfo.pxSizePadded (2) - pxinfo.padSize(2);
            z_start = pxinfo.padSize(3) + 1;
            z_end = pxinfo.pxSizePadded (3) - pxinfo.padSize(3);
            % Crop to original subscripts:
            resizedImage{it} = imgToResize{it}(x_start:x_end,y_start:y_end,z_start:z_end);
            % If currently cropped, pad:
        case false
            % Pad image using padSize variables
            resizedImage{it} = padarray(imgToResize{it},pxinfo.padSize,0,'both');
    end
end

end