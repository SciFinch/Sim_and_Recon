function [ isPaddedFlag ] = CheckPadding(imagesToCheck,pixelInfo)
% Usage: [ isPaddedFlag ] = CheckPadding(imagesToCheck,pixelInfo)
% A simple checker function which determines whether an image is padded
% or not, in line with the FOV padding in the Sim_and_Recon project.
%   True => FOV is padded, and vice-versa

% Find size of input image(s)
if iscell(imagesToCheck)
    if numel(imagesToCheck) > 1
        % Check all images in cell are same size
        leadingSize = size(imagesToCheck{1});
        for it = 2:numel(imagesToCheck)
            thisSize = size(imagesToCheck{it});
            if sum(abs(thisSize - leadingSize)) > 10*eps
                error('Cannot check padding: Not all images in cell array are the same size!');
            end            
        end 
        currSize = leadingSize;
    else
        currSize = size(imagesToCheck{1});
    end
else
    currSize = size(imagesToCheck);
end

% Check whether the [first] image is padded or not, according the the
% definitions stored in the pixelInfo struct:
if sum(abs(currSize(1:3) - pixelInfo.pxSize)) < 10*eps
	isPaddedFlag = false;
elseif sum(abs(currSize(1:3) - pixelInfo.pxSizePadded)) < 10*eps
	isPaddedFlag = true;
else
    error('Unable to determine whether images are padded or not');
end

end