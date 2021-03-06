function [imageArray,imageSize,imagePixelDims] = AutoLoadImage(fullFileName,varName)
% AUTOLOADIMAGE [imageArray,imageSize,imagePixelDims] = AutoLoadImage(fullFileName)
% 	This function will automatically load medical image formats and return a 3D
%	image array, its size in pixels, and the physical dimension of each pixel
%	in millimetres.
%	This is still in development, so currently can only process gipls, niftis,
%	and matlab .mat formats. If mat format is required, an optional variable
%	name can be provided (as a string in varName), which will select the specific
%	variable from the .mat file corresponding to the image.

if isempty(fullFileName)
    error('No image filename has been provided');
end

if nargin == 1
	varName = [];
end

% - identify possible scalefactor sometimes used to preserve decimal places
% note: Must start with highest number, otherwie 10x case would always be
% first one identified. Could try regex instead.
if strfind(fullFileName, 'x10000')
    scalefactor = 10000;
elseif strfind(fullFileName, 'x1000')
    scalefactor = 1000;
elseif strfind(fullFileName, 'x100')
    scalefactor = 100;
elseif strfind(fullFileName, 'x10')
    scalefactor = 10;
else
    scalefactor = 1;
end

% - identify filetype and load
formatStr = [];
if strfind(fullFileName, '.gipl')
    % is in gipl format
    formatStr = 'gipl';
elseif strfind(fullFileName, '.mat')
    % is in mat format
    formatStr = 'mat';
elseif strfind(fullFileName, '.nii') 
    % is in nifti format (not currently supported)
    formatStr = 'nifti';
elseif strfind(fullFileName, '.hv')
	% is in interfile format
	%formatStr = 'interfile';
	error('Support for loading interfile has not yet been implemented')
elseif strfind(fullFileName, '.dat')
	error('Filename provided is .dat format, which is not [yet] supported')
else
    error('Unable to determine the motion field file format');
end

switch formatStr
case 'gipl'
	oldPath = addpath('C:\Users\db12\Dropbox\Work_Docs\MATLAB Scripts\');
	[imageArray, imageSize, imagePixelDims] = giplread(fullFileName);
	imageArray = imageArray/scalefactor;
    path(oldPath); clear oldPath;
case 'nifti'
	error('nifti file loader hasnt been checked yet');	
	addpath('C:\Users\db12\Google Drive\CommonLib\DeedsRegistration\toolbox_matlab_nifti\');
	dataStruct = MRIread(fullFileName);
	imageArray = dataStruct.vol/scalefactor;
	imageSize = dataStruct;
	imagePixelDims = dataStruct;
case 'interfile'
	warning('interfile reading is a stub');
	imageArray = [];
	imageSize = [];
	imagePixelDims = [];
case 'mat'
	warning('mat file loader hasnt been checked yet');
	% warning: this assumes the opposite function was used to save the data
	if ~isempty(varName)
		dataStruct = load(fullFileName,varName);
		imageSize = [];
		imagePixelDims = [];
	else
		dataStruct = load(fullFileName);
	end
case 'dat'
	warning('Read-in for .dat format is currently a stub');
	imageArray = [];
	imageSize = [];
	imagePixelDims = [];
otherwise
	error('Unrecognised format string - something has gone wrong somewhere');
end
	

end