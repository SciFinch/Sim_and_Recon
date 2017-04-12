function [imageArray,imageSize,imagePixelDims] = AutoLoadImage(fullFileName,varName)
% AUTOLOADIMAGE [imageArray,imageSize,imagePixelDims] = AutoLoadImage(fullFileName)
% 	This function will automatically load medical image formats and return a 3D
%	image array, its size in pixels, and the physical dimension of each pixel
%	in millimetres.
%	This is still in development, so currently can only process gipls, niftis,
%	and matlab .mat formats. If mat format is required, an optional variable
%	name can be provided (as a string in varName), which will select the specific
%	variable from the .mat file corresponding to the image.

if nvargin == 1
	varName = [];
end

% - identify possible scalefactor sometimes used to preserve decimal places
scalefactor = 1;
if strfind(fullFileName, 'x10')
    scalefactor = 10;
elseif strfind(fullFileName, 'x100')
    scalefactor = 100;
elseif strfind(fullFileName, 'x1000')
    scalefactor = 1000;
elseif strfind(fullFileName, 'x10000')
    scalefactor = 10000;
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
	formatStr = 'interfile';
	error('Support for loading interfile has not yet been implemented')
elseif strfind(fullFilename, '.dat')
	error('Filename provided is .dat format, which is not [yet] supported')
else
    error('Unable to determine the motion field file format');
end

switch formatStr
case 'gipl'
	addpath('C:\Users\db12\Dropbox\Work_Docs\MATLAB Scripts\giplread.m');
	[imageArray, imageSize, imagePixelDims] = giplread(filename);
	imageArray = imageArray/scalefactor;
case 'nifti'
	error('nifti file loader hasnt been checked yet');	
	addpath('C:\Users\db12\Google Drive\CommonLib\DeedsRegistration\toolbox_matlab_nifti\');
	dataStruct = MRIread(filename);
	imageArray = dataStruct.vol/scalefactor;
	imageSize = dataStruct;
	imagePixelDims = dataStruct;
case 'interfile'
	warning('interfile reading is a stub');
	imageArray = [];
	imageSize = [];
	imagePixelDims = [];
case 'mat'
	error('mat file loader hasnt been checked yet');
	% warning: this assumes the opposite function was used to save the data
	if ~isempty(varName)
		dataStruct = load(filename,varName);
		imageSize = [];
		imagePixelDims = [];
	else
		dataStruct = load(filename);
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