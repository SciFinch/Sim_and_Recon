function [ transImg ] = TransformImage(imgToTrans,pxinfo,dHF,dAP,dRL,interpType,padVal)
%TRANSFORMIMAGE [ transImg ] = TransformImage(emMap,pxinfo,dHF,dAP,dRL,interpType,padVal)
%   This function will return a cell array consisting of generic images
%   transformed according to the motion fields provided.
%   This requires specification of the interpolation method and the values
%   used to pad the outside of the array.
%   Note that this is just a one-image version of TransformMaps, and it
%   might be worth replacing one for the other

% Note:
% If forward transformation, only 1 img is used. For back transformations,
% all images are transformed by their respective transformation. To make
% this compatible for either direction transformation for a group of a
% transformations, fwd transformations are applied to set of copies of the
% image to be transformed [this is, of course, unnecessary for back
% transformations]

%% quick test for suitibility
dR = sqrt(dHF{1}.^2 + dAP{1}.^2 + dRL{1}.^2);
mx_dR = max(dR(:));
if mx_dR > 25.0
    warning('Largest displacements above 25mm - MFs might not be scaled correctly');
end
clear mx_dR dR;

%% Initialise
% Find number of transformations to process:
nTrans = numel(dHF); 

% Number of images to process (must be either 1 or nTrans)
if iscell(imgToTrans)
    nImgs = numel(imgToTrans);
else
    nImgs = 1;
end

toTransCell = cell(1,nTrans);

switch nImgs
% - Case nImgs == 1: transform image into multiple positions:
case 1 
    % Copy image to transform
    for it = 1:nTrans
        toTransCell{it} = imgToTrans;
    end
% - Case nImgs == nTrans: 1-to-1 matching of each image to each transformation:
case nTrans
    % Prep each image for transform
    for it = 1:nTrans
        toTransCell{it} = imgToTrans{it};
    end
% - Otherwise, throw error because conditions not met
otherwise
    error('Number of images provided is inconsistent with the number of tranformations');
end

fprintf('== Transforming %d image(s) with %d motion field(s)\n',nImgs,nTrans)

% Assuming orientation [HF AP RL]
hf_vec = linspace(pxinfo.pxdims(1),pxinfo.pxdims(1)*pxinfo.pxSizePadded(1),pxinfo.pxSizePadded(1));
ap_vec = linspace(pxinfo.pxdims(2),pxinfo.pxdims(2)*pxinfo.pxSizePadded(2),pxinfo.pxSizePadded(2));
rl_vec = linspace(pxinfo.pxdims(3),pxinfo.pxdims(3)*pxinfo.pxSizePadded(3),pxinfo.pxSizePadded(3));

[coords.hf,coords.ap,coords.rl] = ndgrid(hf_vec,ap_vec,rl_vec);
clear *vec;

%% Perform transformations
transImg = cell(1,nTrans);

parfor it = 1:nTrans
    x = coords.hf;
    y = coords.ap;
    z = coords.rl;

    % Still assuming correct orientations...
    xdx = x + dHF{it};
    ydx = y + dAP{it};
    zdz = z + dRL{it};

    transImg{it} = interpn(x,y,z,toTransCell{it},xdx,ydy,zdz,interpType,padVal);
end
clear toTransCell x y z xdx ydy zdz coords;

% If input was not a cell and only one transformation was requred, 
% ensure output is not a cell either:
if (~iscell(imgToTrans) && (nTrans == 1))
    transImg = transImg{1};
end

warning('TransformImgs has not been data-tested yet');

end

