function comprSino = CompressSino(sinoToCompr,currSpan,mashFactor,normFlag)
% comprSino = CompressSino(sinoToCompr,currSpan,mashFactor,normFlag)
%
% Compresses a [span-11] mMR sinogram into a segment-0 sinogram using the
% SSRB algorithm. This is done using a pre-estimated* Michelogram. This can
% also include mashing, where the mashFactor provides how many pixels
% are combined into a single pixel in the output. Note that currently this
% only works if the mashFactor is a divisor of the number of azimuthal
% angles in the sinogram (252 for mMR sinograms).
%
% normFlag is used for normalisation (1 for attenuation, 0 for prompts)
% - use this if sinogram values are to be averaged or summed, respectively
%
% *There are multiple ways to improve this function, including
% (1) Making the Michelogram calculation direct
% (2) Making the Michelogram calc non-empirical (as is at present)
% (3) Intro'ing support for non-span-11 sinograms, inc tying parameters to
% the span itself
% (4) Including a maximal angle of acceptance
% [see ...\Work_Docs\compressSinogram.m for Michelogram calc]
% Written by DB, 14/3/2017

if currSpan ~= 11
    error('Span other than span 11 given; currently unsupported');
end

% Load axial positions for each axial slice
temp = load('C:\Users\db12\Dropbox\Work_Docs\MATLAB_constFiles\CompressSino.mat');
slicePositions = temp.slicePositions;
clear temp;

% Span-11 Sinogram Specifics
maxSeg = 5;
nSegments = 2*maxSeg+1;
segNum = containers.Map(1:currSpan,[0 -1 1 -2 2 -3 3 -4 4 -5 5]);
nAxialBins = containers.Map(-maxSeg:maxSeg,[27,49,71,93,115,127,115,93,71,49,27]);

% Initialise variables
if mod(size(sinoToCompr,2),mashFactor)
    error('mashFactor must be a factor of the initial sinogram size')
end

comprSino = zeros([size(sinoToCompr,1),size(sinoToCompr,2),nAxialBins(0)]);
seg0_posns_rel = slicePositions{1};
normCount = zeros(1,nAxialBins(0)); % counts contributions to each seg0 slice

% Perform compression
for comprSliceNum = 1:nAxialBins(0)
    comprSlicePos = seg0_posns_rel(comprSliceNum);
    z = 1;
    for s = 1:nSegments
        local_z = 0;
        zmin = 0;
        for ss = 1:(s-1)
            zmin = zmin + nAxialBins(segNum(ss));
        end
        for z = (zmin+1):(zmin+nAxialBins(segNum(s)))
            local_z = local_z + 1;
            
            fullSlicePos = slicePositions{s};
            if abs(fullSlicePos(local_z) - comprSlicePos) < 1E-4
                comprSino(:,:,comprSlicePos*2+1) = comprSino(:,:,comprSlicePos*2+1) + sinoToCompr(:,:,z);
                normCount(comprSlicePos*2+1) = normCount(comprSlicePos*2+1) + 1;
            end
        end
    end
end


% Apply normalisation if necessary
if normFlag == 1
    for z = 1:nAxialBins(0)
        if normCount(z) > 0
            comprSino(:,:,z) = comprSino(:,:,z)/normCount(z);
        end
    end
end

% Perform mashing
if mashFactor > 1
    comprSino_mashed = zeros([size(comprSino,1),size(comprSino,2)/mashFactor,size(comprSino,3)]);
    % relies on integer conversion - is there a better way to do this?
    for az = 1:size(comprSino,2)
        comprSino_mashed(:,ceil(az/mashFactor),:) = ...
            comprSino_mashed(:,ceil(az/mashFactor),:) + comprSino(:,az,:);
    end
    if normFlag == 1
        comprSino_mashed = comprSino_mashed/mashFactor;
    end
    comprSino = comprSino_mashed;
end

end