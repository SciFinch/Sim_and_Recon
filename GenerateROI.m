function roiImage = GenerateROI(datasetStr,regionID,pixelInfo)
%GENERATEROI roiImage = GenerateROI(dataset,regionID,pixelInfo);
%   A region of interest is sometimes required to encourage convergence to
%   an appropriate solution in motion estimating MCIR. This is something
%   that can be played with, but generally depends on which region is
%   required (specified by regionID), and will depend on the specific
%   distribution of a given dataset's FDG map. The main function of this
%   code is to tidy this up and make it more repeatable.


switch datasetStr
    case 'simple'
        roiImage = ones(pixelInfo.pxSize);
    case 'stub'
        if strcmp(regionID,'all')
            roiImage = ones(pixelInfo.pxSize);
            
        elseif strcmp(regionID,'RHS')
            roiImage = ones(pixelInfo.pxSize);
            warning('RHS ROI specification is a stub')
        end
    otherwise
        error('Dataset string not recognised')
end

end

