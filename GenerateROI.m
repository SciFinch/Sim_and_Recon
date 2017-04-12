function roiImage = GenerateROI(datasetStr,regionID,pixelInfo)
%GENERATEROI roiImage = GenerateROI(dataset,regionID,pixelInfo);
%   A region of interest is sometimes required to encourage convergence to
%   an appropriate solution in motion estimating MCIR. This is something
%   that can be played with, but generally depends on which region is
%   required (specified by regionID), and will depend on the specific
%   distribution of a given dataset's FDG map. The main function of this
%   code is to tidy this up and make it more repeatable.

if strcmp(datasetStr,'stub')
    
    if strcmp(regionID,'all')
        roiImage = ones(344,344,127);
        
    elseif strcmp(regionID,'RHS')
        roiImage = ones(344,344,127); 
        warning('RHS ROI specification is a stub')
    end
    
else
    error('Dataset string not recognised')
end

end

