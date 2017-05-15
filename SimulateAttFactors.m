function [attFactorSinograms] = SimulateAttFactors(transMuMap,projectorType,pixelInfo)
%SIMULATEATTFACTORS [attFactorSinograms] = SimulateAttFactors(transMuMap,projectorType,pixelInfo)
% Given a set of mu-maps, this function calculates the probability of an average 511keV annihilation photon surviving
% from detector A to detector B through when traversing the field of view, which contains a distribution
% of optically-dense media parameterised per mm.
% NOTE: This function expects the mu-maps to be provided in PER MILLIMETRE scale

% Check variables are compatible:
if iscell(transMuMap)
    nToProject = numel(transMuMap);
else
    nToProject = size(transMuMap,4);
    for it = 1:nToProject
        temp{it} = transMuMap(:,:,:,it);
    end
    transMuMap = temp;
    clear temp;
end

% Very simple check for millimetres rather than centimetres
mx_mu = max(transMuMap{1}(:));

% per cm (most often used), cortical bone is approx 0.15 /cm. Thus, check for 0.015 /mm as a maximum value.
if mx_mu > 0.02
    error('Mu-maps appear to be in per-cm. Please ensure they are provided in per mm.');
else
    fprintf(' - Simulating attenuation factor sinogram(s)\n')
end

% The extent to which the attenuation will be modelled in this function depends on the type of projector
% - APIRL simulates its own AF sinograms
switch projectorType
    % - MATLAB and CECR projectors require manual AF modelling
    % MATLAB projectors:
    case 1
        totalMus = FwdProject(transMuMapprojectorType,pixelInfo);
        for it = 1:nToProject
            attFactorSinograms{it} = exp( -totalMus{it} );
        end
    % CECR projectors
    case 2
        totalMus = FwdProject(transMuMap,projectorType,pixelInfo);
        for it = 1:nToProject
            attFactorSinograms{it} = exp( -totalMus{it} );
        end
    % APIRL calculates ACF sinograms directly
    case 3
        InitialiseApirl;
        for it = 1:nToProject
            ACFsino = PET.ACF(transMuMap{it}); %%%!!!
            attFactorSinograms{it} = ACFsino;
            attFactorSinograms{it}(attFactorSinograms{it}~=0) = 1./attFactorSinograms{it}(attFactorSinograms{it}~=0);
        end
    otherwise
        error('Unrecognised projectorType provided');
end

end


