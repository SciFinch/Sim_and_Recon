function [normFactorSinogram] = SimulateNormFactors(projectorType,pixelInfo,datainfo)
%SIMULATENORMFACTORS [normFactorSinogram] = SimulateNormFactors(projectorType,pixelInfo,datainfo)
%   This function SIMULATES normalisation factors that may affect data acquisition. Currently
%   this only includes data-based NF sims, or 100% efficient detectors. This will eventually include
%   simulation of the individual detector efficiencies and direct calc'ing of the NF sinio
%
%	Check whether norm is being simulated either
% 	(a) from data
% 	(b) from random number generation for each detector
% 	(c) with 100% detector eficiency
%


% check whether normalisation factors will be derived from data
if exist('datainfo.simNorm_pn','var')
    if projectorType == 3
        InitialiseApirl;
        normFactorSinogram = PET.NCF();
    else
        error('Norm sim from file has been opted (naively), but current projector is not APIRL')
    end
    % Otherwise calculate efficiencies and simulate these
else
    warning('TO DO: Sim NFs is currently only returning 100% efficient NFs');
    % If somehow this function knows that it should calculate detector efficiencies,
    %normFactorSinogram = SimulateDetectorEfficiencies(projectorType);
    switch projectorType
        % - MATLAB and CECR projectors require manual modelling
        % MATLAB projectors:
        case 1
            normImage = ones(pixelInfo.pxSize);
            normFactorSinogram = FwdProject(normImage,projectorType,pixelInfo);
        % CECR projectors
        case 2
            normImage = ones(pixelInfo.pxSize);
            normFactorSinogram = FwdProject(normImage,projectorType,pixelInfo);
        % APIRL NF method currently not used
        case 3
            %InitialiseApirl;
            warning('Calculating NF sino for APIRL without data not currently implemented')
            normImage = ones(pixelInfo.pxSize);
            normFactorSinogram = FwdProject(normImage,projectorType,pixelInfo);            
        otherwise
            error('Unrecognised projectorType provided');
    end
    
    
end
