function [AttFactorSinograms] = SimulateAttFactors(transMuMap,projectorType,pixelInfo)
%SIMULATEATTFACTORS [AttFactorSinograms] = SimulateAttFactors(transMuMap,projectorType,pixelInfo)
% Given a set of mu-maps, this function calculates the probability of an average 511keV annihilation photon surviving 
% from detector A to detector B through when traversing the field of view, which contains a distribution
% of optically-dense media parameterised per mm.
% NOTE: This function expects the mu-maps to be provided in PER MILLIMETRE scale

% Very simple check for millimetres rather than centimetres
if iscell(transEmMap)
	mx_mu = max(transMuMap{1}(:));
else
	mx_mu = max(transMuMap(:));
end
% per cm (most often used), cortical bone is approx 0.15 /cm. Thus, check for 0.015 /mm as a maximum value.
if mx_mu > 0.02
	error('Mu-maps appear to be in per-cm. Please ensure they are provided in per mm.');
else
	fprintf('== Simulating attenuation factor sinogram(s)\n')
end

% The extent to which the attenuation will be modelled in this function depends on the type of projector
% - APIRL simulates its own AF sinograms
switch projectorType
% - MATLAB and CECR projectors require manual modelling
	% MATLAB projectors:
	case 1
		for it = 1:nProjs
			totalMus = FwdProject(transMuMap{it},projectorType);
			AttFactorSinograms{it} = exp(-totalMus);
		end
	% CECR projectors
	case 2
		for it = 1:nProjs
			totalMus = FwdProject(transMuMap{it},projectorType);
			AttFactorSinograms{it} = exp(-totalMus);
		end
	end
	% APIRL calculates ACF sinograms directly
	case 3
		InitialiseApirl;
		for it = 1:nProjs
			ACFsino = PET.ACF(transMuMap{it}); %%%!!!
			AttFactorSinograms{it} = ACFsino;
			AttFactorSinograms{it}(AttFactorSinograms{it}~=0) = 1./AttFactorSinograms{it}(AttFactorSinograms{it}~=0);
		end
otherwise
	error('Unrecognised projectorType provided');
end

end


