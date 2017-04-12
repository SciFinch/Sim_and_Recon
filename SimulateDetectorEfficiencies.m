function [normFactorSinogram] = SimulateDetectorEfficiencies(projectorType)
% SIMULATEDETECTOREFFICIENCIES [normFactorSinogram] = SimulateDetectorEfficiencies(projectorType)
% This function calculates the efficiency of each detector ad calculates the normalisation
% factor for each detector pair present in the span-11 sinogram

warning('SimulateDetectorEfficiencies is a stub');
normFactorSinogram = ones([344 252 837]);

minEfficiency = 0.8;
currEfficiency = minEfficiency + randn()*(1 - minEfficiency));
%detctorList = LoadDetectorInfo_span11();

% Something here that generates all possible detectors in the sinogram, then sims an efficiency
% rating for each one. The NF is then the product of these (might need to hard-code something
% like the minimum acceptable efficiency (like 80%), which should be derived from real data)


end