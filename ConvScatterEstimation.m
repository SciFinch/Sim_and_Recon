function [ scatterEstimate ] = ConvScatterEstimation(trueSino,attFactorSino,randSino,mthdStr)
% CONVSCATTERESTIMATION [scatterEstimate ] = ConvScatterEstimation(promptSino,attFactorSino,randSino,mthdStr)
%	Provides a convolution-based scatter estimate based on various methods,
%	selected using mthdStr
%
%	Current methods:
%	- Bailey93 - early convolutional paper, uses spatially-invariant 1/exp(-ar) filter
%	- Bentourkia99 - uses a spatially-dependent filter (not yet implemented)

nDirPlanes = 127;

switch mthdStr
	case 'Bailey93'
		% This is the convolution-subtraction method suggested by Bailey & Meikle in 1993
		% Uses a 1/exp(-ar) kernel, independent of axial position
		
		% Determine source data and find direct-plane sinograms
		unscatteredSino = CompressSino(trueSino.*attFactorSino + randSino); % Bailey's model (99% certain...)
		warning('Correct input for Bailey93 scatterConv estimation not checked');

		% Main parameters (defaults: alpha = 0.081, scatterFrac = 0.5*totCounts)
		alpha = 0.081;
		% scatterFrac = 0.5; %not necessary here
		nIterations = 4;

		% Generate the 2D scatter kernel (viewgram)
		xvec = linspace(-size(unscatteredSino,1)/2,size(unscatteredSino,1)/2,size(unscatteredSino,1));
		yvec = linspace(-size(unscatteredSino,3)/2,size(unscatteredSino,3)/2,size(unscatteredSino,3));
		[x,y] = ndgrid(xvec,yvec);
		clear yvec xvec;
		r = sqrt(x.^2 + y.^2);
		scatterKernel = exp(-alpha*r);
		
		% Generate scatter estimates
		scatterEstimate_dirPlane = zeros(size(unscatteredSino));
		for it = 1:size(dirPrompts,2)
		    currView = squeeze(unscatteredSino(:,it,:));
		    currScatterEst = currView; %0th (initial) estimate
		    for kk = 1:nIterations
		    	% update estimate
		    	% Formula: g_u' = g_0 - k(conv2(g_u,scatterKernel,'same'))
		    	currScatterEst = conv2(currScatterEst,scatterKernel,'same');
		    end
		    scatterEstimate_dirPlane(:,it,:) = currScatterEst;
		end

		scatterEstimate_dirPlane = scatterEstimate_dirPlane/max(scatterEstimate_dirPlane(:))*max(unscatteredSino(:));
%		scatterEstimate = CalcSpan11Sino(scatterEstimate_dirPlane);

		% STUB:
		scatterEstimate = zeros(size(trueSino));
		scatterEstimate(:,:,1:nDirPlanes) = scatterEstimate_dirPlane;
		warning('The end of ConvScatterEstimation is stubbed: CalcSpan11Sino has not been written yet')

	case 'Bentourkia99'
		error('The Bentourkia99 convScatter estimation method has not been implemented yet');
	otherwise
		error('Unrecognised convScatter Estimator requested');
end

