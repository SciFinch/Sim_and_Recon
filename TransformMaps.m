function [transEmMap,transMuMap] = TransformMaps(emMap,muMap,pxinfo,dHF,dAP,dRL,interpType,padVal)
%TRANSFORMMAPS [transEmMap,transMuMap] = TransformMaps(emMap,muMap,pxinfo,dHF,dAP,dRL,interpType,padVal);
%   This function will return two cell arrays, consisting of the provided
%   emMap and muMap transformed according to the motion fields provided.
%   This requires specification of the interpolation method and the values
%   used to pad the outside of the array.

if iscell(emMap)
	if ~(numel(emMap) == numel(muMap))
		error('Number of emission maps does not equal number of mu-maps');
	end
end

[ transEmMap ] = TransformImage(emMap,pxinfo,dHF,dAP,dRL,interpType,padVal);
[ transMuMap ] = TransformImage(muMap,pxinfo,dHF,dAP,dRL,interpType,padVal);

end

