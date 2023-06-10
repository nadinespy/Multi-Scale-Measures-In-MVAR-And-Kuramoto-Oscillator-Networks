function [ phi ] = PhiGMatlab(lS, scan_method)

dim = size(lS,1)/2;
S = lS(1:dim,1:dim);
G = lS(1:dim, (1+dim):end);

if mod(dim,2) ~= 0
  error('Odd dimensions not supported');
end

s = dec2bin(1:(2^(dim-1) - 1), dim);
allBipartitions = cellfun(@str2num, mat2cell(s, ones([1,size(s,1)]), ones([1, size(s,2)])));

if strcmp(scan_method, 'atomic')
  phi = phi_G_Gauss(S, G, S, 1:dim);
elseif strcmp(scan_method, 'bipartitions')
  phi = min(arrayfun(@(i) phi_G_Gauss(S, G, S, allBipartitions(i,:)+1), 1:size(allBipartitions,1)));
elseif strcmp(scan_method, 'even_bipartitions')
  evenBipartitions = allBipartitions(sum(allBipartitions,2) == dim/2,:);
  phi = min(arrayfun(@(i) phi_G_Gauss(S, G, S, allBipartitions(i,:)+1), 1:size(allBipartitions,1)));
end

