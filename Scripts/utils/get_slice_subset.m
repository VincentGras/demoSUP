function subset = get_slice_subset(nslc, nslctot, s)

% pick a subset of all slices
% author : Vincent Gras
% contact : vincent.gras@cea.fr


if (nargin <3)
    s = 1:nslctot;
end

j = mod(s,  round(numel(s)/nslc))<1;
j = find(j);
j = ceil(j - mean(j) + mean(s));
subset = false(1, nslctot);
subset(j) = true;
