function [NRMSE, NRMSEs] = calcNRMSE(fa, target, dim, indx)

% Compute Normalized Root Means Squares Error along dimension dim
% if indx is provided, then slice the input vector fa into blocks and
% compute NRMSE within each block (case of multiple "subjects")
% author : Vincent Gras
% contact : vincent.gras@cea.fr

if (nargin < 3)
    dim = find(size(fa)>1, 1, 'first');    
end

if (nargin < 4)
    indx = shiftdim( ones(size(fa,dim), 1), 1-dim);
end

s = unique (indx);
NRMSEs = cell(size(s));

for j = 1:numel(s)
    
    subset = indx == s(j);
    fa_ = slice(fa, subset, dim);
    NRMSEs{j} = sqrt(sum( abs(fa_/target - 1).^2, dim) ) / sqrt(nnz(subset));
end

NRMSEs = cat(dim, NRMSEs{:});
NRMSE = mean(NRMSEs, dim);

function fa = slice(fa, subset, dim)

perm = 1:ndims(fa);
perm(dim) = 1;
perm(1) = dim;
fa = permute( fa, perm);
sz = size(fa);
fa = fa(subset(:), :);
sz(1) = nnz(subset);
fa = reshape(fa, sz);
fa = ipermute( fa, perm);

