function data = applyMask(Mask, data, fillval)

% apply Mask to a 3-d image or a stack of 3-d images
% author : Vincent Gras
% contact : vincent.gras@cea.fr


if (nargin < 3)
    fillval = NaN;
end



if (~islogical(Mask))
    data = applyMask(Mask>0, data, fillval);
else
    
    
    nvol = numel(data)/numel(Mask);
    n = max(ndims(data), ndims(Mask)+1);
    s = repmat({':'}, 1, n-1);
    
    assert(rem(nvol,1)==0, 'Mask size is not compatible with data size');
    
    for i = 1:nvol
        
        s_i = [s, i];        
        data_ = data(s_i{:});
        data_(~Mask) = fillval;
        data(s_i{:}) = data_;        
        
    end
end


