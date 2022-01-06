function [p, q, r] = dichotomicsearch(arr, x, varargin)

% Computes the arrays (p, q, r) such that for each i, arr(q(i)) is the best approximation of x(i)
%  and [arr(p(i)), arr(r(i))] provides a bounding for x(i)
% Usage :
%    [q, p, r] = dichotomicsearch(t, x)
% Varargin-parameters:
%    sorted.t     [false]    indicate that t is sorted
%    sorted.x     [false]    indicate that x is sorted
% Note :
% - If x(i) = NaN, q(i) is set to NaN (hence, in such a case, t(q) is
%   an invalid statement)
% - If x(i) < min(t), then p(i) = 0; similarly if x(i) > max(t), r(i) =
%   numel(t)+1
% Author : Vincent Gras
% Contact : vincent.gras@cea.fr

assert (~any(isnan(arr)), 'arr contains NaNs !');
assert (numel(arr)>1, 'numel(arr) < 2');
assert(ndims(arr) <= 2 && min(size(arr)) == 1, 'arr is multidimensional');
 
% Number of elements in arr
arrLength = length(arr);
arrMin = arr(1);
arrMax = arr(end);

% Number of elements in x
xLength = length(x);

% Initialisation of q, p, r
p = ones(size(x));
q = ones(size(x));
r = ones(size(x));

o = isnan(x(:));

% handle the case where x contains NaNs
if (any(o))
    
    p(o) = NaN;
    q(o) = NaN;
    r(o) = NaN;
    [p(~o),q(~o),r(~o)] = dichotomicsearch(arr, x(~o));
    
    return;
    
end

% For all first indexes i such that x(i) <= t(1), do q(i) = r(i) = 1

i = find (x >= arrMin, 1, 'first');

if (isempty(i))
    i = 1;
end
  
% Start dichotomy

a = 1;
      
while (i <= xLength && x(i) < arrMax) 
   
   % index 'a' doesn't need to be set to 1 because we know that x(i) >= t(a) (x is sorted)
   
   b = arrLength;  
       
   while (true)
           
       % we have a < b
       
       if (b - a <= 1)
                         
          break;
        
      else
          c = floor((a + b) / 2);    
          if (x(i) <= arr(c))
            b = c;
          else
            a = c;
          end
      end
   end
   
   p(i) = a;
   r(i) = b;
 
   if (abs(x(i) - arr(a)) <= abs(x(i) - arr(b)))
       q(i) = a;
   else
       q(i) = b;
   end
      
   i = i + 1;
     
end  
 
% Here, i is such that x(j) >= arrMax for all j >= i

if (i <= xLength)
    p(i:xLength) = arrLength;
    q(i:xLength) = arrLength;
    r(i:xLength) = arrLength;
end

% Update p and r for the elements x(i) which 
% are not in the interval [min(t), max(t)]

p(x < arrMin) = 0;
r(x < arrMin) = 1;
p(x > arrMax) = arrLength;
r(x > arrMax) = arrLength + 1;

