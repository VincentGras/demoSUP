function [dsq, dtq] = interp1dsignal(dt, ds, dtq)

% 1d signal interpolation adapted to RF pulse manipulation (FA preserving)
% author : Vincent Gras
% Contact vincent.gras@cea.fr


n = numel(dt);
T = sum(dt);
if (isscalar(dtq))
    nq = round(T/dtq);
    dtq = repmat(dtq, 1, nq);
    % adapt last dt to preserve pulse duration
    dtq(end) = T - sum(dtq(1:end-1));
end

%nq=numel(dtq);
tq = cumsum(dtq); % - 0.5 * dtq(1);
%tq(end)-sum(dt)
s = zeros(size(ds, 1), numel(tq));
dsq = zeros(size(ds, 1), numel(tq));

t = [0, cumsum(dt)];
t(end) = t(end)-2*eps;
[i1, i, i2] = pTXUtils.dichotomicsearch(t, tq);

s_ = cumsum(ds, 2);

for j = 1:numel(tq)
    
    k = i1(j);

    if (k > n)
        s(:,j) = s_(:, k-1);
    elseif (k > 1)
        s(:,j) = s_(:, k-1) + (tq(j)-t(k))/dt(k) * ds(:, k);
    else
        s(:,j) = (tq(j)-t(k))/dt(k) * ds(:, k);
    end
    
    if (j == 1)
        dsq(:,j) = s(:,j);
    else
        dsq(:,j) = s(:,j) - s(:,j-1);
    end
end

