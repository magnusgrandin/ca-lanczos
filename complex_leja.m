%%
%% Given a vector of complex numbers x, compute a Leja ordering of x:
%% return a vector y and a permutation array idx such that x(idx) = y
%% is the Leja ordering of x.
%%
%% Note: the Leja ordering of a set of points is not unique:  see e.g.
%% Baglama, Calvetti and Reichel:  "Fast Leja points", ETNA, Vol. 7,
%% 1998, pp. 124-140.
%%
function [y, idx] = complex_leja (x)

    function [y,idx] = complex_leja_no_multiplicities(x, n)
        y = zeros(n,1);
        y(1:n) = x(1:n);
        idx = 1:n;
        yprod = abs(y(1:n));

        [max_val, max_idx] = max(yprod(1:n));
        tmp_val = y(1);
        tmp_idx = idx(1);
        y(1) = y(max_idx);
        idx(1) = idx(max_idx);
        y(max_idx) = tmp_val;
        idx(max_idx) = tmp_idx;

        for k = 2:n
            for cur = k:n
                yprod(cur) = prod(abs(y(cur) - y(1:k-1)));
            end
            % max_idx: index of the y value to pick
            [max_val,max_idx] = max(yprod(k:n)); 
            if (max_val == 0)
                ME = MException('leja:multiple', ...
                                'Multiple shifts require special handling');
                throw(ME);
            end
            tmp_val = y(k);
            tmp_idx = idx(k);
            y(k) = y(max_idx);
            idx(k) = idx(max_idx);
            y(max_idx) = tmp_val;
            idx(max_idx) = tmp_idx;
        end
    end

    % Make sure x is a column vector.
    if (size(x,1) < size(x,2))
        x = x';
    end
    [n,n2] = size(x);
    if (n2 ~= 1)
        error ('x must be a vector');
    end
    clear n2;

    % Repeated values in x require special handling, which we do not do.
    [y,idx] = complex_leja_no_multiplicities(x, n);
end

% -------------------------------------------------------------------------
%  Copyright (2012, 2013)  Magnus Grandin <magnus.grandin@it.uu.se>
% -------------------------------------------------------------------------
