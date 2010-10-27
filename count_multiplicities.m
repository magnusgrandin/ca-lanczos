% x: length n vector of numbers
% n: number of elements in x
function [y, multiplicities, num_unique] = count_multiplicities(x, n)

    % Matlab's "unique", just like Matlab's "sort", sorts complex numbers by
    % their absolute values.  "unique" sorts by default.  ii indexes the
    % first occurrence of each unique value in x.  y == x(ii) and x == y(jj).
    [y, ii, jj] = unique(x, 'first');

    % Array of multiplicies (number of occurrences) of the values in x.  
    % This is indexed by the indices of y, so that multiplicies(k) gives
    % the number of times that y(k) occurs in x.
    num_unique = length(y);

    % Short circuit branch
    if (num_unique == n)
        multiplicities = ones(1,n);
        return;
    end

    % The "if 0" branch is broken when x is not sorted.
    if 0
        multiplicities = zeros(num_unique, 1);
        for k = 1 : num_unique-1
            multiplicities(k) = ii(k+1) - ii(k);
        end
        multiplicities(num_unique) = n - ii(num_unique) + 1;
    else
        % This breaks ii and jj ...
        [val, idx] = sort(x);
        [y, ii, jj] = unique(val, 'first');

        multiplicities = zeros(num_unique, 1);
        for k = 1 : num_unique-1
            multiplicities(k) = ii(k+1) - ii(k);
        end
        multiplicities(num_unique) = n - ii(num_unique) + 1;
    end
end    
    

