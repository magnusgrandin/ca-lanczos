%% function [y, outidx] = nonmodified_leja(x, n, mults)
%
%   Compute the (nonmodified) Leja ordering of the n shifts in x.
%
%   Inputs:
%   x: Row vector of shifts for the Newton Krylov basis. 
%      Elements of x must be unique, but may occur with multiplicities
%      specified by the "mults" parameter (see below).  Elements of x 
%      may be complex and may occur in any order.
%   n: Number of (unique) shifts (length of x)
%   mults: mults(j) is the number of occurrences of x(j) in the original
%          collection of shifts (which this function does not see).
% 
%   Outputs:
%   y: All the elements of x, ordered in the (nonmodified) Leja ordering.
%   outidx: index array such that x(outidx) == y.
%
%   NOTE: The order of the complex shifts matters, not for the input here
%   (as it does for the modified Leja ordering), but definitely for the
%   output here.  That means, do NOT take y', because that computes the
%   conjugate transpose.
function [y, outidx] = nonmodified_leja(x, n, mults)

    function [y, outidx] = leja_start(x, n)
        if (n < 1)
            y = [];
            outidx = [];
        elseif (n == 1)
            y = x(1);
            outidx = [1];
        else
            [val, j] = max(abs(x(1:n)));
            y = [x(j)];
            outidx = [j];
        end
    end

    function [y, outidx, capacity] = leja_helper(x, n, mults, ...
                                                 inidx, outidx, y, ...
                                                 num_points, ...
                                                 capacity)
        if (isempty(inidx))
            if (nargin < 8)
                capacity = 1;
            end
            return;
        else
            % Later we will update "capacity" to be a running estimate of
            % the capacity of the region for which we are finding Leja points.
            if (nargin < 8)
                capacity = 1;
            else
                if (num_points > 1)
                    old_capacity = capacity;

                    % Update the capacity estimate.
                    y_last = y(num_points);
                    for i = 1:num_points-1
                        capacity = prod( abs(y_last - x(outidx(i))) ...
                                     .^ (mults(outidx(i)) * ...
                                         (1.0/num_points)) );
                    end
                    % Unscale all the shifts by the old capacity
                    % estimate, and rescale by the new capacity estimate.
                    % We do this all at once to avoid overflow on
                    % unscaling, if the old capacity estimate is large.
                    x = x ./ (capacity / old_capacity);
                    y = y ./ (capacity / old_capacity);

                    %disp(sprintf('New capacity estimate: %e', capacity));
                end
            end
            zprod = zeros(n,1);
            count = 0;
            for j = inidx
                count = count + 1;
                % Scale each term in the product by the estimated
                % capacity, in order to avoid overflow (which was
                % observed in practice when no scaling was done, for the
                % modified Leja ordering at least).
                for i = 1:num_points
                    zprod(count) = ...
                        prod( (abs(x(j) - x(outidx(i)) ./ capacity) .^ mults(outidx(i))) );
                end
            end
            [max_zprod, k] = max( zprod(1:count) );
            %disp(sprintf('New max(zprod): %e', max_zprod));
            j = inidx(k);
            % A zero maximum indicates either that there are multiple
            % shifts, or that the product underflowed.  An Inf maximum
            % indicates that the product overflowed.  Both of these
            % situations can be fixed by appropriate scaling; I'm not
            % implementing that until I see it's necessary.
            if (max_zprod == 0)
                error(['Product to maximize is zero; either there are ' ...
                       'multiple shifts, or the product underflowed']);
            elseif (max_zprod == Inf)
                error(['Product to maximize is Inf; must have ' ...
                       'overflowed']);
            end
            inidx = setdiff(inidx, [j]);
            outidx = [outidx, j];
            y = [y, x(j)];
            new_num_points = num_points + 1;
            [y, outidx, capacity] = leja_helper( x, n, mults, inidx, outidx, ...
                                                 y, new_num_points, capacity );
        end
    end

    %disp(sprintf('Shifts before Leja ordering:'));
    %x
        
    [y, outidx] = leja_start(x, n);
    inidx = setdiff(1:n, outidx);
    [y, outidx, capacity] = leja_helper(x, n, mults, inidx, outidx, y, ...
                                        length(outidx));

    % Unscale the computed Leja points by the capacity estimate.
    y = y .* capacity;

    %disp(sprintf('Shifts after Leja ordering:'));
    %y
end

% -------------------------------------------------------------------------
%  Contribution by Erin Carson, University of California at Berkeley
% -------------------------------------------------------------------------
