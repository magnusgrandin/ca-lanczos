% Compute the modified Leja ordering of the n shifts in x.
%
% Inputs:
% x: Row vector of shifts for the Newton Krylov basis. 
%    Elements of x must be unique, but may occur with multiplicities
%    specified by the "mults" parameter (see below).  Any complex
%    elements of x must occur in a complex conjugate pair, with the
%    elements of the pair adjacent, and the member of the pair with 
%    positive imaginary part occurring first.
% n: Number of (unique) shifts (length of x)
% mults: mults(j) is the number of occurrences of x(j) in the original
%        collection of shifts (which this function does not see).
%
% Outputs:
% y: All the elements of x, ordered in the modified Leja ordering.
% outidx: index array such that x(outidx) == y.
%
% NOTE: x' means the conjugate transpose, which is NOT what you want, as
% it messes up the order of the complex conjugate pairs in x.  The same
% thing goes for y; do NOT take y'.
function [y, outidx] = modified_leja(x, n, mults)

    function bool = is_conj_pair(a, b)
        bool = 0;
        % Note: this test rejects real a and b (a == b), even though a
        % and b could be called a (trivial) complex conjugate pair in
        % that case.
        if (real(a) == real(b) && imag(a) == -imag(b) && imag(a) ~= 0)
            bool = 1;
%         elseif (abs(real(a) - real(b)) <= eps && ...
%                 abs(imag(a) + imag(b)) <= eps && ...
%                 imag(a) ~= 0 && imag(b) ~= 0)
%             % We accept "approximate conjugates" also, just in case...
%             bool = 1;
        end
    end

    function [y, outidx] = modified_leja_start(x, n)
        if (n < 1)
            y = [];
            outidx = [];
        elseif (n == 1)
            y = x(1);
            outidx = [1];
        else
            [val, j] = max(abs(x(1:n)));
            if (imag(x(j)) == 0) % x(j) is real
                y = [x(j)];
                outidx = [j];
            elseif (j > 1 && is_conj_pair(x(j-1), x(j)))
                if (imag(x(j-1)) < 0)
                    error('Complex conjugate pair out of order at indices %d and %d', ...
                          j-1, j);
                end
                y = [x(j-1), x(j)];
                outidx = [j-1, j];
            elseif (j < n && is_conj_pair(x(j), x(j+1)))
                if (imag(x(j)) < 0)
                    %error('Complex conjugate pair out of order at indices %d and %d', ...
                    %      imag(x(j)), imag(x(j+1)));
                    x(j) = real(x(j));
                    x(j+1) = real(x(j+1));
                    
                end
                y = [x(j), x(j+1)];
                outidx = [j, j+1];
            else % error-handling case
                if (j == 1)
                    error('Complex shift, not in a pair, occurs at beginning of input');
                elseif (j == n)
                    error('Complex shift, not in a pair, occurs at end of input');
                end
            end
        end
    end

    function [y, outidx, capacity] = modified_leja_helper(x, n, mults, ...
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
                    capacity = prod( abs(y_last - x(outidx(1:num_points-1))) ...
                                     .^ (mults(outidx(1:num_points-1)) * ...
                                         (1.0/num_points)) );
                    % Unscale all the shifts by the old capacity
                    % estimate, and rescale by the new capacity estimate.
                    % We do this all at once to avoid overflow on
                    % unscaling, if the old capacity estimate is large.
                    %
                    % Hopefully complex division by a positive real
                    % scalar gets the same answer, whether the 
                    % imaginary part is positive or negative --
                    % otherwise finding complex conjugate pairs will 
                    % need to be "approximate."
                    x = x ./ (capacity / old_capacity);
                    y = y ./ (capacity / old_capacity);

                    disp(sprintf('New capacity estimate: %e', capacity));
                end
            end
            zprod = zeros(n,1);
            count = 0;
            for j = inidx
                count = count + 1;
                % Scale each term in the product by the estimated
                % capacity, in order to avoid overflow (which was
                % observed in practice).
                zprod(count) = ...
                    prod( (abs(x(j) - x(outidx)) ./ capacity) .^ mults(outidx) );
            end
            [max_zprod, k] = max( zprod(1:count) );
            disp(sprintf('New max(zprod): %e', max_zprod));
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
            if (imag(x(j)) == 0) % x(j) is real
                inidx = setdiff(inidx, [j]);
                outidx = [outidx, j];
                y = [y, x(j)];
                new_num_points = num_points + 1;
            elseif (j > 1 && is_conj_pair(x(j-1), x(j)))
                if (imag(x(j-1)) < 0)
                    error('Complex conjugate pair out of order at indices %d and %d', ...
                         j-1, j);
                 
                end
                inidx = setdiff(inidx, [j-1, j]);
                outidx = [outidx, j-1, j];
                y = [y, x(j-1), x(j)];
                new_num_points = num_points + 2;
            elseif (j < n && is_conj_pair(x(j), x(j+1)))
                if (imag(x(j)) < 0)
                    error('Complex conjugate pair out of order at indices %d and %d', ...
                          j, j+1);
                    
                   
                end
                inidx = setdiff(inidx, [j, j+1]);
                outidx = [outidx, j, j+1];
                y = [y, x(j), x(j+1)];
                new_num_points = num_points + 2;
            else % error-handling case
                if (j == 1)
                    error('Complex shift, not in a pair, occurs at beginning of input');
                elseif (j == n)
                    error('Complex shift, not in a pair, occurs at end of input');
                end
            end
            [y, outidx, capacity] = ...
                modified_leja_helper( x, n, mults, inidx, outidx, y, ...
                                      new_num_points, capacity );
        end
    end

    disp(sprintf('Shifts before Leja ordering:'));
    x
        
    [y, outidx] = modified_leja_start(x, n);
    inidx = setdiff(1:n, outidx);
    [y, outidx, capacity] = modified_leja_helper(x, n, mults, inidx, outidx, ...
                                                 y, length(outidx));

    % Unscale the computed modified Leja points by the capacity estimate.
    y = y .* capacity;

    disp(sprintf('Shifts after Leja ordering:'));
    y
end
