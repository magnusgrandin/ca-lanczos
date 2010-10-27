%%
%% Input: a vector of n possibly complex numbers x.  If there are
%% complex numbers in x, they must occur in complex conjugate pairs.
%% (We do preprocessing to group complex conjugate pairs.)
%%
%% LAPACK's DHSEQR (upper Hessenberg eigenvalue solver) returns the
%% eigenvalues in the right order so that it shouldn't require such
%% preprocessing, but Matlab doesn't seem to do this all the time...
%%
%% Output: a modified Leja ordering of x.
%% This ordering is like the Leja ordering, but preserves the adjacency
%% and order of complex conjugate pairs.
%%
%% Note: the Leja ordering of a set of points is not unique, so we cannot
%% expect the modified Leja ordering to be unique either.
%%
function [y,idx] = real_leja (x)

    % Make sure x is a row vector.  Don't take x = x', because if there
    % are complex shifts, the ' will disturb their order.
    if (size(x,1) > size(x,2))
        if (size(x,2) ~= 1)
            error ('x must be a vector');
        end
        % x is a column vector
        n = size(x,1);
        temp = x;
        x = zeros(1,n);
        for k = 1:n
            x(k) = temp(k);
        end
        clear temp;
    else
        if (size(x,1) ~= 1)
            error ('x must be a vector');
        end
        % x is a row vector
        n = size(x,2);
    end
    
    % Count repeated entries in x.  Repeated entries require special
    % handling in the Leja ordering.
    [y, mults, num_unique] = count_multiplicities(x, n);

    % Ensure that complex conjugate pairs always occur adjacently, with
    % the element of the pair with positive real part always occurring
    % first.  We do this as follows:
    %
    % 1. Sort the elements of y according to their real part.  Matlab
    %    uses a stable sort, so this won't mess up anything from before.
    % 2. Apply the sort order to y and multiplicities from above.
    % 3. Go through y in order, checking for complex conjugate pairs
    %    (identical real parts, and imaginary parts same except for
    %    opposite sign).  Swap so the positive-imaginary-part element of
    %    each pair occurs first.
    %
    % What if there are repeated complex conjugate pairs?  There won't
    % be, because "unique" already ate them!

    [tmp_val, tmp_idx] = sort(real(y));
    clear tmp_val;
    y = y(tmp_idx);
    mults = mults(tmp_idx);
    k = 1;
    num_unique
    while k < num_unique
        if (imag(y(k)) ~= 0)
            if (real(y(k)) == real(y(k+1)) && imag(y(k)) == -imag(y(k+1)))
                
                y(k) = real(y(k)) + i * abs(imag(y(k)));
                y(k+1) = real(y(k)) - i * abs(imag(y(k)));
              
                k = k + 2;
            end
        else
            k = k + 1;
        end
    end

    % FIXME (mfh 29,30 Sep 2009) If x contains nonunique entries, then the
    % "idx" output doesn't make sense.  It indexes into the uniquified
    % shifts, not the input shifts.  I'll just return it for now...
    [y, idx] = modified_leja(y, n, mults);
end
