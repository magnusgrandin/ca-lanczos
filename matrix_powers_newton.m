% Return s-step Newton basis with s shifts lambda(1:s) for the given matrix
% A and starting vector v.  Basis contains s+1 vectors.
%
% A: sparse matrix or other data structure suitable for the first
%    argument of the SpMV function.
% v: starting vector for the s-step basis
% s: s-step basis length (returns s+1 vectors)
% lambda: s approximate eigenvalues to be used as shifts
% modifiedp: 1 if the modified Leja ordering and corresponding 
%            modified Newton basis are to be used, else 0.  
%            Defaults to zero.
%
function V = matrix_powers_newton(A, v, s, lambda, modifiedp)
    if (nargin < 5)
        modifiedp = 0;
    end
   

    n = length(v);
    V = zeros(n, s+1);
    V(:, 1) = v;

    if (modifiedp == 0)
        for k = 1:s
            w = SpMV(A, V(:, k));
            V(:, k+1) = w - lambda(k) * V(:, k);
        end
    else
        for k = 1:s
            w = SpMV(A, V(:, k));
            if (imag(lambda(k)) > 0)
                V(:, k+1) = w - real(lambda(k)) * V(:, k);
            elseif (imag(lambda(k)) < 0)
                if (k == 1)
                    error(['k==1, but shift %e has a negative imaginary ' ...
                           'part'], lambda(k));
                end
                V(:, k+1) = w - real(lambda(k)) * V(:, k) + ...
                    imag(lambda(k))^2 * V(:, k-1);
            elseif (isreal(lambda(k)) || imag(lambda(k)) == 0)
                V(:, k+1) = w - real(lambda(k)) * V(:, k);
            else
                error('Should never get here');
            end
        end
    end

    if (0)
        r = rank(V);
        disp(sprintf('Numerical rank of basis (should be %d): %e', s+1, r));
    end
end % function newton_basis
