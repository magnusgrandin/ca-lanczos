%%% Compute the (s+1) x s matrix $\underline{B}$ such that for the
%%% Newton basis $V_{s+1}$, $A V_s = V_{s+1} \underline{B}$.
%%%
%%% lambda:  Newton basis parameters
%%% s:       basis size 
%%% modifiedp: 1 if the modified Leja ordering and corresponding 
%%%            modified Newton basis are to be used, else 0.  
%%%            Defaults to zero.
%%% B_:      $\underline{B}$
%%%
function B_ = newton_basis_matrix (lambda, s, modifiedp)
    if (nargin < 3)
        modifiedp = 0;
    end
  
    B_ = zeros (s+1, s);
    if (modifiedp == 0)
        for k = 1 : s
            B_(k,k) = lambda(k);
            B_(k+1,k) = 1.0; % sigma(k) / sigma(k+1); % scaling factors
        end
    else
        for k = 1 : s
            shift = lambda(k);
            if (imag(shift) > 0)
                if (lambda(k) ~= conj(lambda(k+1)))
                    lambda
                    error('Modified Leja ordering broken at k = %d, %d, by values %e, %e', ...
                          k, k+1, lambda(k), lambda(k+1));
                elseif (k == s)
                    lambda
                    if (imag(shift) ~= 0)
                        error(['Complex shift %e occurs at end (k = %d) of shifts (s = %d), ' ...
                               'without its conjugate'], shift, s, k);
                    end
                end
                B_(k,k) = real(shift);
            elseif (imag(shift) < 0)
                if (k == 1)
                    lambda
                    error (['newton_basis_matrix: imaginary part is negative ' ...
                            'for k = 1, which should be impossible if you ' ...
                            'computed the modified Leja ordering ' ...
                            'correctly']);
                elseif (lambda(k-1) ~= conj(lambda(k)))
                    lambda
                    error('Modified Leja ordering broken at k = %d, %d', ...
                          k-1, k);
                end
                B_(k,k) = real(shift);
                B_(k-1,k) = -imag(shift)^2;
            else
                B_(k,k) = shift;
            end
            B_(k+1,k) = 1.0;
        end
    end
end

