%% function V = matrix_powers_monomial(A,q,s)
% 
%   Compute s matrix-vector multiplications of A and q using 
%   the monomial basis

function V = matrix_powers_monomial(A,q,s)
    V = zeros(length(q),s);
    V(:,1) = A*q;
    for i = 2:s
        V(:,i) = A*V(:,i-1);
    end
end

% -------------------------------------------------------------------------
%  Copyright (2012, 2013)  Magnus Grandin <magnus.grandin@it.uu.se>
% -------------------------------------------------------------------------
