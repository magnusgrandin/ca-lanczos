%% function [Q,R] = tsqr(A)
%
%   Compute a QR factorization of matrix A, and shift
%   the signs of the resulting factorization matrices such that R 
%   only has positive values on its diagonals.

function [Q,R] = tsqr(A)
   [Q,R] = qr(A,0);
   d = sign(diag(R));
   R = diag(d)*R;
   Q = Q*diag(d);
end

% -------------------------------------------------------------------------
%  Copyright (2012, 2013)  Magnus Grandin <magnus.grandin@it.uu.se>
% -------------------------------------------------------------------------
