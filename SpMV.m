% Temporary fix. At the moment this function just computes the product of a
% matrix A and a vector v, but should be extended to support other data
% structures and function references as well.
function Av = SpMV(A,v)

Av = A*v;

% -------------------------------------------------------------------------
%  Copyright (2012, 2013)  Magnus Grandin <magnus.grandin@it.uu.se>
% -------------------------------------------------------------------------
