%% function [Q,R] = cholqr(X)

function [Q,R] = cholqr(X)

   G = X'*X;
   R = chol(G);
   %R = R';
   Q = X/R;
   
% -------------------------------------------------------------------------
%  Copyright (2012, 2013)  Magnus Grandin <magnus.grandin@it.uu.se>
% -------------------------------------------------------------------------
