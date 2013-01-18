function [Q,R] = cholqr(X)

   G = X'*X;
   R = chol(G);
   %R = R';
   Q = X/R;