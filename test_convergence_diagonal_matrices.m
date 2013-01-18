% Test script for explicitly restarted CA-Lanczos, running a predefined 
% number of Lanczos iterations, tracking the convergence and orthogonality
% of the largest eigenpair.
%
% The matrices used here are diagonal, with equally spaced eigenvalues from
% 1..K, where K is the desired condition number of each matrix. The value
% of N is constant for all matrices generated.

N = 500;
niter = 480; % This number must be a multiple of 4,6,8 and 10.
basis = 'newton';
orthos = {'periodic'};

r = ones(N,1);
for orth = orthos
    for k = 2
        K = 10^k;
        a = linspace(1,K,N);
        A = sparse(diag(a));
        test_ca_lanczos(A,r,niter,orth,basis);
    end
end

% -------------------------------------------------------------------------
%  Copyright (2012, 2013)  Magnus Grandin <magnus.grandin@it.uu.se>
% -------------------------------------------------------------------------
