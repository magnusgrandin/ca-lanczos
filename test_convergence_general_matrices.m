% Test script for explicitly restarted CA-Lanczos, running a predefined 
% number of Lanczos iterations, tracking the convergence and orthogonality
% of the largest eigenpair.
%
function test_convergence_general_matrices(matrix,niter,basis,orthos)

% Check input
if nargin < 2
    disp('-- ''niter'' defaults to 120.')
    niter = 120;
end
if nargin < 3
    disp('-- ''basis'' defaults to ''newton''.');
    basis = 'newton';
end
if nargin < 4
    disp('-- ''orthos'' defaults to {''local'',''full'',''selective'',''periodic''}.');
    orthos = {'local','full','selective','periodic'};
end

disp(['Matrix: ', matrix]);
matrix_path_prefix = '../matrices';
load ([matrix_path_prefix,'/',matrix]);

A = Problem.A;
N = size(A,1);
r = ones(N,1);
for orth = orthos
    test_ca_lanczos(A,r,niter,orth,basis);
end

% -------------------------------------------------------------------------
%  Copyright (2012, 2013)  Magnus Grandin <magnus.grandin@it.uu.se>
% -------------------------------------------------------------------------
