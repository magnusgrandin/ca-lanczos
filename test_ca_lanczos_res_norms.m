%% function test_ca_lanczos_res_norms(matrix, steps, orth, basis)

function test_ca_lanczos_res_norms(matrix, steps, orth, basis)

% Load matrix
load(matrix);
Problem
A=Problem.A;
n=size(A,1);

% Total number of Lanczos steps (maximum)
if nargin >= 2 && isempty(steps)
    steps = 0;
end
if nargin >= 2 && steps ~= 0
    m=steps;
else
    m=120;
end

% Starting vector
r0=rand(n,1);

% Orhtogonality option
if nargin < 3 || ~exist('orth','var')
    orth = 'local';
end

% Basis option
if nargin < 4 || ~exist('basis','var')
    basis = 'newton';
end

% Standard Lanczos
time_st=cputime;
[T_st,V_st]=lanczos(A,r0,m,orth);
time_st=cputime-time_st;
T_st = T_st(1:m,1:m);
V_st = V_st(:,1:m);
res_norms(:,1) = compute_res_norms(A,V_st,T_st,0,40);

% CA-Lanczos, s = 4
time_ca_4=cputime;
[T_ca_4,V_ca_4]=ca_lanczos(A,r0,4,m/4,basis,orth);
time_ca_4=cputime-time_ca_4;
res_norms(:,2) = compute_res_norms(A,V_ca_4,T_ca_4,0,40);

% CA-Lanczos, s = 6
time_ca_6=cputime;
[T_ca_6,V_ca_6]=ca_lanczos(A,r0,6,m/6,basis,orth);
time_ca_6=cputime-time_ca_6;
res_norms(:,3) = compute_res_norms(A,V_ca_6,T_ca_6,0,40);

% CA-Lanczos, s = 8
time_ca_8=cputime;
[T_ca_8,V_ca_8]=ca_lanczos(A,r0,8,m/8,basis,orth);
time_ca_8=cputime-time_ca_8;
res_norms(:,4) = compute_res_norms(A,V_ca_8,T_ca_8,0,40);

% CA-Lanczos, s = 10
time_ca_10=cputime;
[T_ca_10,V_ca_10]=ca_lanczos(A,r0,10,m/10,basis,orth);
time_ca_10=cputime-time_ca_10;
res_norms(:,5) = compute_res_norms(A,V_ca_10,T_ca_10,0,40);

figure; semilogy(res_norms,'*-','MarkerSize',8);
legend('Lanczos','CA-Lanczos(4)','CA-Lanczos(6)','CA-Lanczos(8)','CA-Lanczos(10)','CA-Lanczos(20)')


end


function [res_norm] = compute_res_norms(A,Q,T,num_smallest,num_largest)
    m = size(T,1);
    res_norm = zeros(num_smallest+num_largest,1);
    
    [Vp,Dp] = eig(T);
    [Dp,pix] = sort(diag(Dp));
    
    nnorms = 1;
    for i = 1:num_smallest
        d = Dp(i);
        x = Q*Vp(:,pix(i));
        res_norm(nnorms) = norm(A*x-d*x)/norm(d*x);
        nnorms = nnorms+1;
    end
    for i = 1:num_largest
        d = Dp(m-num_largest+i);
        x = Q*Vp(:,pix(m-num_largest+i));
        res_norm(nnorms) = norm(A*x-d*x)/norm(d*x);
        nnorms = nnorms+1;
    end
end

% -------------------------------------------------------------------------
%  Copyright (2012, 2013)  Magnus Grandin <magnus.grandin@it.uu.se>
% -------------------------------------------------------------------------
