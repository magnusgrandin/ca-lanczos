% Test script for explicitly restarted CA-Lanczos, computing the 10 largest
% eigenvalues.
%

matrix = 'mesh2e1'; 

disp(['Matrix: ', matrix]);

matrix_path_prefix = '../matrices';
load ([matrix_path_prefix,'/',matrix]);

A = Problem.A;
%max_norm = max(max(abs(A)));
%A = A/max_norm;
N = size(A,1);
maxvecs = 56;
nweigs = 10;
S = 8 %[1,2,4,6,8,10];
basis = 'newton';
orth = 'full';
tol = 1.0e-8;

if N < 1000
    eref = eig(A);
else
    opts = [];
    opts.tol = tol;
    eref = eigs(A,nweigs,'lm',opts);
end
    
r = ones(N,1);
f = fopen('eigerr.dat','w');
fprintf(f,'#tol: %.4e\n',tol);
fprintf(f,'#info: basis,orth,K,s\n\n');
for s = S
%    fprintf(f,'#%s,%s,%.1e,%d\n',basis,orth,K,s);
    [E,V,nres,rnorms,ortherr] = restarted_ca_lanczos(A,r,maxvecs,nweigs,s,basis,orth,tol);
    nconv = length(E);
    %E = E*max_norm;

    fprintf(f,'%d\n',nconv);
    fprintf(f,'%.16e ',full(E)); fprintf(f,'\n');
    fprintf(f,'%.16e ',eref); fprintf(f,'\n');
    fprintf(f,'%d\n',nres);
    for i = 1:nconv
    %    fprintf(f,'%.6e ',rnorms(1:nres,i)); fprintf(f,'\n');
    end
    fprintf(f,'%.6e ',ortherr(1:nres)); fprintf(f,'\n');
    figure('Position',[1 1 500 500]);
    subplot(2,1,1);
    set(gca,'FontSize',16);
    semilogy(rnorms);
    if nres > 1
        xlim([1 nres]);
    end
    ref_line = tol*ones(nres,1);
    hold on; plot(ref_line, 'k--');
    ylabel('Relative residual norm');
    subplot(2,1,2);
    set(gca,'FontSize',16);
    semilogy(ortherr);
    if nres > 1
        xlim([1 nres]);
    end
    xlabel('Iteration');
    ylabel('Orthogonalization error');
end

    