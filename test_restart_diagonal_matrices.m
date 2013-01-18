% Test script for explicitly restarted CA-Lanczos, computing the 10 largest
% eigenvalues.
%
% The matrices used here are diagonal, with equally spaced eigenvalues from
% 1..K, where K is the desired condition number of each matrix. The value
% of N is constant for all matrices generated.

N = 5000;
maxvecs = 60;
nweigs = 10;
S = 4 %[1,2,4,6,8,10];
basis = 'newton';
orth = 'full';
tol = 1.0e-8;

r = ones(N,1);
f = fopen('eigerr.dat','w');
fprintf(f,'#tol: %.4e\n',tol);
fprintf(f,'#info: basis,orth,K,s\n\n');
for k = 4
    K = 10^k;
    a = linspace(1,K,N);
    A = sparse(diag(a));
    for s = S
        fprintf(f,'#%s,%s,%.1e,%d\n',basis,orth,K,s);
        [E,V,nres,rnorms,ortherr] = restarted_ca_lanczos(A,r,maxvecs,nweigs,s,basis,orth,tol);
        nconv = length(E);
        eref = a(end:-1:end-nconv+1)';
        fprintf(f,'%d\n',nconv);
        fprintf(f,'%.16e ',E); fprintf(f,'\n');
        fprintf(f,'%.16e ',eref); fprintf(f,'\n');
        fprintf(f,'%d\n',nres);
        for i = 1:nconv
            fprintf(f,'%.6e ',rnorms(1:nres,i)); fprintf(f,'\n');
        end
        fprintf(f,'%.6e ',ortherr(1:nres)); fprintf(f,'\n');
        figure('Position',[1 1 500 500]);
        subplot(2,1,1);
        set(gca,'FontSize',16);
        semilogy(rnorms);
        ylabel('Relative residual norm');
        ref_line = tol*ones(nres,1);
        hold on; plot(ref_line, 'k--');
        subplot(2,1,2);
        set(gca,'FontSize',16);
        semilogy(ortherr);
        xlabel('Iteration');
        ylabel('Orthogonalization error');
    end
end

% -------------------------------------------------------------------------
%  Copyright (2012, 2013)  Magnus Grandin <magnus.grandin@it.uu.se>
% -------------------------------------------------------------------------
    