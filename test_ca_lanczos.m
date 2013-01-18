%% function test_ca_lanczos(A, r0, steps, orth, basis, eig_s, eig_l)

function test_ca_lanczos(A, r0, steps, orth, basis, eig_s, eig_l)

% Total number of Lanczos steps (maximum)
if nargin >= 2 && isempty(steps)
    steps = 0;
end
if nargin >= 2 && steps ~= 0
    m=steps;
else
    m=120;
end

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
[T_st,V_st,r_st,o_st]=lanczos(A,r0,m,orth);
time_st=cputime-time_st;

% CA-Lanczos, for several values of s
time_ca_4=cputime;
[T_ca_4,V_ca_4,r_ca_4,o_ca_4]=ca_lanczos(A,r0,4,m,basis,orth);
time_ca_4=cputime-time_ca_4;
time_ca_8=cputime;
[T_ca_8,V_ca_8,r_ca_8,o_ca_8]=ca_lanczos(A,r0,8,m,basis,orth);
time_ca_8=cputime-time_ca_8;
time_ca_12=cputime;
[T_ca_12,V_ca_12,r_ca_12,o_ca_12]=ca_lanczos(A,r0,12,m,basis,orth);
time_ca_12=cputime-time_ca_12;
time_ca_16=cputime;
[T_ca_16,V_ca_16,r_ca_16,o_ca_16]=ca_lanczos(A,r0,16,m,basis,orth);
time_ca_16=cputime-time_ca_16;

% Compute reference smallest and largest eigenvalues
opts.disp=0;
se=0; le=0;
disp('=======');
if ~exist('eig_s','var') || isempty(eig_s)
    [V_,se,fs] = eigs(A,1,'sm',opts);
    if(fs == 0)
        disp(['Smallest eigenvalue: ' num2str(se)]);
    else
        disp('Smallest eigenvalue did not converge using eigs.');
    end
end
if ~exist('eig_l','var') || isempty(eig_l)
    [V_,le,fl] = eigs(A,1,'lm',opts);
    if(fl == 0)
        disp(['Largest eigenvalue:  ' num2str(le)]);
    else
        disp('Largest eigenvalue did not converge using eigs.');
    end
end
disp('=======');
    
% Compute eigenvalues/vectors for standard Lanczos
eig_T_st=eig(T_st);
disp('Lanczos:');
disp(['Execution time:                ', num2str(time_st)]);
disp(['Number of iterations:          ', num2str(size(T_st,1))]);
disp(['Error in smallest eigenvalue:  ', num2str(abs(se-min(eig_T_st))/abs(se))]);
disp(['Error in largest eigenvalue:   ', num2str(abs(le-max(eig_T_st))/abs(le))]);

% Compute eigenvalues/vectors for CA-Lanczos
eig_T_ca_4 = eig(T_ca_4);
disp('CA-Lanczos, s=4');
disp(['Execution time:                ', num2str(time_ca_4)]);
disp(['Number of iterations:          ', num2str(size(T_ca_4,1))]);
disp(['Error in smallest eigenvalue:  ', num2str(abs(se-min(eig_T_ca_4))/abs(se))]);
disp(['Error in largest eigenvalue:   ', num2str(abs(le-max(eig_T_ca_4))/abs(le))]);
eig_T_ca_8 = eig(T_ca_8);
disp('CA-Lanczos, s=8');
disp(['Execution time:                ', num2str(time_ca_8)]);
disp(['Number of iterations:          ', num2str(size(T_ca_8,1))]);
disp(['Error in smallest eigenvalue:  ', num2str(abs(se-min(eig_T_ca_8))/abs(se))]);
disp(['Error in largest eigenvalue:   ', num2str(abs(le-max(eig_T_ca_8))/abs(le))]);
eig_T_ca_12 = eig(T_ca_12);
disp('CA-Lanczos, s=12');
disp(['Execution time:                ', num2str(time_ca_12)]);
disp(['Number of iterations:          ', num2str(size(T_ca_12,1))]);
disp(['Error in smallest eigenvalue:  ', num2str(abs(se-min(eig_T_ca_12))/abs(se))]);
disp(['Error in largest eigenvalue:   ', num2str(abs(le-max(eig_T_ca_12))/abs(le))]);
eig_T_ca_16 = eig(T_ca_16);
disp('CA-Lanczos, s=16');
disp(['Execution time:                ', num2str(time_ca_16)]);
disp(['Number of iterations:          ', num2str(size(T_ca_16,1))]);
disp(['Error in smallest eigenvalue:  ', num2str(abs(se-min(eig_T_ca_16))/abs(se))]);
disp(['Error in largest eigenvalue:   ', num2str(abs(le-max(eig_T_ca_16))/abs(le))]);
disp('=======');

%Define some colors for plotting
blue         = [   0   0   1 ];
mediumgreen  = [   0 .60   0 ];
red          = [   1   0   0 ];
orchid       = [ .85 .44 .84 ];
pink         = [   1 .41 .71 ];
salmon       = [ .98 .50 .45 ];
orange       = [   1 .5   0 ];
royalblue    = [ .25 .41 .88 ];
blueviolet   = [ .58 .40 .64 ];
%orange       = [ .95 .19 .33 ];

% Plot the eigenvalues
figure; hold on;
plot((0:size(T_st,1)-1)/(size(T_st,1)-1),sort(eig_T_st),'*','Color',mediumgreen);
plot((0:size(T_ca_4,1)-1)/(size(T_ca_4,1)-1),sort(eig_T_ca_4),'o','Color',red);
plot((0:size(T_ca_8,1)-1)/(size(T_ca_8,1)-1),sort(eig_T_ca_8),'s','Color',royalblue);
plot((0:size(T_ca_12,1)-1)/(size(T_ca_12,1)-1),sort(eig_T_ca_12),'d','Color',orange);
plot((0:size(T_ca_16,1)-1)/(size(T_ca_16,1)-1),sort(eig_T_ca_16),'^','Color',blueviolet);
title('Spectrum of A');

% % Plot the convergence of the smallest eigenpair
% figure;
% semilogy(1:size(T_st,end),r_st(:,1),'*:','Color',mediumgreen,'MarkerSize',8);  hold on;
% semilogy((1:size(T_ca_4,end)/4)*4,r_ca_4(:,1),'o:','Color',red,'MarkerSize',8,'LineWidth',1.5);
% semilogy((1:size(T_ca_8,end)/8)*8,r_ca_8(:,1),'o:','Color',royalblue,'MarkerSize',8,'LineWidth',1.5);
% semilogy((1:size(T_ca_12,end)/12)*10,r_ca_12(:,1),'o:','Color',orange,'MarkerSize',8,'LineWidth',1.5);
% %semilogy((1:size(T_ca_16,end)/16)*16,r_ca_16(:,1),'o:','Color',salmon,'MarkerSize',8,'LineWidth',1.5);
% title('Convergence of smallest eigenpair');
% ylabel('||Ay - \lambda y||');
% xlabel('Iterations');
% legend('Lanczos','CA-Lanczos(4)','CA-Lanczos(8)','CA-Lanczos(12)','CA-Lanczos(16)');
% if (min(r_st(:,1)) < 1.0e-14) ...
%    || (min(r_ca_4(:,1)) < 1.0e-14) || (min(r_ca_8(:,1)) < 1.0e-14) || ...
%       (min(r_ca_12(:,1)) < 1.0e-14) ...%|| (min(r_ca_16(:,1)) < 1.0e-14)
% 
%     ylim([1.0e-16 1.0e02]);
% end

% Plot the convergence of the largest eigenpair
figure('Position',[1 1 500 500]);
subplot(2,1,1);
set(gca,'FontSize',16)
semilogy(1:size(T_st,1),r_st(:,1),'*:','Color',mediumgreen,'MarkerSize',8);  hold on;
semilogy((1:size(T_ca_4,1)/4)*4,r_ca_4(:,1),'o:','Color',red,'MarkerSize',8,'LineWidth',1);
semilogy((1:size(T_ca_8,1)/8)*8,r_ca_8(:,1),'s:','Color',royalblue,'MarkerSize',8,'LineWidth',1);
semilogy((1:size(T_ca_12,1)/12)*12,r_ca_12(:,1),'d:','Color',orange,'MarkerSize',8,'LineWidth',1);
semilogy((1:size(T_ca_16,1)/16)*16,r_ca_16(:,1),'^:','Color',blueviolet,'MarkerSize',8,'LineWidth',1);
%title('Convergence of largest eigenpair');
%ylabel('||Ay - \lambda y||');
ylabel('Relative residual norm');
%xlabel('Iteration');
legend('Lanczos','s=4','s=8','s=12','s=16','Location','North','Orientation','Horizontal');
ylim([1.0e-16 1.0e03]);
set(gca,'YTick',[1e-15,1e-10,1e-05,1]);
% if (min(r_st(:,2)) < 1.0e-14) ...
%    || (min(r_ca_4(:,2)) < 1.0e-14) || (min(r_ca_6(:,2)) < 1.0e-14) ...
%    || (min(r_ca_8(:,2)) < 1.0e-14) || (min(r_ca_10(:,2)) < 1.0e-14) ...%|| (min(r_ca_16(:,2)) < 1.0e-14)
% 
%     ylim([1.0e-16 1.0e03]);
% end

% Plot the orthogonality 
subplot(2,1,2);
set(gca,'FontSize',16)
semilogy(1:size(T_st,1),o_st,'*:','Color',mediumgreen,'MarkerSize',8); hold on;
semilogy((1:size(T_ca_4,1)/4)*4,o_ca_4,'o:','Color',red,'MarkerSize',8,'LineWidth',1);
semilogy((1:size(T_ca_8,1)/8)*8,o_ca_8,'s:','Color',royalblue,'MarkerSize',8,'LineWidth',1);
semilogy((1:size(T_ca_12,1)/12)*12,o_ca_12,'d:','Color',orange,'MarkerSize',8,'LineWidth',1);
semilogy((1:size(T_ca_16,1)/16)*16,o_ca_16,'^:','Color',blueviolet,'MarkerSize',8,'LineWidth',1);
%title('Orthogonality of basis vectors');
%ylabel('|| I - Q_m^*Q ||_F');
ylabel('Orthogonalization error');
xlabel('Iteration');
%legend('Lanczos','s=4','s=6','s=8','s=10','s=16','Location','North');
ylim([1.0e-16 1e02]);
set(gca,'YTick',[1e-15,1e-10,1e-05,1]);
         
% if (min(o_st) < 1.0e-14) ...
%    || (min(o_ca_4) < 1.0e-14) || (min(o_ca_6) < 1.0e-14) ...
%    || (min(o_ca_8) < 1.0e-14) || (min(o_ca_10) < 1.0e-14) ...%|| (min(o_ca_16) < 1.0e-14)
% 
%     ylim([1.0e-16 1.0e02]);
% end

% -------------------------------------------------------------------------
%  Copyright (2012, 2013)  Magnus Grandin <magnus.grandin@it.uu.se>
% -------------------------------------------------------------------------
