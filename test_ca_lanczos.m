function test_ca_lanczos(matrix, steps, eig_s, eig_l)

% Load matrix
load(matrix);
Problem
A=Problem.A;
n=size(A,1);
%N=1000; A=sparse(1:N, 1:N, rand(N,1), N, N);

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

% Lanczos options
opt.break = 0; %TODO: Remove
opt.orth = 'full';

% Standard Lanczos
time_st=cputime;
[T_st,V_st,r_st,o_st]=lanczos(A,r0,m,opt.break,opt.orth);
time_st=cputime-time_st;

% CA-Lanczos, for several values of s
time_ca_4=cputime;
[T_ca_4,V_ca_4,r_ca_4,o_ca_4]=ca_lanczos(A,r0,4,m/4,'monomial',opt.orth);
time_ca_4=cputime-time_ca_4;
time_ca_6=cputime;
[T_ca_6,V_ca_6,r_ca_6,o_ca_6]=ca_lanczos(A,r0,6,m/6,'monomial',opt.orth);
time_ca_6=cputime-time_ca_6;
time_ca_8=cputime;
[T_ca_8,V_ca_8,r_ca_8,o_ca_8]=ca_lanczos(A,r0,8,m/8,'monomial',opt.orth);
time_ca_8=cputime-time_ca_8;
time_ca_10=cputime;
[T_ca_10,V_ca_10,r_ca_10,o_ca_10]=ca_lanczos(A,r0,10,m/10,'monomial',opt.orth);
time_ca_10=cputime-time_ca_10;
time_ca_20=cputime;
%[T_ca_20,V_ca_20,r_ca_20,o_ca_20]=ca_lanczos(A,r0,20,m/20,'newton',opt.orth);
%time_ca_20=cputime-time_ca_20;

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
eig_T_ca_6 = eig(T_ca_6);
disp('CA-Lanczos, s=6');
disp(['Execution time:                ', num2str(time_ca_6)]);
disp(['Number of iterations:          ', num2str(size(T_ca_6,1))]);
disp(['Error in smallest eigenvalue:  ', num2str(abs(se-min(eig_T_ca_6))/abs(se))]);
disp(['Error in largest eigenvalue:   ', num2str(abs(le-max(eig_T_ca_6))/abs(le))]);
eig_T_ca_8 = eig(T_ca_8);
disp('CA-Lanczos, s=8');
disp(['Execution time:                ', num2str(time_ca_8)]);
disp(['Number of iterations:          ', num2str(size(T_ca_8,1))]);
disp(['Error in smallest eigenvalue:  ', num2str(abs(se-min(eig_T_ca_8))/abs(se))]);
disp(['Error in largest eigenvalue:   ', num2str(abs(le-max(eig_T_ca_8))/abs(le))]);
eig_T_ca_10 = eig(T_ca_10);
disp('CA-Lanczos, s=10');
disp(['Execution time:                ', num2str(time_ca_10)]);
disp(['Number of iterations:          ', num2str(size(T_ca_10,1))]);
disp(['Error in smallest eigenvalue:  ', num2str(abs(se-min(eig_T_ca_10))/abs(se))]);
disp(['Error in largest eigenvalue:   ', num2str(abs(le-max(eig_T_ca_10))/abs(le))]);
%eig_T_ca_20 = eig(T_ca_20);
%disp('CA-Lanczos, s=20');
%disp(['Execution time:                ', num2str(time_ca_20)]);
%disp(['Number of iterations:          ', num2str(size(T_ca_20,1))]);
%disp(['Error in smallest eigenvalue:  ', num2str(abs(se-min(eig_T_ca_20))/abs(se))]);
%disp(['Error in largest eigenvalue:   ', num2str(abs(le-max(eig_T_ca_20))/abs(le))]);
disp('=======');

%Define some colors for plotting
blue         = [   0   0   1 ];
mediumgreen  = [   0 .60   0 ];
red          = [   1   0   0 ];
orchid       = [ .85 .44 .84 ];
pink         = [   1 .41 .71 ];
salmon       = [ .98 .50 .45 ];
orange       = [   1 .55   0 ]; 
royalblue    = [.25  .41 .88 ];

% Plot the eigenvalues
figure; hold on;
plot((0:size(T_st,1)-1)/(size(T_st,1)-1),sort(eig_T_st),'*','Color',mediumgreen);
plot((0:size(T_ca_4,1)-1)/(size(T_ca_4,1)-1),sort(eig_T_ca_4),'o','Color',red);
plot((0:size(T_ca_6,1)-1)/(size(T_ca_6,1)-1),sort(eig_T_ca_6),'o','Color',orchid);
plot((0:size(T_ca_8,1)-1)/(size(T_ca_8,1)-1),sort(eig_T_ca_8),'o','Color',royalblue);
plot((0:size(T_ca_10,1)-1)/(size(T_ca_10,1)-1),sort(eig_T_ca_10),'o','Color',orange);
%plot((0:size(T_ca_20,1)-1)/(size(T_ca_20,1)-1),sort(eig_T_ca_20),'o','Color',salmon);
title('Spectrum of A');

% Plot the convergence of the smallest eigenpair
figure;
semilogy(1:size(T_st,1),r_st(:,1),'*:','Color',mediumgreen,'MarkerSize',8);  hold on;
semilogy((1:size(T_ca_4,1)/4)*4,r_ca_4(:,1),'o:','Color',red,'MarkerSize',8,'LineWidth',1.5);
semilogy((1:size(T_ca_6,1)/6)*6,r_ca_6(:,1),'o:','Color',orchid,'MarkerSize',8,'LineWidth',1.5);
semilogy((1:size(T_ca_8,1)/8)*8,r_ca_8(:,1),'o:','Color',royalblue,'MarkerSize',8,'LineWidth',1.5);
semilogy((1:size(T_ca_10,1)/10)*10,r_ca_10(:,1),'o:','Color',orange,'MarkerSize',8,'LineWidth',1.5);
%semilogy((1:size(T_ca_20,1)/20)*20,r_ca_20(:,1),'o:','Color',salmon,'MarkerSize',8,'LineWidth',1.5);
title('Convergence of smallest eigenpair');
ylabel('||Ay - \lambda y||');
xlabel('Iterations');
legend('Lanczos','CA-Lanczos(4)','CA-Lanczos(6)','CA-Lanczos(8)','CA-Lanczos(10)','CA-Lanczos(20)');
if (min(r_st(:,1)) < 1.0e-14) ...
   || (min(r_ca_4(:,1)) < 1.0e-14) || (min(r_ca_6(:,1)) < 1.0e-14) ...
   || (min(r_ca_8(:,1)) < 1.0e-14) || (min(r_ca_10(:,1)) < 1.0e-14) ...%|| (min(r_ca_20(:,1)) < 1.0e-14)

    ylim([1.0e-16 1.0e02]);
end

% Plot the convergence of the largest eigenpair
figure;
semilogy(1:size(T_st,1),r_st(:,2),'*:','Color',mediumgreen,'MarkerSize',8);  hold on;
semilogy((1:size(T_ca_4,1)/4)*4,r_ca_4(:,2),'o:','Color',red,'MarkerSize',8,'LineWidth',1.5);
semilogy((1:size(T_ca_6,1)/6)*6,r_ca_6(:,2),'o:','Color',orchid,'MarkerSize',8,'LineWidth',1.5);
semilogy((1:size(T_ca_8,1)/8)*8,r_ca_8(:,2),'o:','Color',royalblue,'MarkerSize',8,'LineWidth',1.5);
semilogy((1:size(T_ca_10,1)/10)*10,r_ca_10(:,2),'o:','Color',orange,'MarkerSize',8,'LineWidth',1.5);
%semilogy((1:size(T_ca_20,1)/20)*20,r_ca_20(:,2),'o:','Color',salmon,'MarkerSize',8,'LineWidth',1.5);
title('Convergence of largest eigenpair');
ylabel('||Ay - \lambda y||');
xlabel('Iterations');
legend('Lanczos','CA-Lanczos(4)','CA-Lanczos(6)','CA-Lanczos(8)','CA-Lanczos(10)','CA-Lanczos(20)');
if (min(r_st(:,2)) < 1.0e-14) ...
   || (min(r_ca_4(:,2)) < 1.0e-14) || (min(r_ca_6(:,2)) < 1.0e-14) ...
   || (min(r_ca_8(:,2)) < 1.0e-14) || (min(r_ca_10(:,2)) < 1.0e-14) ...%|| (min(r_ca_20(:,2)) < 1.0e-14)

    ylim([1.0e-16 1.0e02]);
end

% Plot the orthogonality 
figure;
semilogy(1:size(T_st,1),o_st,'*:','Color',mediumgreen,'MarkerSize',8); hold on;
semilogy((1:size(T_ca_4,1)/4)*4,o_ca_4,'o:','Color',red,'MarkerSize',8,'LineWidth',1.5);
semilogy((1:size(T_ca_6,1)/6)*6,o_ca_6,'o:','Color',orchid,'MarkerSize',8,'LineWidth',1.5);
semilogy((1:size(T_ca_8,1)/8)*8,o_ca_8,'o:','Color',royalblue,'MarkerSize',8,'LineWidth',1.5);
semilogy((1:size(T_ca_10,1)/10)*10,o_ca_10,'o:','Color',orange,'MarkerSize',8,'LineWidth',1.5);
%semilogy((1:size(T_ca_20,1)/20)*20,o_ca_20,'o:','Color',salmon,'MarkerSize',8,'LineWidth',1.5);
title('Orthogonality of basis vectors');
ylabel('|| I - Q_m^*Q ||_F');
xlabel('Iterations');
legend('Lanczos','CA-Lanczos(4)','CA-Lanczos(6)','CA-Lanczos(8)','CA-Lanczos(10)','CA-Lanczos(20)');
         
if (min(o_st) < 1.0e-14) ...
   || (min(o_ca_4) < 1.0e-14) || (min(o_ca_6) < 1.0e-14) ...
   || (min(o_ca_8) < 1.0e-14) || (min(o_ca_10) < 1.0e-14) ...%|| (min(o_ca_20) < 1.0e-14)

    ylim([1.0e-16 1.0e02]);
end
