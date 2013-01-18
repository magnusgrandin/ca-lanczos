numTimeSteps = 200;
numLanczosSteps = 24;
numRestarts = 1;
s = 6;
N=512; dt=0.025;
t=dt:dt:numTimeSteps*dt;
e=ones(N,1);
range=[-10 10];
h=1/N*(range(2)-range(1));
x=range(1)+h/2:h:range(2)-h/2;

w=0.5; m=1;
displ = 4;

H = sparse(1:N, [2:N 1], 4*e/3, N, N) - sparse(1:N, [3:N 1 2], e/12, N, N);
H = (H + H') - sparse(1:N, 1:N, 5*e/2, N, N);
H=-H/(2*(h^2));
H=H + sparse(1:N, 1:N, 0.5*x.^2, N, N);

psi = (1/(pi*w^2))^(1/4)*exp(-0.5*(((x-displ)/w).^2));
psi = psi';

tf = numTimeSteps*dt;
x0 = cos(w*tf)*displ;
p = -m*w*sin(w*tf)*displ;
g = -0.5*(m*w*displ^2*cos(w*tf)*sin(w*tf) + w*tf);
sigma = 1/sqrt(m*w);
%psi_ref = (1/(pi*w^2))^(1/4)*exp(-0.5*((x-x0)/sigma).^2).*(cos(p*(x-x0)+g)-1i*sin(p*(x-x0)+g));
%psi_ref = psi_ref';
psi_ref = expm(-1i*tf*H)*psi;

eig_refnc = zeros(numLanczosSteps,1);

smalleig_err_stlan = zeros(numTimeSteps,1);
smalleig_err_calan_newton = zeros(numTimeSteps,1);
smalleig_err_calan_monomial = zeros(numTimeSteps,1);
%smalleig_err_sslan = zeros(numTimeSteps,1);
largeeig_err_stlan = zeros(numTimeSteps,1);
largeeig_err_calan_newton = zeros(numTimeSteps,1);
largeeig_err_calan_monomial = zeros(numTimeSteps,1);
%largeeig_err_sslan = zeros(numTimeSteps,1);

psi_refnc = zeros(N,numTimeSteps+1);
psi_stlan = zeros(N,numTimeSteps+1);
psi_calan_newton = zeros(N,numTimeSteps+1);
psi_calan_monomial = zeros(N,numTimeSteps+1);
%psi_sslan = zeros(N,numTimeSteps+1);
psi_refnc = psi;
psi_stlan = psi;
psi_calan_newton = psi;
psi_calan_monomial = psi;
%psi_sslan = psi;

% Reference solution: Propagated matrix exponential
%disp('Reference solution.');
%time_refnc = cputime;
% for k=2:numTimeSteps+1
%     psi_refnc = expm(-1i*dt*H) * psi_refnc;
%     fprintf('.');
% end
%fprintf('\n');
%time_refnc = cputime-time_refnc;

% Std Lanczos
disp('Standard Lanczos.');
tic;
for k=2:numTimeSteps+1
    [T,Q] = lanczos_prop(H, psi_stlan, numLanczosSteps, dt, 1.0e-10, false);
    nLanczos = size(T,1);
    [V,D]=eig(T);
    %psi_stlan = ...
    %    Q*(V*(expm(-1i*dt*D)*(V'*(eye(nLanczos,1)*norm(psi_stlan)))));
    psi_stlan = ...
        Q*(expm(-1i*dt*T)*(eye(nLanczos,1)*norm(psi_stlan)));
    eig_stlan = diag(D);
    n_stlan = nLanczos;
    T_stlan = T;
    Q_stlan = Q;
    V_stlan = V;
    D_stlan = D;
end
time_stlan = toc;
       
% CA-Lanczos
disp('CA-Lanczos, newton basis.');
tic;
% Take first step with regular Lanczos to get inital eigenvalue estimate
[T,Q] = lanczos_prop(H, psi_calan_newton, numLanczosSteps, dt, 1.0e-10, false);
[V,D] = eig(T);
nLanczos = size(T,1);
psi_calan_newton = ...
    Q*(expm(-1i*dt*T)*(eye(nLanczos,1)*norm(psi_calan_newton)));
eigest = diag(D);
%eigest = [];
for k=3:numTimeSteps+1
    [T,Q] = ca_lanczos_prop(H, psi_calan_newton, s, numLanczosSteps/s, dt, 1.0e-10,'newton',eigest,false);
    [V,D] = eig(T);
    %eigest = diag(D);
    nLanczos = size(T,1);
    %psi_calan_newton = ...
    %    Q*(V*(expm(-1i*dt*D)*(V\(eye(nLanczos,1)*norm(psi_calan_newton)))));
    psi_calan_newton = ...
        Q*(expm(-1i*dt*T)*(eye(nLanczos,1)*norm(psi_calan_newton)));
    eig_calan_newton = diag(D);
    n_calan_newton = nLanczos;
    T_calan_newton = T;
    Q_calan_newton = Q;
    V_calan_newton = V;
    D_calan_newton = D;
end
time_calan_newton = toc;

% CA-Lanczos
disp('CA-Lanczos, monomial basis.');
tic;
for k=2:numTimeSteps+1
    [T,Q] = ca_lanczos_prop(H, psi_calan_monomial, s, numLanczosSteps/s, dt, 1.0e-10,'monomial',eigest,false);
    nLanczos = size(T,1);
    [V,D] = eig(T);
    %psi_calan_monomial = ...
    %    Q*(V*(expm(-1i*dt*D)*(V\(eye(nLanczos,1)*norm(psi_calan_monomial)))));
    psi_calan_monomial = ...
        Q*(expm(-1i*dt*T)*(eye(nLanczos,1)*norm(psi_calan_monomial)));
    eig_calan_monomial = diag(D);
    n_calan_monomial = nLanczos;
    T_calan_monomial = T;
    Q_calan_monomial = Q;
    V_calan_monomial = V;
    D_calan_monomial = D;
end
time_calan_monomial = toc;

% % s-Step Lanczos
% disp('s-step Lanczos.');
% tic;
% for k=2:numTimeSteps+1
%     [T,Q] = sstep_lanczos_prop(H, psi_sslan, s, numLanczosSteps/s, dt, 1.0e-10);
% %    [T,Q] = sStepLanczos(H, psi_sslan(:,k-1), s, numLanczosSteps/s);
% %    nLanczos = numLanczosSteps;
%     [V,D] = eig(T);
%     nLanczos = size(T,1);
%     psi_sslan = ...
%         Q*(V*(expm(-1i*dt*D)*(V\(norm(psi_sslan)*eye(nLanczos,1)))));
%     eig_sslan = diag(D);
%     n_sslan = nLanczos;
%     T_sslan = T;
%     Q_sslan = Q;
%     V_sslan = V;
%     D_sslan = D;
% end
% time_sslan = toc;

% Plot eigenvalues on the interval [0 1]
figure; hold on;

%plot((0:N-1)/(N-1), sort(eig_refnc, 'descend'), 'k-');
plot((0:n_stlan-1)/(n_stlan-1),sort(real(eig_stlan),'descend'),'rs');
plot((0:n_calan_newton-1)/(n_calan_newton-1),sort(real(eig_calan_newton),'descend'),'bo');
plot((0:n_calan_monomial-1)/(n_calan_monomial-1),sort(real(eig_calan_monomial),'descend'),'g^');
%plot((0:n_sslan-1)/(n_sslan-1),sort(real(eig_sslan),'descend'),'m*');

legend('Lanczos','CA-Lanczos (newton)','CA-Lanczos (monomial)'); %'s-Step Lanczos'

% Compute errors in eigenvalues
smalleig_err_stlan(k) = abs(min(abs(eig_refnc)) - min(abs(eig_stlan)));
smalleig_err_calan_newton(k) = abs(min(abs(eig_refnc)) - min(abs(eig_calan_newton)));
smalleig_err_calan_monomial(k) = abs(min(abs(eig_refnc)) - min(abs(eig_calan_monomial)));
%smalleig_err_sslan(k) = abs(min(abs(eig_refnc)) - min(abs(eig_sslan)));
largeeig_err_stlan(k) = abs(max(abs(eig_refnc)) - max(abs(eig_stlan)));
largeeig_err_calan_newton(k) = abs(max(abs(eig_refnc)) - max(abs(eig_calan_newton)));
largeeig_err_calan_monomial(k) = abs(max(abs(eig_refnc)) - max(abs(eig_calan_monomial)));
%largeeig_err_sslan(k) = abs(max(abs(eig_refnc)) - max(abs(eig_sslan)));

%disp('Matrix exponential:');
%disp(['Execution time:                ', num2str(time_refnc)]);
disp('Standard Lanczos:');
disp(['Execution time:                ', num2str(time_stlan)]);
disp(['Number of Lanczos iterations:  ', num2str(n_stlan)]);
disp(['Error in smallest eigenvalue:  ', num2str(max(smalleig_err_stlan))]);
disp(['Error in largest eigenvalue:   ', num2str(max(largeeig_err_stlan))]);
disp('CA-Lanczos, newton basis:');
disp(['Execution time:                ', num2str(time_calan_newton)]);
disp(['Number of Lanczos iterations:  ', num2str(n_calan_newton)]);
disp(['Error in smallest eigenvalue:  ', num2str(max(smalleig_err_calan_newton))]);
disp(['Error in largest eigenvalue:   ', num2str(max(largeeig_err_calan_newton))]);
disp('CA-Lanczos, monomial basis:');
disp(['Execution time:                ', num2str(time_calan_monomial)]);
disp(['Number of Lanczos iterations:  ', num2str(n_calan_monomial)]);
disp(['Error in smallest eigenvalue:  ', num2str(max(smalleig_err_calan_monomial))]);
disp(['Error in largest eigenvalue:   ', num2str(max(largeeig_err_calan_monomial))]);
%disp('s-Step Lanczos:');
%disp(['Execution time:                ', num2str(time_sslan)]);
%disp(['Number of Lanczos iterations:  ', num2str(n_sslan)]);
%disp(['Error in smallest eigenvalue:  ', num2str(max(smalleig_err_sslan))]);
%disp(['Error in largest eigenvalue:   ', num2str(max(largeeig_err_sslan))]);

figure;
subplot(1,3,1);
plot(x, real(psi_stlan), x, imag(psi_stlan), x, abs(psi_stlan));
legend('Re(\Psi_{k})', 'Im(\Psi_{k})', 'abs(\Psi_{k})');
subplot(1,3,2);
plot(x, real(psi_calan_newton), x, imag(psi_calan_newton), x, abs(psi_calan_newton));
legend('Re(\Psi_{k})', 'Im(\Psi_{k})', 'abs(\Psi_{k})');
subplot(1,3,3);
plot(x, real(psi_calan_monomial), x, imag(psi_calan_monomial), x, abs(psi_calan_monomial));
legend('Re(\Psi_{k})', 'Im(\Psi_{k})', 'abs(\Psi_{k})');
%subplot(1,3,3);
%plot(x, real(psi_sslan), x, imag(psi_sslan), x, abs(psi_sslan));
%legend('Re(\Psi_{k})', 'Im(\Psi_{k})', 'abs(\Psi_{k})');

err_stlan = abs(psi_ref-psi_stlan); max(err_stlan)
err_calan_newton = abs(psi_ref-psi_calan_newton); max(err_calan_newton)
err_calan_monomial = abs(psi_ref-psi_calan_monomial); max(err_calan_monomial)
%err_sslan = abs(psi_ref-psi_sslan); max(err_sslan)

%     mass = params.mass[0];
%     omega = params.potential[0][0].arg[0];
%     gaussWaveFunction(Psi, 0, params.displacement, sqrt(1.0/(mass*omega)), params);

%     ohmft = omega*params.timeMax;
%     for(d = 0; d < NUM_DIMENSIONS; d++)
%     {
%         x0[d]    = (cos(ohmft))*params.displacement[d];
%         p[d]     = (-mass*omega*sin(ohmft))*params.displacement[d];
%         gamma[d] = -((mass*omega*pow(params.displacement[d],2.0)*cos(ohmft)*sin(ohmft))+(ohmft))/2.0;
%     }
%     gaussWaveFunctionComplex(PsiRef, 0, x0, sqrt(1.0/(mass*omega)), p, gamma, params);

