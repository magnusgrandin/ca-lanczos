numTimeSteps = 100;
numLanczosSteps = 24;
numRestarts = 1;
s =6;
N=256; dt=0.1;
t=dt:dt:numTimeSteps*dt;
e=ones(N,1);
range=[-10 10];
h=1/N*(range(2)-range(1));
x=range(1)+h/2:h:range(2)-h/2;

w=1.0; m=0;
T=0; D=10.90240 ; alpha=0.3; xe=8;
displ = 2;

H = sparse(1:N, [2:N 1], 4*e/3, N, N) - sparse(1:N, [3:N 1 2], e/12, N, N);
H = (H + H') - sparse(1:N, 1:N, 5*e/2, N, N);
H=-H/(2*(h^2));
H=H + sparse(1:N, 1:N, 0.5*x.^2, N, N);
%H=H + sparse(1:N, 1:N, T+D*(1.0-exp(-alpha*((x-displ)-xe).^2)), N, N);

%psi(:,1) = (1.0/(w*sqrt(2*pi))) * exp(-0.5*((x-m)/w).^2)';
%psi(:,1) = 6*(1.0/(w*sqrt(2*pi))) * exp(-0.5*((x-m)/w).^2)';
psi = (1/(pi*w^2))^(1/4)*exp(-0.5*(((x-displ)/w).^2))';
%psi(:,1) = exp(-(20*(x-m)).^2)';
%psi = psi/norm(psi);

eig_refnc = zeros(numLanczosSteps,1);

smalleig_err_stlan = zeros(numTimeSteps,1);
smalleig_err_calan = zeros(numTimeSteps,1);
smalleig_err_sslan = zeros(numTimeSteps,1);
largeeig_err_stlan = zeros(numTimeSteps,1);
largeeig_err_calan = zeros(numTimeSteps,1);
largeeig_err_sslan = zeros(numTimeSteps,1);

psi_refnc = zeros(N,numTimeSteps+1);
psi_stlan = zeros(N,numTimeSteps+1);
psi_calan = zeros(N,numTimeSteps+1);
psi_sslan = zeros(N,numTimeSteps+1);
psi_refnc(:,1) = psi;
psi_stlan(:,1) = psi;
psi_calan(:,1) = psi;
psi_sslan(:,1) = psi;

%for 2:numTimeSteps+1
    % Reference solution
%    opts.isreal = 0;
%    opts.issym  = 1;
%    opts.disp   = 0;
%    opts.v0     = psi_refnc(:,k-1);
%    opts.maxit  = numLanczosSteps;
%    opts.p      = numLanczosSteps+1;
%    [V,D] = eigs(H,numLanczosSteps,0,opts);
%    psi_refnc(:,k) = V*(expm(-1i*dt*D)*(V'*psi_refnc(:,k-1)));
%    eig_refnc = diag(D);
    %psi_refnc(:,k) = expm(-1i*dt*H) * psi_refnc(:,k-1);
%end

% Reference eigenvalues
eig_refnc = eig(H);

% Reference solution: Matrix exponential
%disp('Reference solution.');
%time_refnc = cputime;
%for k=2:numTimeSteps+1
%    psi_refnc(:,k) = expm(-1i*dt*H) * psi_refnc(:,k-1);
%    fprintf('.');
%end
%fprintf('\n');
%time_refnc = cputime-time_refnc;

% Std Lanczos
disp('Standard Lanczos.');
time_stlan = cputime;
for k=2:numTimeSteps+1
    [T,Q] = lanczos_prop(H, psi_stlan(:,k-1), numLanczosSteps, dt, 1.0e-10, true);
    [V,D]=eig(T);
    nLanczos = size(T,1);
    psi_stlan(:,k) = ...
        Q*(V*(expm(-1i*dt*D)*(V'*(eye(nLanczos,1)*norm(psi_stlan(:,k-1))))));
    eig_stlan = diag(D);
    n_stlan = nLanczos;
    T_stlan = T;
    Q_stlan = Q;
    V_stlan = V;
    D_stlan = D;
end
time_stlan = cputime - time_stlan;
       
% CA-Lanczos
disp('CA-Lanczos.');
time_calan = cputime;
for k=2:numTimeSteps+1
    disp(['timestep ', num2str(k),'.']);
    [T,Q] = ca_lanczos_prop(H, psi_calan(:,k-1), s, numLanczosSteps/s, dt, 1.0e-10,'monomial', true);
    [V,D] = eig(T);
    % Investigate whether or not we should use the inverse or
    % the complex conjugate    |  here
    %                          V  
    nLanczos = size(T,1);
    psi_calan(:,k) = ...
        Q*(V*(expm(-1i*dt*D)*(V\(eye(nLanczos,1)*norm(psi_calan(:,k-1))))));
    eig_calan = diag(D);
    n_calan = nLanczos;
    T_calan = T;
    Q_calan = Q;
    V_calan = V;
    D_calan = D;
end
time_calan = cputime - time_calan;

% s-Step Lanczos
disp('s-step Lanczos.');
time_sslan = cputime;
for k=2:numTimeSteps+1
    [T,Q] = sstep_lanczos_prop(H, psi_sslan(:,k-1), s, numLanczosSteps/s, dt, 1.0e-10);
%    [T,Q] = sStepLanczos(H, psi_sslan(:,k-1), s, numLanczosSteps/s);
%    nLanczos = numLanczosSteps;
    [V,D] = eig(T);
    nLanczos = size(T,1);
    psi_sslan(:,k) = ...
        Q*(V*(expm(-1i*dt*D)*(V\(norm(psi_sslan(:,k-1))*eye(nLanczos,1)))));
    eig_sslan = diag(D);
    n_sslan = nLanczos;
    T_sslan = T;
    Q_sslan = Q;
    V_sslan = V;
    D_sslan = D;
end
time_sslan = cputime - time_sslan;

% Plot eigenvalues on the interval [0 1]
figure(1); hold on;
for k=2:numTimeSteps+1
    plot((0:N-1)/(N-1), sort(eig_refnc, 'descend'), 'k-');
    plot((0:n_stlan-1)/(n_stlan-1),sort(real(eig_stlan),'descend'),'rs');
    plot((0:n_calan-1)/(n_calan-1),sort(real(eig_calan),'descend'),'bo');
    plot((0:n_sslan-1)/(n_sslan-1),sort(real(eig_sslan),'descend'),'g*');
    
    legend('Reference','Lanczos','CA-Lanczos','s-step Lanczos'); %'s-Step Lanczos'
     
     % Compute errors in eigenvalues
     smalleig_err_stlan(k) = abs(min(abs(eig_refnc)) - min(abs(eig_stlan)));
     smalleig_err_calan(k) = abs(min(abs(eig_refnc)) - min(abs(eig_calan)));
     smalleig_err_sslan(k) = abs(min(abs(eig_refnc)) - min(abs(eig_sslan)));
     largeeig_err_stlan(k) = abs(max(abs(eig_refnc)) - max(abs(eig_stlan)));
     largeeig_err_calan(k) = abs(max(abs(eig_refnc)) - max(abs(eig_calan)));
     largeeig_err_sslan(k) = abs(max(abs(eig_refnc)) - max(abs(eig_sslan)));
 end
%disp('Matrix exponential:');
%disp(['Execution time:                ', num2str(time_refnc)]);
disp('Standard Lanczos:');
disp(['Execution time:                ', num2str(time_stlan)]);
disp(['Number of Lanczos iterations:  ', num2str(n_stlan)]);
disp(['Error in smallest eigenvalue:  ', num2str(max(smalleig_err_stlan))]);
disp(['Error in largest eigenvalue:   ', num2str(max(largeeig_err_stlan))]);
disp('CA-Lanczos:');
disp(['Execution time:                ', num2str(time_calan)]);
disp(['Number of Lanczos iterations:  ', num2str(n_calan)]);
disp(['Error in smallest eigenvalue:  ', num2str(max(smalleig_err_calan))]);
disp(['Error in largest eigenvalue:   ', num2str(max(largeeig_err_calan))]);
disp('s-Step Lanczos:');
disp(['Execution time:                ', num2str(time_sslan)]);
disp(['Number of Lanczos iterations:  ', num2str(n_sslan)]);
disp(['Error in smallest eigenvalue:  ', num2str(max(smalleig_err_sslan))]);
disp(['Error in largest eigenvalue:   ', num2str(max(largeeig_err_sslan))]);

figure;
subplot(1,3,1);
plot(x, real(psi_stlan(:,1)), x, imag(psi_stlan(:,1)), ...
     x, real(psi_stlan(:,numTimeSteps+1)), x, imag(psi_stlan(:,numTimeSteps+1)), ...
     x, abs(psi_stlan(:,numTimeSteps+1)));
legend('Re(\Psi_{0})', 'Im(\Psi_{0})', 'Re(\Psi_{k})', 'Im(\Psi_{k})', 'abs(\Psi_{k})');
subplot(1,3,2);
plot(x, real(psi_calan(:,1)), x, imag(psi_calan(:,1)), ...
     x, real(psi_calan(:,numTimeSteps+1)), x, imag(psi_calan(:,numTimeSteps+1)), ...
     x, abs(psi_calan(:,numTimeSteps+1)));
legend('Re(\Psi_{0})', 'Im(\Psi_{0})', 'Re(\Psi_{k})', 'Im(\Psi_{k})', 'abs(\Psi_{k})');
subplot(1,3,3);
plot(x, real(psi_sslan(:,1)), x, imag(psi_sslan(:,1)), ...
     x, real(psi_sslan(:,numTimeSteps+1)), x, imag(psi_sslan(:,numTimeSteps+1)), ...
     x, abs(psi_sslan(:,numTimeSteps+1)));
legend('Re(\Psi_{0})', 'Im(\Psi_{0})', 'Re(\Psi_{k})', 'Im(\Psi_{k})', 'abs(\Psi_{k})');
figure;

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

