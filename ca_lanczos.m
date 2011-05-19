%%
%   function [T,Q,rnorm,orthl] = ca_lanczos(A,r,s,t,basis,stop,orth)
%   
%   Communication avoiding Lanczos algorithm, as described in 
%   M. Hoemmen, Communication-Avoiding Krylov Subspace Methods, PhD thesis,
%   University of California Berkeley (2010).
%   
%   Input:
%     A     - the large sparse symmetric matrix
%     r     - initial vector
%     s     - SpMV kernel step size
%     t     - number of iterations (restart length = s*t)
%     basis - basis to use {'monomial'|'newton'}
%     stop  - flag to tell whether we should stop on convergence of first
%             Ritz-pair (optional) {0|1}
%     orth  - orthogonalization strategy to use (optional)
%             {'local'|'full'|'periodic'|'select'}
%     
%   Output:
%     T     - Lanczos projection matrix [(s*t) x (s*t)]
%     Q     - basis vectors [n x (s*t)]
%     rnorm - vector of the residual norms in each iteration (optional)
%     orthl - vector of level of orthogonality in each iteration (optional)
%
function [T,Q,ritz_rnorm,orth_err] = ca_lanczos(A,r,s,t,basis,orth)

    global g_ca_lanczos_do_compute_ritz_rnorm;
    global g_ca_lanczos_do_compute_orth_err;
    
    % Check input arguments
    if nargin < 6
        orth = 'local';
    else
        if isnumeric(orth)
            orth = num2str(orth);
        end
        if strcmpi(orth,'local')==0 && strcmpi(orth,'full')==0 ...
                && strcmpi(orth,'selective')==0 && strcmpi(orth,'periodic')==0
            disp(['ca_lanczos.m: Invalid option value for orth: ', orth]);
            disp('    expected {''local''|''full''|''periodic''|''selective''}');
            return;
        end
    end
   
    % Check required output arguments
    g_ca_lanczos_do_compute_ritz_rnorm = false;
    g_ca_lanczos_do_compute_orth_err = false;
    if nargout >= 3
        g_ca_lanczos_do_compute_ritz_rnorm = true;
    end
    if nargout >= 4 
        g_ca_lanczos_do_compute_orth_err = true;
    end
 
    if strcmpi(orth,'local')
        disp('Local orthogonalization');
    elseif strcmpi(orth,'full')
        disp('Full orthogonalization');
    elseif strcmpi(orth,'periodic')
        disp('Periodic orthogonalization');
    elseif strcmpi(orth,'selective')
        disp('Selective orthogonalization');
    end

    % Make sure that t (number of outer iterations) is an integer
    t = ceil(t);
    
    % Normalize the starting vector
    q = r/norm(r);
    
    % Fix change-of-basis matrix
    if strcmpi(basis,'monomial')
        I = eye(s+1);
        Bk = I(:,2:s+1);        
    elseif strcmpi(basis,'newton')
        % Run standard Lanczos for 2s steps
        T = lanczos(A,r,2*s,'full');
        basis_eigs = eig(T);
        basis_shifts = leja(basis_eigs,'nonmodified');
        Bk = newton_basis_matrix(basis_shifts, s,1);
    else
        disp(['ERROR: Unknown basis type: ', basis]);
    end
    
    if strcmpi(orth,'local')
        [Q,T,ritz_rnorm,orth_err] = ca_lanczos_basic(A, q, s, t, Bk, basis);
    elseif strcmpi(orth,'full')
        [Q,T,ritz_rnorm,orth_err] = ca_lanczos_basic(A, q, s, t, Bk, basis,'fro');
    elseif strcmpi(orth,'periodic')
        [Q,T,ritz_rnorm,orth_err] = ca_lanczos_periodic(A, q, s, t, Bk, basis);
    elseif strcmpi(orth,'selective')
        [Q,T,ritz_rnorm,orth_err] = ca_lanczos_selective(A, q, s, t, Bk, basis);
    else
        % Do nothing
    end
end

function [ritz_rnorm] = compute_ritz_rnorm(A,Q,Vp,Dp)
    ritz_rnorm = [];
    % Residual norm for smallest eigenpair
    [d_s,i_s] = min(diag(Dp));
    x_s = Q*Vp(:,i_s);
    ritz_rnorm(1) = norm(A*x_s-d_s*x_s)/norm(d_s*x_s);
    % Residual norm for largest eigenpair
    [d_l,i_l] = max(diag(Dp));
    x_l = Q*Vp(:,i_l);
    ritz_rnorm(2) = norm(A*x_l-d_l*x_l)/norm(d_l*x_l);
end

function [orth_err] = compute_orth_err(Q)      
    orth_err = norm(eye(size(Q,2))-Q'*Q,'fro');
end

% Compute matrix powers
function V = matrix_powers(A, q, s, Bk, basis)
    if strcmpi(basis,'monomial')
        V(:,1) = q;
        V(:,2:s+1) = matrix_powers_monomial(A, q, s);
    elseif strcmpi(basis,'newton')
        basis_shifts = diag(Bk);
        V(:,1:s+1) = matrix_powers_newton(A, q, s, basis_shifts, 1);
    end
end

function Q = reorthogonalize(Q)
    M = length(Q);
    for i=1:M
        Q_ = Q{i};
        [Q_,Rk_] = projectAndNormalize(Q(1:i-1), Q_);
        Q{i} = Q_;
    end
end

function y = blockVectorTimesVector(Q, v)
    n = size(Q{1},1);
    y = zeros(n,1);
    cols = 0;
    for l=1:length(Q)-1
        ncols = size(Q{l},2);
        y = y + Q{l}*v(cols+1:cols+ncols);
        cols = cols+ncols;
    end
    ncols = size(Q{end},2);
    y = y + Q{end}(:,1:ncols-1)*v(cols+1:cols+ncols-1);
end

%% Returns a columnvector with n elements, in which all elements are 
%% zero except the last element. 
function vec = eyeshvec(len)
    vec = eye(len,1);
    vec=circshift(vec,len-1);
end
 
%% CA-Lanczos with local or full orthogonalization.
function [Q,T,rnorm,ortherr] = ca_lanczos_basic(A, q, s, t, Bk, basis, orth)
    
    global g_ca_lanczos_do_compute_ritz_rnorm;
    global g_ca_lanczos_do_compute_orth_err;
    
    if nargin < 7
        orth = 'local';
    end
    
    b = zeros(t+1,1);
    rnorm = zeros(t,2);
    ortherr = zeros(t,1);   

    has_converged = false;
    k = 0;
    while (k <= t) && (has_converged == false)

        k = k+1;

        if k > 1
            q = Q(:,(k-1)*s+1);
        end
        
        V = matrix_powers(A,q,s,Bk,basis);
        
        if k == 1
            % QR factorization
            [Q(:,1:s+1),Rk] = normalize(V(:,1:s+1));
            % Compute first part of tridiagonal matrix
            T = Rk*Bk/Rk(1:s,1:s);
            % Compute next beta
            b(k) = T(s+1,s);
        else
            if strcmpi(orth,'local')
                % Orthogonalize against previous block of basis vectors
                [Q_,Rk_] = projectAndNormalize({Q(:,(k-2)*s+1:(k-1)*s+1)},V(:,2:s+1),false);
                Q(:,(k-1)*s+2:k*s+1) = Q_(:,1:s);
                Rkk_s = Rk_{1};
                Rk_s = Rk_{2};
            elseif strcmpi(orth,'fro')
                % Orthogonality against all previous basis vectors
                [Q_,Rk_] = projectAndNormalize({Q(:,(k-2)*s+1:(k-1)*s+1)},V(:,2:s+1),false);
                Rkk_s = Rk_{1};
                Rk_s = Rk_{2};
                [Q_] = projectAndNormalize({Q},Q_,true);
                Q(:,(k-1)*s+2:k*s+1) = Q_(:,1:s);
            end
            
            % Compute Tk (tridiagonal sub-matrix of T)
            Rkk = [zeros(s,1), Rkk_s(1:s,:)];
            Rk  = [eye(s+1,1), [Rkk_s(s+1,1:s);Rk_s]];
            zk    = Rk(1:s,s+1);
            rho   = Rk(s+1,s+1);
            rho_t = Rk(s,s);
            bk = Bk(s+1,s); %should be one with monomial basis
            e1 = eye(s,1);
            es = eyeshvec(s);
            Tk = Rk(1:s,1:s)*Bk(1:s,:)/Rk(1:s,1:s) ...
                + (bk/rho_t)*zk*es' ...
                - b(k-1)*e1*es'*Rkk(1:s,1:s)/Rk(1:s,1:s);

            % Compute the next beta
            b(k) = bk*(rho/rho_t);
            
            % Extend T
            T11 = T(1:s*(k-1),1:s*(k-1));
            T12 = b(k-1)*eyeshvec(s*(k-1))*eye(s,1)';
            T21 = b(k-1)*eye(s,1)*eyeshvec(s*(k-1))';
            T22 = Tk;
            T31 = zeros(1,s*(k-1));
            T32 = b(k)*eyeshvec(s)';
            T = [T11, T12; T21, T22; T31, T32];
 
        end
        
        % Compute the ritz-norm, if it is required
        if g_ca_lanczos_do_compute_ritz_rnorm
            [Vp,Dp] = eig(T(1:s*k,1:s*k));
            rnorm(k,:) = compute_ritz_rnorm(A,Q(:,1:s*k),Vp,Dp);
        end
        
        % Compute the orthogonalization error, if it is required
        if g_ca_lanczos_do_compute_orth_err
            ortherr(k) = compute_orth_err(Q);
        end
    end
    
    % Fix output
    T = T(1:s*(k-1),1:s*(k-1));
    Q = Q(:,1:s*(k-1));
    rnorm = rnorm(1:(k-1),:);
    ortherr = ortherr(1:(k-1));
    
end

function [Q,T,rnorm,ortherr] = ca_lanczos_selective(A, q, s, t, Bk, basis)

    global g_ca_lanczos_do_compute_ritz_rnorm;
    global g_ca_lanczos_do_compute_orth_err;
    
    norm_A = normest(A);
    norm_sqrt_eps = norm_A*sqrt(eps); %(1/s^2)
    n = size(A,1);
    b = zeros(t+1,1);
    nritz = 0;
    QR = zeros(n,s*t);
    rnorm = zeros(t,2);
    ortherr = zeros(t,1);
    omega = [];
    
    has_converged = false;
    k = 0;
    
    while (k <= t) && (has_converged == false)

        k = k+1;
        if k > 1
            q = Q(:,(k-1)*s+1);
        end
        
        V = matrix_powers(A,q,s,Bk,basis);
        
        if k == 1
            % QR factorization
            [Q(:,1:s+1),Rk] = normalize(V(:,1:s+1));
            % Compute first part of tridiagonal matrix
            T = Rk*Bk/Rk(1:s,1:s);
            % Compute next beta
            b(k) = T(s+1,s);
            
        else
            % Orthogonalize against previous block of basis vectors and the
            % already converged ritz vectors
            [Q_,Rk_] = projectAndNormalize({Q(:,(k-2)*s+1:(k-1)*s+1),QR(:,1:nritz)},V(:,2:s+1),false);
            Q(:,(k-1)*s+2:k*s+1) = Q_(:,1:s);
            Rkk_s = Rk_{1};
            Rk_s = Rk_{3};
%             if nritz > 0
%                 Q(:,(k-1)*s+2:k*s+1) = projectAndNormalize({QR(:,1:nritz)},Q(:,(k-1)*s+2:k*s+1),false);
%             end
            
            % Compute Tk (tridiagonal sub-matrix of T)
            Rkk = [zeros(s,1), Rkk_s(1:s,:)];
            Rk = [eye(s+1,1), [Rkk_s(s+1,1:s);Rk_s]];
            zk = Rk(1:s,s+1);
            rho = Rk(s+1,s+1);
            rho_t = Rk(s,s);
            bk = Bk(s+1,s); %should be one with monomial basis
            e1 = eye(s,1);
            es = eyeshvec(s);
            Tk = Rk(1:s,1:s)*Bk(1:s,:)/Rk(1:s,1:s) ...  %% <-- Warnings for bad scaling of matrix
                + (bk/rho_t)*zk*es' ...
                - b(k-1)*e1*es'*Rkk(1:s,1:s)/Rk(1:s,1:s);

            % Compute the next beta
            b(k) = bk*(rho/rho_t);
            
            % Extend T
            T11 = T(1:s*(k-1),1:s*(k-1));
            T12 = b(k-1)*eyeshvec(s*(k-1))*eye(s,1)';
            T21 = b(k-1)*eye(s,1)*eyeshvec(s*(k-1))';
            T22 = Tk;
            T31 = zeros(1,s*(k-1));
            T32 = b(k)*eyeshvec(s)';
            T = [T11, T12; T21, T22; T31, T32];
        end
        
%         Estimate orthogonalization error, reorthogonalize if necessary
%         alpha = diag(T,0);
%         beta  = diag(T,-1);
%         omega = update_omega(omega,alpha,beta,norm_A, s);
%         err = max(max(abs(omega - eye(size(omega)))));
%         if err >= norm_A*sqrt(eps)
%             disp('.')
%             nritz = 0;
%             [Vp,Dp] = eig(T(1:s*k,1:s*k));
%             for i = 1:k*s
%                 if beta(i)*abs(Vp(s*k,i)) < norm_sqrt_eps
%                     nritz = nritz+1;
%                     y = Q(:,1:end-1)*Vp(:,i);
%                     QR(:,nritz) = y;
%                 end
%             end
%             QR(:,1:nritz) = normalize(QR(:,1:nritz));
%             omega = reset_omega(omega, norm_A, s);
%         end
        
        
        beta = diag(T,-1);
        [Vp,Dp] = eig(T(1:s*k,1:s*k));
        
        nritz_new = 0;
        for i = 1:k*s
            if beta(i)*abs(Vp(s*k,i)) < norm_sqrt_eps
                nritz_new = nritz_new+1;
            end
        end
        if nritz_new > nritz
            nritz = 0;
            for i = 1:k*s
                if beta(i)*abs(Vp(s*k,i)) < norm_sqrt_eps
                    nritz = nritz+1;
                    y = Q(:,1:end-1)*Vp(:,i);
                    QR(:,nritz) = y;
                end
            end
            QR(:,1:nritz) = normalize(QR(:,1:nritz));
        end

        disp(['nritz=' num2str(nritz)])
       
%         if nritz > 0
%             Q(:,(k-1)*s+2:k*s+1) = projectAndNormalize({QR(:,1:nritz)},Q(:,(k-1)*s+2:k*s+1),false);
%         end
        
        % Compute the ritz-norm, if it is required
        if g_ca_lanczos_do_compute_ritz_rnorm
            %[Vp,Dp] = eig(T(1:s*k,1:s*k));
            rnorm(k,:) = compute_ritz_rnorm(A,Q(:,1:s*k),Vp,Dp);
        end
        
        % Compute the orthogonalization error, if it is required
        if g_ca_lanczos_do_compute_orth_err
            ortherr(k) = compute_orth_err(Q);
        end
    end
      
    % Fix output
    T = T(1:s*(k-1),1:s*(k-1));
    Q = Q(:,1:s*(k-1));
    rnorm = rnorm(1:(k-1),:);
    ortherr = ortherr(1:(k-1));
end

function [Q,T,rnorm,ortherr] = ca_lanczos_periodic(A, q, s, t, Bk, basis)

    global g_ca_lanczos_do_compute_ritz_rnorm;
    global g_ca_lanczos_do_compute_orth_err;
    
    b = zeros(t+1,1);
    rnorm = zeros(t,2);
    ortherr = zeros(t,1);
    omega = [];
    norm_A = normest(A);
    
    has_converged = false;
    k = 0;

    while (k <= t) && (has_converged == false)

        k = k+1;

        if k > 1
            q = Q(:,(k-1)*s+1);
        end
        
        V = matrix_powers(A,q,s,Bk,basis);
        
        if k == 1
            % QR factorization
            [Q(:,1:s+1),Rk] = normalize(V(:,1:s+1));
            % Compute first part of tridiagonal matrix
            T = Rk*Bk/Rk(1:s,1:s);
            % Compute next beta
            b(k) = T(s+1,s);
            
        else
            % Orthogonalize against previous block of basis vectors
            [Q_,Rk_] = projectAndNormalize({Q(:,(k-2)*s+1:(k-1)*s+1)},V(:,2:s+1),true);
            Q(:,(k-1)*s+2:k*s+1) = Q_(:,1:s);
            Rkk_s = Rk_{1};
            Rk_s = Rk_{2};
            
            % Compute Tk (tridiagonal sub-matrix of T)
            Rkk = [zeros(s,1), Rkk_s(1:s,:)];
            Rk = [eye(s+1,1), [Rkk_s(s+1,1:s);Rk_s]];
            zk = Rk(1:s,s+1);
            rho = Rk(s+1,s+1);
            rho_t = Rk(s,s);
            bk = Bk(s+1,s); %should be one with monomial basis
            e1 = eye(s,1);
            es = eyeshvec(s);
            Tk = Rk(1:s,1:s)*Bk(1:s,:)/Rk(1:s,1:s) ...
                + (bk/rho_t)*zk*es' ...
                - b(k-1)*e1*es'*Rkk(1:s,1:s)/Rk(1:s,1:s);

            % Compute the next beta
            b(k) = bk*(rho/rho_t);
            
            % Extend T
            T11 = T(1:s*(k-1),1:s*(k-1));
            T12 = b(k-1)*eyeshvec(s*(k-1))*eye(s,1)';
            T21 = b(k-1)*eye(s,1)*eyeshvec(s*(k-1))';
            T22 = Tk;
            T31 = zeros(1,s*(k-1));
            T32 = b(k)*eyeshvec(s)';
            T = [T11, T12; T21, T22; T31, T32];
        end
        
        % Compute the ritz-norm, if it is required
        if g_ca_lanczos_do_compute_ritz_rnorm
            [Vp,Dp] = eig(T(1:s*k,1:s*k));
            rnorm(k,:) = compute_ritz_rnorm(A,Q(:,1:s*k),Vp,Dp);
        end
        
        % Compute the orthogonalization error, if it is required
        if g_ca_lanczos_do_compute_orth_err
            ortherr(k) = compute_orth_err(Q);
        end
        
        % Estimate orthogonalization error, reorthogonalize if necessary
        alpha = diag(T,0);
        beta  = diag(T,-1);
        omega = update_omega(omega,alpha,beta,norm_A, s);
        err = max(max(abs(omega - eye(size(omega)))));
        if err >= norm_A*sqrt(eps)
            Q(:,(k-1)*s+1:k*s+1) = projectAndNormalize({Q(:,1:(k-1)*s)},Q(:,(k-1)*s+1:k*s+1),true);
            %R = Q(:,1:(k-1)*s+1)'*Q(:,(k-1)*s+2:k*s+1);
            %Q(:,(k-1)*s+2:k*s+1) = Q(:,(k-1)*s+2:k*s+1) - Q(:,1:(k-1)*s+1)*R;
            omega = reset_omega(omega, norm_A, s);
        end
    end
      
    % Fix output
    T = T(1:s*(k-1),1:s*(k-1));
    Q = Q(:,1:s*(k-1));
    rnorm = rnorm(1:(k-1),:);
    ortherr = ortherr(1:(k-1));

end

function omega = update_omega(omega_in, alpha, beta, anorm, s)

    % Get iteration number and block size
    n = length(alpha);
    m = size(omega_in,1)-1;
    
    % Estimate of contribution to roundoff errors from A*v:  fl(A*v) = A*v + f, 
    T = eps*anorm;
    binv = 1.0/beta(n);
    
    if isempty(omega_in) 
        omega = zeros(s+1,s+1);
        omega(1,1) = 1;
        omega(1,2) = 0;
        omega(2,1) = binv*T;
        omega(2,2) = 1;

        for j = 2:s
            % k == 1, omega(j,k-1) == 0
            omega(j+1,1) = beta(2)*omega(j,2) + (alpha(1)-alpha(j))*omega(j,1) - beta(j)*omega(j-1,1);
            if omega(j+1,1) > 0
                omega(j+1,1) = binv*(omega(j+1,1) + T);
            else
                omega(j+1,1) = binv*(omega(j+1,1) - T);
            end
            
            % Update remaining components.
            for k=2:j-1
                omega(j+1,k) = beta(k+1)*omega(j,k+1) + (alpha(k)-alpha(j))*omega(j,k) + beta(k)*omega(j,k-1) - beta(j)*omega(j-1,k);
                if omega(j+1,k) > 0 
                    omega(j+1,k) = binv*(omega(j+1,k) + T);
                else
                    omega(j+1,k) = binv*(omega(j+1,k) - T);
                end
            end
            
            % k == j, k == j+1
            omega(j+1,j) = binv*T;
            omega(j+1,j+1) = 1;
        end
        
    else
        omega = zeros(n+1,n+1);
        omega(1:m+1,1:m+1) = omega_in;
        
        for j = m+1:m+s
            % k == 1, omega(j,k-1) == 0
            omega(j+1,1) = beta(2)*omega(j,2) + (alpha(1)-alpha(j))*omega(j,1) - beta(j)*omega(j-1,1);
            if omega(j+1,1) > 0
                omega(j+1,1) = binv*(omega(j+1,1) + T);
            else
                omega(j+1,1) = binv*(omega(j+1,1) - T);
            end
            
            % Update remaining components.
            for k=2:j-1
                omega(j+1,k) = beta(k+1)*omega(j,k+1) + (alpha(k)-alpha(j))*omega(j,k) + beta(k)*omega(j,k-1) - beta(j)*omega(j-1,k);
                if omega(j+1,k) > 0 
                    omega(j+1,k) = binv*(omega(j+1,k) + T);
                else
                    omega(j+1,k) = binv*(omega(j+1,k) - T);
                end
            end
            
            % k == j, k == j+1
            omega(j+1,j) = binv*T;
            omega(j+1,j+1) = 1;
        end
    end
end

function omega = reset_omega(omega_in, anorm, s)
    T = eps*anorm;
    m = size(omega_in,1)-s-1;
    omega = omega_in;
    for j = m+1:m+s
        for k = 1:j
            omega(j+1,k) = T;
        end
        omega(j+1,j+1) = 1;
    end   
end
