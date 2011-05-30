function [Q,final_eigs,ritz_rnorm,orth_err] = restarted_ca_lanczos(A, r, max_lanczos, s, n_wanted_eigs, basis, orth, tol)

    % Check input arguments
    if nargin < 4
        s = 6;
    end
    if nargin < 5
        n_wanted_eigs = 10;
    end
    if nargin < 6
        basis = 'newton';
    end
    if nargin < 7
        orth = 'local';
    else
        if isnumeric(orth)
            orth = num2str(orth);
        end
        if strcmpi(orth,'local')==0 && strcmpi(orth,'full')==0 ...
                && strcmpi(orth,'selective')==0 && strcmpi(orth,'periodic')==0
            disp(['lanczos.m: Invalid option value for orth: ', orth]);
            disp('    expected {''local''|''full''|''periodic''|''selective''}');
            return;
        end
    end
    if nargin < 8
        tol = 1.0e-06;
    end
    tol = tol*normest(A);

    % Check required output arguments
    do_ritz_norm = false;
    if nargout >= 3
        do_ritz_norm = true;
    end
    do_orth_err = false;
    if nargout >= 4
        do_orth_err = true;
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

    restart_strategy = 'largest'; % 'largest','smallest','closest_conv','random'

    % Vector length
    n = length(r);
        
    % Normalize the initial vector.
    q = r/norm(r);

    Q = zeros(n,max_lanczos+s);
    conv_eigs = [];
    conv_rnorms = [];
    
    restart = true;
    nconv = 0;
    while(restart)
               
        % Get the number of iterations to do next.
        iters = ceil((max_lanczos-nconv)/s);
        
        Q_conv = Q(:,1:nconv);
               
        if strcmpi(orth,'local')
            [Q_new,T,ritz_rnorm,orth_err] = lanczos_basic(A, Q_conv, q, iters, s, basis, 'local', do_ritz_norm, do_orth_err);
        elseif strcmpi(orth,'full')
            [Q_new,T,ritz_rnorm,orth_err] = lanczos_basic(A, Q_conv, q, iters, s, basis, 'fro', do_ritz_norm, do_orth_err);
        elseif strcmpi(orth,'periodic')
            [Q_new,T,ritz_rnorm,orth_err] = lanczos_periodic(A, Q_conv, q, iters, s, basis, do_ritz_norm, do_orth_err);
        elseif strcmpi(orth,'selective')
            [Q_new,T,ritz_rnorm,orth_err] = lanczos_selective(A, Q_conv, q, iters, s, basis, do_ritz_norm, do_orth_err);
        else
            % Do nothing
        end
 
        % Compute residual norm estimates of all computed ritz-pairs.
        ritz_norms = zeros(s*iters,1);
        [Vp,Dp] = eig(T(1:s*iters,1:s*iters));
        beta = T(s*iters+1,s*iters);
        for i = 1:s*iters
            y = Vp(:,i);
            ritz_norms(i) = beta*abs(y(s*iters));
        end
        
        % Rearrange the eigenvalues such that the converged ones are first.
        k = 0;
        for i = 1:s*iters
            if ritz_norms(i) < tol
                k = k + 1;
                % Push converged eigenvalues to the left.
                Dp_temp = Dp(i,i); Dp(i,i) = Dp(k,k); Dp(k,k) = Dp_temp;
                % Push corresponding eigenvectors to the left.
                Vp_temp = Vp(:,i); Vp(:,i) = Vp(:,k); Vp(:,k) = Vp_temp;
                % Push converged residual ritz norms to the left.
                rn_temp = ritz_norms(i); ritz_norms(i) = ritz_norms(k); ritz_norms(k) = rn_temp;
                ritz_norms(k)
                Dp(k,k)
            end
        end
        
        if strcmpi(restart_strategy,'largest')
            % Generate new starting vector from the largest non-converged basis vector.
            l = k+1;
            for j = k+1:s*iters
                if Dp(j,j) > Dp(l,l)
                    l = j;
                end
            end
            Dp(l,l)
            q = Q_new*Vp(:,l);  
        elseif strcmpi(restart_strategy,'smallest')
            % Generate new starting vector from the largest non-converged basis vector.
            l = k+1;
            for j = k+1:s*iters
                if Dp(j,j) < Dp(l,l)
                    l = j;
                end
            end
            Dp(l,l)
            q = Q_new*Vp(:,l);             
        elseif strcmpi(restart_strategy,'closest_conv')
            % Generate new starting vector from the non-converged basis vector that is 
            % closest to convergence
            min_norm = inf;
            ix = k+1;
            for i = k+2:s*iters
                if ritz_norms(i) < min_norm
                    min_norm = ritz_norms(i);
                    ix = i;
                end
            end
            q = Q_new*Vp(:,ix);
        elseif strcmpi(restart_strategy,'random')
            % Generate new starting vector from a random choice of the non-converged ones.
            ix = (k+1) + round(rand(1)*(s*iters-k-1));
            q = Q_new*Vp(:,ix);
        end
        
        % Normalize the new starting vector
        q = q/norm(q);
                
        for i = 1:k
            Q(:,i+nconv) = Q_new*Vp(:,i);
            conv_eigs = [conv_eigs; Dp(i,i)];
            conv_rnorms = [conv_rnorms; ritz_norms(i)];
        end
        
        %if strcmpi(orth,'local')
            q = projectAndNormalize({Q(:,1:nconv+k)},q,false);
        %end
               
        % Update the count of converged eigenvalues
        nconv = nconv+k;

        % Check if we should continue iterations
        restart = check_wanted_eigs(conv_eigs, diag(Dp(k+1:s*iters,k+1:s*iters)), n_wanted_eigs);
        if ~restart
            sort_eigs = sort(conv_eigs,'descend');
            final_eigs = sort_eigs(1:n_wanted_eigs);
        end
    end
end

function restart = check_wanted_eigs(conv_eigs, eigs, num_wanted_eigs)
    if isempty(eigs)
        % Quick exit, no new eigs
        restart = true;
    else
        max_eig = max(eigs);
        largest_conv_eigs = find(conv_eigs > max_eig);
        if length(largest_conv_eigs) >= num_wanted_eigs
            restart = false;
        else
            restart = true;
        end
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

%% Returns a columnvector with n elements, in which all elements are 
%% zero except the last element. 
function vec = eyeshvec(len)
    vec = eye(len,1);
    vec=circshift(vec,len-1);
end
 
function [Q,T,rnorm,ortherr] = lanczos_basic(A, Q_conv, q, maxiter, s, basis, orth, do_ritz_norm, do_orth_err)

    if nargin < 7
        orth = 'local';
    end
  
    max_lanczos = maxiter*s;
    n = length(q);
    Q = zeros(n,maxiter);
    Q(:,1) = q;
    b = zeros(maxiter+1,1);
    T = [];
    rnorm = zeros(maxiter,2);
    ortherr = zeros(maxiter,1);
    Bk = [];
    
    % Fix change-of-basis matrix
    if strcmpi(basis,'monomial')
        I = eye(s+1);
        Bk = I(:,2:s+1);        
    elseif strcmpi(basis,'newton')
        % Run standard Lanczos for 2s steps
        T = lanczos(A,q,2*s,'full');
        basis_eigs = eig(T);
        basis_shifts = leja(basis_eigs,'nonmodified');
        Bk = newton_basis_matrix(basis_shifts, s,1);
    end
    
    has_converged = false;
    k = 0;
    while (k <= maxiter) && (has_converged == false)

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
                [Q_] = projectAndNormalize({Q_conv},Q_,false);
                Q(:,(k-1)*s+2:k*s+1) = Q_(:,1:s);
            elseif strcmpi(orth,'fro')
                % Orthogonality against all previous basis vectors
                [Q_,Rk_] = projectAndNormalize({Q(:,(k-2)*s+1:(k-1)*s+1)},V(:,2:s+1),false);
                Rkk_s = Rk_{1};
                Rk_s = Rk_{2};
                [Q_] = projectAndNormalize({Q_conv,Q},Q_);
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
        if do_ritz_norm
            [Vp,Dp] = eig(T(1:s*k,1:s*k));
            rnorm(k,:) = compute_ritz_rnorm(A,Q(:,1:s*k),Vp,Dp);
        end
        
        % Compute the orthogonalization error, if it is required
        if do_orth_err
            ortherr(k) = compute_orth_err(Q);
        end
    end
    
    % Fix output
    T = T(1:s*(k-1)+1,1:s*(k-1));
    Q = Q(:,1:s*(k-1));
    rnorm = rnorm(1:(k-1),:);
    ortherr = ortherr(1:(k-1));
end

function [Q,T,rnorm,ortherr] = lanczos_selective(A, Q_conv, q, maxiter, s, basis, do_ritz_norm, do_orth_err)

    global g_lanczos_do_compute_ritz_rnorm;
    global g_lanczos_do_compute_orth_err;
    
    n = length(q);
    Q = zeros(n,maxiter);
    QR = zeros(n,maxiter);
    Q(:,1) = q;
    alpha = zeros(1,maxiter);
    beta = zeros(1,maxiter);
    rnorm = zeros(maxiter,2);
    ortherr = zeros(maxiter,1);
    norm_A = normest(A);
    norm_sqrt_eps = norm_A*sqrt(eps);
    nritz = 0;
    
    has_converged = false;
    j = 1;
    
    while (j <= maxiter) && (has_converged == false)
        r = A*Q(:,j);
        if j > 1
            r=r-beta(j-1)*Q(:,j-1);
        end
        alpha(j) = r'*Q(:,j);
        r = r - alpha(j)*Q(:,j);
        beta(j) = sqrt(r'*r);
        Q(:,j+1) = r/beta(j);
     
        % Orthogonalize against the locked basis vectors
        Q(:,j+1) = orthogonalize(Q_conv,Q(:,j+1));

        T = diag(alpha(1:j)) + diag(beta(1:j-1),1) + diag(beta(1:j-1),-1);
        [Vp,Dp] = eig(T);
        nritz_new = 0;
        for i = 1:j
            if beta(j)*abs(Vp(j,i)) < norm_sqrt_eps
                nritz_new = nritz_new+1;
            end
        end
        if nritz_new > nritz
            nritz = 0;
            for i = 1:j
                if beta(j)*abs(Vp(j,i)) < norm_sqrt_eps
                    nritz = nritz+1;
                    y = Q(:,1:j)*Vp(:,i);
                    QR(:,nritz) = y;
                end
            end
        end       
        if nritz > 0
            Q(:,j+1) = orthogonalize(QR(:,1:nritz),Q(:,j+1));
        end
     
        % Compute the ritz-norm, if it is required
        if g_lanczos_do_compute_ritz_rnorm && j > 2
            T = diag(alpha(1:j-1)) + diag(beta(1:j-2),1) + diag(beta(1:j-2),-1);
            [Vp,Dp] = eig(T);
            rnorm(j,:) = compute_ritz_rnorm(A,Q(:,1:j-1),Vp,Dp);
        end
        
        % Compute the orthogonalization error, if it is required
        if g_lanczos_do_compute_orth_err && j > 2
            ortherr(j) = compute_orth_err(Q(:,1:j-1));
        end
                       
        j = j+1;
    end
    
    T = diag(alpha(1:j-1)) + diag(beta(1:j-2),1) + diag(beta(1:j-2),-1);
    T = [T; zeros(1,j-2) beta(j-1)];
    Q = Q(:,1:j-1);
    rnorm = rnorm(1:j-1,:);
    ortherr = ortherr(1:j-1);
end
   
function [Q,T,rnorm,ortherr] = lanczos_periodic(A, Q_conv, q, maxiter, s, basis, do_ritz_norm, do_orth_err)

    global g_lanczos_do_compute_ritz_rnorm;
    global g_lanczos_do_compute_orth_err;

    n = length(q);
    Q = zeros(n,maxiter);
    Q(:,1) = q;
    alpha = zeros(1,maxiter);
    beta = zeros(1,maxiter);
    rnorm = zeros(maxiter,2);
    ortherr = zeros(maxiter,1);
    omega = [];
    norm_A = normest(A);
    
    has_converged = false;
    j = 1;
    
    while (j <= maxiter) && (has_converged == false)
        r = A*Q(:,j);
        if j > 1
            r=r-beta(j-1)*Q(:,j-1);
        end
        alpha(j) = r'*Q(:,j);
        r = r - alpha(j)*Q(:,j);
        beta(j) = sqrt(r'*r);
        Q(:,j+1) = r/beta(j);

        Q(:,j+1) = project({Q_conv},Q(:,j+1));

        % Compute the ritz-norm, if it is required
        if g_lanczos_do_compute_ritz_rnorm && j > 2
            T = diag(alpha(1:j-1)) + diag(beta(1:j-2),1) + diag(beta(1:j-2),-1);
            [Vp,Dp] = eig(T);
            rnorm(j,:) = compute_ritz_rnorm(A,Q(:,1:j-1),Vp,Dp);
        end
        
        % Compute the orthogonalization error, if it is required
        if g_lanczos_do_compute_orth_err && j > 2
            ortherr(j) = compute_orth_err(Q(:,1:j-1));
        end
                       
        % Estimate orthogonalization error, reorthogonalize if necessary
        omega = update_omega(omega, j, alpha, beta, norm_A);
        err = max(max(abs(omega - eye(size(omega)))));
        if err >= norm_A*sqrt(eps)
            Q(:,j+1) = project({Q(:,1:j)},Q(:,j+1));
            omega = reset_omega(omega, j, norm_A);
        end
        
        j = j+1;
    end
    
    T = diag(alpha(1:j-1)) + diag(beta(1:j-2),1) + diag(beta(1:j-2),-1);
    T = [T; zeros(1,j-2) beta(j-1)];
    Q = Q(:,1:j-1);
    rnorm = rnorm(1:j-1,:);
    ortherr = ortherr(1:j-1);
end

function omega = update_omega(omega_in, n, alpha, beta, anorm)

    % Estimate of contribution to roundoff errors from A*v:  fl(A*v) = A*v + f, 
    T = eps*anorm;
    binv = 1.0/beta(n);
    
    omega = zeros(n+1,n+1);
    
    % In the first iteration, just initiate omega
    if isempty(omega_in) 
        omega(1,1) = 1;
        omega(1,2) = 0;
        omega(2,1) = binv*T;
        omega(2,2) = 1;
    else
        omega(1:n,1:n) = omega_in;
        omega(n+1,1) = beta(2)*omega(n,2) + (alpha(1)-alpha(n))*omega(n,1) - beta(n)*omega(n-1,1);
        if omega(n+1,1) > 0
            omega(n+1,1) = binv*(omega(n+1,1) + T);
        else
            omega(n+1,1) = binv*(omega(n+1,1) - T);
        end
        for k=2:n-1
            omega(n+1,k) = beta(k+1)*omega(n,k+1) + (alpha(k)-alpha(n))*omega(n,k) + beta(k)*omega(n,k-1) - beta(n)*omega(n-1,k);
            if omega(n+1,k) > 0 
                omega(n+1,k) = binv*(omega(n+1,k) + T);
            else
                omega(n+1,k) = binv*(omega(n+1,k) - T);
            end
        end
        omega(n+1,n) = binv*T;
        omega(n+1,n+1) = 1;
    end
end

function omega = reset_omega(omega_in, n, anorm)
    omega = omega_in;
    T = eps*anorm;
    for k = 1:n
        omega(n+1,k) = T;
    end
    omega(n+1,n+1) = 1;
end