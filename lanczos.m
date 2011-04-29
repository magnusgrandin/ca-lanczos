%%
%   function [T,V,rnorm,orthl] = lanczos(A,r0,maxiter,stop,orth)
%
%   The symmetric Lanczos algorithm
% 
%   Input:
%     A       - the large sparse symmetric matrix
%     r0      - initial vector
%     maxiter - mximum order of the Krylov space
%     stop    - break on convergence of any eigenpair [{0} | 1]
%     orth    - reorthogonalize lanczos vectors [{0} | 1]
%
%   Output variable:
%     T       - the tridiagonal lanczos matrix
%     V       - the orthonormal matrix of lanczos vectors
%     rnorm   - the residual norm in each step
%     orthl   - the level of orthogonality in each step (||I-Vj'*Vj||)
%--------------------------------------------------------------------------
function [T,Q,ritz_rnorm,orth_err] = lanczos(A,r,maxiter,orth)
 
    global g_lanczos_do_compute_ritz_rnorm;
    global g_lanczos_do_compute_orth_err;
    
    if nargin < 4
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

    % Check required output arguments
    g_lanczos_do_compute_ritz_rnorm = false;
    g_lanczos_do_compute_orth_err = false;
    if nargout >= 3
        g_lanczos_do_compute_ritz_rnorm = true;
    end
    if nargout >= 4 
        g_lanczos_do_compute_orth_err = true;
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

    % Normalize the starting vector
    q = r/norm(r);

    if strcmpi(orth,'local')
        [Q,T,ritz_rnorm,orth_err] = lanczos_basic(A, q, maxiter);
    elseif strcmpi(orth,'full')
        [Q,T,ritz_rnorm,orth_err] = lanczos_basic(A, q, maxiter,'fro');
    elseif strcmpi(orth,'periodic')
        [Q,T,ritz_rnorm,orth_err] = lanczos_periodic(A, q, maxiter);
    elseif strcmpi(orth,'selective')
        [Q,T,ritz_rnorm,orth_err] = lanczos_selective(A, q, maxiter);
    else
        % Do nothing
    end
end

function V = orthogonalize(V,iter,orth)
    % Reorthogonalization of basis vectors
    Rkk = V(:,1:iter)'*V(:,iter+1);
    V(:,iter+1) = V(:,iter+1) - V(:,1:iter)*Rkk;       
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

function [Q,T,rnorm,ortherr] = lanczos_basic(A,q,maxiter,orth)

    global g_lanczos_do_compute_ritz_rnorm;
    global g_lanczos_do_compute_orth_err;

    if nargin < 4
        orth = 'local';
    end
    
    n = length(q);
    Q = zeros(n,maxiter);
    Q(:,1) = q;
    alpha = zeros(1,maxiter);
    beta = zeros(1,maxiter);
    rnorm = zeros(maxiter,2);
    ortherr = zeros(maxiter,1);

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
        
        if strcmpi(orth,'fro') == 1
            Q = orthogonalize(Q,j,orth);
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
    rnorm = rnorm(1:j-1,:);
    ortherr = ortherr(1:j-1);
end

function [Q,T,rnorm,ortherr] = lanczos_selective(A,q,maxiter)

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
        
        T = diag(alpha(1:j)) + diag(beta(1:j-1),1) + diag(beta(1:j-1),-1);
        [Vp,Dp] = eig(T);
        nritz_new = 0;
        for i = 1:j
            if beta(i)*abs(Vp(j,i)) < norm_sqrt_eps
                nritz_new = nritz_new+1;
            end
        end
        if nritz_new > nritz
            nritz = 0;
            for i = 1:j
                if beta(i)*abs(Vp(j,i)) < norm_sqrt_eps
                    nritz = nritz+1;
                    y = Q(:,1:j)*Vp(:,i);
                    QR(:,nritz) = y;
                end
            end
        end       
        if nritz > 0
            Q(:,j+1) = projectAndNormalize({QR(:,1:nritz)},Q(:,j+1),false);
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
    rnorm = rnorm(1:j-1,:);
    ortherr = ortherr(1:j-1);
end
   
function [Q,T,rnorm,ortherr] = lanczos_periodic(A,q,maxiter)

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
            Q(:,j+1) = projectAndNormalize({Q(:,1:j)},Q(:,j+1),false);
            omega = reset_omega(omega, j, norm_A);
        end
        
        j = j+1;
    end
    
    T = diag(alpha(1:j-1)) + diag(beta(1:j-2),1) + diag(beta(1:j-2),-1);
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