function [Q,final_eigs,ritz_rnorm,orth_err] = restarted_lanczos(A, r, max_lanczos, n_wanted_eigs, orth, tol)

    restart_strategy = 'largest'; % 'largest','smallest','closest_conv','random'
    
    global g_lanczos_do_compute_ritz_rnorm;
    global g_lanczos_do_compute_orth_err;
    
    if nargin < 5
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
    if nargin < 6
        tol = 1.0e-06;
    end
    tol = tol*normest(A);


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

    % Vector length
    n = length(r);
        
    % Normalize the initial vector.
    q = r/norm(r);

    Q = zeros(n,max_lanczos);
    conv_eigs = [];
    conv_rnorms = [];
    
    restart = true;
    nconv = 0;
    while(restart)
               
        % Get the number of iterations to do next.
        iters = max_lanczos - nconv;
        
        Q_conv = Q(:,1:nconv);
               
        if strcmpi(orth,'local')
            [Q_new,T,ritz_rnorm,orth_err] = lanczos_basic(A, Q_conv, q, iters);
        elseif strcmpi(orth,'full')
            [Q_new,T,ritz_rnorm,orth_err] = lanczos_basic(A, Q_conv, q, iters, 'fro');
        elseif strcmpi(orth,'periodic')
            [Q_new,T,ritz_rnorm,orth_err] = lanczos_periodic(A, Q_conv, q, iters);
        elseif strcmpi(orth,'selective')
            [Q_new,T,ritz_rnorm,orth_err] = lanczos_selective(A, Q_conv, q, iters);
        else
            % Do nothing
        end
 
        % Compute residual norm estimates of all computed ritz-pairs.
        ritz_norms = zeros(iters,1);
        [Vp,Dp] = eig(T(1:iters,1:iters));
        beta = T(iters+1,iters);
        for i = 1:iters
            y = Vp(:,i);
            ritz_norms(i) = beta*abs(y(iters));
        end
        
        % Rearrange the eigenvalues such that the converged ones are first.
        k = 0;
        for i = 1:iters
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
            l = min(k+1,iters);
            for j = k+1:iters
                if Dp(j,j) > Dp(l,l)
                    l = j;
                end
            end
            Dp(l,l)
            q = Q_new*Vp(:,l);  
        elseif strcmpi(restart_strategy,'smallest')
            % Generate new starting vector from the largest non-converged basis vector.
            l = min(k+1,iters);
            for j = k+1:iters
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
            ix = min(k+1,iters);
            for i = k+2:iters
                if ritz_norms(i) < min_norm
                    min_norm = ritz_norms(i);
                    ix = i;
                end
            end
            q = Q_new*Vp(:,ix);
        elseif strcmpi(restart_strategy,'random')
            % Generate new starting vector from a random choice of the non-converged ones.
            ix = min((k+1) + round(rand(1)*(iters-k-1)), iters);
            q = Q_new*Vp(:,ix);
        end
        
        % Normalize the new starting vector
        q = q/norm(q);
                
        for i = 1:k
            Q(:,i+nconv) = Q_new*Vp(:,i);
            conv_eigs = [conv_eigs; Dp(i,i)];
            conv_rnorms = [conv_rnorms; ritz_norms(i)];
        end
        
        if strcmpi(orth,'local')
            q = projectAndNormalize({Q(:,1:nconv+k)},q,false);
        end
               
        % Update the count of converged eigenvalues
        nconv = nconv+k;
        
        if nconv == max_lanczos;
            % We have reached max number of Lanczos vectors. No further restarts.
            restart = false;
        else
            % Check if we should continue iterations
            restart = check_wanted_eigs(conv_eigs, diag(Dp(k+1:iters,k+1:iters)), n_wanted_eigs);
        end
        
        if ~restart
            sort_eigs = sort(conv_eigs,'descend');
            final_eigs = sort_eigs(1:n_wanted_eigs);
        end
    end
end

function restart = check_wanted_eigs(conv_eigs, eigs, num_wanted_eigs)
    max_eig = max(eigs);
    largest_conv_eigs = find(conv_eigs > max_eig);
    if length(largest_conv_eigs) >= num_wanted_eigs
        restart = false;
    else
        restart = true;
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

function v = orthogonalize(Q,v)
    if ~isempty(Q)
        Rkk = Q'*v;
        v = v - Q*Rkk;
    end
end

function [orth_err] = compute_orth_err(Q)      
    orth_err = norm(eye(size(Q,2))-Q'*Q,'fro');
end

function [Q,T,rnorm,ortherr] = lanczos_basic(A,Q_conv,q,maxiter,orth)

    global g_lanczos_do_compute_ritz_rnorm;
    global g_lanczos_do_compute_orth_err;

    if nargin < 5
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
            % Orthogonalize against the locked vectors and all previously
            % computed vectors in this Lanczos instance.
            Q(:,j+1) = orthogonalize(Q_conv,Q(:,j+1));
            Q(:,j+1) = orthogonalize(Q(:,1:j),Q(:,j+1));
        else
            % Only orthogonlize against the locked vectors.
            Q(:,j+1) = orthogonalize(Q_conv,Q(:,j+1));
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

function [Q,T,rnorm,ortherr] = lanczos_selective(A,Q_conv,q,maxiter)

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
   
function [Q,T,rnorm,ortherr] = lanczos_periodic(A,Q_conv,q,maxiter)

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