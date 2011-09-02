% TODO (110530): Remove spurious/repeated eigenvalues in local orth. 

function [conv_eigs,Q_conv,num_restarts,conv_rnorms,orth_err] = restarted_lanczos(A, r, max_lanczos, n_wanted_eigs, orth, tol)

    max_restarts = 100;

    % Check input arguments
    if nargin < 3
        disp('Usage:');
        disp('  [E,{V},{r},{o}] = restarted_ca_lanczos(A, r, max_iter, {n_wanted_eigs}, {s}, {basis}, {orth}, {tol})')
    end
    if nargin < 4 || isempty(n_wanted_eigs)
        n_wanted_eigs = 10;
    end
    if nargin < 5 || isempty(orth)
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
    norm_A = normest(A);
    if nargin < 6 || isempty(tol)
        tol = 1.0e-06;
    end
    tol = tol*norm_A;

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

    Q = zeros(n,max_lanczos);
    Q_conv = [];
    conv_eigs = [];
    conv_rnorms = [];
    orth_err = [];
    
    num_restarts = 0;
    restart = true;
    nconv = 0;
    while(restart && num_restarts < max_restarts)
               
        num_restarts = num_restarts + 1;
        
        % Get the number of iterations to do next.
        iters = max_lanczos - nconv;
        
        if strcmpi(orth,'local')
            [Q_new,T] = lanczos_basic(A, Q_conv, q, iters, 'local');
        elseif strcmpi(orth,'full')
            [Q_new,T] = lanczos_basic(A, Q_conv, q, iters, 'fro');
        elseif strcmpi(orth,'periodic')
            [Q_new,T] = lanczos_periodic(A, Q_conv, q, iters);
        elseif strcmpi(orth,'selective')
            [Q_new,T] = lanczos_selective(A, Q_conv, q, iters);
        else
            % Do nothing
        end
 
        % Update the maximum orthogonalization error
        if nargout >= 4
            Q_ = [Q_conv Q_new];
            orth_err = [orth_err; norm(eye(size(Q_,2))-Q_'*Q_,'fro')];
        end

        % Compute residual norm estimates of all computed ritz-pairs.
        ritz_norms = zeros(iters,1);
        [Vp,Dp] = eig(T(1:iters,1:iters));
        beta = T(iters+1,iters);
        for i = 1:iters
            y = Vp(:,i);
            ritz_norms(i) = beta*abs(y(iters)) + eps*norm_A;
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
        
        
        for i = 1:k
            Q(:,i+nconv) = Q_new*Vp(:,i);
            conv_eigs = [conv_eigs; Dp(i,i)];
            conv_rnorms = [conv_rnorms; ritz_norms(i)];
        end
        
        % Update the count of converged eigenvalues
        nconv = nconv+k;
        Q_conv = Q(:,1:nconv);
        
        restart = check_wanted_eigs(conv_eigs, diag(Dp(k+1:iters,k+1:iters)), n_wanted_eigs);

        if restart
            q = generateStartVector(diag(Dp),Vp,Q_new,ritz_norms,k,restart_strategy);                
            q = project({Q(:,1:nconv+k)},q,true);
        end
    end
    
    if ~restart
        [sort_eigs,ixs] = sort(conv_eigs,'descend');
        sort_rnorms = conv_rnorms(ixs);
        conv_eigs = sort_eigs(1:n_wanted_eigs);
        conv_rnorms = sort_rnorms(1:n_wanted_eigs);
        Q_conv = Q_conv(:,ixs);
        Q_conv = Q_conv(:,1:n_wanted_eigs);
        disp(['Converged in ' num2str(num_restarts) ' restarts.']);
        disp(['Max residual norm: ' num2str(max(conv_rnorms))]);
    else       
        disp(['Did not converge.']);
        conv_eigs = [];
        conv_rnorms = [];
        Q_conv = [];
        orth_err = [];
    end
end

function q = generateStartVector(eig_vals,eig_vecs,Q,ritz_norms,k,strategy)
    if nargin < 4
        strategy = 'largest';
    end
    m = length(eig_vals);
    if strcmpi(strategy,'largest')
        % Generate new starting vector from the largest non-converged basis vector.
        l = k+1;
        for j = k+1:m
            if eig_vals(j) > eig_vals(l)
                l = j;
            end
        end
        q = Q*eig_vecs(:,l);
    elseif strcmpi(strategy,'smallest')
        % Generate new starting vector from the largest non-converged basis vector.
        l = k+1;
        for j = k+1:m
            if eig_vals(j) < eig_vals(l,l)
                l = j;
            end
        end
        q = Q*eig_vecs(:,l);
    elseif strcmpi(strategy,'closest_conv')
        % Generate new starting vector from the non-converged basis vector that is
        % closest to convergence
        min_norm = inf;
        ix = k+1;
        for i = k+2:m
            if ritz_norms(i) < min_norm
                min_norm = ritz_norms(i);
                ix = i;
            end
        end
        q = Q*eig_vecs(:,ix);
    elseif strcmpi(strategy,'random')
        % Generate new starting vector from a random choice of the non-converged ones.
        ix = (k+1) + round(rand(1)*(m-k-1));
        q = Q*eig_vecs(:,ix);
    end
    
    % Normalize the new vector
    q = q/norm(q);
end

function restart = check_wanted_eigs(conv_eigs, eigs, num_wanted_eigs)
    if isempty(eigs)
        % Quick exit, no new eigs
        restart = true;
    else
%         max_eig = max(eigs);
%         largest_conv_eigs = find(conv_eigs > max_eig);
%         if length(largest_conv_eigs) >= num_wanted_eigs
%             restart = false;
%         else
%             restart = true;
%         end
        if length(conv_eigs) >= num_wanted_eigs
            restart = false;
        else
            restart = true;
        end
    end
end

function v = orthogonalize(Q,v)
    if ~isempty(Q)
        Rkk = Q'*v;
        v = v - Q*Rkk;
    end
end

function [Q,T] = lanczos_basic(A,Q_conv,q,maxiter,orth)

    if nargin < 5
        orth = 'local';
    end
    
    n = length(q);
    Q = zeros(n,maxiter);
    Q(:,1) = q;
    alpha = zeros(1,maxiter);
    beta = zeros(1,maxiter);

    has_converged = false;
    j = 1;
    
    while (j <= maxiter) && (has_converged == false)
        r = A*Q(:,j);
        if j > 1
            r=r-beta(j-1)*Q(:,j-1);
        end
        if strcmpi(orth,'local')
            [r,R_] = project({Q(:,j),Q_conv},r,true);
        elseif strcmpi(orth,'fro')
            [r,R_] = project({Q(:,j),Q_conv,Q(:,1:j)},r,true);
        end
        alpha(j) = R_{1};
        beta(j) = sqrt(r'*r);
        Q(:,j+1) = r/beta(j);
        j = j+1;
    end
    
    T = diag(alpha(1:j-1)) + diag(beta(1:j-2),1) + diag(beta(1:j-2),-1);
    T = [T; zeros(1,j-2) beta(j-1)];
    Q = Q(:,1:j-1);
end

function [Q,T] = lanczos_selective(A,Q_conv,q,maxiter)
    
    n = length(q);
    Q = zeros(n,maxiter);
    QR = zeros(n,maxiter);
    Q(:,1) = q;
    alpha = zeros(1,maxiter);
    beta = zeros(1,maxiter);
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
        [r,R_] = project({Q(:,j),Q_conv,QR(:,1:nritz)},r,true);
        alpha(j) = R_{1};
        beta(j) = sqrt(r'*r);
        Q(:,j+1) = r/beta(j);    

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
            QR(:,1:nritz) = normalize(QR(:,1:nritz));
        end
     
        j = j+1;
    end
    
    T = diag(alpha(1:j-1)) + diag(beta(1:j-2),1) + diag(beta(1:j-2),-1);
    T = [T; zeros(1,j-2) beta(j-1)];
    Q = Q(:,1:j-1);
end
   
function [Q,T] = lanczos_periodic(A,Q_conv,q,maxiter)

    n = length(q);
    Q = zeros(n,maxiter);
    Q(:,1) = q;
    alpha = zeros(1,maxiter);
    beta = zeros(1,maxiter);
    omega = [];
    norm_A = normest(A);
    
    has_converged = false;
    j = 1;
    
    while (j <= maxiter) && (has_converged == false)
        r = A*Q(:,j);
        if j > 1
            r=r-beta(j-1)*Q(:,j-1);
        end
        [r,R_] = project({Q(:,j),Q_conv},r,true);
        alpha(j) = R_{1};
        beta(j) = sqrt(r'*r);
        Q(:,j+1) = r/beta(j);

        
        % Estimate orthogonalization error, reorthogonalize if necessary
        omega = update_omega(omega, j, alpha, beta, norm_A);
        err = max(max(abs(omega - eye(size(omega)))));
        if err >= norm_A*sqrt(eps)
            Q(:,j:j+1) = projectAndNormalize({Q(:,1:j-1),Q_conv},Q(:,j:j+1),true);
            omega = reset_omega(omega, j, norm_A);
        end
        
        j = j+1;
    end
    
    T = diag(alpha(1:j-1)) + diag(beta(1:j-2),1) + diag(beta(1:j-2),-1);
    T = [T; zeros(1,j-2) beta(j-1)];
    Q = Q(:,1:j-1);
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