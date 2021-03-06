%% function [conv_eigs,Q_conv,num_restarts,rnorms,orth_err] = 
%    restarted_ca_lanczos(A, r, max_lanczos, n_wanted_eigs, s, basis, orth, tol)

function [conv_eigs,Q_conv,num_restarts,rnorms,orth_err] = restarted_ca_lanczos(A, r, max_lanczos, n_wanted_eigs, s, basis, orth, tol)

    max_restarts = 200;

    % Check input arguments
    if nargin < 3
        disp('Usage:');
        disp('  [E,{V},{r},{o}] = restarted_ca_lanczos(A, r, max_iter, {n_wanted_eigs}, {s}, {basis}, {orth}, {tol})')
    end
    if nargin < 4 || isempty(n_wanted_eigs)
        n_wanted_eigs = 10;
    end
    if nargin < 5 || isempty(s)
        s = 6;
    end
    if nargin < 6 || isempty(basis)
        basis = 'newton';
    end
    if nargin < 7 || isempty(orth)
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
    if nargin < 8 || isempty(tol)
        tol = 1.0e-08;
    end
    tol = tol*norm_A;
    %tol = sqrt(tol);

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
    
    % Fix change-of-basis matrix
    Bk = [];
    if strcmpi(basis,'monomial')
        I = eye(s+1);
        Bk = I(:,2:s+1);        
    elseif strcmpi(basis,'newton')
        % Run standard Lanczos for 2s steps
        T = lanczos(A,q,2*s,'local');
        basis_eigs = eig(T);
        basis_shifts = leja(basis_eigs,'nonmodified');
        Bk = newton_basis_matrix(basis_shifts, s,1);
    end

    Q = zeros(n,max_lanczos+n_wanted_eigs);
    Q_conv = [];
    conv_eigs = [];
    conv_rnorms = [];
    rnorms = zeros(max_restarts,n_wanted_eigs);
    orth_err = [];  
    
    num_restarts = 0;
    restart = true;
    nconv = 0;
    while(restart && num_restarts < max_restarts)
               
        num_restarts = num_restarts+1;
        
        % Get the number of iterations to do next.
        %iters = floor((max_lanczos-nconv)/s);
        iters = floor(max_lanczos/s);
        
        if iters == 0
            % Check if we got all the eigenpairs we were after.
            restart = check_wanted_eigs(conv_eigs, diag(Dp(k+1:s*iters,k+1:s*iters)), n_wanted_eigs);
            break;
        end
        
        if strcmpi(orth,'local')
            [Q_new,T] = lanczos_basic(A, Q_conv, q, Bk, iters, s, basis, 'local');
        elseif strcmpi(orth,'full')
            [Q_new,T] = lanczos_basic(A, Q_conv, q, Bk, iters, s, basis, 'fro');
        elseif strcmpi(orth,'periodic')
            [Q_new,T] = lanczos_periodic(A, Q_conv, q, Bk, iters, s, basis);
        elseif strcmpi(orth,'selective')
            [Q_new,T] = lanczos_selective(A, Q_conv, q, Bk, iters, s, basis);
        else
            % Do nothing
        end
 
        % Compute residual norm estimates of all computed ritz-pairs.
        ritz_norms = zeros(s*iters,1);
        [Vp,Dp] = eig(T(1:s*iters,1:s*iters));
        beta = T(s*iters+1,s*iters);
        for i = 1:s*iters
            y = Vp(:,i);
            ritz_norms(i) = beta*abs(y(s*iters));% + eps*norm_A;
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
        
        % Save the converged eigenspace information
        for i = 1:k
            Q(:,nconv+i) = Q_new*Vp(:,i);
            conv_eigs = [conv_eigs; Dp(i,i)];
            conv_rnorms = [conv_rnorms; ritz_norms(i)];
        end
    
        if nargout >= 3    
            % Since the previously converged eigenpairs have not been
            % touched, just copy the corresponding residual norms.
            if num_restarts > 1
                rnorms(num_restarts,1:nconv) = rnorms(num_restarts-1,1:nconv);
            end
            % Save the relative residual norms of the newly converged
            % eigenpairs.
            for i = 1:k
                l = conv_eigs(nconv+i);
                x = Q(:,nconv+i);
                rnorms(num_restarts,nconv+i) = norm(A*x-l*x)/norm(l*x);
            end
            [d,ix] = sort(diag(Dp(k+1:end,k+1:end)),'descend');
            for i = 1:n_wanted_eigs-nconv-k
                l = d(i);
                x = Q_new*Vp(:,ix(i));
                rnorms(num_restarts,nconv+i+k) = norm(A*x-l*x)/norm(l*x);
            end
%            rn_sort = sort(ritz_norms);
%            rnorms(num_restarts,:) = [conv_rnorms',rn_sort(k+1:n_wanted_eigs-nconv)'];
        end
            
        % Update the maximum orthogonalization error
        if nargout >= 4
            Q_ = [Q_conv Q_new];
            orth_err = [orth_err; norm(eye(size(Q_,2))-Q_'*Q_,'fro')];
        end

        % Update the count of converged eigenvalues
        nconv = nconv+k;
        Q_conv = Q(:,1:nconv);
               
        % Check if we should continue iterations
        restart = check_wanted_eigs(conv_eigs, diag(Dp), n_wanted_eigs);

        if restart
            q = generateStartVector(diag(Dp),Vp,Q_new,ritz_norms,k,restart_strategy);   
            %q = projectAndNormalize({Q_conv},q,true);
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
        [sort_eigs,ixs] = sort(conv_eigs,'descend');
        sort_rnorms = conv_rnorms(ixs);
        conv_eigs = sort_eigs(1:nconv);
        conv_rnorms = sort_rnorms(1:nconv);
        Q_conv = Q_conv(:,ixs);
        Q_conv = Q_conv(:,1:nconv);
        disp(['Did not converge in ' num2str(num_restarts) ' restarts.']);
        disp(['Max residual norm: ' num2str(max(conv_rnorms))]);
    end
end

function q = generateStartVector(eig_vals,eig_vecs,Q,ritz_norms,k,strategy)
    if nargin < 6
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
        eig_vals(l)
        q = Q*eig_vecs(:,l);
    elseif strcmpi(strategy,'smallest')
        % Generate new starting vector from the smallest non-converged basis vector.
        l = k+1;
        for j = k+1:m
            if eig_vals(j) < eig_vals(l)
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
%        max_eig = max(eigs);
%        largest_conv_eigs = find(conv_eigs > max_eig);
%        if length(largest_conv_eigs) >= num_wanted_eigs
%            restart = false;
%        else
%            restart = true;
%        end
        if length(conv_eigs) >= num_wanted_eigs
            restart = false;
        else
            restart = true;
        end
    end
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
 
function [Q,T] = lanczos_basic(A, Q_conv, q, Bk, maxiter, s, basis, orth)

    if nargin < 7
        orth = 'local';
    end
  
    n = length(q);
    Q = zeros(n,maxiter*s);
    Q(:,1) = q;
    b = zeros(maxiter+1,1);
    T = [];
    
    k = 0;
    while k <= maxiter

        k = k+1;

        if k > 1
            q = Q(:,(k-1)*s+1);
        end
        
        V = matrix_powers(A,q,s,Bk,basis);
        
        if k == 1
            % Orthogonalize initial basis vectors
            [Q_,Rk] = normalize(V(:,1:s+1));
            % Orthogonalize against already converged basis vectors
            [Q(:,1:s+1),R_] = projectAndNormalize({Q_conv},Q_,true);           
            % Compute first part of tridiagonal matrix
            T = Rk*Bk/Rk(1:s,1:s);
            % Compute next beta
            b(k) = T(s+1,s);
            
        else
            if strcmpi(orth,'local')
                % Orthogonalize against previous block of basis vectors
                [Q_,Rk_] = projectAndNormalize({Q(:,(k-2)*s+1:(k-1)*s+1),Q_conv},V(:,2:s+1),true);
                Q(:,(k-1)*s+2:k*s+1) = Q_(:,1:s);
                Rkk_s = Rk_{1};
                Rk_s = Rk_{3};
            elseif strcmpi(orth,'fro')
                % Orthogonality against all previous basis vectors
                [Q_,Rk_] = projectAndNormalize({Q(:,(k-2)*s+1:(k-1)*s+1)},V(:,2:s+1),true);
                Rkk_s = Rk_{1};
                Rk_s = Rk_{2};
                Q(:,(k-1)*s+2:k*s+1) = projectAndNormalize({Q_conv,Q(:,1:(k-2)*s)},Q_,true);
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
    end
    
    % Fix output
    T = T(1:s*(k-1)+1,1:s*(k-1));
    Q = Q(:,1:s*(k-1));
end

function [Q,T] = lanczos_selective(A, Q_conv, q, Bk, maxiter, s, basis)
    
    n = length(q);
    Q = zeros(n,maxiter*s);
    QR = [];
    b = zeros(maxiter+1,1);
    T = [];
    norm_sqrt_eps = normest(A)*sqrt(eps);
    nritz = 0;
    
    Q(:,1) = q;
    
    k = 0;
    while k <= maxiter

        k = k+1;

        if k > 1
            q = Q(:,(k-1)*s+1);
        end
        
        V = matrix_powers(A,q,s,Bk,basis);
        
        if k == 1
            % Orthogonalize initial basis vectors
            [Q_,Rk] = normalize(V(:,1:s+1));
            % Orthogonalize against already converged basis vectors
            [Q(:,1:s+1),R_] = projectAndNormalize({Q_conv},Q_,true);
            % Compute first part of tridiagonal matrix
            T = Rk*Bk/Rk(1:s,1:s);
            % Compute next beta
            b(k) = T(s+1,s);
            
        else
            % Orthogonalize against previous block of basis vectors and the
            % already converged ritz vectors
            [Q_,Rk_] = projectAndNormalize({Q(:,(k-2)*s+1:(k-1)*s+1),Q_conv,QR(:,1:nritz)},V(:,2:s+1),true);
            Q(:,(k-1)*s+2:k*s+1) = Q_(:,1:s);
            Rkk_s = Rk_{1};
            Rk_s = Rk_{4};
            
            % Compute Tk (tridiagonal sub-matrix of T)
            Rkk = [zeros(s,1), Rkk_s(1:s,:)];
            Rk = [eye(s+1,1), [Rkk_s(s+1,1:s);Rk_s]];
            zk = Rk(1:s,s+1);
            rho = Rk(s+1,s+1);
            rho_t = Rk(s,s);
            bk = Bk(s+1,s);
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

        [Vp,Dp] = eig(T(1:s*k,1:s*k));

        nritz_new = 0;
        for i = 1:k*s
            if b(k)*abs(Vp(s*k,i)) < norm_sqrt_eps
                nritz_new = nritz_new+1;
            end
        end
        if nritz_new > nritz
            nritz = 0;
            for i = 1:k*s
                if b(k)*abs(Vp(s*k,i)) < norm_sqrt_eps
                    nritz = nritz+1;
                    y = Q(:,1:k*s)*Vp(:,i);
                    QR(:,nritz) = y;
                end
            end
            QR(:,1:nritz) = normalize(QR(:,1:nritz));
        end

        disp(['nritz=' num2str(nritz)])
    end
      
    % Fix output
    T = T(1:s*(k-1)+1,1:s*(k-1));
    Q = Q(:,1:s*(k-1));

end
   
function [Q,T] = lanczos_periodic(A, Q_conv, q, Bk, maxiter, s, basis)

    n = length(q);
    Q = zeros(n,maxiter*s);
    b = zeros(maxiter+1,1);
    T = [];
    omega = [];
    norm_A = normest(A);

    Q(:,1) = q;

    k = 0;
    while k <= maxiter

        k = k+1;

        if k > 1
            q = Q(:,(k-1)*s+1);
        end
        
        V = matrix_powers(A,q,s,Bk,basis);
        
        if k == 1
            % Orthogonalize initial basis vectors
            [Q_,Rk] = normalize(V(:,1:s+1));
            % Orthogonalize against already converged basis vectors
            [Q(:,1:s+1),R_] = projectAndNormalize({Q_conv},Q_,true);
            % Compute first part of tridiagonal matrix
            T = Rk*Bk/Rk(1:s,1:s);
            % Compute next beta
            b(k) = T(s+1,s);
            
        else
            % Orthogonalize against previous block of basis vectors and
            % previously converged vectors.
            [Q_,Rk_] = projectAndNormalize({Q(:,(k-2)*s+1:(k-1)*s+1),Q_conv},V(:,2:s+1),true);
            Q(:,(k-1)*s+2:k*s+1) = Q_(:,1:s);
            Rkk_s = Rk_{1};
            Rk_s = Rk_{3};
            
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
        
        % Estimate orthogonalization error, reorthogonalize if necessary
        alpha = diag(T,0);
        beta  = diag(T,-1);
        omega = update_omega(omega,alpha,beta,norm_A, s);
        err = 0;
        for i = 1:s
            omega_row = omega((k-1)*s+i+1,:);
            row_err = max(omega_row(1:i));
            if(row_err > err)
                err = row_err;
            end
        end
        if err >= sqrt(eps/(k*s)) %norm_A
            Q(:,(k-1)*s+1:k*s+1) = projectAndNormalize({Q(:,1:(k-1)*s)},Q(:,(k-1)*s+1:k*s+1),true);
            omega = reset_omega(omega, norm_A, s);
        end
    end
      
    % Fix output
    T = T(1:s*(k-1)+1,1:s*(k-1));
    Q = Q(:,1:s*(k-1));
end

function omega = update_omega(omega_in, alpha, beta, anorm, s)

    % Get iteration number and block size
    n = length(alpha);
    m = size(omega_in,1)-1;
    
    % Estimate of contribution to roundoff errors from A*v:  fl(A*v) = A*v + f, 
    T = eps*anorm;
    
    if isempty(omega_in) 
        omega = zeros(s+1,s+1);
        omega(1,1) = 1;
        omega(1,2) = 0;
        omega(2,1) = T/beta(1);
        omega(2,2) = 1;
        
        for j = 2:s
            binv = 1.0/beta(j);
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
            binv = 1.0/beta(j);
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

% -------------------------------------------------------------------------
%  Copyright (2012, 2013)  Magnus Grandin <magnus.grandin@it.uu.se>
% -------------------------------------------------------------------------
