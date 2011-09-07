function [conv_eigs,Q_conv,num_restarts,conv_rnorms,orth_err] = impl_restarted_ca_lanczos(A, r, max_lanczos, n_wanted_eigs, s, basis, orth, tol)

    max_restarts = 10;

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
    
    % Fix change-of-basis matrix
    Bk = [];
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

    Q = zeros(n,max_lanczos+s);
    Q_conv = [];
    conv_eigs = [];
    conv_rnorms = [];
    orth_err = [];
       
    m = max_lanczos;
    k = n_wanted_eigs;
    p = m-k;

    % Compute initial m-step Lanczos factorization
    if strcmpi(orth,'local')
        [Qm,Tm] = lanczos_basic(A, [], q, Bk, floor(m/s), s, basis, 'local');
    elseif strcmpi(orth,'full')
        [Qm,Tm] = lanczos_basic(A, [], q, Bk, floor(m/s), s, basis, 'fro');
    else
        % Do nothing
    end  
    
    num_restarts = 0;
    restart = true;
    nconv = 0;
    while(restart && num_restarts < max_restarts)
               
        num_restarts = num_restarts+1;
       
        [Vm,Dm] = eig(Tm(1:m,1:m));
        mu = sort(diag(Dm),'descend'); % TODO: This only covers the case of the largest eigenvalues
        mu = mu(1+k:p+k);
        
        Q = eye(m);
        for j = 1:p
            [Qj,Rj] = qr(Tm(1:m,1:m)-mu(j)*eye(m));
            Tm = Qj'*Tm(1:m,1:m)*Qj;
            Q = Q*Qj;
        end

        beta = Tm(k+1,k);
        sigma = Q(m,k);
        rm = Qm(:,m+1);
        rk = Qm(:,k+1)*beta + rm*sigma;
        Vk = Qm(:,1:m)*Q(:,1:k);
        Tk = Tm(1:k,1:k);
                      
        % Complete the m-step Lanczos factorization by applying p more steps
        if strcmpi(orth,'local')
            [Qp,Tp] = lanczos_basic(A, Vk, rk, Bk, floor(p/s), s, basis, 'local');
        elseif strcmpi(orth,'full')
            [Qp,Tp] = lanczos_basic(A, Vk, rk, Bk, floor(p/s), s, basis, 'fro');
        else
            % Do nothing
        end  
        Tm = [Tk, beta*eyeshvec(k)*eye(p,1)'; beta*eye(p,1)*eyeshvec(k)', Tp(1:p,1:p)];
        Qm = Qp;
        
        % Compute residual norm estimates of all computed ritz-pairs.
        %beta = Tm(k+1,k);
        new_conv = true;
        while new_conv
            [Vp,Dp] = eig(Tm(nconv+1:m,nconv+1:m));
            iter_count = 0;
            for i = 1:m-nconv
                y = Vp(:,i);
                ritz_norm = beta*abs(y(m-nconv));
                if ritz_norm <= tol
                    nconv = nconv+1;
                    conv_rnorm(nconv) = ritz_norm;
                    Tm(nconv:m,nconv:m) = deflate(y,Tm(nconv:m,nconv:m));
                    break;
                end
                iter_count = iter_count + 1;
            end
            if iter_count == m-nconv
                new_conv = false;
            end
        end
        
        disp('');

        % 1. Initially, lock Ritz values as they converge until k values have been locked
        % 2. Continue to iterate, lock every new Ritz value that is better than any of 
        %    the already locked ones. Purge all unwanted but yet converged Ritz values.
        % 3. Continue with (2) until the next Ritz value to converge is not better. Replace 
        %    the (k+1)st basis vector with a random vector and orthogonalize it against the 
        %    previous ones (the k locked ones?). Go back to step (2)
        % 4. When step (3) has been executed two consecutive times with no replacement of 
        %    existing locked Ritz values, the iteration is stopped.
        %


%         % Get the number of iterations to do next.
%         iters = floor((max_lanczos-nconv)/s);
%         if iters == 0
%             % Check if we got all the eigenpairs we were after.
%             restart = check_wanted_eigs(conv_eigs, diag(Dp(k+1:s*iters,k+1:s*iters)), n_wanted_eigs);
%             break;
%         end
%         
%         if strcmpi(orth,'local')
%             [Q_new,T] = lanczos_basic(A, Q_conv, q, Bk, iters, s, basis, 'local');
%         elseif strcmpi(orth,'full')
%             [Q_new,T] = lanczos_basic(A, Q_conv, q, Bk, iters, s, basis, 'fro');
%         elseif strcmpi(orth,'periodic')
%             [Q_new,T] = lanczos_periodic(A, Q_conv, q, Bk, iters, s, basis);
%         elseif strcmpi(orth,'selective')
%             [Q_new,T] = lanczos_selective(A, Q_conv, q, Bk, iters, s, basis);
%         else
%             % Do nothing
%         end
%  
%         % Update the maximum orthogonalization error
%         if nargout >= 4
%             Q_ = [Q_conv Q_new];
%             orth_err = [orth_err; norm(eye(size(Q_,2))-Q_'*Q_,'fro')];
%         end
% 
%         % Compute residual norm estimates of all computed ritz-pairs.
%         ritz_norms = zeros(s*iters,1);
%         [Vp,Dp] = eig(T(1:s*iters,1:s*iters));
%         beta = T(s*iters+1,s*iters);
%         for i = 1:s*iters
%             y = Vp(:,i);
%             ritz_norms(i) = beta*abs(y(s*iters)) + eps*norm_A;
%         end
%                
%         % Rearrange the eigenvalues such that the converged ones are first.
%         k = 0;
%         for i = 1:s*iters
%             if ritz_norms(i) < tol
%                 k = k + 1;
%                 % Push converged eigenvalues to the left.
%                 Dp_temp = Dp(i,i); Dp(i,i) = Dp(k,k); Dp(k,k) = Dp_temp;
%                 % Push corresponding eigenvectors to the left.
%                 Vp_temp = Vp(:,i); Vp(:,i) = Vp(:,k); Vp(:,k) = Vp_temp;
%                 % Push converged residual ritz norms to the left.
%                 rn_temp = ritz_norms(i); ritz_norms(i) = ritz_norms(k); ritz_norms(k) = rn_temp;
%                 ritz_norms(k)
%                 Dp(k,k)
%             end
%         end
%         
%         
%         for i = 1:k
%             Q(:,i+nconv) = Q_new*Vp(:,i);
%             conv_eigs = [conv_eigs; Dp(i,i)];
%             conv_rnorms = [conv_rnorms; ritz_norms(i)];
%         end
%         % Update the count of converged eigenvalues
%         nconv = nconv+k;
%         Q_conv = Q(:,1:nconv);
% 
%         % Check if we should continue iterations
%         restart = check_wanted_eigs(conv_eigs, diag(Dp), n_wanted_eigs);
% 
%         if restart
%             q = generateStartVector(diag(Dp),Vp,Q_new,ritz_norms,k,restart_strategy);   
%             q = projectAndNormalize({Q_conv},q,true);
%         end
    end

%    if ~restart
        [sort_eigs,ixs] = sort(conv_eigs,'descend');
        sort_rnorms = conv_rnorms(ixs);
        conv_eigs = sort_eigs(1:n_wanted_eigs);
        conv_rnorms = sort_rnorms(1:n_wanted_eigs);
        Q_conv = Q_conv(:,ixs);
        Q_conv = Q_conv(:,1:n_wanted_eigs);
        disp(['Converged in ' num2str(num_restarts) ' restarts.']);
        disp(['Max residual norm: ' num2str(max(conv_rnorms))]);
%     else       
%         disp(['Did not converge.']);
%         conv_eigs = [];
%         conv_rnorms = [];
%         Q_conv = [];
%         orth_err = [];
%     end
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
 
function [Q,T] = lanczos_basic(A, Vk, q, Bk, maxiter, s, basis, orth)

    if nargin < 7
        orth = 'local';
    end
  
    n = length(q);
    Q = zeros(n,maxiter*s);
    k = 0;
    if isempty(Vk)
        Q(:,1) = q;
    else
        k = size(Vk,2);
        Q(:,1:k) = Vk;
        Q(:,k+1) = q;
    end
    b = zeros(maxiter+1,1);
    T = [];
    
    iter = 0;
    while iter <= maxiter

        iter = iter+1;

        if iter > 1
            q = Q(:,(iter-1)*s+1+k);
        end
        
        V = matrix_powers(A,q,s,Bk,basis);
        
        if iter == 1
            % Orthogonalize initial basis vectors
            [Q(:,1+k:s+1+k),Rk] = normalize(V(:,1:s+1));
            % Compute first part of tridiagonal matrix
            T = Rk*Bk/Rk(1:s,1:s);
            % Compute next beta
            b(iter) = T(s+1,s);
            
        else
            if strcmpi(orth,'local')
                % Orthogonalize against previous block of basis vectors
                [Q_,Rk_] = projectAndNormalize({Q(:,(iter-2)*s+1:(iter-1)*s+1)},V(:,2:s+1),true);
                Q(:,(iter-1)*s+2+k:iter*s+1+k) = Q_(:,1:s);
                Rkk_s = Rk_{1};
                Rk_s = Rk_{3};
            elseif strcmpi(orth,'fro')
                % Orthogonality against all previous basis vectors
                [Q_,Rk_] = projectAndNormalize({Q(:,1:(iter-2)*s+k),Q(:,(iter-2)*s+1+k:(iter-1)*s+1+k)},V(:,2:s+1),true);
                Rkk_s = Rk_{2};
                Rk_s = Rk_{3};
                Q(:,(iter-1)*s+2+k:iter*s+1+k) = Q_;
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
                - b(iter-1)*e1*es'*Rkk(1:s,1:s)/Rk(1:s,1:s);

            % Compute the next beta
            b(iter) = bk*(rho/rho_t);
            
            % Extend T
            T11 = T(1:s*(iter-1),1:s*(iter-1));
            T12 = b(iter-1)*eyeshvec(s*(iter-1))*eye(s,1)';
            T21 = b(iter-1)*eye(s,1)*eyeshvec(s*(iter-1))';
            T22 = Tk;
            T31 = zeros(1,s*(iter-1));
            T32 = b(iter)*eyeshvec(s)';
            T = [T11, T12; T21, T22; T31, T32];
 
        end
    end
    
    % Fix output
    T = T(1:s*(iter-1)+1,1:s*(iter-1));
    Q = Q(:,1:s*(iter-1)+k+1);
end

function T_new = deflate(y,T)
    k = length(y);
    Q = zeros(k,k);
    Q(:,1) = y;
    s = y(1)^2;
    t0 = abs(y(1));
    for j = 2:k
        s = s + y(j)^2;
        t = sqrt(s);
        if t0 ~= 0
            g = (y(j)/t)/t0;
            Q(1:j-1,j) = -y(1:j-1)*g;
            Q(j,j) = t0/t;
        else
            Q(j-1,j) = 1;
        end
        t0 = t;
    end
    T_new = Q'*T*Q;
end

