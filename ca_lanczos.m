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
function [T,Q,ritz_rnorm,orth_err] = ca_lanczos(A,r,s,iter,basis,orth)

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
    t = ceil(iter/s);
    
    % Normalize the starting vector
    q = r/sqrt(r'*r);
    
    if ~strcmpi(basis,'monomial') && ~strcmpi(basis,'newton')
        disp(['ERROR: Unknown basis type: ', basis]);
    end
    
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

    if strcmpi(orth,'local')
        [Q,T,ritz_rnorm,orth_err] = ca_lanczos_basic(A, q, Bk, t, s, basis, 'local');
    elseif strcmpi(orth,'full')
        [Q,T,ritz_rnorm,orth_err] = ca_lanczos_basic(A, q, Bk, t, s, basis, 'fro');
    elseif strcmpi(orth,'periodic')
        [Q,T,ritz_rnorm,orth_err] = ca_lanczos_periodic(A, q, Bk, t, s, basis);
    elseif strcmpi(orth,'selective')
        [Q,T,ritz_rnorm,orth_err] = ca_lanczos_selective(A, q, Bk, t, s, basis);
    else
        % Do nothing
    end
    
end

function [ritz_rnorm] = compute_ritz_rnorm(A,Q,Vp,Dp)
    m = size(Vp,1);
    ritz_rnorm = zeros(1,m);
    [d,ix] = sort(diag(Dp),'descend');
    for i = 1:m
        l = d(i);
        x = Q*Vp(:,ix(i));
        ritz_rnorm(i) = norm(A*x-l*x)/norm(l*x);
    end
end

function [orth_err] = compute_orth_err(Q,s)
%    orth_err = norm(eye(size(Q,2))-Q'*Q,'fro');
    j = size(Q,2);
    if j > s+1
        orth_err = max(max(abs(Q(:,1:j-s-1)'*Q(:,j-s:j))));
    else
        orth_err = max(max(abs(Q(:,1:j)'*Q(:,1:j) - eye(s+1))));
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
function [Q,T,rnorm,ortherr] = ca_lanczos_basic(A, q, Bk, iter, s, basis, orth)

    if nargin < 7
        orth = 'local';
    end
  
    rnorm = zeros(iter,iter*s);
    ortherr = zeros(iter,1);   

    n = length(q);
    Q = zeros(n,iter*s);
    Q(:,1) = q;
    b = zeros(iter+1,1);
    T = [];
    
    k = 0;
    while k < iter

        k = k+1;

        if k > 1
            q = Q(:,(k-1)*s+1);
        end
        
        V = matrix_powers(A,q,s,Bk,basis);
        
        if k == 1
            % Orthogonalize initial basis vectors
            [Q(:,1:s+1),Rk] = normalize(V(:,1:s+1));
            % Compute first part of tridiagonal matrix
            T = Rk*Bk/Rk(1:s,1:s);
            % Compute next beta
            b(k) = T(s+1,s);
            
        else
            if strcmpi(orth,'local')
                % Orthogonalize against previous block of basis vectors
                [Q_,Rk_] = projectAndNormalize({Q(:,(k-2)*s+1:(k-1)*s+1)},V(:,2:s+1),true);
                Q(:,(k-1)*s+2:k*s+1) = Q_(:,1:s);
                Rkk_s = Rk_{1};
                Rk_s = Rk_{2};
            elseif strcmpi(orth,'fro')
                % Orthogonality against all previous basis vectors                
                [Q_,Rk_] = projectAndNormalize({Q(:,(k-2)*s+1:(k-1)*s+1)},V(:,2:s+1),true);
                Rkk_s = Rk_{1};
                Rk_s = Rk_{2};
                Q(:,(k-1)*s+2:k*s+1) = Q_;
                Q(:,(k-1)*s+2:k*s+1) = projectAndNormalize({Q(:,1:(k-1)*s+1)},Q(:,(k-1)*s+2:k*s+1));
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
        if nargout >= 3
            [Vp,Dp] = eig(T(1:s*k,1:s*k));
            rnorm(k,1:s*k) = compute_ritz_rnorm(A,Q(:,1:s*k),Vp,Dp);
        end
        
        % Compute the orthogonalization error
        if nargout >= 4
            ortherr(k) = compute_orth_err(Q(:,1:s*k+1),s);
        end

    end
    
    % Fix output
    T = T(1:s*k,1:s*k);
    Q = Q(:,1:s*k);
    rnorm = rnorm(1:k,:);
    ortherr = ortherr(1:k);
end


function [Q,T,rnorm,ortherr] = ca_lanczos_selective(A, q, Bk, iter, s, basis)

    rnorm = zeros(iter,iter*s);
    ortherr = zeros(iter,1);

    n = length(q);
    Q = zeros(n,iter*s);
    QR = [];
    b = zeros(iter+1,1);
    T = [];
    norm_sqrt_eps = normest(A)*sqrt(eps);
    nritz = 0;
    nbreaks = 0;
    
    Q(:,1) = q;
    
    k = 0;
    while k < iter

        k = k+1;

        if k > 1
            q = Q(:,(k-1)*s+1);
        end
        
        V = matrix_powers(A,q,s,Bk,basis);
        
        if k == 1
            % Orthogonalize initial basis vectors
            [Q(:,1:s+1),Rk] = normalize(V(:,1:s+1));
            % Compute first part of tridiagonal matrix
            T = Rk*Bk/Rk(1:s,1:s);
            % Compute next beta
            b(k) = T(s+1,s);
            
        else
            % Orthogonalize against previous block of basis vectors and the
            % already converged ritz vectors
            [Q_,Rk_] = projectAndNormalize({Q(:,(k-2)*s+1:(k-1)*s+1),QR(:,1:nritz)},V(:,2:s+1),true);
            Q(:,(k-1)*s+2:k*s+1) = Q_(:,1:s);
            Rkk_s = Rk_{1};
            Rk_s = Rk_{3};
            
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
            nbreaks = nbreaks + 1;
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

        % Compute the ritz-norm, if it is required
        if nargout >= 3
            [Vp,Dp] = eig(T(1:s*k,1:s*k));
            rnorm(k,1:s*k) = compute_ritz_rnorm(A,Q(:,1:s*k),Vp,Dp);
        end
        
        % Compute the orthogonalization error
        if nargout >= 4
            ortherr(k) = compute_orth_err(Q(:,1:s*k+1),s);
        end
    end
      
    % Fix output
    T = T(1:s*k,1:s*k);
    Q = Q(:,1:s*k);
    rnorm = rnorm(1:k,:);
    ortherr = ortherr(1:k);
    
    disp(['--- Number of orthogonalization breaks: ' num2str(nbreaks)]);
end
         

function [Q,T,rnorm,ortherr] = ca_lanczos_periodic(A, q, Bk, iter, s, basis)

    rnorm = zeros(iter,iter*s);
    ortherr = zeros(iter,1);

    n = length(q);
    Q = zeros(n,iter*s);
    b = zeros(iter+1,1);
    T = [];
    omega = [];
    norm_A = normest(A);
    nbreaks = 0;

    Q(:,1) = q;

    k = 0;
    while k < iter

        k = k+1;

        if k > 1
            q = Q(:,(k-1)*s+1);
        end
        
        V = matrix_powers(A,q,s,Bk,basis);
        
        if k == 1
            % Orthogonalize initial basis vectors
            [Q(:,1:s+1),Rk] = normalize(V(:,1:s+1));
            % Compute first part of tridiagonal matrix
            T = Rk*Bk/Rk(1:s,1:s);
            % Compute next beta
            b(k) = T(s+1,s);
            
        else
            % Orthogonalize against previous block of basis vectors and
            % previously converged vectors.
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
        
        % Estimate orthogonalization error, reorthogonalize if necessary
        alpha = diag(T,0);
        beta  = diag(T,-1);
        omega = update_omega(omega,alpha,beta,norm_A, s);
        err = 0;
        for i = 1:s
            omega_row = omega((k-1)*s+i+1,:);
            row_err = max(abs(omega_row(1:(k-1)*s+i)));
            if(row_err > err)
                err = row_err;
            end
        end
        if err >= sqrt(eps) %norm_A*
            nbreaks = nbreaks + 1;
            Q(:,(k-1)*s+1:k*s+1) = projectAndNormalize({Q(:,1:(k-1)*s)},Q(:,(k-1)*s+1:k*s+1),true);
            omega = reset_omega(omega, norm_A, s);
        end
        
        % Compute the ritz-norm, if it is required
        if nargout >= 3
            [Vp,Dp] = eig(T(1:s*k,1:s*k));
            rnorm(k,1:s*k) = compute_ritz_rnorm(A,Q(:,1:s*k),Vp,Dp);
        end
        
        % Compute the orthogonalization error
        if nargout >= 4
            ortherr(k) = compute_orth_err(Q(:,1:s*k+1),s);
        end
    end
      
    % Fix output
    T = T(1:s*k,1:s*k);
    Q = Q(:,1:s*k);
    rnorm = rnorm(1:k,:);
    ortherr = ortherr(1:k);
    
    disp(['--- Number of orthogonalization breaks: ' num2str(nbreaks)]);
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
