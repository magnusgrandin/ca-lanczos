function [conv_eigs,Q_conv,num_restarts] = impl_restarted_ca_lanczos(A, r, max_lanczos, n_wanted_eigs, s, basis, orth, tol)

    max_restarts = 40;

    % Check input arguments
    if nargin < 3
        disp('ERROR: Wrong number of input arguments.')
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
    if nargin < 8 || isempty(tol)
        tol = 1.0e-06;
    end
    % Adjust the tolerance according to the norm of A.
    norm_A = normest(A);
    tol = tol*norm_A;

    if strcmpi(orth,'local')
        disp('Local orthogonalization');
    elseif strcmpi(orth,'full')
        disp('Full orthogonalization');
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

    Vk = zeros(n,n_wanted_eigs);
    conv_eigs = [];
    conv_rnorms = [];
    orth_err = [];

    if mod(n_wanted_eigs,s) ~= 0
        % TODO: Fix!!
        disp('Warning: Number of wanted eigs is not a multiple of s.');
        return;
    end
    
    k = n_wanted_eigs;
    p = s*floor((max_lanczos-k)/s);
    m = k+p;
    
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
    while(restart && (num_restarts < max_restarts))
               
        num_restarts = num_restarts+1;
              
        [Vm,Dm] = eig(Tm(nconv+1:m,nconv+1:m));
        eigsSorted = sort(diag(Dm),'descend'); % TODO: This only covers the case of the largest eigenvalues
     
        beta = Tm(m+1,m);
        x = [];
        new_conv = 0;
        for i = nconv+1:k
            y = Vm(:,i);
            y(m-nconv)
            if abs(beta*y(m-nconv)) < tol
               x = y;
               nconv = nconv + 1;
               break;
            end
        end
        Tm = deflate(x,Tm);
        if nconv >= k
            restart = false;
            break;
        end
        
        % Check stopping criteria
%         hasConverged = true;
%         beta = Tm(m+1,m);
%         for i = 1:m
%             d = Dm(i,i);
%             y = Vm(:,i);
%             if abs(beta*y(m))/abs(d) > tol
%                 hasConverged = false;
%                 break;
%             end
%         end
%         if hasConverged
%             restart = false;
%             disp('Converged.');
%             break;
%         end

%         beta = Tm(k+1,k);
%         if abs(beta*Vm(m,k)) < tol
%             restart = false;
%             disp('Converged.');
%             break;
%         end

        u = eigsSorted(k-nconv+1:k-nconv+p);

        %%%%%%%%%%%%%
        I = eye(m-nconv);
        for j = 1:p
            [Qj,Rj] = qr(Tm(nconv+1:m,nconv+1:m)-u(j)*I);
            Qm(:,nconv+1:m) = Qm(:,nconv+1:m)*Qj;
            Tm(nconv+1:m,nconv+1:m) = Qj'*Tm(nconv+1:m,nconv+1:m)*Qj;
%            Tm(nconv+1:m,nconv+1:m) = Rj*Qj'+ u*I;
        end
        %%%%%%%%%%%%%
        
        b = Tm(k+1,k);
        Vk(:,1:k) = Qm(:,1:k);
        rk = b*Qm(:,k+1) + Qj(m-nconv,k);
        Tk = Tm(1:k,1:k);%; zeros(1,k-1), b];
                      
        b_ = sqrt(rk'*rk);
        % Complete the m-step Lanczos factorization by applying p more steps
        if strcmpi(orth,'local')
            [Qp,Tp] = lanczos_basic(A, Vk, rk, Bk, floor(p/s), s, basis, 'local');
        elseif strcmpi(orth,'full')
            [Qp,Tp] = lanczos_basic(A, Vk, rk, Bk, floor(p/s), s, basis, 'fro');
        else
            % Do nothing
        end
        
        Tm = [Tk(1:k,1:k), b_*eyeshvec(k)*eye(p,1)'; b_*eye(p,1)*eyeshvec(k)', Tp(1:p,1:p); Tp(p+1,p)*eyeshvec(k+p)'];
        Qm = Qp;

        % For deflation:
        % 1. Initially, lock Ritz values as they converge until k values have been locked
        % 2. Continue to iterate, lock every new Ritz value that is better than any of 
        %    the already locked ones. Purge all unwanted but yet converged Ritz values.
        % 3. Continue with (2) until the next Ritz value to converge is not better. Replace 
        %    the (k+1)st basis vector with a random vector and orthogonalize it against the 
        %    previous ones (the k locked ones?). Go back to step (2)
        % 4. When step (3) has been executed two consecutive times with no replacement of 
        %    existing locked Ritz values, the iteration is stopped.
        %

        eig(Tm(1:k,1:k))
        
    end

    [Vk,Dk] = eig(Tm(1:k,1:k));
    conv_eigs = sort(diag(Dk),'descend');
    Q_conv = zeros(n,k);
    for i = 1:k
        Q_conv(:,i) = Qm(:,1:k)*Vk(:,i);
    end
    
    if ~restart
        disp(['Converged in ' num2str(num_restarts) ' restarts.']);
    else       
        disp(['Did not converge.']);
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
    T_new = T;
    if isempty(y)
        return;
    else
        k = size(y,1);
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
        T_new(1:k,1:k) = Q'*T(1:k,1:k)*Q;
    end
end

% function [V,H] = deflate(V, H, X, j)
%     
%     % Quick exit if there are no Ritz-pairs to deflate.
%     if isempty(X)
%         return;
%     end
%     
%     m = size(H,2);
%     i = size(X,2);
%     [Q,R] = qr(X);
%     H(j+1:m,j+1:m) = Q'*H(j+1:m,j+1:m)*Q;
%     H(1:j,j+1:m) = H(1:j,j+1:m)*Q;
%     V(:,j+1:m) = V(:,j+1:m)*Q;
%     [P,S] = qr(H(1+j+i:m,1+j+i:m)); %% Not sure about this (step 3 in alg 6.2)
%     %[H(1+j+i:m,1+j+i:m),P] = Hessred(H(1+j+i:m,1+j+i:m));
%     H(1+j+i:m,1+j+i:m) = P'*H(1+j+i:m,1+j+i:m)*P;
%     H(1:j+i,1+j+i:m) = H(1:j+i,1+j+i:m)*P;
%     V(:,1+j+i:m) = V(:,1+j+i:m)*P;   
% end


function [A, P] = Hessred(A)
% Reduction of the square matrix A to the upper
% Hessenberg form using Householder reflectors.
% The reflection vectors are stored in columns of
% the matrix V. Matrix A is overwritten with its
% upper Hessenberg form.
[m,n] =size(A);
if A == triu(A,-1)
   V = eye(m);
   return
end
V = [];
for k=1:m-2
   x = A(k+1:m,k);
   v = Housv(x);
   A(k+1:m,k:m) = A(k+1:m,k:m) - 2*v*(v'*A(k+1:m,k:m));
   A(1:m,k+1:m) = A(1:m,k+1:m) - 2*(A(1:m,k+1:m)*v)*v';
   v = [zeros(k,1);v];
   V = [V v];
end
P = eye(m)-2*V*V';
end

function u = Housv(x)
% Householder reflection unit vector u from the vector x.
m = max(abs(x));
u = x/m;
if  u(1) == 0
   su = 1;
else
   su = sign(u(1));
end
u(1) = u(1)+su*norm(u);
u = u/norm(u);
u = u(:);
end
% function getLejaPoints(a,b,theta_k,theta_p,z_in)
%     p = length(z_in);
%     z = [z_in; zeros(p,1)];
%     a = min(a,theta_k);
%     b = max(b,theta_k);
%     k = p;
%     while k < 2*p
%         if k == 0
%             z0 = b;
%         else
%             w = abs(z(k)-theta_k);
%             for i = 0 
%         end
%         k = k+1;
%     end
% end
