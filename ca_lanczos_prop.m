function [T,Q,lsteps] = ca_lanczos_prop(A,r0,s,m,dt,tol,basis,adaptive)
    
    if nargin < 6
        tol = 1.0e-10;
    end
    if nargin < 7
        basis = 'newton';
    end
    if nargin < 8
        adaptive = false;
    end
        
    n = length(r0);
    b = zeros(m+1,1);
    V = zeros(n,s+1);
    Q = {};

    nrm = norm(r0);
    q = (1/nrm)*r0;
    V(:,1) = q;
    
    % Fix basis vectors
    if strcmpi(basis,'monomial')
        I = eye(s+1);
        Bk = I(:,2:s+1);        
    elseif strcmpi(basis,'newton')
        % Run standard Lanczos for 2s steps
        T = lanczos(A,r0,2*s,tol);
        basis_eigs = eig(T);
        basis_shifts = leja(real(basis_eigs),'nonmodified');
        Bk = newton_basis_matrix(basis_shifts,s,1);
    else
        disp(['ERROR: Unknown basis', basis]);
    end
    
    hasConverged = false;
    k = 0;

    while (k < m) && (hasConverged == false)
        
        k = k+1;

        if k > 1
            q = Q{k-1}(:,s+1);
        end
        
        % Compute matrix powers
        if strcmpi(basis,'monomial')
            V(:,1) = q;
            V(:,2:s+1) = matrix_powers(A, q, s);
        elseif strcmpi(basis,'newton')
            V(:,1:s+1) = newton_basis(A, q, s, basis_shifts,1);
        end
        
        if k == 1
            % QR factorization
            [Q{1},Rk] = cholqr(V(:,1:s+1));
            
            % Compute first part of tridiagonal matrix
            T = Rk*Bk/Rk(1:s,1:s);
            
            % Compute next beta
            b(k) = T(s+1,s);
            
        else
            % Orthogonality against previous block of basis vectors
            Rkk_s = Q{k-1}'*V(:,2:s+1);
            V(:,2:s+1) = V(:,2:s+1) - Q{k-1}*Rkk_s;
            [Q_,Rk_s] = cholqr(V(:,2:s+1));
            Q{k} = [Q{k-1}(:,s+1) Q_(:,1:s)];
            Q{k-1} = Q{k-1}(:,1:s);
            
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
        
        % Compute the residual and stop iterations if tolerance,
        % fulfilled if adaptive flag is set.
        if (k*s >= 3) && (adaptive == true)
            [Vp,D] = eig(T(1:k*s,1:k*s));
            matexp = Vp*expm(-1i*dt*D)*Vp';
            residual = abs(dt*b(k)*matexp(k*s,1)*nrm);
            if residual < tol
                disp(num2str(residual));
                hasConverged = true;
                break;
            end
        end
    end
    lsteps = k*s;
    T = real(T(1:lsteps,1:lsteps));
    Q_ = Q{k};
    Q = [Q{1:k-1} Q_(:,1:s)];
end

%% Returns a columnvector with n elements, in which all elements are 
%% zero except the last element. 
function vec = eyeshvec(len)
    vec = eye(len,1);
    vec=circshift(vec,len-1);
end

%% Compute s matrix-vector multiplications of A and q
function V = matrix_powers(A,q,s)
    V = zeros(length(q),s);
    V(:,1) = A*q;
    for i = 2:s
        V(:,i) = A*V(:,i-1);
    end
end


