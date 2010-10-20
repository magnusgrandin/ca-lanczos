function [T,Q,lsteps] = ca_lanczos_prop(A,r0,s,m,dt,tol)
    
    n = length(r0);
    
    b = zeros(s*m+2,1);
    V = zeros(n,s*m+1);
    Q = zeros(n,s*m+1);
    
    b(1) = norm(r0);
    Q(:,1) = (1/b(1))*r0;
    V(:,1) = Q(:,1);

    for k = 0:m-1
        
        % Fix basis vectors (first attempt: Monomial basis)
        I = eye(s+1);
        Bk = I(:,2:s+1);        
               
        % Compute matrix powers
        V(:,s*k+2:s*(k+1)+1) = matrix_powers(A, Q(:,s*k+1), s);
        
        if k == 0
            % QR factorization
            [Q(:,1:s+1),Rk] = pqr(V(:,1:s+1));
            
            % Compute first part of tridiagonal matrix
            T = Rk*Bk/Rk(1:s,1:s);

            % Compute next beta
            b(s+1) = T(s+1,s);
            
        else
            % BGS update
            Rkk_s = Q(:,s*(k-1)+1:s*k+1)' * V(:,s*k+2:s*(k+1)+1);
            V(:,s*k+2:s*(k+1)+1) = ...
                V(:,s*k+2:s*(k+1)+1) - Q(:,s*(k-1)+1:s*k+1)*Rkk_s;
 
            % QR factorization
            [Q(:,s*k+2:s*(k+1)+1),Rk_s] = pqr(V(:,s*k+2:s*(k+1)+1));
            
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
                - b(s*k+1)*e1*es'*Rkk(1:s,1:s)/Rk(1:s,1:s);

            % Compute the next beta
            b(s*(k+1)+1) = bk*(rho/rho_t);
            
            % Extend T
            T11 = T(1:s*k,1:s*k);
            T12 = b(s*k+1)*eyeshvec(s*k)*eye(s,1)';
            T21 = b(s*k+1)*eye(s,1)*eyeshvec(s*k)';
            T22 = Tk;
            T31 = zeros(1,s*k);
            T32 = b(s*(k+1)+1)*eyeshvec(s)';
            T = [T11, T12; T21, T22; T31, T32];
        end
        
        % Compute the residual and stop iterations if tolerance fulfilled
        if(s*(k+1) > 2)
            [Vp,D] = eig(T(1:s*(k+1),1:s*(k+1)));
            matexp = Vp*expm(-1i*dt*D)*Vp';
            residual = abs(dt*T(s*(k+1)+1,s*(k+1))*matexp(s*(k+1),1)*b(1));
            if residual < tol
                break;
            end
        end
        
    end
    lsteps = (k+1)*s;
    T = T(1:lsteps,1:lsteps);
    Q = Q(:,1:lsteps);
end

%% Returns a columnvector with n elements, in which all elements are 
%% zero except the last element. 
function vec = eyeshvec(len)
    vec = eye(len,1);
    vec=circshift(vec,len-1);
end

%% A temporary fix.
%% This function does a QR factorization of matrix A, and shifts
%% the signs of the resulting factorization matrices such that R 
%% only has positive values on its diagonals.
function [Q,R] = pqr(A)
    [Q,R] = qr(A,0);
    d = sign(diag(R));
    R = diag(d)*R;
    Q = Q*diag(d);
end
 
%% Compute s matrix-vector multiplications of A and q
function V = matrix_powers(A,q,s)
    V = zeros(length(q),s);
    V(:,1) = A*q;
    for i = 2:s
        V(:,i) = A*V(:,i-1);
    end
end


