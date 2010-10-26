function [T,Q,rnorm,orth] = ca_lanczos(A,r,s,t,opt)
    
    % Set default values for options
    may_break  = 0;
    do_reorth  = 0;
    do_restart = 0;
    
    % Check options and set values accordingly
    if nargin == 5
        if isfield(opt,'break') && (opt.break == 0 || opt.break == 1)
            may_break = opt.break;
        end
        if isfield(opt,'reorth') && (opt.reorth == 0 || opt.reorth == 1)
            % currently not implemented for CA-Lanczos
            do_reorth = opt.reorth;
        end
        if isfield(opt,'restart') && (opt.restart == 0 || opt.restart == 1) 
            % currently not implemented
            do_restart = opt.restart;
        end
    end

    if (do_reorth == 1) || (do_restart == 1) || (may_break == 1) || (nargout > 2)
        solve_eigs_every_step = 1;
    else
        solve_eigs_every_step = 0;
    end

    t = ceil(t);             % To make sure that t is an integer
    maxlanczos = s*t;

    n = length(r);
    b = zeros(s+1,1);
    V = zeros(n,maxlanczos+1);
    Q = zeros(n,maxlanczos+1);
    rnorm = zeros(t,2);
    orth = zeros(t,1);

    b0 = norm(r);
    Q(:,1) = (1/b0)*r;
    V(:,1) = Q(:,1);
    
    for k = 0:t-1
    
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
            b(k+1) = T(s+1,s);
            
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
                - b(k)*e1*es'*Rkk(1:s,1:s)/Rk(1:s,1:s);

            % Compute the next beta
            b(k+1) = bk*(rho/rho_t);
            
            % Extend T
            T11 = T(1:s*k,1:s*k);
            T12 = b(k)*eyeshvec(s*k)*eye(s,1)';
            T21 = b(k)*eye(s,1)*eyeshvec(s*k)';
            T22 = Tk;
            T31 = zeros(1,s*k);
            T32 = b(k+1)*eyeshvec(s)';
            T = [T11, T12; T21, T22; T31, T32];
 
        end   
        
        if do_reorth == 1
            for i = 1:s
                Rkk_ = Q(:,1:s*k+(i-1))'*Q(:,s*k+i);
                Q(:,s*k+i) = Q(:,s*k+i) - Q(:,1:s*k+(i-1))*Rkk_;
            end
        end
        
        if solve_eigs_every_step == 1
            % Find all eigenpairs of T
            [Vp,Dp] = eig(T(1:s*(k+1),1:s*(k+1)));
        end
        
        if (do_reorth == 1) || (do_restart == 1) || (may_break == 1) || (nargout >= 3)
            % Residual norm for smallest eigenpair
            [d_s,i_s] = min(diag(Dp));
            s_s = Vp(s*(k+1),i_s);
            rnorm(k+1,1) = b(k+1)*abs(s_s);
            
            % Residual norm for largest eigenpair
            [d_l,i_l] = max(diag(Dp));
            s_l = Vp(s*(k+1),i_l);
            rnorm(k+1,2) = b(k+1)*abs(s_l);
        end
        
        if nargout >= 4
            % Level of orthogonality
            orth(k+1) = norm(eye(k*s+1)-Q(:,1:k*s+1)'*Q(:,1:k*s+1),'fro');
        end
        
        if (may_break == 1) && (min([rnorm(k+1,1) rnorm(k+1,2)]) < norm(T)*sqrt(eps))
            break;
        end

    end
    
    T = T(1:s*(k+1),:);
    Q = Q(:,1:s*(k+1));
    
    if nargout >= 3
        rnorm = rnorm(1:k+1,:);
    end
    if nargout >= 4
        orth = orth(1:k+1);
    end
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

% function [Q,R] = BGSreorth(V,mk)
%     [n,m] = size(V);
%     M = mod(V,mk)+1;
%     
%     for k = 1:M
%             
%         Vk = V(:,mk*(k-1)+1:mk*k);
%         Qk = Q(:,1:mk*(k-1));
%         % First pass of block CGS
%         if k == 1
%             Yk = Vk;
%         else
%             Rk = Qk*Vk;
%             Yk = Vk - Qk*Vk;
%         end
%         
%         % Second pass of block CGS
%         if k == 1
%             Zk = Yk;
%         else
%             Rk = Q(:,1:mk*(k-1))*Yk;
%             Yk = V(:,mk*(k-1)+1:mk*k) - Q(:,1:mk*(k-1))*V(:,mk*(k-1)+1:mk*k);
%         end
%         
%         % Compute column norms of Yk and Zk mark the orthogonalization as
%         % good if the norm of each column of Zk is no smaller than half the
%         % norm of the corresponding column of Yk
%         good_cols = zeros(1,mk);
%         for i = 1:mk
%             colnorm_Yk = norm(Yk(:,i));
%             colnorm_Zk = norm(Zk(:,1));
%             if colnorm_Zk >= colnorm_Yk/2
%                 good_cols(i) = 1;
%             end
%         end
%         
%         % Were the column norms good?
%         colnorms_good = min(good_cols);
% 
%         if colnorms_good == false
%             
%         end
%         
%         
%         Rkk = Q(:,1)'*V(:,m*k+1:m*(k+1));
%         Yk  = V(:,m*k+1:m*(k+1)) - Q(:,1:m*(k-2)+1:m*(k-1))*Rkk;
%         [Q(:,m*(k-1)+1:m*k),R(m*(k-1)+1:m*k,m*(k-1)+1:m*k)] = qr(Yk,0);
%     end
%  
% end
% function [Q,R] = blockCGS(V,mk)
%      [n,m] = size(V);
%      M = mod(V,mk);
%      
%     for j = 1:M
%          Rkk = Q(:,1:)
%          
%          Rkk = Q(:,1)'*V(:,m*k+1:m*(k+1));
%          Yk  = V(:,m*k+1:m*(k+1)) - Q(:,1:m*(k-2)+1:m*(k-1))*Rkk;
%          [Q(:,m*(k-1)+1:m*k),R(m*(k-1)+1:m*k,m*(k-1)+1:m*k)] = qr(Yk,0);
%      end
%  end
