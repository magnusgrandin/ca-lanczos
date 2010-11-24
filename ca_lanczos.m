function [T,Q,rnorm,orth] = ca_lanczos(A,r,s,t,basis,may_break,reorth)

    if nargin < 6
        may_break = 0;
    end
    if nargin < 7
        reorth = 0;
    end
    
    if (nargout >= 3) || (may_break == 1) || (reorth == 1)
        solve_eigs_every_step = 1;
    else
        solve_eigs_every_step = 0;
    end

    t = ceil(t);             % To make sure that t is an integer
    maxlanczos = s*t;
    
    n = length(r);
    b = zeros(s+1,1);
    V = zeros(n,maxlanczos+1);
    
    % Allocate Q as a cell array, one cell for each block, allowing
    % the last block to have fewer than s columns.
    Q = cell(1,ceil(maxlanczos/s)+1);
    for i = 1:length(Q)
        Q{i} = zeros(n,s);
    end
    Q{end} = zeros(n,1);
    
    rnorm = zeros(t,2);
    orth = zeros(t,1);
    
    b0 = norm(r);
    Q{1}(:,1) = (1/b0)*r;
    V(:,1) = Q{1}(:,1);
    
    has_converged = false;
    k = 1;
    
    % Fix basis vectors
    if strcmpi(basis,'monomial')
        I = eye(s+1);
        Bk = I(:,2:s+1);        
    elseif strcmpi(basis,'newton')
        % Run standard Lanczos for 2s steps
        T = lanczos(A,r,2*s);
        basis_eigs = eig(T);
        basis_shifts = leja(basis_eigs,'modified');
        Bk = newton_basis_matrix(basis_shifts, s,1);
    else
        disp(['ERROR: Unknown basis', basis]);
    end
    
  
    while (k <= t) && (has_converged == false)
    
        % Compute matrix powers
        if strcmpi(basis,'monomial')
            V(:,1) = Q{k}(:,1);
            V(:,2:s+1) = matrix_powers(A, Q{k}(:,1), s);
        elseif strcmpi(basis,'newton')
            V(:,1:s+1) = newton_basis(A, Q{k}(:,1), s, basis_shifts,1);
        end
        
        if k == 1
            % QR factorization
            [Q_,Rk] = tsqr(V(:,1:s+1));
            Q{1}(:,1:s) = Q_(:,1:s);
            Q{2}(:,1) = Q_(:,s+1);
            
            % Compute first part of tridiagonal matrix
            T = Rk*Bk/Rk(1:s,1:s);

            % Compute next beta
            b(k) = T(s+1,s);
            
        else

%             % BGS update
%             Q_ = [Q{k-1}(:,1:s),Q{k}(:,1)];
%             Rkk_s = Q_(:,1:s+1)' * V(:,2:s+1);
%             V(:,2:s+1) = V(:,2:s+1) - Q_(:,1:s+1)*Rkk_s;
% 
%             % QR factorization
%             Q_ = zeros(n,s+1);
%             [Q_(:,2:s+1),Rk_s] = tsqr(V(:,2:s+1));
%             Q{k}(:,2:s) = Q_(:,2:s);
%             Q{k+1}(:,1) = Q_(:,s+1);
            
            Q_ = {[Q{k-1}(:,1:s),Q{k}(:,1)], V(:,2:s+1)};
            [Q_,Rk_] = rr_tsqr_bgs(Q_);
            Q{k}(:,2:s) = Q_{2}(:,1:s-1);
            Q{k+1}(:,1) = Q_{2}(:,s);
            Rkk_s = Rk_(1:s+1,end-s+1:end);
            Rk_s = Rk_(end-s+1:end,end-s+1:end);
            
            
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
        
        if solve_eigs_every_step == 1
            [Vp,Dp] = eig(T(1:s*k,1:s*k));
        end

        rnormest = 0;
        if (nargout >= 3) || (reorth == 1) || (may_break == 1)
            % Residual norm for smallest eigenpair
            [d_s,i_s] = min(diag(Dp));
            s_s = Vp(s*k,i_s);
            rnormest(1) = b(k)*abs(s_s);
            %rnorm(k+1,1) = rnormest(1);
            x_s = zeros(n,1);
            for i = 1:k
                x_s = x_s+Q{i}*Vp(s*(i-1)+1:s*i,i_s);
            end
            rnorm(k+1,1) = norm(A*x_s-d_s*x_s)/norm(d_s*x_s);
            
            % Residual norm for largest eigenpair
            [d_l,i_l] = max(diag(Dp));
            s_l = Vp(s*k,i_l);
            rnormest(2) = b(k)*abs(s_l);
            %rnorm(k+1,2) = rnormest(2);
            x_l = zeros(n,1);
            for i = 1:k
                x_l = x_l+Q{i}*Vp(s*(i-1)+1:s*i,i_l);
            end
            rnorm(k+1,2) = norm(A*x_l-d_l*x_l)/norm(d_l*x_l);
        end
        
        % Check stopping criteria 
        if may_break == 1
            min_rnorm = min([rnormest(1), rnormest(2)]);
            if min_rnorm < norm(T)*sqrt(eps);
                has_converged = true;
            end
        end

        if reorth == 1
%            Rkk_ = Q(:,1:s*k+1)'*Q(:,s*k+2:s*(k+1)+1);
%            Q_(:,s*k+2:s*(k+1)+1) = ...
%                Q(:,s*k+2:s*(k+1)+1) - Q(:,1:s*k+1)*Rkk_;
%            [Q(:,s*k+2:s*(k+1)+1),Rkk_] = tsqr(Q_(:,s*k+2:s*(k+1)+1));
            Q_ = Q(1:k);
            Q_{k} = [Q{k} Q{k+1}(:,1)];
            [Q_,Rk_] = rr_tsqr_bgs(Q_);
            for i = 1:k-1
                Q{i} = Q_{i};
            end
            Q{k} = Q_{k}(:,1:s);
            Q{k+1}(:,1) = Q_{k}(:,s+1);
        end
     
        % Level of orthogonality
        if nargout >= 4
            Q_ = [];
            for i = 1:k
                Q_ = [Q_ Q{i}];
            end
            orth(k+1) = norm(eye(s*k)-Q_(:,1:s*k)'*Q_(:,1:s*k),'fro');
        end
        

        k = k+1;
    end
    
    T = T(1:s*(k-1),:);
    Q_ = [];
    for i = 1:k-1
        Q_ = [Q_ Q{i}];
    end
    Q = Q_;
    
    if nargout >= 3
        rnorm = rnorm(1:(k-1),:);
    end
    if nargout >= 4
        orth = orth(1:(k-1));
    end
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
