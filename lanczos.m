%--------------------------------------------------------------------------
% The symmetric Lanczos algorithm
% 
% Input variables:
% A       - the large sparse matrix
% r0      - starting vector
% maxiter - maximum order of the Krylov space
% opt     - structure of options
%           opt.break   : break on convergence of any eigenpair [{0} | 1]
%           opt.reorth  : reorthogonalize lanczos vectors [{0} | 1]
%           opt.restart : use restarting [{0} | 1]
%
% Output variable:
% T     - the tridiagonal lanczos matrix
% V     - the orthonormal matrix of lanczos vectors
% rnorm - the residual norm in each step
% orth  - the level of orthogonality in each step (||I-Vj'*Vj||)
%--------------------------------------------------------------------------
function [T,V,rnorm,orth] = lanczos(A,r0,maxiter,may_break,reorth)
    
    if nargin < 4
        may_break = 0;
    end
    if nargin < 5
        reorth = 0;
    end
    
    if (nargout >= 3) || (may_break == 1) || (reorth == 1)
        solve_eigs_every_step = 1;
    else
        solve_eigs_every_step = 0;
    end
    
    n = length(r0);
    V = zeros(n,maxiter);
    V(:,1) = r0/(sqrt(r0'*r0));
    alpha = zeros(1,maxiter);
    beta = zeros(1,maxiter);
    rnorm = zeros(maxiter,2);
    orth = zeros(maxiter,1);

    has_converged = false;
    j = 1;
    
    while (j <= maxiter) && (has_converged == false)
        r = A*V(:,j);
        if j > 1
            r=r-beta(j-1)*V(:,j-1);
        end
        alpha(j) = r'*V(:,j);
        r = r - alpha(j)*V(:,j);
        beta(j) = sqrt(r'*r);
        V(:,j+1) = r/beta(j);
        
        V = orthogonalize(V,j,reorth);
        
        if solve_eigs_every_step == 1
            T = diag(alpha(1:j)) + diag(beta(1:j-1),1) + diag(beta(1:j-1),-1);
            [Vp,Dp] = eig(T);
        end
        
        rnormest = 0;
        if (nargout >= 3) || (reorth == 1) || (may_break == 1)
            % Residual norm for smallest eigenpair
            [d_s,i_s] = min(diag(Dp));
            s_s = Vp(j,i_s);
            rnormest(1) = beta(j)*abs(s_s);
            %rnorm(j,1) = rnormest(1);
            x_s = V(:,1:j)*Vp(:,i_s);
            rnorm(j,1) = norm(A*x_s-d_s*x_s)/norm(d_s*x_s);
           
            % Residual norm for largest eigenpair
            [d_l,i_l] = max(diag(Dp));
            s_l = Vp(j,i_l);
            rnormest(2) = beta(j)*abs(s_l);
            %rnorm(j,2) = rnormest(2);
            x_l = V(:,1:j)*Vp(:,i_l);
            rnorm(j,2) = norm(A*x_l-d_l*x_l)/norm(d_l*x_l);
        end
        
        if nargout >= 4
            % Level of orthogonality
            orth(j) = norm(eye(j)-V(:,1:j)'*V(:,1:j),'fro');
        end
        
        % Check stopping criteria 
        if may_break == 1
            min_rnorm = min([rnormest(1), rnormest(2)]);
            if min_rnorm < norm(T)*sqrt(eps);
                has_converged = true;
            end
        end
        
        if reorth == 1
        end

        
        j = j+1;
    end
    
    if solve_eigs_every_step == 0
        T = diag(alpha(1:j-1)) + diag(beta(1:j-2),1) + diag(beta(1:j-2),-1);
    end

    if nargout >= 3
        rnorm = rnorm(1:j-1,:);
    end
    if nargout >= 4
        orth = orth(1:j-1);
    end
end

function V = orthogonalize(V,iter,orth)

    % Reorthogonalization of basis vectors
    % TODO: RR-TSQR-BGS?
    %[Q(:,1:m),R(1:m,1:m)] = qr(V(:,1:m),0);
    %for k = 2:j/m
    %    Rkk = Q(:,m*(k-2)+1:m*(k-1))'*V(:,m*(k-1)+1:m*k);
    %    Yk  = V(:,m*(k-1)+1:m*k) - Q(:,m*(k-2)+1:m*(k-1))*Rkk;
    %    [Q(:,m*(k-1)+1:m*k),R(m*(k-1)+1:m*k,m*(k-1)+1:m*k)] = qr(Yk,0);
    %end
    %Rkk = V(:,1:j)'*V(:,j+1);
    %V(:,j+1) = V(:,j+1) - V(:,1:j)*Rkk;
    
    % Copy V to temporary cell (block) array
    b = 4; %TODO: take this as option.
    V_ = cell(1,ceil((iter+1)/b));
    for i = 1:ceil((iter+1)/b)
        cols = b*(i-1)+1:min(b*i,iter+1);
        V_{i} = V(:,cols);
    end
    % Do reorthogonalization on blocks
    [V_,R_] = rr_tsqr_bgs(V_);
    % Copy reorthogonalized blocks back to V
    for i = 1:ceil((iter+1)/b)
        cols = b*(i-1)+1:min(b*i,iter+1);
        V(:,cols) = V_{i};
    end

end


    