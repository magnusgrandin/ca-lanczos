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
function [T,V,rnorm,orth] = lanczos(H,r0,maxiter,opt)

    % Reorthogonalization interval;
    m = 5;
    
    % Set default values for options
    may_break  = 0;
    do_reorth  = 0;
    do_restart = 0;
    
    % Check options and set values accordingly
    if nargin == 4
        if isfield(opt,'break') && (opt.break == 0 || opt.break == 1)
            may_break = opt.break;
        end
        if isfield(opt,'reorth') && (opt.reorth == 0 || opt.reorth == 1)
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
    
    n = length(r0);
    V = zeros(n,maxiter);
    V(:,1) = r0/(sqrt(r0'*r0));
    alpha = zeros(1,maxiter);
    beta = zeros(1,maxiter);
    rnorm = zeros(maxiter,2);
    orth = zeros(maxiter,1);

    for j = 1:maxiter
        r = H*V(:,j);
        if j > 1
            r=r-beta(j-1)*V(:,j-1);
        end
        alpha(j) = r'*V(:,j);
        r = r - alpha(j)*V(:,j);
        beta(j) = sqrt(r'*r);
        V(:,j+1) = r/beta(j);  
        
        if do_reorth == 1
            % Reorthogonalization of basis vectors
            % TODO: RR-TSQR-BGS?
            %[Q(:,1:m),R(1:m,1:m)] = qr(V(:,1:m),0);
            %for k = 2:j/m
            %    Rkk = Q(:,m*(k-2)+1:m*(k-1))'*V(:,m*(k-1)+1:m*k);
            %    Yk  = V(:,m*(k-1)+1:m*k) - Q(:,m*(k-2)+1:m*(k-1))*Rkk;
            %    [Q(:,m*(k-1)+1:m*k),R(m*(k-1)+1:m*k,m*(k-1)+1:m*k)] = qr(Yk,0);
            %end
            Rkk = V(:,1:j)'*V(:,j+1);
            V(:,j+1) = V(:,j+1) - V(:,1:j)*Rkk;
        end
        if solve_eigs_every_step == 1
            T = diag(alpha(1:j)) + diag(beta(1:j-1),1) + diag(beta(1:j-1),-1);
            [Vp,Dp] = eig(T);
        end
        
        
        if (do_reorth == 1) || (do_restart == 1) || (may_break == 1) || (nargout >= 3)
            % Residual norm for smallest eigenpair
            [d_s,i_s] = min(diag(Dp));
            s_s = Vp(j,i_s);
            rnorm(j,1) = beta(j)*abs(s_s);
            
            % Residual norm for largest eigenpair
            [d_l,i_l] = max(diag(Dp));
            s_l = Vp(j,i_l);
            rnorm(j,2) = beta(j)*abs(s_l);
        end
        
        if nargout >= 4
            % Level of orthogonality
            orth(j) = norm(eye(j)-V(:,1:j)'*V(:,1:j),'fro');
        end
        
        if (may_break == 1) && (min([rnorm(j,1) rnorm(j,2)]) < norm(T)*sqrt(eps))
            break;
        end
    end
    
    if solve_eigs_every_step == 0
        T = diag(alpha(1:j)) + diag(beta(1:j-1),1) + diag(beta(1:j-1),-1);
    end

    if nargout >= 3
        rnorm = rnorm(1:j,:);
    end
    if nargout >= 4
        orth = orth(1:j);
    end
end

    