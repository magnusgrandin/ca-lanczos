%%
%   function [T,V,rnorm,orthl] = lanczos(A,r0,maxiter,stop,orth)
%
%   The symmetric Lanczos algorithm
% 
%   Input:
%     A       - the large sparse symmetric matrix
%     r0      - initial vector
%     maxiter - mximum order of the Krylov space
%     stop    - break on convergence of any eigenpair [{0} | 1]
%     orth    - reorthogonalize lanczos vectors [{0} | 1]
%
%   Output variable:
%     T       - the tridiagonal lanczos matrix
%     V       - the orthonormal matrix of lanczos vectors
%     rnorm   - the residual norm in each step
%     orthl   - the level of orthogonality in each step (||I-Vj'*Vj||)
%--------------------------------------------------------------------------
function [T,V,rnorm,orthl] = lanczos(A,r0,maxiter,stop,orth)
    
    if nargin < 4
        stop = 0;
    end
    if nargin < 5
        orth = 'local';
    else
        if isnumeric(orth)
            orth = num2str(orth);
        end
        if strcmpi(orth,'local')==0 && strcmpi(orth,'full')==0 && strcmpi(orth,'selective')==0 ...
                && strcmpi(orth,'periodic')==0 && strcmpi(orth,'partial')==0 
            disp(['lanczos.m: Invalid option value for orth: ', orth]);
            disp('    expected {''local''|''full''|''periodic''|''partial''|''selective''}');
            return;
        end
    end
    
    if (nargout >= 3) || (stop == 1) || strcmpi(orth,'local')==0
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
    orthl = zeros(maxiter,1);

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
        
        V = orthogonalize(V,j,orth);
        
        if solve_eigs_every_step == 1
            T = diag(alpha(1:j)) + diag(beta(1:j-1),1) + diag(beta(1:j-1),-1);
            [Vp,Dp] = eig(T);
        end
        
        rnormest = 0;
        if (nargout >= 3) || (stop == 1) || strcmpi(orth,'local')==0 
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
            orthl(j) = norm(eye(j)-V(:,1:j)'*V(:,1:j),'fro');
        end
        
        % Check stopping criteria 
        if stop == 1
            min_rnorm = min([rnormest(1), rnormest(2)]);
            if min_rnorm < norm(T)*sqrt(eps);
                has_converged = true;
            end
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
        orthl = orthl(1:j-1);
    end
end

function V = orthogonalize(V,iter,orth)

    if strcmpi(orth,'full') == 1
        % Reorthogonalization of basis vectors
        Rkk = V(:,1:iter)'*V(:,iter+1);
        V(:,iter+1) = V(:,iter+1) - V(:,1:iter)*Rkk;
        
        % Copy V to temporary cell (block) array
        %b = 4; %TODO: take this as option.
        %V_ = cell(1,ceil((iter+1)/b));
        %for i = 1:ceil((iter+1)/b)
        %    cols = b*(i-1)+1:min(b*i,iter+1);
        %    V_{i} = V(:,cols);
        %end
        % Do reorthogonalization on blocks
        %[V_,R_] = projectAndNormalize(V_);
        % Copy reorthogonalized blocks back to V
        %for i = 1:ceil((iter+1)/b)
        %    cols = b*(i-1)+1:min(b*i,iter+1);
        %    V(:,cols) = V_{i};
        %end
    end   
end


    