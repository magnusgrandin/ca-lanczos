%% function [T,Q,nLanczos] = lanczos_prop(H,r0,maxiter,dt,tol,adaptive)

%   The Lanczos algorithm for time integration
%   
%   Input variables:
%     H        - hamiltonian
%     psi      - starting vector for the lanczos iteration
%     maxiter  - maximum size of the Krylov space
%     dt       - time step size
%     
%   Output variables:
%     T   - tridiagonal Lanczos matrix
%     Q   - matrix of Krylov basis vectors

function [T,Q,nLanczos] = lanczos_prop(H,r0,maxiter,dt,tol,adaptive)

    if nargin < 5
        tol = 1.0e-10;
    end
    if nargin < 6
        adaptive = false;
    end
    
    N = length(r0);
    Q = zeros(N,maxiter);
    nrm = norm(r0);
    Q(:,1) = r0/(nrm);
    alpha = zeros(1,maxiter);
    beta = zeros(1,maxiter);
    
    for j=1:maxiter
        if j == 1      
               r = H*Q(:,j);
        else
            Q(:,j) = r/beta(j-1);  
            r = H*Q(:,j);
            r=r-beta(j-1)*Q(:,j-1);
        end
        alpha(j) = r'*Q(:,j);
        r = r - alpha(j)*Q(:,j);
        beta(j) = sqrt(r'*r);
            
        % Compute the residual and stop iterations if tolerance fulfilled
        if j >= 3 
            T = diag(alpha(1:j)) + diag(beta(1:j-1),1) + diag(beta(1:j-1),-1);
            T = [T; [zeros(1,j-1), beta(j)]];
            [V,D] = eig(T(1:j,1:j));
            matexp = V*expm(-1i*dt*D)*V';
            residual = abs(dt*beta(j)*matexp(j,1)*nrm);
            %disp(['residual(', num2str(j), '): ', num2str(residual)]);
            if (residual < tol) && (adaptive == true)
                disp(num2str(residual));
                break;
            end
        end
    end
    disp(num2str(residual));
    nLanczos = j;
    T = T(1:nLanczos,:);
    Q = Q(:,1:nLanczos);
end
   
% -------------------------------------------------------------------------
%  Copyright (2012, 2013)  Magnus Grandin <magnus.grandin@it.uu.se>
% -------------------------------------------------------------------------
