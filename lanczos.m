%--------------------------------------------------------------------------
% Time-integrator the Lanczos algorithm for computing the matrix 
% exponential
% 
% Input variables:
% H    = Hamiltonian
% psi  = wave function to be proceeded (at time time)
% dt   = time step
% iter = order of the Krylov space
%
% Output variable:
% res  = approximation of the wave function at time time+dt
%--------------------------------------------------------------------------
function [T,v] = lanczos(H,psi,iter)
    N = length(psi);
    v = zeros(N,iter);
    v(:,1) = psi/(sqrt(psi'*psi));
    alpha = zeros(1,iter);
    beta = zeros(1,iter-1);

    for j = 1:iter
        if j == 1      
               r = H*v(:,j);
        else
            v(:,j) = r/beta(j-1);  
            r = H*v(:,j);
            r=r-beta(j-1)*v(:,j-1);
        end
        alpha(j) = r'*v(:,j);
        r = r - alpha(j)*v(:,j);
        if j ~= iter
            beta(j) = sqrt(r'*r);
        end
    end
    
    T = diag(alpha) + diag(beta,1) + diag(beta,-1);
end