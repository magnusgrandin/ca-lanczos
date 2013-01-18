% Function [T,Q] = sstep_lanczos_prop(H,r0,s,m,dt)
% 
% Input:
%    H    - Hamiltonian matrix
%    psi  - Wavefunction to propagate
%    dt   - Time step
%    s    - Number of simultaneous Lanczos steps
%    m    - Number of outer iterations
%
% Output:
%    T    - Tri-diagonal Lanczos matrix
%    v    - Krylov subspace

function [T,Q,lsteps]= sstep_lanczos_prop(H,r0,s,m,dt,tol)
    N = length(r0);
    
    % Initialize some arrays
    E  = initMatrixBlockArray(m,s,s);     %gamma (super-diagonal blocks)
    F  = initMatrixBlockArray(m,s,s);     %      (sub-diagonal blocks)
    G  = initMatrixBlockArray(m+1,s,s);   %alpha (diagonal blocks)
    V = initMatrixBlockArray(m+1,N,s);

    W = zeros(s,s);
    %vav = zeros(s,s);
    v1av1 = zeros(s,s);
    t  = zeros(s,s);
    Q  = zeros(N,s*m);
    P  = zeros(N,s+1);
    T  = zeros(s*m);
    dotP = zeros(2*s,1);
    b  = zeros(s,s);
    c  = zeros(s,s);
    d  = zeros(s,s);

    % Compute norm of starting vector
    nrm = norm(r0);
    
    % Use the normalized wave function as starting vector
    P(:,1) = r0/(nrm);
    
    % Apply hamiltonian s times
    %P = hams(P, H, s);
    for i=2:s+1
        P(:,i) = H*P(:,i-1);
    end
    
    V{2}(:,1:s) = P(:,1:s);
    
    % Compute 2*s dot-products
    %dotP = dotp2s(P, s);
    for i=1:s
        ix = 2*(i-1)+1;
        dotP(ix) = P(:,i)'*P(:,i);
        dotP(ix+1) = P(:,i+1)'*P(:,i);
    end
    
    % Compute inner products (V{k}(:,1), A^(i+j-1)*V{k+1}(:,1))
    %v1av1 = inner_prods(dotP, zeros(s,s), s, V, P, 1);
    for j=1:s
        for i=1:s
            if i+j-s > 0
                v1av1(i,j) = dotP(i+j-s);
            end
        end
    end
    

    %%% Loop starts here
    for k=2:m+1

        % Decompose c and solve for gamma (i.e. E), using W{k-1}
        %c = comp_c(dotP, s);
        c(1:s-1,:) = 0.0;
        c(s,:) = dotP(1:s);
        if k > 2
            for j = 1:s
                E{k-1}(:,j) = W\c(:,j);
            end
        end
        
        % Decompose W (i.e. W{k})
        %W = comp_W(dotP, v1av1, t, s, V, k);
        for j=1:s
            for i=j:s
                W(i,j) = dotP(i+j-1);
                %r = 2*(s+1)-(i+j);
                r = s+2-j;
                for l=r:s
                    W(i,j) = W(i,j) - t(l,i)*v1av1(l,j-1);
                end
                W(j,i) = W(i,j);
            end
        end
        
        % Decompose d and solve for alpha (i.e. G)
        %d = comp_d(W, dotP, v1av1, t, c, s, V, H, k);
        for j=1:s-1
            for i=j:s
                d(i,j) = W(i,j+1)-t(s,j)*c(s,i);
                d(j,i) = d(i,j);
            end
        end
        d(s,s) = dotP(2*s)-t(s,s)*c(s,s);
        for i=1:s
            d(s,s) = d(s,s)-t(i,s)*v1av1(i,s);
        end
        for j = 1:s
            G{k}(:,j)=W\d(:,j);
        end
        
        % Compute v(k+1,1)
        P(:,1) = H*V{k}(:,s) - V{k-1}*E{k-1}(:,s) - V{k}*G{k}(:,s);
        
        %Set nonzero element of F to the residual norm (is this the way?)
        F{k}(1,s) = 1; %sqrt(P(:,1)'*P(:,1));
        
        % Compute the residual and stop iterations if tolerance fulfilled
        if(s*(k-1) > 2)
            % Assemble block tridiagonal Lanczos matrix
            T = assemble_T(E,F,G,m,s);
            [Vp,D] = eig(T(1:s*(k-1),1:s*(k-1)));
            matexp = Vp*expm(-1i*dt*D)/Vp;
            residual = abs(dt*matexp(s*(k-1),1)*norm(P(:,1))*nrm);
            %if residual < tol
            %    break;
            %end
        end
        
        if k < m+1
            % Compute matrix-vector products
            %P = hams(P, H, s);
            for i=2:s+1
                P(:,i) = H*P(:,i-1);
            end
            
            % Compute 2*s dot-products
            %dotP = dotp2s(P, s);
            for i=1:s
                ix = 2*(i-1)+1;
                dotP(ix) = P(:,i)'*P(:,i);
                dotP(ix+1) = P(:,i+1)'*P(:,i);
            end
            
            % Compute inner products (V1{k}, A^(i+j-1)*V1{k+1})
            %v1av1 = inner_prods(dotP, G{k}, s, V, P, k);
            for j=1:s
                for i=1:s
                    if i+j-s > 0
                        v1av1(i,j) = dotP(i+j-s);
                    end
                    r = 2*(s+1)-(i+j);
                    for l=r:s
                        v1av1(i,j) = v1av1(i,j) + G{k}(l,s)*v1av1(l,(i+j)-(s+1));
                    end
                end
            end
            
            %Assemble b and solve for t
            %b = comp_b(v1av1, s, V, P, k);
            for j=2:s
                for i=s-j+2:s
                    b(i,j) = v1av1(i,j-1);
                end
            end
            for j = 1:s
                t(:,j) = W\b(:,j);
            end
            %t(:,2) = E{k-1}(:,1);
            
            % Compute orthogonal vectors
            V{k+1}(:,1) = P(:,1);
            for j=2:s
                V{k+1}(:,j) = P(:,j) - V{k}*t(:,j);
            end
        end
                               
    end %% Main loop
    
    %Assemble Q
    for j=1:k-1
        kx = (j-1)*s;
        Q(:,kx+1:kx+s) = V{j+1};
    end
    
    lsteps = s*(k-1);
    T = T(1:lsteps,1:lsteps);
    Q = Q(:,1:lsteps);
    
end


function M = initMatrixBlockArray(l,s1,s2)
    M = cell(1,l);
    for i=1:length(M)
        M{i} = zeros(s1,s2);
    end
end

function dotp = dotp2s(V, s)
    dp = zeros(1,2*s);
    for i=1:s
        ix = 2*(i-1)+1;
        dp(ix) = V(:,i)'*V(:,i);
        dp(ix+1) = V(:,i+1)'*V(:,i);
    end
    dotp = dp;
end

function hV = hams(V, H, s)
    for i=2:s+1
        V(:,i) = H*V(:,i-1);
    end
    hV = V;
end

function d = comp_d(W, dotP, v1av1, t, c, s, V, H, k)
    d = zeros(s,s);
    for j=1:s-1
        for i=j:s
            d(i,j) = W(i,j+1)-t(s,j)*c(s,i);
            d(j,i) = d(i,j);
        end
    end
    d(s,s) = dotP(2*s)-t(s,s)*c(s,s);
    for i=1:s
        d(s,s) = d(s,s)-t(i,s)*v1av1(i,s);
    end
end

function b = comp_b(v1av1, s, V, P, k)
    b = zeros(s,s);
    for j=2:s
        for i=s-j+2:s
            b(i,j) = v1av1(i,j-1);
        end
    end
end

function c = comp_c(dotP, s)
    c = zeros(s,s);
    c(1:s-1,:) = 0.0;
    c(s,:) = dotP(1:s);
end

function W = comp_W(dotP, v1av1, t, s, V, k)
    W = zeros(s,s);
    for j=1:s
        for i=j:s
            W(i,j) = dotP(i+j-1);
            %r = 2*(s+1)-(i+j);
            r = s+2-j;
            for l=r:s
                W(i,j) = W(i,j) - t(l,i)*v1av1(l,j-1);
            end
            W(j,i) = W(i,j);
        end
    end
end

function ip = inner_prods(dotP, G, s, V, P, k)
    ip = zeros(s,s);
    for j=1:s
        for i=1:s
            if i+j-s > 0
                ip(i,j) = dotP(i+j-s);
            end
            r = 2*(s+1)-(i+j);
            for l=r:s
                ip(i,j) = ip(i,j) + G(l,s)*ip(l,(i+j)-(s+1));
            end
        end
    end
end

function T = assemble_T(E,F,G,m,s)
    %Assemble T (matrix blocks for k = 1 are all zero)
    for k=1:m
        ix = (k-1)*s;
        T(ix+1:ix+s, ix+1:ix+s) = G{k+1};
        if k < m
            T(ix+(s+1):ix+2*s, ix+1:ix+s) = F{k+1};
            T(ix+1:ix+s, ix+(s+1):ix+2*s) = E{k+1};
        end
    end
end
