function [Q,H] = arnoldi(A,Q,H,q,maxvecs,prevvecs,orth)

    if nargin < 4
        orth = 'local';
    end
    
    n = size(Q,1);
    
    % If necessary, allocate enough space in Q and H
    if size(Q,2) < maxvecs+1
        Q = [Q, zeros(n,maxvecs+1-size(Q,2))];
    end
    [Hrows,Hcols] = size(H);
    if Hrows < maxvecs+1 || Hcols < maxvecs
        H_ = H;
        H = zeros(maxvecs+1,maxvecs);
        H(1:Hrows,1:Hcols) = H_;
    end
    
    if prevvecs > 0
        b = sqrt(Q(:,prevvecs+1)'*Q(:,prevvecs+1));
    else
        b = H(Hrows,Hcols);
    end
    
    nvecs = prevvecs+1;
    iter = 0;
    while (nvecs < maxvecs+1) 
        iter = iter+1;
        Q(:,nvecs+1) = A*Q(:,nvecs);
        h = Q(:,1:nvecs)'*Q(:,nvecs+1);
        Q(:,nvecs+1) = Q(:,nvecs+1) - Q(:,1:nvecs)*h;
        g = sqrt(Q(:,nvecs+1)'*Q(:,nvecs+1));
        Q(:,nvecs+1) = (1/g)*Q(:,nvecs+1);
        H(1:nvecs,nvecs) = h;
        nvecs = nvecs+1;
    end
    
    % Fix output
    Q = Q(:,1:nvecs);

end
