function [Q,R,rank] = normalize(X,tol,opt)

    if nargin < 2
        tol = 1.0e-12;
    end
    if nargin < 3
        opt = 'none';
    end
    
    nrows = size(X,1);
    ncols = size(X,2);
    [Q,R] = tsqr(X);
    [U,S,W] = svd(R);
    S_diag = diag(S);
    rank = ncols;
    for i = 1:ncols
        if S_diag(i) <= eps
            rank = i-1;
            break;
        end
    end
    if(rank == ncols)
        % Q,R and rank already computed.
        return;
    elseif strcmpi(opt,'randomizeNullSpace')
        R = S*W';
        Q = Q*U;
        Q = randomizeNullSpace(Q,rank);
    end
end
