function [Q,R,rank] = normalize(X,opt,tol)

    if nargin < 2
        opt = 'None';
    end
    if nargin < 3
        tol = 1.0e-14;
    end
    
    nrows = size(X,1);
    ncols = size(X,2);
    [Q,R] = tsqr(X);
    [U,S,W] = svd(R);
    S_diag = diag(S);
    abs_tol = tol*S_diag(1);
    rank = ncols;
    for i = 1:ncols
        if S_diag(i) <= abs_tol
            rank = i-1;
            break;
        end
    end
    if(rank == ncols)
        % Q,R and rank already computed.
        return;
    elseif strcmpi(opt,'randomizeNullSpace')
        R = S*W';
        Q = Q*U;%, zeros(nrows,ncols-rank)];
        Q = randomizeNullSpace(Q,rank);
    end
end

function Q = randomizeNullSpace(Q,rank)

    disp('Randomize null space.');
    disp(['Rank ' num2str(rank)]);
    nrows = size(Q,1);
    ncols = size(Q,2);
    numNullSpaceCols = ncols - rank;
    fullRankColIndices = 1:rank;
    nullSpaceColIndices = rank+1:ncols;
    Q(:,nullSpaceColIndices) = rand(nrows,numNullSpaceCols);
    [Q(:,nullSpaceColIndices),R_] = project({Q(:,fullRankColIndices)},Q(:,nullSpaceColIndices));
    [Q(:,nullSpaceColIndices),R_] = tsqr(Q(:,nullSpaceColIndices));
    
end
