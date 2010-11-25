%%
%   Rank revealing TSQR with block Gram-Schmidt orthogonalization
%
%   Input:
%     V   - Cell array with blocks of vectors to be orthogonalized
%     tol - The tolerance in rank deficiency (defalt: 1.0e-10)
%
%   Output:
%     Q   - Cell array with blocks of orthogonalized vectors
%     R   - The R factor in V = QR (dense matrix)
%
%
function [Q,R,rank] = rr_tsqr_bgs(V,tol)
    
    if nargin > 1
        eps = tol;
    else
        eps = 1.0e-10;
    end
    
    M = length(V);
    nrows = length(V{1}(:,1));
    ncols = 0;
    for i = 1:M
        ncols = ncols + size(V{i},2);
    end
    R = zeros(ncols,ncols);
    Q = cell(1,M);
    for i = 1:M
        Q{i} = zeros(size(V{i}));
    end
    
    ccols = 0;
    
    for k = 1:M
        
        blockNumCols = size(V{k},2);
        [Q{k},Rk] = projectAndNormalize(Q(1:k-1),V{k});

        crows = 0;
        updateCols = ccols+1:ccols+blockNumCols;
        for i = 1:length(Rk)
            blockNumRows = size(Rk{i},1);
            updateRows = crows+1:crows+blockNumRows;
            R(updateRows,updateCols)=Rk{i};
            crows = crows+blockNumRows;
        end
        ccols = ccols+blockNumCols;
        
    end

end

%  Compute the projection of X onto Q
%
%     C = Q'*X;  X = X - Q_*C;
%
function [X,R] = project(Q,X)
    % project() only works when Q is a cell array and X is not.
    if ~iscell(Q)
        disp('Input Q (arg 1) to project() must be cell (block) array.');
        return;
    end
    if iscell(X)
        disp('Input X (arg 2) project() must be a column matrix.');
        return;
    end    
    % Quick exit if empty array or no blocks.
    if isempty(Q) || isempty(Q{1})
        R = {};
        return;
    end
    numBlocks = length(Q);
    R = cell(1,numBlocks);
    for i = 1:numBlocks
        R{i} = Q{i}'*X;
        X = X - Q{i}*R{i};
    end
end
    
function [Q,R,rank] = normalize(X,tol)

    nrows = size(X,1);
    ncols = size(X,2);
    [Q,R] = tsqr(X);
    [U,S,W] = svd(R);
    S_diag = diag(S);
    rank = ncols;
    for i = 1:ncols
        if S_diag(i) <= eps
            rank = i;
            break;
        end
    end
    if(rank == ncols)
        % Q,R and rank already computed.
        return;
    else
        R = S*W';
        Q = Q*U;
        Q = randomizeNullSpace(Q,rank);
    end
end

function [Q,R] = projectAndNormalize(Q,X)

    % TODO: tune this
    tol = 1/10;
    
    ncols = size(X,2);
    numBlocksQ = length(Q);
    R = cell(1,numBlocksQ+1);
    
    % Compute norms of each column of X before first orthogonalization
    normsBeforeFirst = zeros(ncols,1);
    for i = 1:ncols
        normsBeforeFirst(i) = sqrt(sum(X(:,i).^2));
    end
    
    % First block orthogonalization pass
    [Y,R] = project(Q,X);
    [QY,RY,rank] = normalize(Y);
    % If R is not full rank, the last ncols-rank columns have been
    % randomized and orthogonalized within QY. Orthogonalize those 
    % columns of QY against the previous Q-blocks.
    if rank < ncols
        nullSpaceCols = rank+1:ncols;
        [QY(:,nullSpaceCols),R_] = project(Q,QY(:,nullSpaceCols));
        [QY(:,nullSpaceCols),R_] = tsqr(QY(:,nullSpaceCols));
    end
    
    % Compute norms of each column of X after first orthogonalization
    normsAfterFirst = zeros(ncols,1);
    for i = 1:ncols
        normsAfterFirst(i) = sqrt(sum(RY(:,i).^2));
    end
    
    % If any column norm drops too much, do second
    % pass of orthogonalization
    if max(abs(normsBeforeFirst-normsAfterFirst)./normsBeforeFirst) > tol
        reorthogonalize = true;
    else
        reorthogonalize = false;
    end
    
    if reorthogonalize == false
        Z  = Y;
        QZ = QY;
        RZ = RY;
    else
        disp('second');
        % Get a copy of the previous coefficients
        [Z,R_] = project(Q,Y); 
        % Add the second pass coefficients to the
        % previous ones.
        for i = 1:numBlocksQ;
            R{i} = R{i} + R_{i};
        end
        [QZ,RZ,rank] = normalize(Z);
        % If R is not full rank, the last ncols-rank columns have been
        % randomized and orthogonalized within QZ. Orthogonalize those 
        % columns of QZ against the previous Q-blocks, don't need to 
        % keep the coefficients.
        if rank < ncols
            nullSpaceCols = rank+1:ncols;
            [QZ(:,nullSpaceCols),R_] = project(Q,QZ(:,nullSpaceCols));
            [QZ(:,nullSpaceCols),R_] = tsqr(QZ(:,nullSpaceCols));
        end
    end
    
    Q = QZ;
    R{numBlocksQ+1} = RZ;
end

function Q = randomizeNullSpace(Q,rank)

    nrows = size(Q,1);
    ncols = size(Q,2);
    numNullSpaceCols = ncols - rank;
    nullSpaceColIndices = rank+1:ncols;
    Q(:,nullSpaceColIndices) = rand(nrows,numNullSpaceCols);
    [Q(:,nullSpaceColIndices),R_] = project(Q,Q(:,nullSpaceColIndices));
    [Q(:,nullSpaceColIndices),R_] = tsqr(Q(:,nullSpaceColIndices));
    
end
