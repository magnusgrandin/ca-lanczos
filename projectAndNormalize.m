function [Q,R] = projectAndNormalize(Q,X,doreorth)

    if nargin < 3
        doreorth = true;
    end
    
    % TODO: tune this
    tol = .5;
    
    nrows = size(X,1);
    ncols = size(X,2);
    numBlocksQ = length(Q);
    R = cell(1,numBlocksQ+1);
    
    % Compute norms of each column of X before first orthogonalization
    if doreorth == true
        normsBeforeFirst = zeros(ncols,1);
        for i = 1:ncols
            normsBeforeFirst(i) = sqrt(sum(X(:,i).^2));
        end
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
    
    if doreorth == true
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
            [Z,R_] = project(Q,Y); 
            [QZ,RZ,rank] = normalize(Z);
            % Update coefficients after second pass of orthogonalization
%            RZ = RZ+RY;
%            for i = 1:numBlocksQ
%                R{i} = R{i} + R_{i}*RY;
%            end
            for i = 1:numBlocksQ
                R{i} = R{i} + R_{i};
            end
            % If R is not full rank, the last ncols-rank columns have been
            % randomized and orthogonalized within QZ. Orthogonalize those 
            % columns of QZ against the previous Q-blocks, don't need to 
            % keep the coefficients.
            if rank < ncols
                disp('Rank deficient');
                nullSpaceCols = rank+1:ncols;
                [QZ(:,nullSpaceCols),R_] = project(Q,QZ(:,nullSpaceCols));
                [QZ(:,nullSpaceCols),R_] = tsqr(QZ(:,nullSpaceCols));
            end
        end
    end

    if doreorth == true
        Q = QZ;
        R{numBlocksQ+1} = RZ;
    else
        Q = QY;
        R{numBlocksQ+1} = RY;
    end
end
