function [QZ,RZ] = projectAndNormalize(Q,X,doreorth)

    if nargin < 3
        doreorth = true;
    end
    
    % TODO: tune this
    tol = .5;
    
    nrows = size(X,1);
    ncols = size(X,2);
    numBlocksQ = length(Q);
    
    % Compute norms of each column of X before first orthogonalization
    normsBeforeFirst = zeros(ncols,1);
    if doreorth == true
        for i = 1:ncols
            normsBeforeFirst(i) = sqrt(sum(X(:,i).^2));
        end
    end
    
    % First block orthogonalization pass
    [Y,RY] = project(Q,X,false);
    [QY,R_,rank] = normalize(Y);
    RY{numBlocksQ+1} = R_;
    % If R is not full rank, the last ncols-rank columns have been
    % randomized and orthogonalized within QY. Orthogonalize those 
    % columns of QY against the previous Q-blocks.
%     if rank < ncols
%         disp('Rank deficient');
%         nullSpaceCols = rank+1:ncols;
% %        [QY(:,nullSpaceCols),R_] = project(Q,QY(:,nullSpaceCols),false);
% %        [QY(:,nullSpaceCols),R1] = tsqr(QY(:,nullSpaceCols));
% %        RY(nullSpaceCols,nullSpaceCols) = R1;
%         [QY(:,nullSpaceCols),R_] = project(Q,QY(:,nullSpaceCols),false);
%         [QY(:,nullSpaceCols),R_] = tsqr(QY(:,nullSpaceCols));
%         %[QY,R_] = project(Q,QY,false);
%         %[QY,R_] = tsqr(QY);
%     end
    
    if doreorth == true
        % Compute norms of each column of X after first orthogonalization
        normsAfterFirst = zeros(ncols,1);
        for i = 1:ncols
            normsAfterFirst(i) = sqrt(sum(RY{numBlocksQ+1}(:,i).^2));
        end
        
        % If any column norm drops too much, do second
        % pass of orthogonalization
        if max(abs(normsBeforeFirst-normsAfterFirst)./normsBeforeFirst) > tol
            reorthogonalize = true;
        else
            reorthogonalize = false;
        end
        
        if reorthogonalize == false
            QZ = QY;
            RZ = RY;
        else
            disp('second');
            [Z,RZ] = project(Q,Y,false); 
            [QZ,R_,rank] = normalize(Z,'randomizeNullSpace');
            RZ{numBlocksQ+1} = R_;
            % Update coefficients after second pass of orthogonalization
%            RZ = RZ+RY;
%            for i = 1:numBlocksQ
%                R{i} = R{i} + R_{i}*RY;
%            end
            for i = 1:numBlocksQ
                RZ{i} = RZ{i} + RY{i};
            end
            
            % If R is not full rank, the last ncols-rank columns have been
            % randomized and orthogonalized within QZ. Orthogonalize those 
            % columns of QZ against the previous Q-blocks, don't need to 
            % keep the coefficients.
            if rank < ncols
                disp('Rank deficient');
                nullSpaceCols = rank+1:ncols;
                [QZ(:,nullSpaceCols),R_] = project(Q,QZ(:,nullSpaceCols),false);
                [QZ(:,nullSpaceCols),R_] = tsqr(QZ(:,nullSpaceCols));
            end
        end
    else
        QZ = QY;
        RZ = RY;        
    end
end
