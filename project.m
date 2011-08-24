%  Compute the projection of X onto Q
%
%     C = Q'*X;  X = X - Q_*C;
%
function [X,R] = project(Q,X,doreorth)
    if nargin < 3
        doreorth = false;
    end
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
    numVectors = size(X,2);
    R = cell(1,numBlocks);
    normBefore = zeros(numVectors,1);
    if doreorth
        normBefore = multiVecNorm2(X);
    end
    for i = 1:numBlocks
        if ~isempty(Q{i})
            R{i} = Q{i}'*X;
            X = X - Q{i}*R{i};
        else
            R{i} = [];
        end
    end
    if doreorth
        disp('project(): reorthogonalize');
        R2 = cell(1,numBlocks);
        normAfter = multiVecNorm2(X);
        rho = 0.5;%1/sqrt(2);
        normDiff = rho*normBefore-normAfter;
        if max(normDiff) > 0
            for i = 1:numBlocks
                if ~isempty(Q{i})
                    R2{i} = Q{i}'*X;
                    X = X - Q{i}*R2{i};
                else
                    R2{i} = [];
                end
                R{i} = R{i}+R2{i};
            end
        end
    end
end

function X_norm = multiVecNorm2(X)
    numVectors = size(X,2);
    for j = 1:numVectors
        X_norm(j) = norm(X(:,j));
    end
end
    