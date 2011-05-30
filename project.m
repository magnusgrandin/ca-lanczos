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
        if ~isempty(Q{i})
            R{i} = Q{i}'*X;
            X = X - Q{i}*R{i};
        else
            R{i} = [];
        end
    end
end
 