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

cum_cols = 0;

for k = 1:M
    
    mk = size(V{k},2);
    
    % Get the columns of Q{1..k}
    Q_ = [];
    for j = 1:k-1
        Q_ = [Q_ Q{j}];
    end

    % Compute norms of each column of Vk
    norms_Vk = zeros(mk,1);
    for i = 1:mk
        norms_Vk(i) = sqrt(sum(V{k}(:,i).^2));
    end
    
    % First block orthogonalization pass
    if k == 1
        Yk = V{k};
    else
        Rkk_ = Q_'*V{k};
        Yk = V{k} - Q_*Rkk_;
        R(1:cum_cols,cum_cols+1:cum_cols+mk) = Rkk_;
    end
    [QYk,RYk] = tsqr(Yk);

    norms_Yk = zeros(mk,1);
    for i = 1:mk
        norms_Yk(i) = sqrt(sum(RYk(:,i).^2));
    end
    
    % If any column norm drops too much, do second
    % pass of orthogonalization
    if max(abs(norms_Vk-norms_Yk)./norms_Vk) > 1/10 %TODO: tune this
        disp('second');
        if(k == 1)
            Zk = Yk;
        else
            Rkk_ = Q_'*Yk;
            Zk = Yk - Q_*Rkk_;
            R(1:cum_cols,cum_cols+1:cum_cols+mk) = R(1:cum_cols,cum_cols+1:cum_cols+mk) + Rkk_;
        end
        [QZk,RZk] = tsqr(Zk);
    else
        Zk  = Yk;
        QZk = QYk;
        RZk = RYk;
    end
    
    % Compute SVD of RZk and get its rank
   [U,S,W] = svd(RZk);
   S_diag = diag(S);
   rk = 0;
   for i = 1:mk
       if S_diag(i) > eps
           rk = rk + 1;
       end
   end
   if rk == mk
        Q{k} = QZk;
        R(cum_cols+1:cum_cols+mk,cum_cols+1:cum_cols+mk) = RZk;
    else
        disp('rank deficient');
        numNullSpaceCols = mk - rk;
        nullSpaceColIndices = rk+1:rk+numNullSpaceCols;
        Rk = S*W';
        Q{k} = Q{k}*U;
        Q{k}(:,nullSpaceColIndices) = rand(nrows,numNullSpaceCols);
        R_ = Q{k}'*Q{k}(:,nullSpaceColIndices);
        Q{k}(:,nullSpaceColIndices) = Q{k}(:,nullSpaceColIndices) - Q{k}*R_;
        [Q{k}(:,nullSpaceColIndices),R_] = tsqr(Q{k}(:,nullSpaceColIndices));
        R(cum_cols+1:cum_cols+mk,cum_cols+1:cum_cols+mk) = RZk;
    end
        
    cum_cols = cum_cols+mk;
    
end
R = R(1:cum_cols,1:cum_cols);
end

function [Q,R,rank_Y] = normalize(Y,tol)

nrows = size(Y,1);
ncols = size(Y,2);
[Q,R] = tsqr(Y);
[U,S,W] = svd(RZk);
S_diag = diag(S);
for i = 1:mk
   if S_diag(i) <= eps
           rk = i;
           break;
       end
   end
   
rank_Y = rank(R);
if(rank_Y < ncols)
    numNullSpaceCols = ncols - rank_Y;
    nullSpaceColIndices = rank_Y+1:rank_Y+numNullSpaceCols;
    Q(:,nullSpaceColIndices) = rand(nrows,numNullSpaceCols);
    R_ = Q'*Q(:,nullSpaceColIndices);
    Q(:,nullSpaceColIndices) = Q(:,nullSpaceColIndices)-Q*R_;
    [Q,R_] = tsqr(Q);
end
end
