%%
%   function [Q,R,rank] = rr_tsqr_bgs(V,tol)
%
%   Rank revealing TSQR with block Gram-Schmidt orthogonalization
%
%   Input:
%     opt - 'project','normalize' or 'projectAndNormalize'
%     Q   - cell (block) array with orthogonal vectors
%     V   - vectors to be orthogonalized against the vector blocks in Q
%     tol - the tolerance in rank deficiency (defalt: 1.0e-10)
%
%   Output:
%     Q   - cell array with blocks of orthogonalized vectors
%     R   - the R factor in V = QR (dense matrix)
%
%
function [Q,R,rank] = rr_tsqr_bgs(opt,Q,V,tol)
    
    
    nrows = length(V{1}(:,1));
    Vcols = size(V{1},2);
    nblocks = length(Q);
    R = {};
    if strcmpi(opt,'normalize') == 1
        Qcols = size(Q{1},2)
        R{1} = zeros(Qcols,Vcols);
    elseif strcmpi(opt,'norma')
        for k = 1:M
        Qcols = size(Q{k},2);
        R{k} = zeros(Qcols,Vcols);
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

   


