function Q = randomizeNullSpace(Q,rank)

    disp('Randomize null space.');
    disp(['Rank ' num2str(rank)]);
    nrows = size(Q,1);
    ncols = size(Q,2);
    numNullSpaceCols = ncols - rank;
    nullSpaceColIndices = rank+1:ncols;
    Q(:,nullSpaceColIndices) = rand(nrows,numNullSpaceCols);
    [Q(:,nullSpaceColIndices),R_] = project({Q},Q(:,nullSpaceColIndices));
    [Q(:,nullSpaceColIndices),R_] = tsqr(Q(:,nullSpaceColIndices));
    
end
