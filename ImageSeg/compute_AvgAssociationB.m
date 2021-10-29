function A = compute_AvgAssociationB(W,normlze)

% syntax: [A] = compute_AvgAssociationB(W,normlze);

disp('In compute_AvgAssociationB : ')

old_way=0; % flag to do or not do things using full matrices.

if(normlze) % use the Normalized Average Association (normalizing helps).
    
    A  = spalloc( size(W,1), size(W,2), numel(find(W)) );
    degA = sum(W);           % Degree of incidences into each vertex

    % New Way : Using Batches to be more memory efficient for large
    % matrices.  The idea here is to figure out how far out the
    % non-zero elements reach from the diagonal and then to calculate
    % the normalization part only for those areas where the numerator
    % is not going to be zero.  Do same thing for Normalized Graph
    % Laplacian certainly.  Maybe for Modularity too?..

    disp('New Way to Compute Normalized Avg Assoc.  In Batches')
    tic
    %
    [x,y] = find(W);
    z = max(abs(x-y)); % how far out from diagonal of W non-zero elements reach.
    %
    batchBeg = [1:z:size(W,1)-z];
    batchEnd = [2*z:z:size(W,1)];
    if ( numel(batchEnd) < numel(batchBeg) )
        batchEnd = [batchEnd,size(W,1)];
    end
    %
    disp(['Batch Iteration # out of ',num2str(numel(batchBeg))])
    for i = 1:numel(batchBeg)

        fprintf('%s',[num2str(i),' '])

        st = batchBeg(i);
        nd = batchEnd(i);

        A(st:nd,st:nd) = W(st:nd,st:nd) ./ ( sqrt(degA(st:nd)'*degA(st:nd)) ) ;
    end
    %
    toc


    % Old Way to Compute : Using Full Dense Matrices
    if(old_way)
        disp('Old Way to Compute Normalized Avg Assoc.  All at once.')
        tic
        D = W./sqrt(degA'*degA); % 
        A2 = D;                   % 
        toc

        figure, imagesc(A2-A), title('Difference between batch and full computations of AAnrm')

        keyboard

    end
    
    

else
    
    disp('Unnormalized Avg Association : Full Substitution')
    tic
    A = W;
    toc
end
