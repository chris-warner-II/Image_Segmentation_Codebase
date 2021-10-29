function L = compute_LaplacianB(W,normlze,neg)

% syntax: L = compute_LaplacianB(W,normlze,neg);
%
%   where: L is the output Graph Laplacian matrix L = D - W
%          D is the Incidence Matrix (a weighted diagonal matrix)
%          W is the Adjacency (or Weights) Matrix calculated in calc_weights
%
%   note: (1). normlze is a flag [0,1] to use the Normalized Graph
%              Laplacian ( D^-(1/2)*L*D^-(1/2) )
%         (2). neg is a flag [0,1] to use the Negative Graph Laplacian (-L)
%         (3). Note: by using both flags, you can implement the Normalized
%              Negative Graph Laplacian
%
%   The Graph Laplacian provides a good image segmentation if you find the
%   2nd smallest eigenvector (or the smallest nonzero eigenvector).  The
%   zero eigenvector corresponds to the uninteresting 1-vector.
%
%   TO ADD:  WHAT ABOUT HILLAR'S OTHER KIND OF NORMALIZATION?? (PUT IN)


disp(' Inside compute LaplacianB: ')


%% New Way (yields exactly the same results as old way - all 3 methods)
if(normlze)
    
    % old method to Normalize (Interestingly, this is the Fastest. Set 1 below to see.)
    A = compute_AvgAssociationB(W,normlze);    % off-diagonal elements for NGL is from just avg association.
    L = eye(size(W)) - A;                     % Normalized Graph Laplacian (ones on diagonal)
        
else
    
    D = diag(sum(W));  % Incidence Matrix (incidence vector on diagonal of matrix)
    L = D - W;         % Graph Laplacian

end



if(neg)
    L = -L; % Negative (Normalized or Non.) Graph Laplacian
end



% show some plots and go into debug...
if(0)
    
    figure, 
    subplot(121), imagesc(W), colorbar('SouthOutside'), title('W')
    subplot(122), imagesc(L), colorbar('SouthOutside'), title('L')

    keyboard

end


