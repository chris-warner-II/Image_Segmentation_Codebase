function L = compute_Laplacian(W,normlze,neg)

% syntax: L = compute_Laplacian(W,normlze,project);
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


%% New Way (yields exactly the same results as old way - all 3 methods)
if(normlze)
    
    % old method to Normalize (Interestingly, this is the Fastest. Set 1 below to see.)
%     tic
    degW = sum(W,1);              % degree of incidence into each vertex
    OD = -W./sqrt(degW'*degW);    % off-diagonal elements for NGL
    L = OD;                       % Normalized Graph Laplacian
    for i = 1:size(W)
        L(i,i)=1;                 % ones on the diagonal
    end
%     disp('Old Method to Normalize:')
%     toc
    
%     if(0)
%         % New Method 1 to Normalize
%         tic
%         o1=1./sqrt(D);
%         o1(o1==inf)=0;
%         L1=o1*L*o1';
%         disp('New Method 1 to Normalize:')
%         toc
% 
%         % New Method 2 to Normalize
%         tic
%         o2=full(D)^(-1/2);
%         o2(o2==inf)=0;
%         L2=o1*L*o1';
%         disp('New Method 2 to Normalize:')
%         toc
% 
%         % New Method 3 to Normalize
%         tic
%         o3=full(D).^(-1/2);
%         o3(o3==inf)=0;
%         L3=o1*L*o1';
%         disp('New Method 3 to Normalize:')
%         toc
% 
%         disp('Comparing Outputs from Different Methods - Should all be zero.')
%         a = L1 - L2;
%         b = L1 - L3;       % these should be all 0 matrices up to machine precision.
%         c = L1 - full(L);
%         max(a(:))
%         max(b(:))
%         max(c(:))
% 
%         % Plot old and new methods side by side alongside adjacency matrix
%         figure, subplot(311),imagesc(W),colorbar,title('Adjacency')
%         subplot(312),imagesc(L),colorbar,title('Normalized Laplacian Old')
%         subplot(313),imagesc(L1),colorbar,title('Normalized Laplacian New')
%     end
    
else
    
    Dvec = sum(W);  % The Incidence vector (row sums)
    D = diag(Dvec); % Incidence Matrix (incidence vector on diagonal of matrix)
    L = D - W;      % Graph Laplacian

end

% % Normalize by dividing by maximum value in AA Matrix
% L = L./max(abs(L(:)));

if(neg)
    L = -L; % Negative (Normalized or Non.) Graph Laplacian
end

% if(0) % For Debugging (best to do this with a small random image input
%     figure, 
%     subplot(131), imagesc(W), colorbar('SouthOutside'), title('Adjacency W')
%     subplot(132), imagesc(D), colorbar('SouthOutside'), title('Degree D')
%     subplot(133), imagesc(L), colorbar('SouthOutside'), title('Laplacian L')
%     %
%     keyboard
% end


