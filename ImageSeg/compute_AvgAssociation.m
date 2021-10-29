function A = compute_AvgAssociation(W,normlze)

% syntax: [A] = compute_AvgAssociation(W,normlze);

if(normlze) % use the Normalized Graph Laplacian
    degA = sum(W);           % Degree of incidences into each vertex
    D = W./sqrt(degA'*degA); % off-diagonal elements for NGL (THIS IS THE REAL ONE!!)
    A = D;                   % Normalized Graph Laplacian
%     for i = 1:size(W)
%         A(i,i)=1;                 % ones on the diagonal (THIS IS DONE IN NGL)
%     end

%     keyboard
    
else
    A = W;
end

% % Normalize by dividing by maximum value in AA Matrix
% A = A./max(abs(A(:)));

% keyboard