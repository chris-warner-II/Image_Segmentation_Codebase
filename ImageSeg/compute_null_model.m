function [ NM ] = compute_null_model(W,W_conn,Mask,topo,maskFlg)

% syntax: [ NM ] = compute_null_model(W,W_conn)
%
% where NM is the output null model (expected connectivity of nodes in graph).
%       W is the weights (adjacency) matrix that describes actual connectivity in graph.
%       W_conn is a 2nd input weights matrix calculated from a uniform
%       image. It is used for normalization to turn numbers into
%       probabilities.  Neutral (uniform) images will have largest possible
%       connectivity / weights between all pixels in image.  Neutral image
%       weights matrix tells us about the connectivity of the "hardware"
%       not based at all on the input image.
%       TOPO is a flag to use topographical (di*dj*bij) version, otherwise 
%       use original Newman & Girvan (di*dj) version.

% Allocate Sparse Matrices for data
imsz = size(W,1); % number of pixels in image.
pctNN = max(sum(logical(W_conn)))/(imsz-1);   % ratio of number of connections by number of possible connections
NM = spalloc(imsz,imsz,ceil(pctNN*numel(W)));

% Compute row sums or node incidences
D = sum(W);           % incidence into a node

% Null Model is just NM_{ij} = D_i*D_j*B_{ij} because D & B are normalized
% by W_neutral.  Also, the NM is constrained to have same topology as
% W_neutral because it uses the same "Hardware".

%% NOTE: I may not need the d variable coming out of PIV_vs_Dist functions.
if(topo) % TOPOGRAPHICAL
    
    switch(maskFlg)
        case(0) % Use only Diagonal Means
            tag='SKH TopographicalD';
            [B] = PIV_vs_DiagDist(W,W_conn);
            tag2 = 'DiagD';
%             for i=1:numel(B)
%                 Bm = Bm + diag(B(i)*ones(imsz+1-i,1),i-1) + diag(B(i)*ones(imsz+1-i,1),-(i-1));
%             end
            
        case(1)  % Use Euclidian Distance in Image vs PIV means with Mask
            tag='SKH TopographicalD';
            [B] = PIV_vs_MaskDist(W,Mask,W_conn);
            tag2 = 'MaskD';
%             for i=1:numel(B)
%                 Bm = Bm + B(i).*Mask(i).mask;                  
%             end
        case(2) % Use different masks for different distance & orientation (vert vs horiz) 
            tag='SKH TopographicalDO';
            [B] = PIV_vs_MaskDist(W,Mask,W_conn);
            tag2 = 'MaskDO';
%             for i=1:numel(B)
%                 Bm = Bm + B(i).*Mask(i).mask;
%             end
    end
    D = D./sum(W_conn);       % Mmmmmmm... Seems important.
    NM = (D'*D).*B.*W_conn;
    
else % NONTOPO       

    tag='NG Nontopographical';
    NM = D'*D./sum(W(:)); % probability or pct incidence from total weight of graph
    NM = NM - diag(diag(NM)); % put zeros on diagonal for no self-weights
    tag2 = '';

end


% disp([tag,' ',tag2])

if(0)
% Plot the Null Model
    figure, imagesc(NM),colorbar, title([tag,' Null Model ',tag2],'FontSize',20,'FontWeight','Bold')
    keyboard
    
end

