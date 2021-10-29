function [ Bm ] = PIV_vs_MaskDist(W,Mask,W_conn)

% syntax: [ B ] = PIV_vs_Dist(W,W_conn)
%
% INPUT: [W & W_conn]
% This function takes in an Adjacency Matrix and a "Hardware" Connectivity
% Matrix (which is equivalent to an Adjacency Matrix calculated on a 
% uniform/null image).  
%
% OUTPUT: [B]
% It outputs a vector B which is the Expected Weight Value (averaged over
% the entire image) of pixels separated by a distance and direction given
% by the distance from the Identity Diagonal of the matrix. The first entry
% in the B vector should be zero because the Adjacency Matrix is
% constrained to have no self-connectivity and therefore the diagonal is
% full of zero weights.  Each B entry is calculated by the ratio of the sum 
% of entries in a diagonal of the Adjacency Matrix and the sum of entries
% in the corresponding diagonal of the Connectivity Matrix.  Each B entry
% take a continuous value between 0 and 1.  The B vector is   

%% Allocate Sparse Matrices for data
imsz = size(W,1); % number of pixels in image.
pctNN = max(sum(logical(W_conn)))/(imsz-1);   % ratio of number of connections by number of possible connections
Bm = spalloc(imsz,imsz,ceil(pctNN*numel(W)));


for i = 1:numel(Mask)
%     figure, imagesc(Mask(i).mask),colorbar,title(['Distance = ',num2str(Mask(i).distance)])
%     d(i) = Mask(i).distance;
    B(i) = mean(mean(W(Mask(i).mask)));
end

for i=1:numel(B)
    Bm = Bm + B(i).*Mask(i).mask;                  
end

% %Plot B (Average Weight) vs. d (distance between pixels in image)
% figure, plot(d,full(B),'b','LineWidth',2)
% title('Gradient Box Pixel Intensity Distance Dependence','FontSize',20,'FontWeight','Bold')
% xlabel('Distance (in Pixels) between 2 Pixels','FontSize',18,'FontWeight','Bold')
% ylabel('Average Weight Value for all Pixels in Image Separated by Distance','FontSize',18,'FontWeight','Bold')