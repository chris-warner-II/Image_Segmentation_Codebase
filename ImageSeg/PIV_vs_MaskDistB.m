function [ Bm ] = PIV_vs_MaskDistB(W,Wdist,Mask)

% syntax: [ Bm ] = PIV_vs_MaskDistB(W,Wdist,Mask)
%
% This function takes in the Adjacency matix (W), the distance-dependent
% connectivity (Wdist) and the topographic mask structure (Mask) and
% outputs the distance dependent expectation of an edge-weight between two
% nodes, which is based on the average edge-weight between all pixels in
% the image (or all nodes in the network) separated by the same distance.
%
% INPUT: (W,Mask,Wdist)
%   W     - N^2 x N^2 Adjacency Matrix using both distance-dependent connectivity constraints and pixel differences.
%   Wdist - N^2 x N^2 "Hardware" Connectivity constraints imposed by programmer - Distance dependence of network.
%   Mask  - structure which contains a mask of all pixel pairs in adjacency separated by a each distance less than rmax.
%            [ Mask(i).distance = the distance (scalar value).                                                           ].
%            [ Mask(i).mask     = the binary mask (is N^2 x N^2 like W matrix).                                          ]. 
%
% OUTPUT: [Bm]
%   Bm    - N^2 x N^2 matrix with entries at each location being the average edge-weight of all nodes in the adjacency matrix (W)
%           separated by the same distance which separates the nodes referenced by that location in the adjacency matrix.



%% Pre allocate Sparse Matrix for distance dependence of edge-weights
imsz = size(W,1);                             % number of pixels in image.
pctNN = max(sum(logical(Wdist)))/(imsz-1);   % ratio of number of connections to number of possible connections (useful for sparse allocation)
B = zeros(1,numel(Mask));                     % preallocate B vector to hold average edge-weight at each distance.
Bm = spalloc(imsz,imsz,ceil(pctNN*numel(W))); % preallocate Bm matrix to replicate average edge-weight at each location in Adjacency matrix.


plt_flg = 0; % flag to plot B vs d. Set to 0 if you dont want to see.


%% Loop through Mask distances. For each distance, compute average edge weight of all pixels separated by that distance.
for i = 1:numel(Mask)
    
    if(plt_flg)
        figure, imagesc(Mask(i).mask),colorbar,
        title(['Mask at Distance = ',num2str(Mask(i).distance)])
        d(i) = Mask(i).distance;
    end

    ind = find(Mask(i).mask);
    B(i) = mean(W(ind));
    
end




%% Build up the Bm matrix to output from this function by multiplying the mask for a given distance by the average edge-weight at that distance.
for i=1:numel(B)
    Bm = Bm + B(i).*Mask(i).mask;                  
end






%% Plot B (Average Weight) vs. d (distance between pixels in image)
if(plt_flg)
    figure, plot(d,full(B),'b','LineWidth',2)
    title('Gradient Box Pixel Intensity Distance Dependence','FontSize',20,'FontWeight','Bold')
    xlabel('Distance (in Pixels) between 2 Pixels','FontSize',18,'FontWeight','Bold')
    ylabel('Average Weight Value for all Pixels in Image Separated by Distance','FontSize',18,'FontWeight','Bold')
    
    disp('Make Sense?')
    keyboard
    
end

