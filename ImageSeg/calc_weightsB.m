function [W, Wdist, Mask] = calc_weightsB(imagin,sigpix,sigdist,r_range,maskFlg,topo)

% syntax: [W, Wdist, Mask] = calc_weightsB(imagin,sigpix,sigdist,r_range,maskFlg,topo)
%
% This function takes in an input image and some parameters and calculates
% and outputs an Adjancency Matrix for pixels within that image. 
%
% INPUT: (imagin,sigpix,sigdist,maskFlg,rmax,maskFlg,topo)
%
%    imagin  - N x N input image is greyscale with pixel values (Z) between 0 (black) and 1 (white).
%
%    sigpix  - sigma value on gaussian function that converts pixel difference abs(Zi - Zj) into edge weight contribution between 0 and 1
%
%    sigdist - sigma value on gaussian function that converts pixel distance in the image dist(i,j) into edge weight contribution between 0 and 1
%              [ NOTE: We use sigdist=inf so the distance dependence is binary. Wdist=1 inside rmax and Wdist=0 outside. ]
%              [       This may be important for how we compute null models. Try it as a first pass. Set sigdist=inf.    ]
%
%    r_range - either [rmin, rmax] or just rmax (then rmin defaults to 0)
%    where rmax is maximum distance across which pixels/nodes in network can be connected. Implemented by heaviside/step function.
%
%    maskFlg - a flag to choose which mask to compute. (0=No Mask. 1=Mask of pairwise pixel distances in image. 2=Not Working/Dont Use.)
%
%    topo    - a flag for network topography (0=none, 1=1D linear, 2=2D lattice, 3=2.5D dist & ori.)
%
% OUTPUT: [W, Wdist, Mask]
%
%   W     - N^2 x N^2 Adjacency Matrix using both distance-dependent connectivity constraints and pixel differences.
%
%   Wdist - N^2 x N^2 "Hardware" Connectivity constraints imposed by programmer - Distance dependence of network.
%
%   Mask  - structure which contains a mask of all pixel pairs in adjacency separated by a each distance less than rmax.
%            [ Mask(i).distance = the distance (scalar value).                                                           ].
%            [ Mask(i).mask     = the binary mask (is N^2 x N^2 like W matrix).                                          ].

disp(['In calc_weightsB : Constructing network from image pixel values xxx:'])



% convert r_range into rmin & rmax, whether it is 1 number of 2 numbers.
if( numel(r_range)==2 )
    rmin = r_range(1);
    rmax = r_range(2);
elseif( numel(r_range)==1 )
    rmin = 0;
    rmax = r_range;
else
    disp('r_range variable does not make sense.')
    keyboard
end




%% Preallocate Sparse Weights Matrices (will be mostly zeros because of nearest neighbor connectivity).
%
if isinf(rmax)                      % nNeighbors is the number of direct couplings that each oscillator will have.
    nNeighbors = numel(imagin);
else
    nNeighbors = sum(2*(1:rmax-rmin)*4); % this is an overestimate right now (takes everything in a box)
    % pctNN = max(sum(logical(W_conn)))/(ximg*yimg-1);   % ratio of number of connections by number of possible connections
	% ceil(pctNN*numel(W_conn))
    % TODO: For more efficient code, we can pare it down from 80 to 48 for Rmax = 4.
end
%
Wdist = spalloc( numel(imagin) , numel(imagin), nNeighbors*numel(imagin) );
W = spalloc( numel(imagin) , numel(imagin), nNeighbors*numel(imagin) );

N = numel(imagin);



%% Compute adjacency matrix in batches. Can handle larger image patches.
%
%  A smarter way to do this - to construct Wdist - without using large, full matrices.
%  Nested for loops may run slower, but they should not run into memory issues.


disp('New / Nested For Loop Way') % can handle larger size networks!
tic


% (1). Pixel Locations - to compute Wdist - vectorized to speed up computation.
y_vec = [1:size(imagin,1)]';                                 % (y = 1 = i)
y_mat = repmat(y_vec,1,size(imagin,2));
%
x_vec = [1:size(imagin,2)];                                  % (x = 2 = j)
x_mat = repmat(x_vec,size(imagin,1),1);

% (2). Matrix of distances from pixel (0,0) to each pixel in the image.
rB2 = sqrt( (y_mat-0).^2 + (x_mat-0).^2 );
hood2 = find(unique(rB2)<=rmax & unique(rB2)>rmin); % how many discrete distance values in the grid layout between 0 and rmax.
% hood3 = find(unique(rB2));                     % how many discrete distance values in the 2D grid (ignore rmax).

% (3). Preallocate memory for Mask Structure that contains a binary masks that overlays the adjacency matrix W and indicate all 
%      entries/edges that connect pairs of oscillators separated by a given distance, for each distance between 0 and rmax. 
%
%      Will populate it in the for loops below if the network has topographic organization.
if(topo)
        
    switch(maskFlg)

        case(0) % Case 1: SKH_Adj or NG or AA or GL, working

            % DO 1D Line Topography.
            % SKH_Adj dist WORKING (This would be used for 1D Topographic Modularity!)
            for k = (rmin+1):(rmax+1)
                Mask(k).distance = 0;
                Mask(k).mask = spalloc( N , N, nNeighbors*N ); 
            end
            
        case(1) % Case 2: SKH_Euc or NG or AA or GL, working

            % SKH_Euc dist WORKING (This is used for 2D Topographic Modularity!)
            for k = 1:numel(hood2) % could be hood3 if u want all distances (uncomment above).
                Mask(k).distance = 0;
                Mask(k).mask = spalloc( N , N , nNeighbors*N ); 
            end
            
        case(2) % Case 3: SKH_EUC D&O NOT WORKING YET (IGNORE FOR NOW !)

            % DO 2.5D Lattice Topography (Direction & Orientation)
            for k = 1:numel(hood2) % could be hood3 if u want all distances (uncomment above).
                Mask(k).distance = 0;
                Mask(k).mask = spalloc( N , N , nNeighbors*N ); 
            end
            
    end
    
else            % Case 4: Modularity N&G - Non-topographic.

    Mask(1).distance = 0;
    Mask(1).mask = spalloc( N , N , nNeighbors*N ); 
    
end



% Loop through ... WORKING HERE TO IMPLEMENT RMIN!!
disp(['Iteration # out of ',num2str(size(imagin,1))])
for i = 1:size(imagin,1)     % step through y dimension of image for pixel location #1.
    
    fprintf('%s',[num2str(i),' '])
    
    for j = 1:size(imagin,2) % step through x dimension of image for pixel location #2.
        
        %fprintf('%s',[num2str(j),' '])
        
        rB2 = sqrt( (y_mat-i).^2 + (x_mat-j).^2 ); % distance of all pixels from pixel (i,j)
        hood2 = find(rB2<=rmax & rB2>0);           % all pixels that are between 0 and rmax distance from pixel (i,j)
        
        % This part can produce plots to show what pixels are in neighborhood for error checking.
        if(0)
            h=figure; imagesc(rB2), hold on
            [hood2y,hood2x] = find(rB2<=rmax & rB2>rmin);
            scatter(hood2x,hood2y,'wo'),  % plot neighbors of pixel j,i
            scatter(j,i,'wx'),            % plot pixel j,i
            axis([1 numel(x_vec) 1 numel(y_vec)]), axis ij
            keyboard
            close(h)
        end
        
        
        Wdist( size(imagin,1)*(j-1)+i, hood2 ) = exp( -( rB2(hood2) - 1).^2 ./ (2*sigdist.^2) );
            W( size(imagin,1)*(j-1)+i, hood2 ) = exp( -( rB2(hood2) - 1).^2 ./ (2*sigdist.^2) ) .* ...
                                                 exp( -( ( (imagin(i,j) - imagin(hood2)).^2 ) )./ (2*sigpix^2) );
            % Note: The -1 in the Wdist gaussian exp( -(r-1)^2 / (2*s^2) ) is not standard for a gaussian.  I use it to impose that
            % pixels separated by a distance r=1 have a Wdist component of their weight equal to 1, rather than exp( -1/(2s^2) ).
            
            
            
        % Create a Series of Masking Matrices that Index into weights matrix to get all pixel pairs that are 
        % separated by a given distance in the image (equivalence classes).
        % PROBABLY WANT TO REINSTATE THESE DIFFERENCES AT SOME POINT. PROCESS THEM ALL THE SAME RIGHT NOW.
        % JUST FOCUS ON CASE(1) - THAT IS TOPOGRAPHIC MODULARITY.
        if(topo)
        
            switch(maskFlg)

                case(0) % SKH_Adj or NG or AA or GL, working

                    % DO 1D Line Topography.
                    
                    
                    % i loops thru y dimension (outer for loop)
                    % j loops thru x dimension (inner for loop)
                    % % % y=diag((ones(1,N-i)),i) + diag((ones(1,N-i)),-i);
                    %
                    Mask(i).distance = (i-1);
                    Mask(i).mask = diag((ones(1,N-(i-1))),(i-1)) + diag((ones(1,N-(i-1))),-(i-1));
                    

                    % code to plot off diagonals for explaining 1D topography.
                    if(0)
                        figure, cmap = colormap(bone); cmap(1,:) = [1,1,1];
                        colors4plot = [0,0,0 ; 1,0,0 ; 0,1,0 ; 0,0,1];
                        for x = 1:4
                            y=diag((ones(1,100-x)),x) + diag((ones(1,100-x)),-x)
                            cmap(end,:) = colors4plot(x,:);
                            subplot(1,4,x), 
                            imagesc(y), colormap(cmap), title(['r=',num2str(x)],...
                            'fontsize',20,'fontweight','bold')
                            axis square
                            set(gca,'Xtick',[],'Ytick',[])
                            freezeColors
                        end
                    end


                    
                    
                    
                    
                    
                case(1) % SKH_Euc dist WORKING (This is used for Topographic Modularity!)

                    % DO 2D Lattice Topography.
                    rd = unique(rB2);             % all possible discretized distances between pixel (i,j) and other pixels in the 2D grid of image.
                    rd_ind = find( rd<=rmax  & rd>=rmin );  % all discretized distances less than rmax
                    rd = rd(rd_ind);

                    for k=1:numel(rd)             % for each distance < rmax, fill in an the entries of that mask pertaining to the pixel pairs  
                        %
                        if(i==1 & j==1)
                            Mask(k).mask = zeros(N);
                        end
                        %
                        hood3 = find(rB2==rd(k)); % where one is pixel (i,j) and the other is another pixel a distance d from pixel (i,j) .
                        Mask(k).distance = rd(k);
                        Mask(k).mask( size(imagin,1)*(j-1)+i, hood3 ) = 1;
                    end

                    
                    
                    
                    
                    
                case(2) % SKH_EUC D&O NOT WORKING YET (IGNORE FOR NOW !)

                    % DO 2.5D Lattice Topography (Direction & Orientation)
                    xd = unique(x_dist);
                    yd = unique(y_dist);
                    k=0;
                    for i=1:numel(xd)
                        for j=1:numel(yd)
                            k=k+1;
                            Mask(k).distance = sqrt(xd(i).^2 + yd(j).^2);
                            Mask(k).distX = xd(i);
                            Mask(k).distY = yd(j);
                            Mask(k).mask = ( (x_dist==xd(i)) & (y_dist==yd(j)) );
                        end
                    end
            end % switch mask
            
%             % Loop through each Mask entry and make sure it is NxN.
%             for h = 1:k
%                 
%                 
%             end
            
        
        end % if topo
        
        % Note: Later in PIV_vs_MaskDist, calculate "Diagonal Means" or Distance Dependence as: B(d) = mean(mean(W(Mask(d).mask)))

    end % Loop over x-dimension of image.
    
    % figure, imagesc(y), colorbar

end % Loop over y-dimension of image.

fprintf('\n')

toc



%% Making a figure for Mod_SKHEuc Null Model Masking. 
if(0) % Set to 1 to visualize what masks look like. Rmax must be >5.

    H=figure;
    subplot(2,3,1), colormap(bone), imagesc(1-Mask(1).mask), set(gca,'XTick',[],'YTick',[]), axis square, title(['d=',num2str(Mask(1).distance)],'FontSize',20,'FontWeight','Bold')
    subplot(2,3,2), colormap(bone), imagesc(1-Mask(2).mask), set(gca,'XTick',[],'YTick',[]), axis square, title(['d=',num2str(Mask(2).distance)],'FontSize',20,'FontWeight','Bold')
    subplot(2,3,3), colormap(bone), imagesc(1-Mask(4).mask), set(gca,'XTick',[],'YTick',[]), axis square, title(['d=',num2str(Mask(4).distance)],'FontSize',20,'FontWeight','Bold')
    subplot(2,3,4), colormap(bone), imagesc(1-Mask(7).mask), set(gca,'XTick',[],'YTick',[]), axis square, title(['d=',num2str(Mask(7).distance)],'FontSize',20,'FontWeight','Bold')
    subplot(2,3,5), colormap(bone), imagesc(1-Mask(10).mask), set(gca,'XTick',[],'YTick',[]), axis square, title(['d=',num2str(Mask(10).distance)],'FontSize',20,'FontWeight','Bold')
    subplot(2,3,6), colormap(bone), imagesc(1-Mask(14).mask), set(gca,'XTick',[],'YTick',[]), axis square, title(['d=',num2str(Mask(14).distance)],'FontSize',20,'FontWeight','Bold')
    
    saveGoodImg(H,['2D_Bij_Not_Diagonal.jpg'],[0 0 1 1])
    close(H)    
        
end


disp('Memory usage in calc_weightsB : ')
% whos
check_memory_usage

for k = 1:numel(Mask)
    if ~all(size(Mask(k).mask)==[N,N])
        keyboard
    end
end



