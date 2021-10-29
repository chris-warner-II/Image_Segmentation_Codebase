function [ NM ] = compute_null_modelB( W, Wdist, Mask, imagin, topo, maskFlg )

% syntax: [ NM ] = compute_null_modelB(W,Wdist,Mask,topo,maskFlg)
%
% This function computes the null model matrix(NM) for a network given that
% network's adjacency matrix (W) and prior knowledge about the arrangement
% of nodes in the network. The null model is a simplified, or homogenized,
% or randomized network that maintains some statistics (closely) with the
% Adjacency matrix that it was built with.
%
% INPUT: (W,Wdist,Mask,topo,maskFlg)
%   W       - N^2 x N^2 Adjacency Matrix using both distance-dependent connectivity constraints and pixel differences.
%   Wdist   - N^2 x N^2 "Hardware" Connectivity constraints imposed by programmer - Distance dependence of network.
%   Mask    - structure which contains a mask of all pixel pairs in adjacency separated by a each distance less than rmax.
%            [ Mask(i).distance = the distance (scalar value).                                                           ].
%            [ Mask(i).mask     = the binary mask (is N^2 x N^2 like W matrix).                                          ].
% imagin - input image to compute MaskFull when   
%   topo    - flag whether to use topographic null model (di*dj*bij) or nontopographic null model (di*dj).
%             [0=original nontopographic Newman&Girvan, 1=topographic]
%   maskFlg - a flag to choose which mask to compute. (0=No Mask. 1=Mask of pairwise pixel distances in image. 2=Not Working/Dont Use.)
%
%
% OUTPUT: [NM]
%   NM      - N^2 x N^2 matrix which contains expectation of edge-weight in a "homogenized" (or otherwise random) 
%             network that maintains some statistics measured from the adjacency matrix. The nontopographic null 
%             model (topo=0) is constructed based on node degree (D). The topographic null model (topo=1) is 
%             constructed based on node degree and average edge-weight between nodes separated by a given distance (B).
%



batchSize = 500;

imPix = size(W,1); %total number of pixels in image.

if(topo) % TOPOGRAPHICAL null model
    
    pctNN = max(sum(logical(Wdist)))/(imPix-1);   % ratio of number of connections by number of possible connections
    NM = spalloc(imPix,imPix,ceil(pctNN*numel(W))); % Preallocate sparse matrix for topographical null model.
    
    % Compute edge-weight dependence on distance.
    switch(maskFlg)
        case(0) % Use only Diagonal Means
            tag='Mod topo 1D';
            [B] = PIV_vs_DiagDistB(W,Wdist);
            tag2 = 'DiagD';
            
        case(1)  % Use Euclidian Distance in Image vs PIV means with Mask
            tag='Mod topo 2D';
            [B] = PIV_vs_MaskDistB(W,Wdist,Mask);
            tag2 = 'MaskD';

        case(2) % Use different masks for different distance & orientation (vert vs horiz) 
            tag='SKH TopographicalDO';
            [B] = PIV_vs_MaskDistB(W,Wdist,Mask);
            tag2 = 'MaskDO';
    end
    
    % Compute row sums or node degree (Here, percent degree because dividing by Wdist)
    D = sum(W) ./ sum(Wdist);
    
    disp(['Back out in compute_Null_ModelB: Putting final Null Model Matrix Together'])
    disp('I expect this to run slower. Current Memory Usage:')
    % Took ~5 hours for 154k x 154k matrix which is the full 481x321 image.
    
%     whos
    check_memory_usage
    
    tic
    % Piece Null Model together from D,B & Wdist in batches of 1000 x 154k (like in KuramotoB)
    batchBeg = [1:batchSize:imPix];
    batchEnd = [batchSize:batchSize:imPix];
    if ( numel(batchEnd) < numel(batchBeg) )
        batchEnd = [batchEnd,imPix];
    end
    %
    disp(['Back out in compute_Null_ModelB: Putting final Null Model Matrix Together NM = D*D*B*W_conn'])
    disp(['Batch Iteration # out of ',num2str(numel(batchBeg))])
    for i = 1:numel(batchBeg)
        
        fprintf('%s',[num2str(i),' '])
        
        st = batchBeg(i);
        nd = batchEnd(i);
        num = nd-st+1;

        % Topographic Null Model is just NM_{ij} = D_i*D_j*B_{ij} because D & B are normalized by Wdist.  
        NM(st:nd,:) = ( repmat(D(st:nd)',1,imPix).*repmat(D,num,1) ).*B(st:nd,:); 
        
    end
    
    % Normalize the Null Model properly so its overall weight matches overall weight of Adjacency matrix
    c = sum(W(:))./sum(NM(:));
    NM = c.*NM;

    fprintf('\n')
    
    toc
    
else % NONTOPOGRAPHIC Null Model       

    NM = zeros(imPix,imPix); % preallocate Nontopographical Null Model Matrix
    
    % Compute row sums or node degree
    D = sum(W);

    disp('NOTE: Non-topographic Null Model Is a Dense Matrix')
    tag='Mod NG';
    tag2 = '';
    
    % NEW WAY - USING BATCHES TO HANDLE BIG MATRICES.
    disp('New Way for big matrices using batches.')
    tic
    % Piece Null Model together working in batches to get whole large matrix.
    batchBeg = [1:batchSize:imPix];
    batchEnd = [batchSize:batchSize:imPix];
    if ( numel(batchEnd) < numel(batchBeg) )
        batchEnd = [batchEnd,imPix];
    end
    %
    disp(['Back out in compute_Null_ModelB: Putting final Null Model Matrix Together NM = D*D*B*W_conn'])
    disp(['Batch Iteration # out of ',num2str(numel(batchBeg))])
    for i = 1:numel(batchBeg)
        
        fprintf('%s',[num2str(i),' '])
        
        st = batchBeg(i);
        nd = batchEnd(i);
        num = nd-st+1;

        NM(st:nd,:) = ( repmat(D(st:nd)',1,imPix).*repmat(D,num,1) ); % NM_ij = D_i*D_j
    end
    
    NM = NM ./ sum(W(:)); % normalize by total amount of edge-weights in Adjacency matrix
    
    fprintf('\n')
    
    toc
    
end




    


%% Null Model Consistency.  Featured in the Motivate Modularity writeup. 
%  For NG, SKHAdj and SKHEuc null models, shows how well each one matches 
%  degree distributions and distance dependence for a real image patch.
if(1)
    

    %     NOTE: HAVE TO RECALCULATE MASK FOR ALL PIXEL DISTANCES IN THE IMAGE. TO
    %     CONSTRUCT THIS FIGURE FOR THE PAPER COMPARING EACH MODULARITY DISTANCE
    %     AND DEGREE DEPENDENCE FOR (N&G, 1D-TM and 2D-TM).
    N = size(imagin,1); % assumes square input image patch.
    
    % (1). Pixel Locations - to compute Wdist - vectorized to speed up computation.
    y_vec = [1:size(imagin,1)]';                                 % (y = 1 = i)
    y_mat = repmat(y_vec,1,size(imagin,2));
    %
    x_vec = [1:size(imagin,2)];                                  % (x = 2 = j)
    x_mat = repmat(x_vec,size(imagin,1),1);

    % (2). Matrix of distances from pixel (0,0) to each pixel in the image.
    rB2 = sqrt( (y_mat-0).^2 + (x_mat-0).^2 );
    hood2 = find(unique(rB2));                     % how many discrete distance values in the 2D grid (ignore rmax).
    %
    for k = 1:numel(hood2) % Preallocate memory to construct MaskFull using all distances all distances 
        MaskFull(k).distance = 0;
        MaskFull(k).mask = spalloc( N , N , numel(imagin)*N ); 
    end
    %
    % Fill in MaskFull to plot below.
    disp(['Iteration # out of ',num2str(size(imagin,1))])
    for i = 1:size(imagin,1)     % step through y dimension of image for pixel location #1.
        fprintf('%s',[num2str(i),' '])

        for j = 1:size(imagin,2) % step through x dimension of image for pixel location #2.
            %fprintf('%s',[num2str(j),' '])

            rB2 = sqrt( (y_mat-i).^2 + (x_mat-j).^2 ); % distance of all pixels from pixel (i,j)

            % 2D Lattice Topography. - SKH_Euc dist WORKING (This is used for Topographic Modularity!)
            rB2 = sqrt( (y_mat-i).^2 + (x_mat-j).^2 ); % distance of all pixels from pixel (i,j)
            rd = unique(rB2);             % all possible discretized distances between pixel (i,j) and other pixels in the 2D grid of image.
            for k=1:numel(rd)             % for each distance < rmax, fill in an the entries of that mask pertaining to the pixel pairs  
                hood3 = find(rB2==rd(k)); % where one is pixel (i,j) and the other is another pixel a distance d from pixel (i,j) .
                MaskFull(k).distance = rd(k);
                MaskFull(k).mask( size(imagin,1)*(j-1)+i, hood3 ) = 1;
            end

        end % Loop over x-dimension of image.

    end % Loop over y-dimension of image.
   

    distW = zeros(1,numel(MaskFull));
    distNM = zeros(1,numel(Mask));
    
%     if(~topo) % Mod_N&G
%         
%         % HOW I WAS DOING IT.  REINSTATE LATER. DOING THIS NOW FOR PLOTTING BELOW        
%         for i = 1:size(W)
%            distW(i) = mean(diag(W,i));
%            distNM(i) = mean(diag(NM,i));
%         end
% 
%     elseif(topo & ~maskFlg) % Mod_SKHAdj
% 
%         % HOW I WAS DOING IT.  REINSTATE LATER. DOING THIS NOW FOR PLOTTING BELOW
%         for i = 1:size(W)
%            distW(i) = mean(diag(W,i));
%            distNM(i) = mean(diag(NM,i));
%         end
%        
%     elseif(topo & maskFlg==1) % Mod_SKHEuc
% 
%         for i = 1:numel(Mask)
%            distW(i) = mean(mean(W.*Mask(i).mask));
%            distNM(i) = mean(mean(NM.*Mask(i).mask));
%         end
%         
%         
%     else % Should not get here.  Problem if so.
%         
%         disp('Should not get here.  Problem if so. Inside compute null modelB')
%         keyboard
%         
%     end


    for i = 1:numel(MaskFull)
        ind = find(MaskFull(i).mask);
        distW(i) = mean(W(ind));
        distNM(i) = mean(NM(ind));
    end


    % Plot Consistency - How well degree and distance dependence of null model matches adjacency matrix.
    H=figure; 

    subplot(321), imagesc(W),colorbar, title(['Adjacency ',tag2],'FontSize',20,'FontWeight','Bold'), axis square
    set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','Bold'), 
    xlabel(['$\sum_{ij} A_{ij} = $',num2str(sum(W(:)),5)],'FontSize',18,'FontWeight','Bold','Interpreter','LaTex')

    subplot(322), imagesc(NM),colorbar, title([tag,' Null Model ',tag2],'FontSize',20,'FontWeight','Bold'), axis square
    set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','Bold'), 
    xlabel(['$\sum_{ij} N_{ij} = $',num2str(sum(NM(:)),5)],'FontSize',18,'FontWeight','Bold','Interpreter','LaTex')

    x = mean(W); % sort nodes based on degree and plot sorted nodes (easier to see)
    y = mean(NM);
    [X,I] = sort(x);
    Y = y(I);
    %
    subplot(312),hold on, plot(Y,'b','LineWidth',2), plot(X,'r--','LineWidth',2), 
    xlabel(['Sorted Node #'],'FontSize',18,'FontWeight','Bold')
    ylabel(['node degree'],'FontSize',18,'FontWeight','Bold')
    legend('Null Model','Adjacency','Location','SouthEast')
    set(gca,'FontSize',16,'FontWeight','Bold')

%     distMax = sqrt(2*imPix); % cludge: NOTE THIS ASSUMES SQUARE IMAGE PATCH AND IS WRONG IF RECTANGULAR.
%                                             % pythagorean rule: dist diag  across image C = sqrt(x^2 + y^2)
    
    subplot(313),hold on, 
    plot([MaskFull.distance],distNM,'b','LineWidth',2), 
    plot([MaskFull.distance],distW,'r--','LineWidth',2),
    xlabel(['Distance (in pixels) between nodes'],'FontSize',18,'FontWeight','Bold')
    ylabel(['avg. edge-weight'],'FontSize',18,'FontWeight','Bold')
    legend('Null Model','Adjacency','Location','SouthEast')
    set(gca,'FontSize',16,'FontWeight','Bold')

    tagB = tag;
    tagB(tagB==' ')='_';
    
    
    save( ['./Consistency_AvsNM_',tagB,'.mat'], 'imagin', 'W', 'NM', 'distW', 'distNM', 'MaskFull' )

    saveGoodImg(H,['./Consistency_AvsNM_',tagB,'.jpg'],[0 0 1 1])
    
    disp('xxx')
    
end

