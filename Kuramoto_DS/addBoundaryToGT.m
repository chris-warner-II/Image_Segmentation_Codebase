[dirPre,sizeGoodIm] = onCluster;
imgsDir = [dirPre,'images/BSDS_patch/101x101_ds1/'];

files = dir([imgsDir,'*.mat']);


for F = 1:numel(files)

    
    load([imgsDir,files(F).name])
    
    % If the boundary information is not in the data structure, make it and save it.
    %if ~exist('MC','var')
        
        disp(['Adding boundary ground truth to image patch info. ',num2str(F),' / ',num2str(numel(files))])
        
        % Boundary Ground Truth on Full Image
        bDfull = gTfull{1}.Boundaries;
        for i = 2:numel(gTfull)
            bDfull = bDfull + gTfull{i}.Boundaries;
        end

        % Boundary Ground Truth on Image Patch
        bD = gTfull{1}.Boundaries(pach.ypbeg:pach.ypfin,pach.xpbeg:pach.xpfin);
        for i = 2:numel(gTfull)
            bD = bD + gTfull{i}.Boundaries(pach.ypbeg:pach.ypfin,pach.xpbeg:pach.xpfin);
        end


        
        % Blur ground truth boundary lines for analysis
        MC.bDc_blurs_info = {'(bD>0)','(bD>0)blur_binom2','(bD>0)blur_binom3','(bD>1)','(bD>1)blur_binom2','(bD>1)blur_binom3','bD'};
        bDc_blurs(:,:,1) = single(bD>0);
        bDc_blurs(:,:,2) = single(logical( blur( single(bD>0),1,'binom2' ) ) );
        bDc_blurs(:,:,3) = single(logical( blur( single(bD>0),1,'binom3' ) ) );
        bDc_blurs(:,:,4) = single(bD>1);
        bDc_blurs(:,:,5) = single(logical( blur( single(bD>1),1,'binom2' ) ) );
        bDc_blurs(:,:,6) = single(logical( blur( single(bD>1),1,'binom3' ) ) ); % note: 'power' increases with increased blurring
        bDc_blurs(:,:,7) = bD;

        
        % Compute boundary based Gradient Metric using image pixels
        [MC.F,MC.M,MC.S,MC.D] = compute_BoundaryGradientMetric(im,bDc_blurs,numel(gT),0,'Im Pix');
        [MC.Fp,MC.Mp,MC.Sp,MC.Dp] = compute_BoundaryGradientMetric(pi.*im,bDc_blurs,numel(gT),0,'Im Pix');
        
        save([imgsDir,files(F).name],'yimg', 'ximg', 'yFimg', 'xFimg', 'ds_fctr', 'gT', 'gTfull', 'im', 'imFull', 'bD', 'bDfull', 'pach','MC') 

    %end

    
    clear MC


    % Make some plots
    if(0)
        figure,
        subplot(231), imagesc(im), axis square off, colormap('bone'), title('Image Patch')
        subplot(232), imagesc(bD), axis square off, colorbar, title('Boundaries')
        subplot(233), imagesc(bD>1), axis square off, title('Consensus')
        %
        subplot(234), imagesc(imFull), axis square off, colormap('bone'), title('Full Image')
        subplot(235), imagesc(bDfull), axis square off, colorbar, title('Boundaries')
        subplot(236), imagesc(bDfull>1), axis square off, title('Consensus')

        pause
    end
    
    
    
    %
    
    

end



