function SegmentMethod(InputImages,imgNums2proc,imdim,AvgAssoc,GraphLap,normlze,... 
                        Modularity,maxEnt,topo,maskFlg,NormCut,meanThresh,IsoDiff,...
                        shiftEvals,maxMEiter,sigpix,sigdist,Rmax,save_Evecs_Mat,dataFileChk,...
                        DoKurProcessing,DoEigProcessing,SigW) % Kscale, SigW, TiScale)
                    
                    
                   
                    
                    
     % rng(1234567); % Random Number Generator Seed.  Rand used in Ising Model.


     
     % NOTE: Kscale & TiScale that are input are not used anymore.  THey
     % are defined below inside. Kscale is computed from rowsums and # of
     % neighbors M. TiScale is set to be 1.
     

    %% Make a directories to put output images and data in
    [dirPre,sizeGoodIm] = onCluster;



    %% Strings to tell user what specific method is was used in calculation
    str = '';

    if(shiftEvals) 
        str=[str,' Shft'];
    end
    if(normlze) 
        str=[str,' Norm'];
    end
    str = [str,' '];


    if GraphLap
        neg = 1;        % flag {0,1} to use negative matrix caluclated by method (Q -> -Q)
        proj = 1;       % flag {0,1} to project out the Largest Eigenvector (for Power Method)
    else
        neg = 0;        
        proj = 0;
    end



    blur_image = 0; % flag whether to blur image with a gaussian filter kernel before doing network computation.
    DoG_image = 0;  % flag to use difference of gaussian filters to "blur" image before network computaiton.

    doPlots = 0; % make a single plot for each:  eigenvectors & kuramoto
    doMCA   = 1; % do metaClustering analysis (both for Eig & Kur)

    % compareCircVsLin = 1; % flag to directly compare linear vs circular divisive margin on strawman & Evec1.

    bracket_Kscale = 1;



    %% Plot distance dependence gaussian and pixel itensity difference gaussians (for development & troubleshooting)
    if(0) 
        h=figure; hold on,
        colour = ['rgkcmbrgkcmbrgkcmbrgkcmbrgkcmbrgkcmbrgkcmbrgkcmb'];

        % (1). plot Weight vs. Distance between pixels
        for r = 1:numel(Rmax)
            rad = Rmax(r);
            Xd = linspace(1,rad);
            subplot(211), hold on
            for i=1:numel(sigdist)
                Wd = exp( -( Xd - 1 ) ./ (2*sigdist(i).^2) );
        %         Wd(Xd>Rmax)=0; % pixels beyond maximum radius have zero weight.
                plot(Xd,Wd,colour(i),'LineWidth',2)
                leg{i} = ['\sigma = ',num2str(sigdist(i))];
            end  
        end
        legend(leg)
        xlabel('Pixel Separation, Euclidian Distance','FontSize',16,'FontWeight','Bold')
        ylabel('Weight Strength','FontSize',16,'FontWeight','Bold')
        title('Pixel Distance','FontSize',20,'FontWeight','Bold')
        text(4,0.8,['rmax = ',num2str(Rmax)],'FontSize',16,'FontWeight','Bold')
        grid on

        % (2). plot Weight vs. Difference in pixel value
        Xd2 = linspace(0,1);
        subplot(212),hold on
        for i=1:numel(sigpix)
            Wd2 = exp( -( Xd2 ) ./ (2*sigpix(i).^2) );
            plot(Xd2,Wd2,colour(i),'LineWidth',2)
            leg2{i} = ['sigpix = ',num2str(sigpix(i))];
        end
        legend(leg2)
        xlabel('Difference in Pixel Intensity Values','FontSize',16,'FontWeight','Bold')
        ylabel('Weight Strength','FontSize',16,'FontWeight','Bold')
        title('Pixel Similarity','FontSize',20,'FontWeight','Bold')
        grid on
        % set(h,'Position',[1 1 1280 700]);
        % saveas(h,['./output/ImgSeg/simpleExamples/',dirIn,'/Adjacency_Matrix_Gaussians'],'jpg');

        % (3). plot a color image in 2D of distance vs. pixel diff vs. weight
        figure, imagesc(Xd2,Xd,Wd'*Wd2),
        %
        hcb=colorbar; 
        colorTitleHandle = get(hcb,'Title');
        titleString = 'Graph Weight';
        set(colorTitleHandle ,'String',titleString,'FontSize',16,'FontWeight','Bold');
        %
        ylabel('Distance between Pixels','FontSize',16,'FontWeight','Bold')
        xlabel('Pixel Value Difference','FontSize',16,'FontWeight','Bold')
        title(['Weight Strength - - - (\sigma_p=',num2str(sigpix(end),2),' , \sigma_d=',num2str(sigdist(end),2),')'],'FontSize',20,'FontWeight','Bold')

        keyboard

    end







    %% Build up method string based on input flags to this function.
    if(GraphLap)
        method = ['GL'];
    end
    %
    if(Modularity)
        method = ['Mod '];
        if(maxEnt)
            if(topo) % if u want to use full (not sparse) MaxEnt null model
                method = [method,'SKH ME']; % ,num2str(NMiter),'iter'
            else
                method = [method,'N&G ME']; % ,num2str(NMiter),'iter'
            end
        else
            if(topo)
                method = [method,'SKH'];
                switch(maskFlg)
                case(0)
                    method = [method,'Adj']; % diagonal in Adjacency
                case(1)
                    method = [method,'Euc']; % mask Euclidian dist
                case(2)
                    method = [method,'D&O']; % mask distance & orientation
                end
            else
                method = [method,'N&G'];
            end
        end
    end
    %
    if(AvgAssoc)
        method = ['AA'];
    end
    %
    if(meanThresh)
        method = ['THmean'];
    end
    %
    if(normlze)
        method = [method,'nrm'];
    end

    if(IsoDiff)
        method = ['IsoDiff'];
    end


    method(method==' ')='_';
    ds_fctr = 1; % hard coding this right now.  It will be replaced if u need to downsample..




    % Make directories to put mat files in with results from Eigenvector segmentation for each method dont want spaces in directory names.
    % NOTE : DOING THIS NOW BELOW INSIDE FOR LOOPS WITH DIRMATKUR DIRMATEIG DIRIMGKUR & DIRIMGEIG
%     if ~exist([dirPre,'output/Kuramoto/NetsFromImgs/',InputImages,'_',num2str(imdim(1)),'x',num2str(imdim(2)),'_ds',num2str(ds_fctr),'/data/spectral/',method],'dir')
%         mkdir([dirPre,'output/Kuramoto/NetsFromImgs/',InputImages,'_',num2str(imdim(1)),'x',num2str(imdim(2)),'_ds',num2str(ds_fctr),'/data/spectral/',method]);
%     end
%     %
%     if ~exist([dirPre,'output/Kuramoto/NetsFromImgs/',InputImages,'_',num2str(imdim(1)),'x',num2str(imdim(2)),'_ds',num2str(ds_fctr),'/imgs/spectral/',method],'dir')
%         mkdir([dirPre,'output/Kuramoto/NetsFromImgs/',InputImages,'_',num2str(imdim(1)),'x',num2str(imdim(2)),'_ds',num2str(ds_fctr),'/imgs/spectral/',method]);
%     end
%     %
%     if ~exist([dirPre,'output/Kuramoto/NetsFromImgs/',InputImages,'_',num2str(imdim(1)),'x',num2str(imdim(2)),'_ds',num2str(ds_fctr),'/data/Kur_PIF_Fourier1/',method],'dir')
%         mkdir([dirPre,'output/Kuramoto/NetsFromImgs/',InputImages,'_',num2str(imdim(1)),'x',num2str(imdim(2)),'_ds',num2str(ds_fctr),'/data/Kur_PIF_Fourier1/',method]);
%     end
%     %
%     if ~exist([dirPre,'output/Kuramoto/NetsFromImgs/',InputImages,'_',num2str(imdim(1)),'x',num2str(imdim(2)),'_ds',num2str(ds_fctr),'/imgs/Kur_PIF_Fourier1/',method],'dir')
%         mkdir([dirPre,'output/Kuramoto/NetsFromImgs/',InputImages,'_',num2str(imdim(1)),'x',num2str(imdim(2)),'_ds',num2str(ds_fctr),'/imgs/Kur_PIF_Fourier1/',method]);
%     end
        





    %% Loop thru preprocessed image MAT files in a directory...
    if ~ischar(imdim)
        dataDir = [dirPre,'images/',InputImages,'/',num2str(imdim(1)),'x',num2str(imdim(2)),'_ds',num2str(ds_fctr),'/'];
    else
        dataDir = [dirPre,'images/',InputImages,'/ds',num2str(ds_fctr),'/'];
    end


    files = dir([dataDir,'*.mat']);

    if ~exist('imgNums2proc','var')
        imgNums2proc = [2000,4000]; % [1,numel(files)]
    end


    %matlabpool()
    
    % Set up gaussian kernel to blur image if blur_image flag==1.
    if(blur_image)
        sigB = 1; % sigma std on optimal gaussian blur filter.
        [kern] = construct_gaussian_kernel(sigB);
    elseif(DoG_image)
        sigC = 3.4; % 1; %
        sigS = 16; % 8; %
        % Krat=0.01; 
        Krat = (sigC/sigS)^2; % this will make the integrated volume under center and surround gaussians equal.
        [kern] = construct_DoG_kernel(sigC,sigS,Krat);
    else
        sigB = 0; % this means no blurring was applied. (Raw Pixels)
    end




    for iii = imgNums2proc(1):imgNums2proc(2)  % <-- for BSDS img patches now. % 1:numel(files) % loop over files in directory
        iid = files(iii);
        fname = iid.name;
        disp(['Image # ',num2str(iii),'/',num2str(numel(files)),' - ',fname])

        load([dataDir,iid.name]);


        if(blur_image | DoG_image) 
            im = imfilter(im,kern,'symmetric');            
        end
        
        
        
        if(0)
            explore_BSDS_patches % look at patch.
            keyboard
        end


        %% Do Network Computation (eigenvector or kuramoto on GL,AA,Mod,etc) to Segment Input Image
        for jjj = 1:numel(sigpix) % loop through sigma pix values - fall off for pixel similarity dependence
        for kkk = 1:numel(sigdist) % loop through sigma dist values - fall off for distance dependence
        for lll = 1:numel(Rmax) % loop through values of maximum radius for distance dependence



            sigP = sigpix(jjj);
            sigD = sigdist(kkk);
            rmax = Rmax(lll);

            kurParams_setup

            % make appropriate directories for saving if they do not exist already 
            sD = num2str(sigD);    sD(sD=='.')='p';
            sP = num2str(sigP);    sP(sP=='.')='p'; % turning sigma and rmax parameters into strings for filenames
            rM = num2str(rmax);    rM(rM=='.')='p';
            sW = num2str(SigW);    sW(sW=='.')='p';
            if(blur_image)
                sB = num2str(sigB);    sB(sB=='.')='p';
            elseif(DoG_image)
                sB = ['C',num2str(sigC),'_S',num2str(sigS),'_Kr',num2str(Krat)];    sB(sB=='.')='p';
            else
                sB = num2str(sigB);    sB(sB=='.')='p';
            end
            
            
            
            if(strcmp(method,'IsoDiff'))
                paramsK = ['rM',rM,'_NF_60_',sW];  % for Kuramoto % ,'_kscale',Ks,'_tscale',Ts,'_runs',num2str(kurParams.runs)
                paramsE = ['rM',rM];  % for Eigen
            else
                paramsK = ['rM',rM,'_sD',sD,'_sP',sP,'_NF_60_',sW];  % for Kuramoto % ,'_kscale',Ks,'_tscale',Ts,'_runs',num2str(kurParams.runs)
                paramsE = ['rM',rM,'_sD',sD,'_sP',sP];  % for Eigen
            end





            % DO DATA FILE CHECK HERE, SO I DONT HAVE TO COMPUTE NETWORK IF FILE ALREADY EXISTS!
            disp('Checking before doing anything if Eig and Kur files already exist.  If so, break out.')
            
            if(blur_image || DoG_image)
                dirMatEig = [dirPre,'output/Kuramoto/NetsFromImgs/',InputImages,'_',num2str(ximg),'x',num2str(yimg),'_ds',num2str(ds_fctr),'_blur_sig',sB,'/data/spectral/',method,'/'];
                dirMatKur = [dirPre,'output/Kuramoto/NetsFromImgs/',InputImages,'_',num2str(ximg),'x',num2str(yimg),'_ds',num2str(ds_fctr),'_blur_sig',sB,'/data/Kur_PIF_Fourier1/',method,'/'];
                %
                dirImgEig = [dirPre,'output/Kuramoto/NetsFromImgs/',InputImages,'_',num2str(ximg),'x',num2str(yimg),'_ds',num2str(ds_fctr),'_blur_sig',sB,'/data/spectral/',method,'/'];
                dirImgKur = [dirPre,'output/Kuramoto/NetsFromImgs/',InputImages,'_',num2str(ximg),'x',num2str(yimg),'_ds',num2str(ds_fctr),'_blur_sig',sB,'/data/Kur_PIF_Fourier1/',method,'/'];
            else
                dirMatEig = [dirPre,'output/Kuramoto/NetsFromImgs/',InputImages,'_',num2str(ximg),'x',num2str(yimg),'_ds',num2str(ds_fctr),'/data/spectral/',method,'/'];
                dirMatKur = [dirPre,'output/Kuramoto/NetsFromImgs/',InputImages,'_',num2str(ximg),'x',num2str(yimg),'_ds',num2str(ds_fctr),'/data/Kur_PIF_Fourier1/',method,'/'];
                %
                dirImgEig = [dirPre,'output/Kuramoto/NetsFromImgs/',InputImages,'_',num2str(ximg),'x',num2str(yimg),'_ds',num2str(ds_fctr),'/data/spectral/',method,'/'];
                dirImgKur = [dirPre,'output/Kuramoto/NetsFromImgs/',InputImages,'_',num2str(ximg),'x',num2str(yimg),'_ds',num2str(ds_fctr),'/data/Kur_PIF_Fourier1/',method,'/'];
            end
            %
            EigMatOut = [dirMatEig,'Evecs_',fname(1:end-4),'_',paramsE,'.mat'];
            KurMatOut = [dirMatKur,'KurMC_',fname(1:end-4),'_',paramsK];
            %
            if( dataFileChk & ...
                    ( ~DoEigProcessing | (DoEigProcessing & exist([EigMatOut],'file') ) ) & ...
                    ( ~DoKurProcessing | (DoKurProcessing & ... 
                    ( ~bracket_Kscale & exist([KurMatOut,'_ksmid.mat'],'file') ) | ...
                    ( bracket_Kscale & exist([KurMatOut,'_kssml.mat'],'file') & exist([KurMatOut,'_ksmid.mat'],'file') & exist([KurMatOut,'_kslrg.mat'],'file')  ) ...
             ) ) )

                tic
                disp(['Both Kur & Eig Data Files already exists for Image Patch: ',fname(1:end-4),' ',paramsK])
                disp('Not constructing network... Next Image patch')
                continue
                toc
                
            else
                
                % create directories for data & img Kur & Eig (doing it in here because BSDS_full requires you have ximg & yimg)
                if ~exist(dirImgEig,'dir')
                    mkdir(dirImgEig)
                end
                %
                if ~exist(dirImgKur,'dir')
                    mkdir(dirImgKur)
                end
                %
                if ~exist(dirMatEig,'dir')
                    mkdir(dirMatEig)
                end
                %
                if ~exist(dirMatKur,'dir')
                    mkdir(dirMatKur)
                end
                
            end



            disp(['Constructing Network From Image: ',fname(1:end-4)])
            %tic

            % Calculate weights from Pixel Intensity Values and Closeness of Pixels
    %                 disp(['Calculating Weights with sigma pix = ',num2str(sigP),'; sigma dist = ',num2str(sigD),'; rmax = ',num2str(rmax)])
            
        [W, Wdist, Mask] = calc_weightsB(im,sigP,sigD,rmax,maskFlg,topo);

            
            
            
            % Plot distributions of W and its degree (rowsums of W) to see
            % if they are Gaussian enough to justify gaussian assumptions
            % in null model calculations
            if(0)

                H=figure;
                subplot(221), imagesc(im), colormap('bone'), freezeColors, axis square, title(['Img: ',fname(1:end-4)],'FontSize',18,'FontWeight','Bold')
                subplot(222), imagesc(W), colormap('jet'), freezeColors, axis square, title('Adjacency','FontSize',18,'FontWeight','Bold')
                subplot(223), hist(W(find(W))), axis square, title('p(A_{ij})','FontSize',18,'FontWeight','Bold')
                subplot(224), hist(sum(W)), axis square, title('p(d_{i})','FontSize',18,'FontWeight','Bold')
                saveGoodImg(H,['../Documentation/Fritz_Null_Model_Statistics_WriteUps/AreDistributionsGaussian/',fname(1:end-4),'_',paramsE],sizeGoodIm)
                close(H)
                break

            end
            
           
            
            %
            % Method: Graph Laplacian (L = D - W).
            if(GraphLap)
                %figure, imagesc(im), colormap('bone')
                Q = compute_LaplacianB(W,normlze,neg);
            end

            %
            % Method: Modularity (Q = W - NM) where NM = null model or expected connectivity.
            if(Modularity)
                Q = compute_ModularityB(W,Wdist,Mask,im,topo,maskFlg,rmax);
            end

            %
            % Method: Average Association (W) - Just the similarity, adjacency, or weights matrix
            if(AvgAssoc)
                Q = compute_AvgAssociationB(W,normlze);   
            end

            %
            % Method: Normalized Cut (a la Shi & Malik)
            if(NormCut)
                disp('That code from Shi''s website isnt working yet, fix it.')
            end

            % Method: Threshold at Average Pixel Value
            if(meanThresh)
                seg{1} = im>mean(im(:));
            end
            
            
            
            % Isotropic Diffusion (Q-matrix contains no information about the image)
            if(IsoDiff)
                
                if exist([dirMatEig,'IsotropicDiffusion_rM',num2str(rM),'.mat'],'file')
                    
                    disp('Loading Existing File for Isotropic Diffusion Matrix & Eigen-Stuffs.')
                    load([dirMatEig,'IsotropicDiffusion_rM',num2str(rM),'.mat'])
                    
                else
   
                    disp('Cant find Existing File for Isotropic Diffusion. Making One.')
                    Q = Wdist;
                    [EVecsML,EValsML] = eig(Q); % MAY NEED TO RERUN IMAGE SEG BECAUSE IT WAS EIGS.
                    
                    % VIP!!! - TO USE THE EIGS MATLAB FUNCTION, Q MUST BE "square and should be large and sparse".
                    % THINK ABOUT FOR WHAT METHODS THESE DO NOT HOLD TRUE! FOR THOSE, we should use eig()
                    
                    EValsML = diag(EValsML);
                    %
                    save([dirMatEig,'IsotropicDiffusion_rM',num2str(rM),'.mat'],'Q','EVecsML','EValsML','imdim')
                
                end
            end


            disp(['We have created the Connectivity Matrix from the input image...'])
            disp(['Back out in SegmentMethod : Below we do metaClusterAnalysis & Kuramoto Simulation'])



            % SINCE IT TAKES SO LONG TO MAKE IT NOW.  MAYBE I SHOULD SAVE Q IN A MAT FILE FOR THIS FULL IMAGE?


            Qsparseness = numel(find(Q(:)))./numel(Q(:))
            

            % Set up netflags data structure to hold onto information about network construction
            netflags_setup_IMG    % script to fill in netflags data structure.
            netParams_setup_IMG

            netflags.method = method;
            netflags.sD = sD;
            netflags.sP = sP;
            netflags.rM = rM;



            % Shift eigenvalues of Q matrix to make them all positive because power 
            % method pulls out eigenvector with largest absolute value eigenvalue.
            if(shiftEvals && ~meanThresh)
                Q = Q + 1.1.*abs(min(Q(:))).*eye(size(Q)); % add C*(identity matrix) 
            end
            % THIS IS TOO LARGE FOR FULL BSDS (400x300) IMAGES.
            
            
            
            
            
            
            
            % For Isotropic Diffusion Case, I want Q (from unnormalized AA with sP=inf) to be identical to Wdist. This if statement checks that.
            if( 0 & ~isempty(find(Q - Wdist,1)) )
                disp('Isotropic Diffusion Confusion : Q should be same as Wdist. It isnt.')
                keyboard
            end
            
            
            
            
            
            


            


            %% Look at distribution of weights in Graph constructed
            if(0)  

                % NOTE: WILL HAVE TO CHANGE GNDTRUTH INDEXING TO MATCH IMAGE PATCH AND SEGMENTATION FOR BSDS.
                x = repmat(reshape(gT{1},1,numel(gT{1})), numel(gT{1}), 1 );
                y = repmat(reshape(gT{1},1,numel(gT{1}))', 1, numel(gT{1}) );
                clusterTruth = (x==y); % This is an NxN matrix that has a 1 when 2 nodes are in same cluster and 0 otherwise.

                [histW, histWX] = hist(Q(:),100);                         % Histogram of all weights in Q Graph/Network.
                [histWin, histWinX] = hist(Q(find(clusterTruth)),100);    % When 2 pixels are within a cluster (Rmax not considered so will have more zeros)
                [histWout, histWoutX] = hist(Q(find(~clusterTruth)),100); % When 2 pixels are in different cluster
                %
                [histWin2, histWin2X] = hist(Q(find(clusterTruth.*Wdist)),100);      % When 2 pixels are within a cluster (Rmax taken into consideration)
                [histWout2, histWout2X] = hist(Q(find(~clusterTruth.*Wdist)),100);   % When 2 pixels are different a cluster 



                h=figure;
                %
                subplot(351)
                imagesc(im), %colorbar
                title({'Example Image Input',InputImages})
                axis off
                %
                subplot(352)
                imagesc(Q), %colorbar
                title({'Graph Coupling Matrix Q',['Method = ',method]})
                axis off
                %
                subplot(353), 
                imagesc(clusterTruth), %colorbar
                title('Ground Truth for Coupling (+ means same cluster)')
                axis off
                %
                subplot(354), 
                imagesc(Wdist), %colorbar
                title(['Distance < ',num2str(Rmax)])
                axis off
                %
                subplot(355), 
                imagesc(Wdist.*clusterTruth), %colorbar
                title('Distance & Connectivity')
                axis off
                %
                subplot(312), hold on
                plot(histWinX,histWin./sum(histWin),'r','LineWidth',2)
                plot(histWoutX,histWout./sum(histWout),'b','LineWidth',2)
                title(['Graph Weight Distribution ignoring Rmax = ',num2str(Rmax)],'FontSize',20,'FontWeight','Bold')
                legend('Within Cluster','Across Cluster')
                set(gca,'FontSize',16,'FontWeight','Bold')

                %
                subplot(313), hold on
                plot(histWin2X,histWin2./sum(histWin2),'m','LineWidth',2)
                plot(histWout2X,histWout2./sum(histWout2),'c','LineWidth',2)
                title(['Graph Weight Distribution inside Rmax = ',num2str(Rmax)],'FontSize',20,'FontWeight','Bold')
                legend('Within Cluster','Across Cluster')
                xlabel('Weight Between Vertices','FontSize',18,'FontWeight','Bold')
                ylabel('% of Occurances','FontSize',18,'FontWeight','Bold')
                set(gca,'FontSize',16,'FontWeight','Bold')

                keyboard

                %

    %                     saveGoodImg(h,[dirPre,'output/ImgSeg/simpleExamples/',InputImages,'/pics_from_segmentation/evecs/'...
    %                             ,OutFname,'_GraphWeightDist'],sizeGoodIm)
    %                     close(h);

            end



            %% Do MetaClustering Analysis on Eigenvectors
            if(DoEigProcessing)


                % If the mat file exists already, move on to next without calculating any eigenvector or eigenvalues.
                if( dataFileChk &  exist(EigMatOut,'file') )
                    disp('Eig Data File already exists:')
                    EigMatOut
                    disp('Next...')
                else


                    % Matlab's EIGS function to compute first 6 eigenvectors & eigenvalues
                    if(1 && ~meanThresh)


                        % Calculate Dominant Eigenvector using MATLAB's eig function
                        if( ~strcmp(method,'IsoDiff') )
                            
            %             try
                            disp(['Running Matlab''s EIGs function on ',num2str(size(Q,1)),'x',num2str(size(Q,2)),' Q matrix.'])
                            tic
                            [EVecsML,EValsML] = eigs(Q);
                            EValsML = diag(EValsML);
                            toc
            %             catch
            %                 % if eig does not converge.
            %                 disp('Eig didnt converge.  I think that is the problem.  Filling Eig with zeros.')
            %                 EigVec = zeros(size(Q));
            %                 EigVal = zeros(size(Q));
            %                 % Note: Not sure filling these with zeros will not lead to
            %                 % errors later.  We will just have to cross that bridge
            %                 % when we come to it.
            %             end
            
                        end

                        



                        % Visualize top 6 Eigenvectors and original image.
                        if(doPlots)

                            % Nonlinearity applied to Eigenvectors to visualize them.
                            vizNonlin=1e-16;

                            hEig=figure; colormap('bone')
                            %
                            subplot(1,2,1), imagesc(im), set(gca,'Xtick',[],'Ytick',[]), axis square, title(['Image Pixels'],'FontSize',20,'FontWeight','Bold')
                            xlabel(['VizNonLin=',num2str(vizNonlin),' --> '],'FontSize',18,'FontWeight','Bold')
                            %
                            % Evecs 1-6
                            subplot(4,6,4),imagesc(reshape( EVecsML(:,1), yimg,ximg)), set(gca,'Xtick',[],'Ytick',[]), axis square, title('Evec1','FontSize',14,'FontWeight','Bold')
                            subplot(4,6,5),imagesc(reshape( EVecsML(:,2), yimg,ximg)), set(gca,'Xtick',[],'Ytick',[]), axis square, title('2','FontSize',14,'FontWeight','Bold')
                            subplot(4,6,6),imagesc(reshape( EVecsML(:,3), yimg,ximg)), set(gca,'Xtick',[],'Ytick',[]), axis square, title('3','FontSize',14,'FontWeight','Bold')
                            subplot(4,6,10),imagesc(reshape( EVecsML(:,4), yimg,ximg)), set(gca,'Xtick',[],'Ytick',[]), axis square, title('4','FontSize',14,'FontWeight','Bold')
                            subplot(4,6,11),imagesc(reshape( EVecsML(:,5), yimg,ximg)), set(gca,'Xtick',[],'Ytick',[]), axis square, title('5','FontSize',14,'FontWeight','Bold')
                            subplot(4,6,12),imagesc(reshape( EVecsML(:,6), yimg,ximg)), set(gca,'Xtick',[],'Ytick',[]), axis square, title('6','FontSize',14,'FontWeight','Bold')
                            %
                            % Evecs 1-6 with nonlinear visualization
                            subplot(4,6,16),imagesc(reshape( EvecVizF(EVecsML(:,1),vizNonlin), yimg,ximg)),  set(gca,'Xtick',[],'Ytick',[]), axis square, title('Evec1v','FontSize',14,'FontWeight','Bold')
                            subplot(4,6,17),imagesc(reshape( EvecVizF(EVecsML(:,2),vizNonlin), yimg,ximg)),  set(gca,'Xtick',[],'Ytick',[]), axis square, title('2v','FontSize',14,'FontWeight','Bold')
                            subplot(4,6,18),imagesc(reshape( EvecVizF(EVecsML(:,3),vizNonlin), yimg,ximg)),  set(gca,'Xtick',[],'Ytick',[]), axis square, title('3v','FontSize',14,'FontWeight','Bold')
                            subplot(4,6,22),imagesc(reshape( EvecVizF(EVecsML(:,4),vizNonlin), yimg,ximg)), set(gca,'Xtick',[],'Ytick',[]), axis square, title('4v','FontSize',14,'FontWeight','Bold')
                            subplot(4,6,23),imagesc(reshape( EvecVizF(EVecsML(:,5),vizNonlin), yimg,ximg)), set(gca,'Xtick',[],'Ytick',[]), axis square, title('5v','FontSize',14,'FontWeight','Bold')
                            subplot(4,6,24),imagesc(reshape( EvecVizF(EVecsML(:,6),vizNonlin), yimg,ximg)), set(gca,'Xtick',[],'Ytick',[]), axis square, title('6v','FontSize',14,'FontWeight','Bold')


                            % Save this figure as a jpeg so I can compare different methods of network construction & Eigenvector vs. Kuramoto
                            EigPicOut = EigMatOut;
                            st = strfind(EigPicOut,'data');
                            nd = st+3;
                            EigPicOut(st:nd) = 'imgs';      % redirect 'data' dir into 'imgs' dir
                            EigPicOut = EigPicOut(1:end-4); % get rid of .mat
                            %
                            saveGoodImg(hEig,EigPicOut,sizeGoodIm)
                            close(hEig);


                        end



                    else % if MATLAB's eig was not used or method was thresholding img pixels at mean.

                        EVecsML = 'MATLAB eigs not used. Segmentation of Image at Mean Pixel Value.';
                        EValsML = 'MATLAB eigs not used. Segmentation of Image at Mean Pixel Value.';

                    end % If statement to calculate eigenvectors using MATLAB eig function and not using mean threshold method.





                    % Power Method to compute dominant (or 2nd) eigenvector of L and Q and W
                    if(0 && ~meanThresh)
                        disp('Power Method Eigenvector Calc')
                        tic

            %                     if(proj)
            %                         EvecTrue = EvecBig2;
            %                         ev = '2';
            %                     else
            %                         EvecTrue = EvecBig1;
            %                         ev = '1';
            %                     end

                        % Works fine now... Even with Normalized Graph Laplacian & Average Association
                        threshPM = 1e-7; % if not small enough, Evector may not converge.
                        maxPMiter = 100000; % if too low, Evector doesnt converge and does not look right.
                        [EVecPM, ConvergencePM] = power_method(Q, ximg, yimg, EvecML, 'random', method, proj, threshPM, maxPMiter, 0,0,1,0);
                        toc

                        % Final Plot Dominant Eigenvector as found from Power Method
                        figure, imagesc(EVecPM); colorbar('FontSize',16,'FontWeight','Bold'); 
                        title(['Eigenvector #',ev,' PM ',method],'FontSize',20,'FontWeight','Bold')

                    else

                        EVecPM = 'Power Method was not run.';

                    end




                    if(doMCA)
%                         kurParams_setup % note:  I need these for the metaClusterAnalysis function, but they dont get 
%                                         % used inside because I am running it in Eigenvector Mode (not phase mode).  


                        % Compute matrix of pairwise hit probability and false alarm probability by setting optimal thresholds.
                        disp('Computing Boundary Discriminability for Eigenvectors for {im,ev1,ev2o,ev3o,ev2,ev3,ev2w,ev3w}')



                        im = netParams.im;
                        ev1 = reshape(EVecsML(:,1), netParams.Ndims);
                        ev2o = reshape(EVecsML(:,2), netParams.Ndims);
                        ev3o = reshape(EVecsML(:,3), netParams.Ndims);

                        ev2 = reshape(sum(EVecsML(:,1:2),2), netParams.Ndims);
                        ev3 = reshape(sum(EVecsML(:,1:3),2), netParams.Ndims);
                        ev2w = reshape(EVecsML(:,1).*EValsML(1) + EVecsML(:,2).*EValsML(2), netParams.Ndims);
                        ev3w = reshape(EVecsML(:,1).*EValsML(1) + EVecsML(:,2).*EValsML(2) + EVecsML(:,3).*EValsML(3), netParams.Ndims);

                        
                        % Compute boundary based Gradient Metric
                        [MC.im.F] = compute_Spatial_Gradient(im, 0);
                        [MC.ev1.F] = compute_Spatial_Gradient(ev1, 0);
                        [MC.ev2o.F] = compute_Spatial_Gradient(ev2o, 0);
                        [MC.ev3o.F] = compute_Spatial_Gradient(ev3o, 0);
                        [MC.ev2.F] = compute_Spatial_Gradient(ev2, 0);
                        [MC.ev3.F] = compute_Spatial_Gradient(ev3, 0);
                        [MC.ev2w.F] = compute_Spatial_Gradient(ev2w, 0);
                        [MC.ev3w.F] = compute_Spatial_Gradient(ev3w, 0);
                        
                        
                        % Plot and check d' metric values for image and various eigenvectors to make sure they make sense.
                        if(0)
                            
                            Hc=figure;
                            ha = tight_subplot(2, 8, [0.05 0.01], 0.05, 0.01);

                            % Image Pixels
                            axes(ha(1)), imagesc(im), colormap(bone), freezeColors, axis square, 
                            title(['\color{black}Img Pix'],'FontSize',20,'FontWeight','Bold') 
                            set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),  % colorbar, cbfreeze
                            %xlabel(['\color{black}d''= ',num2str(MC.im.D(1),2)],'FontSize',18,'FontWeight','Bold')
                            %
                            % Image Pixel Gradients
                            axes(ha(9)), imagesc(MC.im.F), colormap(jet), freezeColors, axis square, 
                            %title(['\color{black}d''= ',num2str(MC.im.D(1),2)],'FontSize',20,'FontWeight','Bold') 
                            set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
                            colorbar('SouthOutside'), cbfreeze

                            % Eigenvector 1
                            axes(ha(2)), imagesc(ev1), colormap(jet), freezeColors, axis square, 
                            title(['\color{black}Evec 1'],'FontSize',20,'FontWeight','Bold') 
                            set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
                            %xlabel(['\color{black}d''= ',num2str(MC.ev1.D(1),2)],'FontSize',18,'FontWeight','Bold')
                            colorbar('SouthOutside'), cbfreeze
                            %
                            % Eigenvector 1 Gradients
                            axes(ha(10)), imagesc(MC.ev1.F), colormap(jet), freezeColors, axis square, 
                            %title(['\color{black}d''= ',num2str(MC.ev1.D(1),2)],'FontSize',20,'FontWeight','Bold') 
                            set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
                            colorbar('SouthOutside'), cbfreeze

                            % Eigenvector 2
                            axes(ha(3)), imagesc(ev2o), colormap(jet), freezeColors, axis square, 
                            title(['\color{black}Evec 2'],'FontSize',20,'FontWeight','Bold') 
                            set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
                            %xlabel(['\color{black}d''= ',num2str(MC.ev2o.D(1),2)],'FontSize',18,'FontWeight','Bold')
                            colorbar('SouthOutside'), cbfreeze
                            %
                            % Eigenvector 2 Gradients
                            axes(ha(11)), imagesc(MC.ev2o.F), colormap(jet), freezeColors, axis square, 
                            %title(['\color{black}d''= ',num2str(MC.ev2o.D(1),2)],'FontSize',20,'FontWeight','Bold') 
                            set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
                            colorbar('SouthOutside'), cbfreeze

                            % Eigenvector 3
                            axes(ha(4)), imagesc(ev3o), colormap(jet), freezeColors, axis square, 
                            title(['\color{black}Evec 3'],'FontSize',20,'FontWeight','Bold') 
                            set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
                            %xlabel(['\color{black}d''= ',num2str(MC.ev3o.D(1),2)],'FontSize',18,'FontWeight','Bold')
                            colorbar('SouthOutside'), cbfreeze
                            %
                            % Eigenvector 3 Gradients
                            axes(ha(12)), imagesc(MC.ev3o.F), colormap(jet), freezeColors, axis square, 
                            %title(['\color{black}d''= ',num2str(MC.ev3o.D(1),2)],'FontSize',20,'FontWeight','Bold') 
                            set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
                            colorbar('SouthOutside'), cbfreeze

                            % Eigenvector 1&2
                            axes(ha(5)), imagesc(ev2), colormap(jet), freezeColors, axis square, 
                            title(['\color{black}Evec 1-2'],'FontSize',20,'FontWeight','Bold') 
                            set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
                            %xlabel(['\color{black}d''= ',num2str(MC.ev1.D(1),2)],'FontSize',18,'FontWeight','Bold')
                            colorbar('SouthOutside'), cbfreeze
                            %
                            % Eigenvector 1&2 Gradients
                            axes(ha(13)), imagesc(MC.ev2.F), colormap(jet), freezeColors, axis square, 
                            %title(['\color{black}d''= ',num2str(MC.ev2.D(1),2)],'FontSize',20,'FontWeight','Bold') 
                            set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
                            colorbar('SouthOutside'), cbfreeze

                            % Eigenvector 1-3
                            axes(ha(6)), imagesc(ev3), colormap(jet), freezeColors, axis square, 
                            title(['\color{black}Evec 1-3'],'FontSize',20,'FontWeight','Bold') 
                            set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
                            %xlabel(['\color{black}d''= ',num2str(MC.ev1.D(1),2)],'FontSize',18,'FontWeight','Bold')
                            colorbar('SouthOutside'), cbfreeze
                            %
                            % Eigenvector 1-3 Gradients
                            axes(ha(14)), imagesc(MC.ev3.F), colormap(jet), freezeColors, axis square, 
                            %title(['\color{black}d''= ',num2str(MC.ev3.D(1),2)],'FontSize',20,'FontWeight','Bold') 
                            set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
                            colorbar('SouthOutside'), cbfreeze
                            
                            % Eigenvector 1&2 Weighted
                            axes(ha(7)), imagesc(ev2), colormap(jet), freezeColors, axis square, 
                            title(['\color{black}Evec 1-2'],'FontSize',20,'FontWeight','Bold') 
                            set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
                            %xlabel(['\color{black}d''= ',num2str(MC.ev1.D(1),2)],'FontSize',18,'FontWeight','Bold')
                            colorbar('SouthOutside'), cbfreeze
                            %
                            % Eigenvector 1&2 Weighted Gradients
                            axes(ha(15)), imagesc(MC.ev2w.F), colormap(jet), freezeColors, axis square, 
                            %title(['\color{black}d''= ',num2str(MC.ev2w.D(1),2)],'FontSize',20,'FontWeight','Bold') 
                            set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
                            colorbar('SouthOutside'), cbfreeze

                            % Eigenvector 1-3 Weighted
                            axes(ha(8)), imagesc(ev3), colormap(jet), freezeColors, axis square, 
                            title(['\color{black}Evec 1-3'],'FontSize',20,'FontWeight','Bold') 
                            set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
                            %xlabel(['\color{black}d''= ',num2str(MC.ev1.D(1),2)],'FontSize',18,'FontWeight','Bold')
                            colorbar('SouthOutside'), cbfreeze
                            %
                            % Eigenvector 1-3 Weighted Gradients
                            axes(ha(16)), imagesc(MC.ev3w.F), colormap(jet), freezeColors, axis square, 
                            %title(['\color{black}d''= ',num2str(MC.ev3w.D(1),2)],'FontSize',20,'FontWeight','Bold') 
                            set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
                            colorbar('SouthOutside'), cbfreeze 
                            
                            keyboard
                            
                        end
                        


                    else

                        disp('NOT DOING META CLUSTER ANALYSIS ON EIGENVECTORS!! CHECK THIS!!')
                        MC = 'meh, didnt do it.';
                        

                    end


%                     clear kurParams % delete kurParams just to be clean here.  I made them but didnt use them.  
%                                     % Will make them again if i need them below.


                    % save an output mat file that I will then pipe through BSDS benchmark code.
                    if(save_Evecs_Mat)  % & ~meanThresh
                        
                        disp(['Saving Eigenvector MetaCluster Analysis Results in ',EigMatOut])
                        save([EigMatOut],'EVecsML','EValsML','EVecPM','netflags','netParams','MC')
                        
                        clear EVecsML EValsML EVecPM MC
                        
                    end


                end % If checking data file and the file exists or doesnt.

            end % if Do Eig Processing



            %% Use Kuramoto Coupled Oscillator Simulation for Segmentation
            if(DoKurProcessing)

                % If the mat file exists already, move on to next without calculating any eigenvector or eigenvalues.
                if( dataFileChk &  exist(KurMatOut,'file') )
                    disp('Kur Data File already exists:')
                    KurMatOut
                    disp('Next...')
                else

                    kurflags_setup  % runflags from the kuramoto simulation code (what plots/movies to make, etc.)
                    kurParams_setup % parameters from kuramoto simulation (phase interaction function, time, runs, NF distribution, etc.)
                    
                    
                    % Compute Kscale
                    M = sum(Q~=0);
                    rowsum = sum(Q);
                    
                    kurParams.sigW = SigW;
                    
                    TiScale = 1;                 % hardcode: It is best to use the whole phase range.
                    kurParams.TiScale = TiScale; % scale determines how to pack filter activations into phase space on phase init.
                    %
                    [kurflags] = kuramoto_output_dirs_setup_IMG(netParams, kurParams, netflags, kurflags, fname);
                    kurflags.dataFileChk = dataFileChk;
                    
                    
                    % Middling Kscale Value
                    Kscale = 15.*max(rowsum./M); % note: 15 comes from 60hz / 4 left over from (pi/2)/2*pi.
                    kurParams.Kscale = Kscale; % scale the weights in the coupling matrix by this scalar value.
                    %
                    [MC, metaCluster] = main_Kuramoto(netParams, kurParams, netflags, kurflags, Kscale.*Q, doMCA,'mid');
                    
                    
                    if(bracket_Kscale)
                        
                        % Small Kscale Value
                        Kscale = 1.5.*max(rowsum./M); % note: 15 comes from 60hz / 4 left over from (pi/2)/2*pi.
                        kurParams.Kscale = Kscale; % scale the weights in the coupling matrix by this scalar value.
                        %
                        [MC, metaCluster] = main_Kuramoto(netParams, kurParams, netflags, kurflags, Kscale.*Q, doMCA,'sml');
                    
                        % Large Kscale Value
                        Kscale = 150.*max(rowsum./M); % note: 15 comes from 60hz / 4 left over from (pi/2)/2*pi.
                        kurParams.Kscale = Kscale; % scale the weights in the coupling matrix by this scalar value.
                        %
                        [MC, metaCluster] = main_Kuramoto(netParams, kurParams, netflags, kurflags, Kscale.*Q, doMCA,'lrg');
                        
%                         % Huge Kscale Value
%                         Kscale = 1500.*max(rowsum./M); % note: 15 comes from 60hz / 4 left over from (pi/2)/2*pi.
%                         kurParams.Kscale = Kscale; % scale the weights in the coupling matrix by this scalar value.
%                         %
%                         [MC, metaCluster] = main_Kuramoto(netParams, kurParams, netflags, kurflags, Kscale.*Q, doMCA,'hug');
                    
                        
                    end
                    
                    


                    % Plot Kuramoto phase results, beginning phase embedding & ground truth
                    if(doPlots)

                        imageOrig = netParams.im;
                        [one, two] = find(imageOrig == min(imageOrig(:))); % find brightest pixel.
                        pixMax(1) = one(1);
                        pixMax(2) = two(1);

                        phaseFinal = reshape(metaCluster.phaseAtClk(:,end),netParams.Ndims(1),netParams.Ndims(2));
                        phaseOffset = pi - phaseFinal(pixMax(1),pixMax(2));           % global phase to add to all phases to make phase at brightest pixel = pi
                        phaseFinal = mod( phaseFinal + phaseOffset , 2*pi);           % rotate all phases by phase offset
                        phaseFinal = abs(phaseFinal - pi) ./ pi;                      % fold circular variabls down to [0 1] number line - ie. pixel space

                        phaseInit = reshape(metaCluster.phaseAtClk(:,1),netParams.Ndims(1),netParams.Ndims(2));
                        phaseOffset = pi - phaseInit(pixMax(1),pixMax(2));           % global phase to add to all phases to make phase at brightest pixel = pi
                        phaseInit = mod( phaseInit + phaseOffset , 2*pi);           % rotate all phases by phase offset
                        phaseInit = abs(phaseInit - pi) ./ pi;                      % fold circular variabls down to [0 1] number line - ie. pixel space


                        % Plot these things: 
                        if iscell(MC)
                            numVertSubs = 4;
                        else
                            numVertSubs = 3;
                        end


                        hKur=figure; colormap('bone')
                        subplot(numVertSubs,2,[1,3]), imagesc(imageOrig), set(gca,'Xtick',[],'Ytick',[]), title('Image Pix & Phase Init','FontSize',20,'FontWeight','Bold'), freezeColors
                        subplot(numVertSubs,2,[2,4]), imagesc(phaseFinal), set(gca,'Xtick',[],'Ytick',[]), title(['Phase after Sim (t=0.3s) '],'FontSize',20,'FontWeight','Bold'), freezeColors
                        %
                        for g = 1:numel(netParams.gT)
                            subplot(numVertSubs,numel(netParams.gT),2*numel(netParams.gT)+g), imagesc(netParams.gT{g}), colormap('jet'), set(gca,'Xtick',[],'Ytick',[]), freezeColors
                        end
                        subplot(numVertSubs,numel(netParams.gT),2*numel(netParams.gT)+1), ylabel('gnd truth','FontSize',18,'FontWeight','Bold')
                        %
                        if iscell(MC)
                            disp('Note: Im gonna need DivMarg')
                            subplot(numVertSubs,1,4),plot(kurParams.tau*kurParams.spp*(0:size(DivMarg,1)-1) , DivMarg)
                            xlabel('Time (sec)','FontSize',18,'FontWeight','Bold')
                            ylabel('Divisive Margin','FontSize',18,'FontWeight','Bold')
                            axis([0 kurParams.Tsec 0 1])
                        end

                        % Save this figure as a jpeg so I can compare different methods of network construction & Eigenvector vs. Kuramoto
                        KurPicOut = KurMatOut;
                        st = strfind(KurPicOut,'data');
                        nd = st+3;
                        KurPicOut(st:nd) = 'imgs';      % redirect 'data' dir into 'imgs' dir
                        KurPicOut = KurPicOut(1:end-4); % get rid of .mat
                        %
                        saveGoodImg(hKur,KurPicOut,sizeGoodIm)
                        %close(hKur);

                    end

                end % if dataFileCheck & mat file exists or not
                
                
                clear MC metaCluster

            end % if do Kur Processing


            close all

            clear netParams netflags


            % end % loop over images in image Ensemble (patches or different random draws from pixel distribution)

        end % loop over max radius for distance dependence (weighting between pixels separated by more than rmax = 0)    
        end % loop over threshold and sigma pairs on distance dependence
        end % loop over different sigma values on pixel similarity weighting

        % display how much memory matlab is currently using
    %     s = whos;
    %     disp(['Matlab is using ~',num2str( (sum([s.bytes]) + 500e6) / 1e9 ),' GB of memory'])
        % s.bytes tells memory for each variable and matlab uses about 500MB on startup.

    end % loop over files in preprocessed directory.

    disp('I HAVE SUCCESSFULLY FINISHED THIS BATCH OF IMAGES AND PARAMETER VALUES! at')
    clock

end



% Renormalize Eigenvector to be between 0 & 1 (same as I did for image
% patch) so that I know its range for further processing.
function [x] = Normlze(x)

    x = x - min(x(:));
    x = x ./max(x(:));

end

%% Code from Saurabh meeting to do benchmark stuff...
% load('/Users/world7one/Desktop/Grad_School/Berkeley/Work/Fritz_Work/Projects/images/BSDS_images/BSR/bench/data/ucm2/2018.mat')
% reg = bwlabel(ucm2 < 0.5); % <--- BIGGER NUMBER MEANS BIGGER SEGMENTS!
% reg = reg(2:2:end, 2:2:end);
% figure, imagesc(reg)


