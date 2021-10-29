function gen_ROC_AUC_scatter(fileType, fileSize, methodType, rMx,sDx,sPx,sWx,KSstr) %

    % syntax:   gen_ROC_AUC_scatter('BSDS_patch','51x51_ds1', ...
    % {'Mod_SKHAdj','GLnrm','AAnrm','Mod_N&G'},[1:10,inf],inf,[0.1,0.2,0.3,0.4],0,1,300)
    %
    % This function will loop through different parameter setting values
    % and will, for each unique setting of those parameter values, run
    % through all the pairs of (Kur,Eig) mat files for all 1500 image
    % patches and will gather up the AUC_ROC data into long vectors that we
    % can scatter plot and calculate statistics (mean, std) on.  The code
    % will also construct matching long vectors that contain the ratio of
    % sizes of the clusters in the pair, the image patch ID from which the
    % pair was taken, and the ground truth # as well as the cluster ID#
    % within the ground truth from with clusters were taken.  It computes
    % the metric, which is simply (AUCmeth - AUCim), and  builds up mean
    % and std for each method in AUC_diff_MethVsImg.  The following are
    % saved in the mat file:
    %
    %     AUC_diff_MethVsImg 
    %     im_ROC_AUC_1D_accum 
    %     kur_ROC_AUC_1D_accum 
    %     kur_ROC_AUC_GS_accum
    %     ev1_ROC_AUC_1D_accum 
    %     ev2o_ROC_AUC_1D_accum 
    %     ev3o_ROC_AUC_1D_accum 
    %     ev2_ROC_AUC_1D_accum 
    %     ev3_ROC_AUC_1D_accum 
    %     ev2w_ROC_AUC_1D_accum 
    %     ev3w_ROC_AUC_1D_accum 
    %     cluster_size_ratio 
    %     img_ptch_ID 
    %     grnd_truth_ID
    %
    % And scatter plot images are saved if the flag is set to 1 below.



    %% Load in mat file and set up directories to different data.
    [dirPre,sizeGoodIm] = onCluster;




    for K = 1:numel(methodType)
    

        % output directory to save the AUC_ROC_results mat files into.
        outDatDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/AUC_ROC_results/',methodType{K},'/'];
        
        fileType
        fileSize
        methodType{K}
        outDatDir
        
        if ~exist(outDatDir)
            mkdir(outDatDir)
        end
        
        % output directory to save the AUC_ROC_results plots and images into.
        outImgDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/imgs/AUC_ROC_results/',methodType{K},'/'];
        if ~exist(outImgDir)
            mkdir(outImgDir)
        end



%         % directory to image patches extracted from BSDS.
%         patchesDir = [dirPre,'images/',fileType,'/',fileSize,'/'];

        % directory to KurMC files and to Kur_metaSummary files (think I need the KurMC ones)
        dirKurMats = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/Kur_PIF_Fourier1/',methodType{K},'/'];

        % directory to the Evecs files
        dirEigMats = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/spectral/',methodType{K},'/'];


        % Making Labels and Legends for plots.
        for i = 1:numel(rMx)
            rM_lgnd{i} =  num2str(rMx(i));
        end
        %
        for j = 1:numel(sPx)
            sP_lgnd{j} =  num2str(sPx(j));
        end



         % Preallocate memory for array of mean values.
         segMethods = {'kur';'im';'ev1';'ev2o';'ev3o'}; % ;'ev2';'ev3';'ev2w';'ev3w'};




        %% Loop through parameter values (preset by user) to extract mean (across all image patches processed) Area Under ROC Curve.
        for j = 1:numel(sPx)

            for i = 1:numel(rMx)


                % Turn parameter values into strings  
                rM = num2str(rMx(i)); % '1';
                sD = num2str(sDx); % 'Inf';
                sP = num2str(sPx(j)); % '0p2';
                sP(sP=='.')='p';
                sW = num2str(sWx); % '0';

                if strcmp(methodType{K},'IsoDiff')
                    fileSpecEig = ['_rM',rM];
                else
                    fileSpecEig = ['_rM',rM,'_sD',sD,'_sP',sP];
                end
                fileSpecKur = [fileSpecEig,'_NF_60_',sW,'_ks',KSstr]; 
                % NEED TO ADD KS PART HERE!! (kslrg, ksmid, kssml)
                filesEig = dir([dirEigMats,'*',fileSpecEig,'.mat']);
                % filesKur = dir([dirKurMats,'*',fileSpecKur,'.mat']);
                
                
                im_ROC_AUC_1D_accum = [];
                kur_ROC_AUC_1D_accum = [];
                kur_ROC_AUC_GS_accum = [];
                ev1_ROC_AUC_1D_accum = [];
                ev2o_ROC_AUC_1D_accum = [];
                ev3o_ROC_AUC_1D_accum = [];
                ev2_ROC_AUC_1D_accum = [];
                ev3_ROC_AUC_1D_accum = [];
                ev2w_ROC_AUC_1D_accum = [];
                ev3w_ROC_AUC_1D_accum = [];
                

                count_kur_1D_biggerthan_GS = 0;
                count_isnan_kur_GS = 0;
                
                
                ImgPatchID = {'start'};
                
                ClusterPairID = [Inf;Inf];
                ClusterSizes = [Inf;Inf];
                GndTruthID = Inf;
                
                
                for k = 1:numel(filesEig)
                  
                    try
                        eig = load([dirEigMats,filesEig(k).name]);                  % load eigen mat file with AUC data. 
                        fileKur = ['KurMC_',eig.netflags.fname,fileSpecKur,'.mat']; % the matching kuramoto file.
                        kur = load([dirKurMats,fileKur]);                           % load kuramoto mat file with AUC data.
                    catch
                        [dirEigMats,filesEig(k).name]
                        [dirKurMats,fileKur]
                        disp('meh!')
                        continue
                    end
                        


                    for G = 1:numel(eig.netParams.gT) % loop over ground truths


                        
                        % NOTE: There is some bug that makes sizeCluster array too long sometimes so that its sum is > # of pixels total.  
                        % Might have to put a fix in here for that at some point... (if cluster_size_ratio is not same size as all the 
                        % ##ROC_AUC_1D_accum vectors)
                        
                        
                        
                        % compute cluster size ratio.
                        c=0;
                        for a = 1:numel(kur.sizeCluster{G})
                            for b = a+1:numel(kur.sizeCluster{G})

                                ClusterPairID(:,end+1) = [eig.idCluster{G}(a); eig.idCluster{G}(b)];
                                
                                ImgPatchID{end+1} = kur.netflags.fname;

                                ClusterSizes(:,end+1) = [kur.sizeCluster{G}(a); kur.sizeCluster{G}(b)];
                                GndTruthID(end+1) = G;
                                
                                c=c+1;
                                
                                if(numel(kur.sizeCluster{G})==2)
                                    
                                    temp_im = eig.AUC_ROC_1D.im{G};
                                    temp_kur1D = kur.AUC_ROC_1D.kur{G};
                                    temp_kurGS = kur.AUC_ROC_GS.kur{G};
                                    temp_ev1 = eig.AUC_ROC_1D.ev1{G};
                                    temp_ev2o = eig.AUC_ROC_1D.ev2o{G};
                                    temp_ev3o = eig.AUC_ROC_1D.ev3o{G};
                                    temp_ev2 = eig.AUC_ROC_1D.ev2{G};
                                    temp_ev3 = eig.AUC_ROC_1D.ev3{G};
                                    temp_ev2w = eig.AUC_ROC_1D.ev2w{G};
                                    temp_ev3w = eig.AUC_ROC_1D.ev3w{G}; 
                                    
                                else
                                
                                    temp_im(c) = eig.AUC_ROC_1D.im{G}(a,b);
                                    temp_kur1D(c) = kur.AUC_ROC_1D.kur{G}(a,b);
                                    temp_kurGS(c) = kur.AUC_ROC_GS.kur{G}(a,b);
                                    temp_ev1(c) = eig.AUC_ROC_1D.ev1{G}(a,b);
                                    temp_ev2o(c) = eig.AUC_ROC_1D.ev2o{G}(a,b);
                                    temp_ev3o(c) = eig.AUC_ROC_1D.ev3o{G}(a,b);
                                    temp_ev2(c) = eig.AUC_ROC_1D.ev2{G}(a,b);
                                    temp_ev3(c) = eig.AUC_ROC_1D.ev3{G}(a,b);
                                    temp_ev2w(c) = eig.AUC_ROC_1D.ev2w{G}(a,b);
                                    temp_ev3w(c) = eig.AUC_ROC_1D.ev3w{G}(a,b);  
                                
                                end
                                
                            end % loop over 1st possible cluster in pairs for 1 gT
                        end % loop over 2nd possible cluster in pairs for 1 gT
                        
                        
                        
                        
                        % NOTE: These could conceivably be NaN's (or at least the kur ones can be) because
                        % I put a catch in there that if the convex hull computation fails, just set its AUC_ROC value to NaN.
                        if(any(isnan(temp_kurGS)))
                            disp('Found a NAN in the kur Grid Search Upper Bound.  Was probably a problem finding convex hull.')
                            [dirKurMats,fileKur]
                            disp(['Ground Truth #',num2str(G)])
                            disp('Replacing Grid Search with 1D value.')
                            count_isnan_kur_GS = count_isnan_kur_GS + numel(find(isnan(temp_kurGS)))
                            temp_kurGS
                            temp_kurGS(isnan(temp_kurGS)) = temp_kur1D(isnan(temp_kurGS))
                        end
                        
                        if any(temp_kur1D - temp_kurGS > 1e-6)
                            [temp_kur1D;temp_kurGS]
                            count_kur_1D_biggerthan_GS = count_kur_1D_biggerthan_GS + numel(find(temp_kur1D > temp_kurGS))
                            filesEig(k).name
                            fileKur
                            %keyboard
                        end
                        
                        
                        im_ROC_AUC_1D_accum = [im_ROC_AUC_1D_accum, temp_im];
                        kur_ROC_AUC_1D_accum = [kur_ROC_AUC_1D_accum, temp_kur1D];
                        kur_ROC_AUC_GS_accum = [kur_ROC_AUC_GS_accum, temp_kurGS];
                        ev1_ROC_AUC_1D_accum = [ev1_ROC_AUC_1D_accum, temp_ev1];
                        ev2o_ROC_AUC_1D_accum = [ev2o_ROC_AUC_1D_accum, temp_ev2o];
                        ev3o_ROC_AUC_1D_accum = [ev3o_ROC_AUC_1D_accum, temp_ev3o];
                        ev2_ROC_AUC_1D_accum = [ev2_ROC_AUC_1D_accum, temp_ev2];
                        ev3_ROC_AUC_1D_accum = [ev3_ROC_AUC_1D_accum, temp_ev3];
                        ev2w_ROC_AUC_1D_accum = [ev2w_ROC_AUC_1D_accum, temp_ev2w];
                        ev3w_ROC_AUC_1D_accum = [ev3w_ROC_AUC_1D_accum, temp_ev3w]; 
                        
                        
                        if(0)
                            im_ROC_AUC_1D_accum
                            ClusterSizes
                            GndTruthID
                            ClusterPairID
                            ImgPatchID
                            keyboard
                        end
                        

                        
                        
                        clear temp_im temp_kur1D temp_kurGS temp_ev1 temp_ev2o temp_ev3o temp_ev2 temp_ev3 temp_ev2w temp_ev3w
                        

                    end % Loop over Ground Truths
                    
                    clear kur eig

                    k    

                end % Loop over mat files in the Eigen Directory.
                
                
                
                
                
                
                % Compute single number metric.  How much better (mean & std) is AUC after applying method vs AUC using image pixels.
                AUC_diff_MethVsImg.Kur1D(1) = mean(kur_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                AUC_diff_MethVsImg.Kur1D(2) = std(kur_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                %
                AUC_diff_MethVsImg.KurGS(1) = mean(kur_ROC_AUC_GS_accum - im_ROC_AUC_1D_accum);
                AUC_diff_MethVsImg.KurGS(2) = std(kur_ROC_AUC_GS_accum - im_ROC_AUC_1D_accum);
                %
                AUC_diff_MethVsImg.Ev1(1) = mean(ev1_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                AUC_diff_MethVsImg.Ev1(2) = std(ev1_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                %
                AUC_diff_MethVsImg.Ev2o(1) = mean(ev2o_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                AUC_diff_MethVsImg.Ev2o(2) = std(ev2o_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                %
                AUC_diff_MethVsImg.Ev3o(1) = mean(ev3o_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                AUC_diff_MethVsImg.Ev3o(2) = std(ev3o_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                %
                AUC_diff_MethVsImg.Ev2(1) = mean(ev2_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                AUC_diff_MethVsImg.Ev2(2) = std(ev2_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                %
                AUC_diff_MethVsImg.Ev3(1) = mean(ev3_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                AUC_diff_MethVsImg.Ev3(2) = std(ev3_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                %
                AUC_diff_MethVsImg.Ev2w(1) = mean(ev2w_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                AUC_diff_MethVsImg.Ev2w(2) = std(ev2w_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                %
                AUC_diff_MethVsImg.Ev3w(1) = mean(ev3w_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                AUC_diff_MethVsImg.Ev3w(2) = std(ev3w_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                
                
                
                
                
                % Get rid of 1st entry on these data-accumulation arrays. Was just used so I could use (end+1) notation. 
                ClusterSizes      = ClusterSizes(:,2:end);
                GndTruthID        = GndTruthID(2:end);
                ClusterPairID     = ClusterPairID(:,2:end);
                ImgPatchID        = ImgPatchID(2:end);
                
                                
                
                % Save a mat file with AUC values for each method.  TO DO..
                save([outDatDir,'AUCdata_',methodType{K},fileSpecKur,'.mat'],'AUC_diff_MethVsImg', 'im_ROC_AUC_1D_accum',...
                    'kur_ROC_AUC_1D_accum', 'kur_ROC_AUC_GS_accum', 'ev1_ROC_AUC_1D_accum', 'ev2o_ROC_AUC_1D_accum', ...
                    'ev3o_ROC_AUC_1D_accum', 'ev2_ROC_AUC_1D_accum', 'ev3_ROC_AUC_1D_accum', 'ev2w_ROC_AUC_1D_accum', ...
                    'ev3w_ROC_AUC_1D_accum', 'ClusterSizes', 'GndTruthID', 'ClusterPairID', 'ImgPatchID', ...
                    'count_kur_1D_biggerthan_GS','count_isnan_kur_GS')
                
                
                % Scatter plots of different methods compared with image pixels performance by AUC measure.
                if(1)
                    % (1). Scatter plot im AUC vs Kuramoto AUC. 
                    h1=figure; hold on % figure for Scatter plot im AUC vs method AUC.
                    axis([0.5 1 0.5 1])
                    xlabel('AUC (im)','FontSize',18,'FontWeight','Bold')
                    ylabel('AUC (kur)','FontSize',18,'FontWeight','Bold')
                    title([methodType{K},' Kuramoto Area Under ROC Curve'],'FontSize',20,'FontWeight','Bold')
                    set(gca,'FontSize',16,'FontWeight','Bold')
                    scatter(im_ROC_AUC_1D_accum, kur_ROC_AUC_1D_accum, 'b.')
                    scatter(im_ROC_AUC_1D_accum, kur_ROC_AUC_GS_accum, 'r.')
                    plot([0.5 1],[0.5 1],'k--')
                    % legend({'Kur 1D','Kur U.B.'},'Location','Best')
                    % axis square
                    % 
                    met1 = mean(kur_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                    met2 = std(kur_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                    met3 = mean(kur_ROC_AUC_GS_accum - im_ROC_AUC_1D_accum);
                    met4 = std(kur_ROC_AUC_GS_accum - im_ROC_AUC_1D_accum);
                    %
                    text(0.85,0.56,['\color{blue}{\{1D\}=',num2str(met1,2),'\pm',num2str(met2,2),'}'],'BackgroundColor','white','FontSize',16,'FontWeight','Bold')
                    text(0.85,0.53,['\color{red}{\{UB\}=',num2str(met3,2),'\pm',num2str(met4,2),'}'],'BackgroundColor','white','FontSize',16,'FontWeight','Bold')
                    saveGoodImg(h1,[outImgDir,'Kur_',methodType{K},fileSpecKur],sizeGoodIm)
                    close(h1)

                    % (2). Scatter plot im AUC vs Evec1 AUC. 
                    h2=figure; hold on % figure for Scatter plot im AUC vs method AUC.
                    axis([0.5 1 0.5 1])
                    xlabel('AUC (im)','FontSize',18,'FontWeight','Bold')
                    ylabel('AUC (ev1)','FontSize',18,'FontWeight','Bold')
                    title([methodType{K},' Eigenvector 1 Area Under ROC Curve'],'FontSize',20,'FontWeight','Bold')
                    set(gca,'FontSize',16,'FontWeight','Bold')
                    scatter(im_ROC_AUC_1D_accum, ev1_ROC_AUC_1D_accum, 'b.')
                    plot([0.5 1],[0.5 1],'k--')
                    %axis square
                    % 
                    met1 = mean(ev1_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                    met2 = std(ev1_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                    %
                    text(0.85,0.56,['\color{blue}{\{1D\}=',num2str(met1,2),'\pm',num2str(met2,2),'}'],'BackgroundColor','white','FontSize',16,'FontWeight','Bold')
                    saveGoodImg(h2,[outImgDir,'Ev1_',methodType{K},fileSpecKur],sizeGoodIm)
                    close(h2)

                    % (3). Scatter plot im AUC vs Evec2 AUC. 
                    h3=figure; hold on % figure for Scatter plot im AUC vs method AUC.
                    axis([0.5 1 0.5 1])
                    xlabel('AUC (im)','FontSize',18,'FontWeight','Bold')
                    ylabel('AUC (ev2)','FontSize',18,'FontWeight','Bold')
                    title([methodType{K},' Eigenvector 2 Area Under ROC Curve'],'FontSize',20,'FontWeight','Bold')
                    set(gca,'FontSize',16,'FontWeight','Bold')
                    scatter(im_ROC_AUC_1D_accum, ev2o_ROC_AUC_1D_accum, 'b.')
                    plot([0.5 1],[0.5 1],'k--')
                    %axis square
                    % 
                    met1 = mean(ev2o_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                    met2 = std(ev2o_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                    %
                    text(0.85,0.56,['\color{blue}{\{1D\}=',num2str(met1,2),'\pm',num2str(met2,2),'}'],'BackgroundColor','white','FontSize',16,'FontWeight','Bold')
                    saveGoodImg(h3,[outImgDir,'Ev2o_',methodType{K},fileSpecKur],sizeGoodIm)
                    close(h3)

                    % (4). Scatter plot im AUC vs Evec3 AUC. 
                    h4=figure; hold on % figure for Scatter plot im AUC vs method AUC.
                    axis([0.5 1 0.5 1])
                    xlabel('AUC (im)','FontSize',18,'FontWeight','Bold')
                    ylabel('AUC (ev3)','FontSize',18,'FontWeight','Bold')
                    title([methodType{K},' Eigenvector 3 Area Under ROC Curve'],'FontSize',20,'FontWeight','Bold')
                    set(gca,'FontSize',16,'FontWeight','Bold')
                    scatter(im_ROC_AUC_1D_accum, ev3o_ROC_AUC_1D_accum, 'b.')
                    plot([0.5 1],[0.5 1],'k--')
                    %axis square
                    % 
                    met1 = mean(ev3o_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                    met2 = std(ev3o_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                    %
                    text(0.85,0.56,['\color{blue}{\{1D\}=',num2str(met1,2),'\pm',num2str(met2,2),'}'],'BackgroundColor','white','FontSize',16,'FontWeight','Bold')
                    saveGoodImg(h4,[outImgDir,'Ev3o_',methodType{K},fileSpecKur],sizeGoodIm)
                    close(h4)

                    % (5). Scatter plot im AUC vs Evec1-2 AUC. 
                    h5=figure; hold on % figure for Scatter plot im AUC vs method AUC.
                    axis([0.5 1 0.5 1])
                    xlabel('AUC (im)','FontSize',18,'FontWeight','Bold')
                    ylabel('AUC (ev1-2)','FontSize',18,'FontWeight','Bold')
                    title([methodType{K},' Eigenvector 1-2 Area Under ROC Curve'],'FontSize',20,'FontWeight','Bold')
                    set(gca,'FontSize',16,'FontWeight','Bold')
                    scatter(im_ROC_AUC_1D_accum, ev2_ROC_AUC_1D_accum, 'b.')
                    plot([0.5 1],[0.5 1],'k--')
                    %axis square
                    % 
                    met1 = mean(ev2_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                    met2 = std(ev2_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                    %
                    text(0.85,0.56,['\color{blue}{\{1D\}=',num2str(met1,2),'\pm',num2str(met2,2),'}'],'BackgroundColor','white','FontSize',16,'FontWeight','Bold')
                    saveGoodImg(h5,[outImgDir,'Ev2_',methodType{K},fileSpecKur],sizeGoodIm)
                    close(h5)

                    % (6). Scatter plot im AUC vs Evec1-3 AUC. 
                    h6=figure; hold on % figure for Scatter plot im AUC vs method AUC.
                    axis([0.5 1 0.5 1])
                    xlabel('AUC (im)','FontSize',18,'FontWeight','Bold')
                    ylabel('AUC (ev1-3)','FontSize',18,'FontWeight','Bold')
                    title([methodType{K},' Eigenvector 1-3 Area Under ROC Curve'],'FontSize',20,'FontWeight','Bold')
                    set(gca,'FontSize',16,'FontWeight','Bold')
                    scatter(im_ROC_AUC_1D_accum, ev3_ROC_AUC_1D_accum, 'b.')
                    plot([0.5 1],[0.5 1],'k--')
                    %axis square
                    % 
                    met1 = mean(ev3_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                    met2 = std(ev3_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                    %
                    text(0.85,0.56,['\color{blue}{\{1D\}=',num2str(met1,2),'\pm',num2str(met2,2),'}'],'BackgroundColor','white','FontSize',16,'FontWeight','Bold')
                    saveGoodImg(h6,[outImgDir,'Ev3_',methodType{K},fileSpecKur],sizeGoodIm)
                    close(h6)

                    % (7). Scatter plot im AUC vs Evec1-2w AUC. 
                    h7=figure; hold on % figure for Scatter plot im AUC vs method AUC.
                    axis([0.5 1 0.5 1])
                    xlabel('AUC (im)','FontSize',18,'FontWeight','Bold')
                    ylabel('AUC (ev1-2w)','FontSize',18,'FontWeight','Bold')
                    title([methodType{K},' Eigenvector 1-2w Area Under ROC Curve'],'FontSize',20,'FontWeight','Bold')
                    set(gca,'FontSize',16,'FontWeight','Bold')
                    scatter(im_ROC_AUC_1D_accum, ev2w_ROC_AUC_1D_accum, 'b.')
                    plot([0.5 1],[0.5 1],'k--')
                    %axis square
                    % 
                    met1 = mean(ev2w_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                    met2 = std(ev2w_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                    %
                    text(0.85,0.56,['\color{blue}{\{1D\}=',num2str(met1,2),'\pm',num2str(met2,2),'}'],'BackgroundColor','white','FontSize',16,'FontWeight','Bold')
                    saveGoodImg(h7,[outImgDir,'Ev2w_',methodType{K},fileSpecKur],sizeGoodIm)
                    close(h7)
                    

                    % (8). Scatter plot im AUC vs Evec1-3w AUC. 
                    h8=figure; hold on % figure for Scatter plot im AUC vs method AUC.
                    axis([0.5 1 0.5 1])
                    xlabel('AUC (im)','FontSize',18,'FontWeight','Bold')
                    ylabel('AUC (ev1-3w)','FontSize',18,'FontWeight','Bold')
                    title([methodType{K},' Eigenvector 1-3w Area Under ROC Curve'],'FontSize',20,'FontWeight','Bold')
                    set(gca,'FontSize',16,'FontWeight','Bold')
                    scatter(im_ROC_AUC_1D_accum, ev3w_ROC_AUC_1D_accum, 'b.')
                    plot([0.5 1],[0.5 1],'k--')
                    %axis square
                    % 
                    met1 = mean(ev3w_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                    met2 = std(ev3w_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum);
                    %
                    text(0.85,0.56,['\color{blue}{\{1D\}=',num2str(met1,2),'\pm',num2str(met2,2),'}'],'BackgroundColor','white','FontSize',16,'FontWeight','Bold')
                    saveGoodImg(h8,[outImgDir,'Ev3w_',methodType{K},fileSpecKur],sizeGoodIm)
                    close(h8)
                end
                
                
                % Does doing ROC cluster analysis on a pair of clusters
                % with unequal sizes bias results for Area under ROC curve?
                % According to this, it seems like not.
                if(0)
                    
                    % Looks like many pairs of clusters are not of equal size. 
                    figure, hist(cluster_size_ratio)
                    
                    % Plot AUC vs cluster size ratios (but I dont really see a trend)
                    figure, scatter(cluster_size_ratio,im_ROC_AUC_1D_accum)
                    figure, scatter(cluster_size_ratio,ev1_ROC_AUC_1D_accum)
                    figure, scatter(cluster_size_ratio,ev2o_ROC_AUC_1D_accum)
                    figure, scatter(cluster_size_ratio,ev3o_ROC_AUC_1D_accum)
                    figure, scatter(cluster_size_ratio,ev3_ROC_AUC_1D_accum)
                    figure, scatter(cluster_size_ratio,ev3w_ROC_AUC_1D_accum)
                    figure, scatter(cluster_size_ratio,ev2_ROC_AUC_1D_accum)
                    
                    % Plot AUC metric (method - img) vs cluster size ratio. (still no trend)
                    figure, scatter(cluster_size_ratio,kur_ROC_AUC_GS_accum - im_ROC_AUC_1D_accum)
                    figure, scatter(cluster_size_ratio,ev1_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum)
                    figure, scatter(cluster_size_ratio,ev2o_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum)
                    figure, scatter(cluster_size_ratio,ev3o_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum)
                    
                end
                
                pause(1)

                clear AUC_diff_MethVsImg im_ROC_AUC_1D_accum kur_ROC_AUC_1D_accum kur_ROC_AUC_GS_accum  ...
                    ev1_ROC_AUC_1D_accum ev2o_ROC_AUC_1D_accum ev3o_ROC_AUC_1D_accum ev2_ROC_AUC_1D_accum  ...
                    ev3_ROC_AUC_1D_accum ev2w_ROC_AUC_1D_accum ev3w_ROC_AUC_1D_accum cluster_size_ratio ...
                    img_ptch_ID grnd_truth_ID


            end  % loop over rM

        end % loop over sP


    
    
    end % loop over net method type.
    

    disp('Function Successfully Completed.')
    clock

end % end main function
