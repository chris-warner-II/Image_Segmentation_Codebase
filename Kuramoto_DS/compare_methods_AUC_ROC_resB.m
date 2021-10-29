function compare_methods_AUC_ROC_resB(fileType, fileSize)

    % syntax:  compare_methods_AUC_ROC_resB(fileType, fileSize)
    %
    %  This function will spit out plots of average Area under ROC curve
    %  with STD error bars for the 4 different netMethods and for 4
    %  different segMethods also.
    %
    % fileType = 'BSDS_patch';
    % fileSize = '51x51_ds1';




    %% Load in mat file and set up directories to different data.
    [dirPre,sizeGoodIm] = onCluster;

    matFilesDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/AUC_ROC_results/'];
    
    imgFilesDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/imgs/AUC_ROC_results/'];
    
    do_AA = 1;
    do_GL = 1;
    do_SK = 1;
    do_NG = 1;
    do_Iso = 1;
    
    
    

    %% Plot mean & std AUC metric for 4 different NetMethods & 4 different SegMethods for set {rM, sP} parameters
    if(1)
        rM = {'1','3','10'}; %,'5'
        sP = {'0p2','0p2','0p2','0p2'};
        KS = {'sml','mid','lrg'};

        for j = 1:numel(KS)
            for i = 1:numel(rM)

                % LOAD IN MAT FILE FOR DIFFERENT NETWORK METHODS YOU WANT TO COMPARE HERE.
                if(do_AA)
                    AA1 = load([matFilesDir,'AAnrm/AUCdata_AAnrm_rM',rM{i},'_sDInf_sP',sP{i},'_NF_60_0_ks',KS{j},'.mat']);
                end
                if(do_GL)
                    GL1 = load([matFilesDir,'GLnrm/AUCdata_GLnrm_rM',rM{i},'_sDInf_sP',sP{i},'_NF_60_0_ks',KS{j},'.mat']);
                end
                if(do_SK)
                    SK1 = load([matFilesDir,'Mod_SKHAdj/AUCdata_Mod_SKHAdj_rM',rM{i},'_sDInf_sP',sP{i},'_NF_60_0_ks',KS{j},'.mat']);
                end
                if(do_NG)
                    NG1 = load([matFilesDir,'Mod_N&G/AUCdata_Mod_N&G_rM',rM{i},'_sDInf_sP',sP{i},'_NF_60_0_ks',KS{j},'.mat']);
                end
                if(do_Iso)
                    Iso = load([matFilesDir,'IsoDiff/AUCdata_IsoDiff_rM',rM{i},'_NF_60_0_ks',KS{j},'.mat']);
                end

                % CHECK TO SEE IF NUMBER OF IMAGE PATCHES USING EACH METHOD IS THE SAME (SHOULD BE)
    %             x1 = any(AA1.ClusterSizes(:)~=GL1.ClusterSizes(:));
    % %            x2 = any(AA1.ClusterSizes(:)~=NG1.ClusterSizes(:));
    %             x3 = any(AA1.ClusterSizes(:)~=SK1.ClusterSizes(:));
    %             x4 = any(AA1.ClusterSizes(:)~=Iso.ClusterSizes(:)); % first make sure they all have same cluster sizes (they should)
    %             if any([x1,x3,x4])
    % % x2,                
    %                 disp('Something is Wrong. Cluster Sizes not same for different Net Methods.')
    %                 break
    %             end
    
                % DIFFERENT FILTERS TO LOOK AT ONLY SUBSETS OF IMAGE PATCHES.

                if(do_AA)
                    indFilter_AA1 = 1:numel(AA1.GndTruthID);                                             % No Filter: Use all the image patches.
                end
                if(do_GL)
                    indFilter_GL1 = 1:numel(GL1.GndTruthID);
                end
                if(do_SK)
                    indFilter_SK1 = 1:numel(SK1.GndTruthID);
                end
                if(do_NG)
                    indFilter_NG1 = 1:numel(NG1.GndTruthID);
                end
                if(do_Iso)
                    indFilter_Iso = 1:numel(Iso.GndTruthID);
                end


                %indFilter = ( (AA1.ClusterSizes(1,:)>100) &(AA1.ClusterSizes(2,:)>100) );           % both clusters > 100 pixels

                %indFilter = (AA1.im_ROC_AUC_1D_accum<0.7);                                          % Clusters well separated by pixels only.
                %indFilter = ( (AA1.im_ROC_AUC_1D_accum>=0.7) & (AA1.im_ROC_AUC_1D_accum<=0.9) );    % Clusters with intermediate performance on pixels only.
                %indFilter =  (AA1.im_ROC_AUC_1D_accum>0.9);                                         % Clusters well separated by pixels only.


                if(do_AA)
                    AA_diff_KurGS = AA1.kur_ROC_AUC_GS_accum(indFilter_AA1) - AA1.im_ROC_AUC_1D_accum(indFilter_AA1);
                    AA_diff_Kur1D = AA1.kur_ROC_AUC_1D_accum(indFilter_AA1) - AA1.im_ROC_AUC_1D_accum(indFilter_AA1);
                    AA_diff_Ev1 =   AA1.ev1_ROC_AUC_1D_accum(indFilter_AA1) - AA1.im_ROC_AUC_1D_accum(indFilter_AA1);
                    AA_diff_Ev2 =  AA1.ev2o_ROC_AUC_1D_accum(indFilter_AA1) - AA1.im_ROC_AUC_1D_accum(indFilter_AA1);
                    AA_diff_Ev3 =  AA1.ev3o_ROC_AUC_1D_accum(indFilter_AA1) - AA1.im_ROC_AUC_1D_accum(indFilter_AA1);
                end
                %
                if(do_GL)
                    GL_diff_KurGS = GL1.kur_ROC_AUC_GS_accum(indFilter_GL1) - GL1.im_ROC_AUC_1D_accum(indFilter_GL1);
                    GL_diff_Kur1D = GL1.kur_ROC_AUC_1D_accum(indFilter_GL1) - GL1.im_ROC_AUC_1D_accum(indFilter_GL1);
                    GL_diff_Ev1 =   GL1.ev1_ROC_AUC_1D_accum(indFilter_GL1) - GL1.im_ROC_AUC_1D_accum(indFilter_GL1);
                    GL_diff_Ev2 =  GL1.ev2o_ROC_AUC_1D_accum(indFilter_GL1) - GL1.im_ROC_AUC_1D_accum(indFilter_GL1);
                    GL_diff_Ev3 =  GL1.ev3o_ROC_AUC_1D_accum(indFilter_GL1) - GL1.im_ROC_AUC_1D_accum(indFilter_GL1);
                end
                %
                if(do_SK)
                    SK_diff_KurGS = SK1.kur_ROC_AUC_GS_accum(indFilter_SK1) - SK1.im_ROC_AUC_1D_accum(indFilter_SK1);
                    SK_diff_Kur1D = SK1.kur_ROC_AUC_1D_accum(indFilter_SK1) - SK1.im_ROC_AUC_1D_accum(indFilter_SK1);
                    SK_diff_Ev1 =   SK1.ev1_ROC_AUC_1D_accum(indFilter_SK1) - SK1.im_ROC_AUC_1D_accum(indFilter_SK1);
                    SK_diff_Ev2 =  SK1.ev2o_ROC_AUC_1D_accum(indFilter_SK1) - SK1.im_ROC_AUC_1D_accum(indFilter_SK1);
                    SK_diff_Ev3 =  SK1.ev3o_ROC_AUC_1D_accum(indFilter_SK1) - SK1.im_ROC_AUC_1D_accum(indFilter_SK1);
                end
                %
                if(do_NG)
                    NG_diff_KurGS = NG1.kur_ROC_AUC_GS_accum(indFilter_NG1) - NG1.im_ROC_AUC_1D_accum(indFilter_NG1);
                    NG_diff_Kur1D = NG1.kur_ROC_AUC_1D_accum(indFilter_NG1) - NG1.im_ROC_AUC_1D_accum(indFilter_NG1);
                    NG_diff_Ev1 =   NG1.ev1_ROC_AUC_1D_accum(indFilter_NG1) - NG1.im_ROC_AUC_1D_accum(indFilter_NG1);
                    NG_diff_Ev2 =  NG1.ev2o_ROC_AUC_1D_accum(indFilter_NG1) - NG1.im_ROC_AUC_1D_accum(indFilter_NG1);
                    NG_diff_Ev3 =  NG1.ev3o_ROC_AUC_1D_accum(indFilter_NG1) - NG1.im_ROC_AUC_1D_accum(indFilter_NG1);
                end
                %
                if(do_Iso)
                    ID_diff_KurGS = Iso.kur_ROC_AUC_GS_accum(indFilter_Iso) - Iso.im_ROC_AUC_1D_accum(indFilter_Iso);
                    ID_diff_Kur1D = Iso.kur_ROC_AUC_1D_accum(indFilter_Iso) - Iso.im_ROC_AUC_1D_accum(indFilter_Iso);
                    ID_diff_Ev1 =   Iso.ev1_ROC_AUC_1D_accum(indFilter_Iso) - Iso.im_ROC_AUC_1D_accum(indFilter_Iso);
                    ID_diff_Ev2 =  Iso.ev2o_ROC_AUC_1D_accum(indFilter_Iso) - Iso.im_ROC_AUC_1D_accum(indFilter_Iso);
                    ID_diff_Ev3 =  Iso.ev3o_ROC_AUC_1D_accum(indFilter_Iso) - Iso.im_ROC_AUC_1D_accum(indFilter_Iso);
                end
                
                
                
                
                
                
                


%                 % Compute Histograms for each of these diff vectors.
%                 edges = [-1:0.05:1];
%                 %
%                 if(do_AA)
%                     N_AA_Kur_GS = histc(AA_diff_KurGS,edges);
%                     N_AA_Kur_1D = histc(AA_diff_Kur1D,edges);
%                     N_AA_Ev1 = histc(AA_diff_Ev1,edges);
%                     N_AA_Ev2 = histc(AA_diff_Ev2,edges);
%                     N_AA_Ev3 = histc(AA_diff_Ev3,edges);
%                 end
%                 %
%                 if(do_GL)
%                     N_GL_Kur_GS = histc(GL_diff_KurGS,edges);
%                     N_GL_Kur_1D = histc(GL_diff_Kur1D,edges);
%                     N_GL_Ev1 = histc(GL_diff_Ev1,edges);
%                     N_GL_Ev2 = histc(GL_diff_Ev2,edges);
%                     N_GL_Ev3 = histc(GL_diff_Ev3,edges);
%                 end
%                 %
%                 if(do_NG)
%                     N_NG_Kur_GS = histc(NG_diff_KurGS,edges);
%                     N_NG_Kur_1D = histc(NG_diff_Kur1D,edges);
%                     N_NG_Ev1 = histc(NG_diff_Ev1,edges);
%                     N_NG_Ev2 = histc(NG_diff_Ev2,edges);
%                     N_NG_Ev3 = histc(NG_diff_Ev3,edges);
%                 end
%                 %
%                 if(do_SK)
%                     N_SK_Kur_GS = histc(SK_diff_KurGS,edges);
%                     N_SK_Kur_1D = histc(SK_diff_Kur1D,edges);
%                     N_SK_Ev1 = histc(SK_diff_Ev1,edges);
%                     N_SK_Ev2 = histc(SK_diff_Ev2,edges);
%                     N_SK_Ev3 = histc(SK_diff_Ev3,edges);
%                 end
%                 %
%                 if(do_Iso)
%                     N_ID_Kur_GS = histc(ID_diff_KurGS,edges);
%                     N_ID_Kur_1D = histc(ID_diff_Kur1D,edges);
%                     N_ID_Ev1 = histc(ID_diff_Ev1,edges);
%                     N_ID_Ev2 = histc(ID_diff_Ev2,edges);
%                     N_ID_Ev3 = histc(ID_diff_Ev3,edges);
%                 end
% 
% 
%                 % Plot Histograms Side By Side (All have same mass. Better to have mass shifted right.)
%                 h0=figure; 
%                 subplot(511),hold on, 
%                 stem(edges+0.00,N_AA_Kur_GS,'k','LineWidth',2)
%                 stem(edges+0.005,N_GL_Kur_GS,'g','LineWidth',2)
%                 stem(edges+0.01,N_SK_Kur_GS,'r','LineWidth',2)
%     %             stem(edges+0.015,N_NG_Kur_GS,'b','LineWidth',2)
%                 stem(edges+0.02,N_ID_Kur_GS,'c','LineWidth',2)
%                 %legend({'AA','GL','SK','NG','ID'})
%                 title(['rM = ',rM{i}],'FontSize',20,'FontWeight','Bold')
%                 %xlabel('\Delta AUC (Mthd - ImPix)')
%                 ylabel({'Counts','Kur GS'},'FontSize',18,'FontWeight','Bold')
%                 set(gca,'FontSize',16,'FontWeight','Bold','YTick',[])
%                 %
%                 subplot(512),hold on, 
%                 stem(edges+0.00,N_AA_Kur_1D,'k','LineWidth',2)
%                 stem(edges+0.005,N_GL_Kur_1D,'g','LineWidth',2)
%                 stem(edges+0.01,N_SK_Kur_1D,'r','LineWidth',2)
%     %             stem(edges+0.015,N_NG_Kur_1D,'b','LineWidth',2)
%                 stem(edges+0.02,N_ID_Kur_1D,'c','LineWidth',2)
%                 %legend({'AA','GL','SK','NG','ID'})
%                 %title([''],'FontSize',20,'FontWeight','Bold')% : rM = ',rM{i}])
%                 %xlabel('\Delta AUC (Mthd - ImPix)')
%                 ylabel('Kur 1D','FontSize',18,'FontWeight','Bold')
%                 set(gca,'FontSize',16,'FontWeight','Bold','YTick',[])
%                 %
%                 subplot(513),hold on, 
%                 stem(edges+0.00,N_AA_Ev1,'k','LineWidth',2)
%                 stem(edges+0.005,N_GL_Ev1,'g','LineWidth',2)
%                 stem(edges+0.01,N_SK_Ev1,'r','LineWidth',2)
%     %             stem(edges+0.015,N_NG_Ev1,'b','LineWidth',2)
%                 stem(edges+0.02,N_ID_Ev1,'c','LineWidth',2)
%                 %legend({'AA','GL','SK','NG','ID'})
%                 %title([''],'FontSize',20,'FontWeight','Bold')% : rM = ',rM{i}])
%                 %xlabel('\Delta AUC (Mthd - ImPix)')
%                 ylabel('Eig 1','FontSize',18,'FontWeight','Bold')
%                 set(gca,'FontSize',16,'FontWeight','Bold','YTick',[])
%                 %
%                 subplot(514),hold on, 
%                 stem(edges+0.00,N_AA_Ev2,'k','LineWidth',2)
%                 stem(edges+0.005,N_GL_Ev2,'g','LineWidth',2)
%                 stem(edges+0.01,N_SK_Ev2,'r','LineWidth',2)
%     %             stem(edges+0.015,N_NG_Ev2,'b','LineWidth',2)
%                 stem(edges+0.02,N_ID_Ev2,'c','LineWidth',2)
%                 %legend({'AA','GL','SK','NG','ID'})
%                 %title([''],'FontSize',20,'FontWeight','Bold')% : rM = ',rM{i}])
%                 %xlabel('\Delta AUC (Mthd - ImPix)')
%                 ylabel('Eig2','FontSize',18,'FontWeight','Bold')
%                 set(gca,'FontSize',16,'FontWeight','Bold','YTick',[])
%                 %
%                 subplot(515),hold on, 
%                 stem(edges+0.00,N_AA_Ev3,'k','LineWidth',2)
%                 stem(edges+0.005,N_GL_Ev3,'g','LineWidth',2)
%                 stem(edges+0.01,N_SK_Ev3,'r','LineWidth',2)
%     %             stem(edges+0.015,N_NG_Ev3,'b','LineWidth',2)
%                 stem(edges+0.02,N_ID_Ev3,'c','LineWidth',2)
%                 legend({'AA','GL','SK','ID'})
%     % 'NG',            
%                 %title([''],'FontSize',20,'FontWeight','Bold')% : rM = ',rM{i}])
%                 xlabel('\Delta AUC (Mthd - ImPix)','FontSize',18,'FontWeight','Bold')
%                 ylabel('Eig 3','FontSize',18,'FontWeight','Bold')
%                 set(gca,'FontSize',16,'FontWeight','Bold','YTick',[])
% 
%                 saveGoodImg(h0,[imgFilesDir,'Comp_NetSeg_Methods_rM',rM{i},'_sP',sP{i},'_whole_dists_AUCi_sml'],sizeGoodIm)
%                 close(h0)











                if(do_AA)
                    
                    % calculate mean & std of results just computed
                    AUC_mnB(1,1) = mean(AA_diff_KurGS);
                    AUC_mnB(1,2) = mean(AA_diff_Kur1D);
                    AUC_mnB(1,3) = mean(AA_diff_Ev1);
                    AUC_mnB(1,4) = mean(AA_diff_Ev2);
                    AUC_mnB(1,5) = mean(AA_diff_Ev3);
                    %
                    AUC_stdB(1,1) = std(AA_diff_KurGS);
                    AUC_stdB(1,2) = std(AA_diff_Kur1D);
                    AUC_stdB(1,3) = std(AA_diff_Ev1);
                    AUC_stdB(1,4) = std(AA_diff_Ev2);
                    AUC_stdB(1,5) = std(AA_diff_Ev3);  
                    
                    % calculate mean & std of results computed in gen_ROC_AUC
                    AUC_mn(1,1) = AA1.AUC_diff_MethVsImg.KurGS(1);
                    AUC_mn(1,2) = AA1.AUC_diff_MethVsImg.Kur1D(1);
                    AUC_mn(1,3) = AA1.AUC_diff_MethVsImg.Ev1(1);
                    AUC_mn(1,4) = AA1.AUC_diff_MethVsImg.Ev2o(1);
                    AUC_mn(1,5) = AA1.AUC_diff_MethVsImg.Ev3o(1);
                    %
                    AUC_std(1,1) = AA1.AUC_diff_MethVsImg.KurGS(2);
                    AUC_std(1,2) = AA1.AUC_diff_MethVsImg.Kur1D(2);
                    AUC_std(1,3) = AA1.AUC_diff_MethVsImg.Ev1(2);
                    AUC_std(1,4) = AA1.AUC_diff_MethVsImg.Ev2o(2);
                    AUC_std(1,5) = AA1.AUC_diff_MethVsImg.Ev3o(2);
                    
                end
                %
                if(do_GL)
                    
                    AUC_mnB(2,1) = mean(GL_diff_KurGS);
                    AUC_mnB(2,2) = mean(GL_diff_Kur1D);
                    AUC_mnB(2,3) = mean(GL_diff_Ev1);
                    AUC_mnB(2,4) = mean(GL_diff_Ev2);
                    AUC_mnB(2,5) = mean(GL_diff_Ev3);
                    %
                    AUC_stdB(2,1) = std(GL_diff_KurGS);
                    AUC_stdB(2,2) = std(GL_diff_Kur1D);
                    AUC_stdB(2,3) = std(GL_diff_Ev1);
                    AUC_stdB(2,4) = std(GL_diff_Ev2);
                    AUC_stdB(2,5) = std(GL_diff_Ev3);
                    
                    
                    AUC_mn(2,1) = GL1.AUC_diff_MethVsImg.KurGS(1);
                    AUC_mn(2,2) = GL1.AUC_diff_MethVsImg.Kur1D(1);
                    AUC_mn(2,3) = GL1.AUC_diff_MethVsImg.Ev1(1);
                    AUC_mn(2,4) = GL1.AUC_diff_MethVsImg.Ev2o(1);
                    AUC_mn(2,5) = GL1.AUC_diff_MethVsImg.Ev3o(1);
                    %
                    AUC_std(2,1) = GL1.AUC_diff_MethVsImg.KurGS(2);
                    AUC_std(2,2) = GL1.AUC_diff_MethVsImg.Kur1D(2);
                    AUC_std(2,3) = GL1.AUC_diff_MethVsImg.Ev1(2);
                    AUC_std(2,4) = GL1.AUC_diff_MethVsImg.Ev2o(2);
                    AUC_std(2,5) = GL1.AUC_diff_MethVsImg.Ev3o(2);
                    
                end
                %
                if(do_SK)
                    
                    AUC_mnB(3,1) = mean(SK_diff_KurGS);
                    AUC_mnB(3,2) = mean(SK_diff_Kur1D);
                    AUC_mnB(3,3) = mean(SK_diff_Ev1);
                    AUC_mnB(3,4) = mean(SK_diff_Ev2);
                    AUC_mnB(3,5) = mean(SK_diff_Ev3);
                    %
                    AUC_stdB(3,1) = std(SK_diff_KurGS);
                    AUC_stdB(3,2) = std(SK_diff_Kur1D);
                    AUC_stdB(3,3) = std(SK_diff_Ev1);
                    AUC_stdB(3,4) = std(SK_diff_Ev2);
                    AUC_stdB(3,5) = std(SK_diff_Ev3);
                    
                    
                    AUC_mn(3,1) = SK1.AUC_diff_MethVsImg.KurGS(1);
                    AUC_mn(3,2) = SK1.AUC_diff_MethVsImg.Kur1D(1);
                    AUC_mn(3,3) = SK1.AUC_diff_MethVsImg.Ev1(1);
                    AUC_mn(3,4) = SK1.AUC_diff_MethVsImg.Ev2o(1);
                    AUC_mn(3,5) = SK1.AUC_diff_MethVsImg.Ev3o(1);
                    %
                    AUC_std(3,1) = SK1.AUC_diff_MethVsImg.KurGS(2);
                    AUC_std(3,2) = SK1.AUC_diff_MethVsImg.Kur1D(2);
                    AUC_std(3,3) = SK1.AUC_diff_MethVsImg.Ev1(2);
                    AUC_std(3,4) = SK1.AUC_diff_MethVsImg.Ev2o(2);
                    AUC_std(3,5) = SK1.AUC_diff_MethVsImg.Ev3o(2);
                    
                end
                %
                if(do_NG)
                    
                    AUC_mnB(4,1) = mean(NG_diff_KurGS);
                    AUC_mnB(4,2) = mean(NG_diff_Kur1D);
                    AUC_mnB(4,3) = mean(NG_diff_Ev1);
                    AUC_mnB(4,4) = mean(NG_diff_Ev2);
                    AUC_mnB(4,5) = mean(NG_diff_Ev3);
                    %
                    AUC_stdB(4,1) = std(NG_diff_KurGS);
                    AUC_stdB(4,2) = std(NG_diff_Kur1D);
                    AUC_stdB(4,3) = std(NG_diff_Ev1);
                    AUC_stdB(4,4) = std(NG_diff_Ev2);
                    AUC_stdB(4,5) = std(NG_diff_Ev3);
                    
                    
                    AUC_mn(4,1) = NG1.AUC_diff_MethVsImg.KurGS(1);
                    AUC_mn(4,2) = NG1.AUC_diff_MethVsImg.Kur1D(1);
                    AUC_mn(4,3) = NG1.AUC_diff_MethVsImg.Ev1(1);
                    AUC_mn(4,4) = NG1.AUC_diff_MethVsImg.Ev2o(1);
                    AUC_mn(4,5) = NG1.AUC_diff_MethVsImg.Ev3o(1);
                    %
                    AUC_std(4,1) = NG1.AUC_diff_MethVsImg.KurGS(2);
                    AUC_std(4,2) = NG1.AUC_diff_MethVsImg.Kur1D(2);
                    AUC_std(4,3) = NG1.AUC_diff_MethVsImg.Ev1(2);
                    AUC_std(4,4) = NG1.AUC_diff_MethVsImg.Ev2o(2);
                    AUC_std(4,5) = NG1.AUC_diff_MethVsImg.Ev3o(2);
                
                end
                %
                if(do_Iso)
                    
                    AUC_mnB(5,1) = mean(ID_diff_KurGS);
                    AUC_mnB(5,2) = mean(ID_diff_Kur1D);
                    AUC_mnB(5,3) = mean(ID_diff_Ev1);
                    AUC_mnB(5,4) = mean(ID_diff_Ev2);
                    AUC_mnB(5,5) = mean(ID_diff_Ev3);
                    %
                    AUC_stdB(5,1) = std(ID_diff_KurGS);
                    AUC_stdB(5,2) = std(ID_diff_Kur1D);
                    AUC_stdB(5,3) = std(ID_diff_Ev1);
                    AUC_stdB(5,4) = std(ID_diff_Ev2);
                    AUC_stdB(5,5) = std(ID_diff_Ev3);
                    
                    
                    AUC_mn(5,1) = Iso.AUC_diff_MethVsImg.KurGS(1);
                    AUC_mn(5,2) = Iso.AUC_diff_MethVsImg.Kur1D(1);
                    AUC_mn(5,3) = Iso.AUC_diff_MethVsImg.Ev1(1);
                    AUC_mn(5,4) = Iso.AUC_diff_MethVsImg.Ev2o(1);
                    AUC_mn(5,5) = Iso.AUC_diff_MethVsImg.Ev3o(1);
                    %
                    AUC_std(5,1) = Iso.AUC_diff_MethVsImg.KurGS(2);
                    AUC_std(5,2) = Iso.AUC_diff_MethVsImg.Kur1D(2);
                    AUC_std(5,3) = Iso.AUC_diff_MethVsImg.Ev1(2);
                    AUC_std(5,4) = Iso.AUC_diff_MethVsImg.Ev2o(2);
                    AUC_std(5,5) = Iso.AUC_diff_MethVsImg.Ev3o(2);
                    
                end




                
                
                %  
                h1=figure; hold on

                errorbar([1:5]+0.0,AUC_mn(1,:),AUC_std(1,:),'k.','LineWidth',2,'MarkerSize',50)
                text(1.0, 1.1*(AUC_mn(1,1)+AUC_std(1,1)), num2str(AA1.count_isnan_kur_GS))
                errorbar([1:5]+0.1,AUC_mn(2,:),AUC_std(2,:),'g.','LineWidth',2,'MarkerSize',50)
                text(1.1, 1.1*(AUC_mn(2,1)+AUC_std(2,1)), num2str(GL1.count_isnan_kur_GS))
                errorbar([1:5]+0.2,AUC_mn(3,:),AUC_std(3,:),'r.','LineWidth',2,'MarkerSize',50)
                text(1.2, 1.1*(AUC_mn(3,1)+AUC_std(3,1)), num2str(SK1.count_isnan_kur_GS))
                errorbar([1:5]+0.3,AUC_mn(4,:),AUC_std(4,:),'b.','LineWidth',2,'MarkerSize',50)
                text(1.3, 1.1*(AUC_mn(4,1)+AUC_std(4,1)), num2str(NG1.count_isnan_kur_GS))
                errorbar([1:5]+0.4,AUC_mn(5,:),AUC_std(5,:),'c.','LineWidth',2,'MarkerSize',50)
                text(1.4, 1.1*(AUC_mn(5,1)+AUC_std(5,1)), num2str(Iso.count_isnan_kur_GS))
                %
                plot([1 5],[0 0], 'k--')
                legend({'AA','GL','SKH','N&G','IsoDiff'})
                xlabel('Seg Method','FontSize',18,'FontWeight','Bold')
                ylabel('<(AUC mthd) - (AUC img)>','FontSize',18,'FontWeight','Bold')
                title(['rM = ',rM{i},' ; sP = ',sP{i}],'FontSize',20,'FontWeight','Bold')
                ylim([-0.4 0.4])

                set(gca,'XTick',[1,2,3,4,5],'XTickLabel',{'Kur UB','Kur 1D','Ev1','Ev2','Ev3'},'FontSize',16,'FontWeight','Bold')
                saveGoodImg(h1,[imgFilesDir,'Comp_NetSeg_Methods_rM',rM{i},'_sP',sP{i},'_KS',KS{j}],sizeGoodIm)
                close(h1)


                %  
                h2=figure; hold on

                errorbar([1:5]+0.0,AUC_mnB(1,:),AUC_stdB(1,:),'k.','LineWidth',2,'MarkerSize',50)
                text(1.0, 1.1*(AUC_mnB(1,1)+AUC_stdB(1,1)), num2str(AA1.count_isnan_kur_GS))
                errorbar([1:5]+0.1,AUC_mnB(2,:),AUC_stdB(2,:),'g.','LineWidth',2,'MarkerSize',50)
                text(1.1, 1.1*(AUC_mnB(2,1)+AUC_stdB(2,1)), num2str(GL1.count_isnan_kur_GS))
                errorbar([1:5]+0.2,AUC_mnB(3,:),AUC_stdB(3,:),'r.','LineWidth',2,'MarkerSize',50)
                text(1.2, 1.1*(AUC_mnB(3,1)+AUC_stdB(3,1)), num2str(SK1.count_isnan_kur_GS))
                errorbar([1:5]+0.3,AUC_mnB(4,:),AUC_stdB(4,:),'b.','LineWidth',2,'MarkerSize',50)
                text(1.3, 1.1*(AUC_mnB(4,1)+AUC_stdB(4,1)), num2str(NG1.count_isnan_kur_GS))
                errorbar([1:5]+0.4,AUC_mnB(5,:),AUC_stdB(5,:),'c.','LineWidth',2,'MarkerSize',50)
                text(1.4, 1.1*(AUC_mnB(5,1)+AUC_stdB(5,1)), num2str(Iso.count_isnan_kur_GS))
                %
                plot([1 5],[0 0], 'k--')
                legend({'AA','GL','SKH','N&G','IsoDiff'})            
                xlabel('Seg Method','FontSize',18,'FontWeight','Bold')
                ylabel('<(AUC mthd) - (AUC img)>','FontSize',18,'FontWeight','Bold')
                title(['rM = ',rM{i},' ; sP = ',sP{i}],'FontSize',20,'FontWeight','Bold')
                ylim([-0.4 0.4])

                set(gca,'XTick',[1,2,3,4,5],'XTickLabel',{'Kur UB','Kur 1D','Ev1','Ev2','Ev3'},'FontSize',16,'FontWeight','Bold')
                saveGoodImg(h2,[imgFilesDir,'Comp_NetSeg_Methods_rM',rM{i},'_sP',sP{i},'_KS',KS{j},'_B'],sizeGoodIm)
                close(h2)


            end
        end

    end
    
    
 
    %% Plot Results for Isotropic Diffusion with Rmax at different values.
    if(0)
        
        rM = {'1','3','5','10'}; % 

        for i = 1:numel(rM)

            Iso = load([matFilesDir,'IsoDiff/AUCdata_IsoDiff_rM',rM{i},'_NF_60_0.mat']);


            % Pack all this info into a matrix to make a surface plot of it later.
            AUC_mn(i,1) = Iso.AUC_diff_MethVsImg.KurGS(1);
            AUC_mn(i,2) = Iso.AUC_diff_MethVsImg.Kur1D(1);
            AUC_mn(i,3) = Iso.AUC_diff_MethVsImg.Ev1(1);
            AUC_mn(i,4) = Iso.AUC_diff_MethVsImg.Ev2o(1);
            AUC_mn(i,5) = Iso.AUC_diff_MethVsImg.Ev3o(1);
            %
            %
            AUC_std(i,1) = Iso.AUC_diff_MethVsImg.KurGS(2);
            AUC_std(i,2) = Iso.AUC_diff_MethVsImg.Kur1D(2);
            AUC_std(i,3) = Iso.AUC_diff_MethVsImg.Ev1(2);
            AUC_std(i,4) = Iso.AUC_diff_MethVsImg.Ev2o(2);
            AUC_std(i,5) = Iso.AUC_diff_MethVsImg.Ev3o(2);
            
            count_isnan(i) = Iso.count_isnan_kur_GS;

        end


        h2=figure; hold on

        errorbar([1:5]+0.0,AUC_mn(1,:),AUC_std(1,:),'k.','LineWidth',2,'MarkerSize',50)
        errorbar([1:5]+0.1,AUC_mn(2,:),AUC_std(2,:),'g.','LineWidth',2,'MarkerSize',50)
        errorbar([1:5]+0.2,AUC_mn(3,:),AUC_std(3,:),'r.','LineWidth',2,'MarkerSize',50)
        errorbar([1:5]+0.3,AUC_mn(4,:),AUC_std(4,:),'b.','LineWidth',2,'MarkerSize',50)
        % errorbar([1:5]+0.4,AUC_mn(5,:),AUC_std(5,:),'c.','LineWidth',2,'MarkerSize',50)
        %
        for i = 1:numel(rM)
            text(1+0.1*(i-1), 1.1*(AUC_mn(i,1)+AUC_std(i,1)), num2str(count_isnan(i)))
        end
        %
        plot([1 5],[0 0], 'k--')
        legend({'rM=1','3','5','10'}) % ,'\infty'
        xlabel('Seg Method','FontSize',18,'FontWeight','Bold')
        ylabel('<(AUC mthd) - (AUC img)>','FontSize',18,'FontWeight','Bold')
        title(['Isotropic Diffusion Only'],'FontSize',20,'FontWeight','Bold')
        ylim([-0.4 0.4])

        set(gca,'XTick',[1,2,3,4,5],'XTickLabel',{'Kur UB','Kur 1D','Ev1','Ev2','Ev3'},'FontSize',16,'FontWeight','Bold')
        saveGoodImg(h2,[imgFilesDir,'Isotropic_Diffusion_Only'],sizeGoodIm)
        close(h2)

    end
    
    
    disp('Function Successfully Completed.')
    clock
    
    keyboard

end % end main function
