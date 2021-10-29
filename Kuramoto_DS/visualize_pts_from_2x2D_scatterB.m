function visualize_pts_from_2x2D_scatterB(fileType, fileSize, methodType, rM,sD,sP,sW,Ts,Ks) %, phaseMovie)

% 'DivMarg_Scatter_2D_tscale0p1_sP0p1_rM1_sD0p25_sW0_Kscale300.mat'

% This script/function will take in a mat file saved from
% exp_SepVPar_avgOverImgPatches.  That file contains the Divisive Margin
% results for Kuramoto, Eigenvector & Strawman segmentations. We devise a
% metric on the scatter points and sort points in order of this metric.
% Then we visualize the image patch and different segmentationsvis
% corresponding to certain points.
%
% fileType = 'BSDS_patch';
% fileSize = '51x51_ds1';
% methodType = 'Mod_SKHAdj';
% fileSpecifier = '_sP0p2_rM4_sDInf_sW0_Ks300_Ts1' or something like.
%
% phaseMovie = 1; flag to make avi movie of oscillator phase evolution during Kuramoto simulation




% 'KurDivMarg','EigDivMarg','SMDivMarg','ImgPtchID','ImgGtID','fileType','fileSize','Tscale','sigP','Rmax','sigD','sigW','Kscale'



%% Load in mat file and set up directories to different data.
[dirPre,sizeGoodIm] = onCluster;



% % for now, delete this. Saving these variables in the explore_compare_DivMarg code.
% sP = '0p2';
% rM = '4';
% sD = 'Inf';
% Ks = '300';
% Ts = '1';
% sW = '0';



% Specify these with strings or with ''.
rM = num2str(rM); % '1';
sD = num2str(sD); % 'Inf';
sP = num2str(sP); % '0p2';
sP(sP=='.')='p';
Ts = num2str(Ts); % '1';
Ks = num2str(Ks); % '300';
sW = num2str(sW); % '0';


fileSpecifier = ['_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts];



paramsStr = fileSpecifier;
paramsStr(paramsStr=='_')=' ';

fileTypeStr = fileType;
fileTypeStr(fileTypeStr=='_')=' ';

fileSizeStr = fileSize;
fileSizeStr(fileSizeStr=='_')=' ';

% mat file saved from exp_SepVPar_avgOverImgPatches.
matFilesDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/DistPairs/',methodType,'/'];
matFiles = dir([matFilesDir,'DistsPW_Results_allPatches*',fileSpecifier,'*.mat']);


% directory to save output plots into.
scatVizJpgsDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/imgs/DistPairs/',methodType,'/'];
if ~exist(scatVizJpgsDir,'dir')
    mkdir(scatVizJpgsDir)
end





% directory to save output image files into (includes a 2D Histogram of the DivMarg Scatter Plots for Kur_Vs_Straw & Eig1vis_Vs_Straw)
imgKurDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/imgs/Kur_PIF_Fourier1/',methodType,'/'];
if ~exist(imgKurDir,'dir')
    mkdir(imgKurDir)
end
%
imgEigDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/imgs/spectral/',methodType,'/'];
if ~exist(imgEigDir,'dir')
    mkdir(imgEigDir)
end




% directory to image patches extracted from BSDS.
patchesDir = [dirPre,'images/',fileType,'/',fileSize,'/'];

% directory to KurMC files and to Kur_metaSummary files (think I need the KurMC ones)
dirKurMats = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/Kur_PIF_Fourier1/',methodType,'/'];

% directory to the Evecs files
dirEigMats = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/spectral/',methodType,'/'];










% ADD SOMETHING HERE TO LOOK AT AVERAGE PERFORMANCE AND ONLY PLOT THINGS WITH POSITIVE AVERAGE PERFORMANCE.
bestPerf = load([matFilesDir,'best_perf_values_rM_vs_sP.mat']);


    % This used to be a for loop, but it is not any more...
    j = numel(matFiles); % anymore, since I am inputing all parameters, I expect this loop to go thru 1 file only the one with largest # of files in it, or last one.

    disp(['File # ',num2str(j),' / ',num2str(numel(matFiles))])
    matFiles(j).name
    load([matFilesDir,matFiles(j).name]);
    
    
    
    
    
    
    
%     % Need to introduce a canonical ordering for image patches using a
%     % fixed input file (fixed parameter settings) because the ordering of
%     % Divisive Margine computed only on strawman image patch can change
%     % based on rM and on whether you use Linear or Circular
%     Ism_fixed = 1;
% This may be when DistsPW_Results_allPatches file doesnt have all 1500
% files. Check it...





    
    vertPatches = 6;  % How many image patches to look at in plot64patches.m figure
    horzPatches = 12; % How many image patches to look at in plot64patches.m figure
    

    % Kuramoto
    if( bestPerf.maxDMM.Kur(1)>0 | bestPerf.maxDSM.Kur(1)>0 )
        [H, Pmtric, DivMarg] =  scatter_plot_DistPairs(KurDistsPW,SMDistsPW_L,dirKurMats,patchesDir,'BSDS_patch','51x51_ds1','Kur',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,GtNumSegs,1);
        
        
        % [H, Pmtric.Eig3o, DivMarg.Eig3o] =  scatter_plot_DistPairs(Eig3oDistsPW(:,goode),SMDistsPW_L(:,goode),RangeUsedEvecs(3,goode),dirEig,dirImg,fileGeneral,fileSize,'EigVec3',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,GtNumSegs,plot2x2Dscatter)
        
        if(H~=0)
            saveGoodImg(H,[scatVizJpgsDir,'Scat2D_Kur',fileSpecifier],sizeGoodIm)
            close(H)
        end
        %
        [H] =  plot64patches(KurDistsPW,SMDistsPW_L,SMDistsPW_L,dirKurMats,patchesDir,'BSDS_patch','51x51_ds1','Kur',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,vertPatches,horzPatches);
        if(H~=0)
            saveGoodImg(H,[scatVizJpgsDir,'VizPatches_Kur_',num2str(vertPatches),'x',num2str(horzPatches),fileSpecifier],sizeGoodIm)
            close(H)
        end 
        
        
        
        [H, Pmtric, DivMarg] =  scatter_plot_DistPairs(KurDistsPW,SMDistsPW_C,dirKurMats,patchesDir,'BSDS_patch','51x51_ds1','Kur',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,GtNumSegs,1);
        if(H~=0)
            saveGoodImg(H,[scatVizJpgsDir,'Scat2D_KurC',fileSpecifier],sizeGoodIm)
            close(H)
        end
        %
        [H] =  plot64patches(KurDistsPW,SMDistsPW_C,SMDistsPW_L,dirKurMats,patchesDir,'BSDS_patch','51x51_ds1','Kur',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,vertPatches,horzPatches);
        if(H~=0)
            saveGoodImg(H,[scatVizJpgsDir,'VizPatches_KurC',num2str(vertPatches),'x',num2str(horzPatches),fileSpecifier],sizeGoodIm)
            close(H)
        end 
 
    end
    
    
    

    % Eigenvector1
    if( bestPerf.maxDMM.Eig1(1)>0 | bestPerf.maxDSM.Eig1(1)>0 )
        [H, Pmtric, DivMarg] =  scatter_plot_DistPairs(Eig1DistsPW_L,SMDistsPW_L,dirEigMats,patchesDir,'BSDS_patch','51x51_ds1','Eig1',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,GtNumSegs,1);
        if(H~=0)
            saveGoodImg(H,[scatVizJpgsDir,'Scat2D_Eig1',fileSpecifier],sizeGoodIm)
            close(H)
        end
        %
        [H] =  plot64patches(Eig1DistsPW_L,SMDistsPW_L,SMDistsPW_L,dirEigMats,patchesDir,'BSDS_patch','51x51_ds1','Eig1',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,vertPatches,horzPatches);
        if(H~=0)
            saveGoodImg(H,[scatVizJpgsDir,'VizPatches_Eig1_',num2str(vertPatches),'x',num2str(horzPatches),fileSpecifier],sizeGoodIm)
            close(H)
        end 
    end

    
    
    % Eigenvector1 with Vizualiztion Nonlinearity
    if( bestPerf.maxDMM.Eig1v(1)>0 | bestPerf.maxDSM.Eig1v(1)>0 )
        [H, Pmtric, DivMarg] =  scatter_plot_DistPairs(Eig1vDistsPW_L,SMDistsPW_L,dirEigMats,patchesDir,'BSDS_patch','51x51_ds1','Eig1v',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,GtNumSegs,1);
        if(H~=0)
            saveGoodImg(H,[scatVizJpgsDir,'Scat2D_Eig1v',fileSpecifier],sizeGoodIm)
            close(H)
        end
        %
        [H] =  plot64patches(Eig1vDistsPW_L,SMDistsPW_L,SMDistsPW_L,dirEigMats,patchesDir,'BSDS_patch','51x51_ds1','Eig1v',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,vertPatches,horzPatches);
        if(H~=0)
            saveGoodImg(H,[scatVizJpgsDir,'VizPatches_Eig1v_',num2str(vertPatches),'x',num2str(horzPatches),fileSpecifier],sizeGoodIm)
            close(H)
        end
    end

    
    
    
    % Eigenvector2 (only)
    if( bestPerf.maxDMM.Eig2o(1)>0 | bestPerf.maxDSM.Eig2o(1)>0 )
        [H, Pmtric, DivMarg] =  scatter_plot_DistPairs(Eig2oDistsPW,SMDistsPW_L,dirEigMats,patchesDir,'BSDS_patch','51x51_ds1','Eig2',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,GtNumSegs,1);
        if(H~=0)
            saveGoodImg(H,[scatVizJpgsDir,'Scat2D_Eig2',fileSpecifier],sizeGoodIm)
            close(H)
        end
        %
        [H] =  plot64patches(Eig2oDistsPW,SMDistsPW_L,SMDistsPW_L,dirEigMats,patchesDir,'BSDS_patch','51x51_ds1','Eig2',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,vertPatches,horzPatches);
        if(H~=0)
            saveGoodImg(H,[scatVizJpgsDir,'VizPatches_Eig2_',num2str(vertPatches),'x',num2str(horzPatches),fileSpecifier],sizeGoodIm)
            close(H)
        end
    end

    
    
    % Eigenvector2 (only) with Vizualiztion Nonlinearity
    if( bestPerf.maxDMM.Eig2ov(1)>0 | bestPerf.maxDSM.Eig2ov(1)>0 )
        [H, Pmtric, DivMarg] =  scatter_plot_DistPairs(Eig2ovDistsPW,SMDistsPW_L,dirEigMats,patchesDir,'BSDS_patch','51x51_ds1','Eig2v',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,GtNumSegs,1);
        if(H~=0)
            saveGoodImg(H,[scatVizJpgsDir,'Scat2D_Eig2v',fileSpecifier],sizeGoodIm)
            close(H)
        end
        %
        [H] =  plot64patches(Eig2ovDistsPW,SMDistsPW_L,SMDistsPW_L,dirEigMats,patchesDir,'BSDS_patch','51x51_ds1','Eig2v',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,vertPatches,horzPatches);
        if(H~=0)
            saveGoodImg(H,[scatVizJpgsDir,'VizPatches_Eig2v_',num2str(vertPatches),'x',num2str(horzPatches),fileSpecifier],sizeGoodIm)
            close(H)
        end
    end

    
    
    
    % Eigenvector3 (only)
    if( bestPerf.maxDMM.Eig3o(1)>0 | bestPerf.maxDSM.Eig3o(1)>0 )
        [H, Pmtric, DivMarg] =  scatter_plot_DistPairs(Eig3oDistsPW,SMDistsPW_L,dirEigMats,patchesDir,'BSDS_patch','51x51_ds1','Eig3',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,GtNumSegs,1);
        if(H~=0)
            saveGoodImg(H,[scatVizJpgsDir,'Scat2D_Eig3',fileSpecifier],sizeGoodIm)
            close(H)
        end
        %
        [H] =  plot64patches(Eig3oDistsPW,SMDistsPW_L,SMDistsPW_L,dirEigMats,patchesDir,'BSDS_patch','51x51_ds1','Eig3',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,vertPatches,horzPatches);
        if(H~=0)
            saveGoodImg(H,[scatVizJpgsDir,'VizPatches_Eig3_',num2str(vertPatches),'x',num2str(horzPatches),fileSpecifier],sizeGoodIm)
            close(H)
        end
    end
    
    
    

    % Eigenvector3 (only) with Vizualiztion Nonlinearity
    if( bestPerf.maxDMM.Eig3ov(1)>0 | bestPerf.maxDSM.Eig3ov(1)>0 )
        [H, Pmtric, DivMarg] =  scatter_plot_DistPairs(Eig3ovDistsPW,SMDistsPW_L,dirEigMats,patchesDir,'BSDS_patch','51x51_ds1','Eig3v',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,GtNumSegs,1);
        if(H~=0)
            saveGoodImg(H,[scatVizJpgsDir,'Scat2D_Eig3v',fileSpecifier],sizeGoodIm)
            close(H)
        end
        %
        [H] =  plot64patches(Eig3ovDistsPW,SMDistsPW_L,SMDistsPW_L,dirEigMats,patchesDir,'BSDS_patch','51x51_ds1','Eig3v',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,vertPatches,horzPatches);
        if(H~=0)
            saveGoodImg(H,[scatVizJpgsDir,'VizPatches_Eig3v_',num2str(vertPatches),'x',num2str(horzPatches),fileSpecifier],sizeGoodIm)
            close(H)
        end
    end
    
    
    
    
    
    % Eigenvectors 1-2
    if( bestPerf.maxDMM.Eig2(1)>0 | bestPerf.maxDSM.Eig2(1)>0 )
        [H, Pmtric, DivMarg] =  scatter_plot_DistPairs(Eig2DistsPW,SMDistsPW_L,dirEigMats,patchesDir,'BSDS_patch','51x51_ds1','Eig1-2',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,GtNumSegs,1);
        if(H~=0)
            saveGoodImg(H,[scatVizJpgsDir,'Scat2D_Eig1-2',fileSpecifier],sizeGoodIm)
            close(H)
        end
        %
        [H] =  plot64patches(Eig2DistsPW,SMDistsPW_L,SMDistsPW_L,dirEigMats,patchesDir,'BSDS_patch','51x51_ds1','Eig1-2',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,vertPatches,horzPatches);
        if(H~=0)
            saveGoodImg(H,[scatVizJpgsDir,'VizPatches_Eig1-2_',num2str(vertPatches),'x',num2str(horzPatches),fileSpecifier],sizeGoodIm)
            close(H)
        end
    end
    
    
    
    
    % Eigenvectors 1-2 with Vizualiztion Nonlinearity
    if( bestPerf.maxDMM.Eig2v(1)>0 | bestPerf.maxDSM.Eig2v(1)>0 )
        [H, Pmtric, DivMarg] =  scatter_plot_DistPairs(Eig2vDistsPW,SMDistsPW_L,dirEigMats,patchesDir,'BSDS_patch','51x51_ds1','Eig1-2v',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,GtNumSegs,1);
        if(H~=0)
            saveGoodImg(H,[scatVizJpgsDir,'Scat2D_Eig1-2v',fileSpecifier],sizeGoodIm)
            close(H)
        end
        %
        [H] =  plot64patches(Eig2vDistsPW,SMDistsPW_L,SMDistsPW_L,dirEigMats,patchesDir,'BSDS_patch','51x51_ds1','Eig1-2v',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,vertPatches,horzPatches);
        if(H~=0)
            saveGoodImg(H,[scatVizJpgsDir,'VizPatches_Eig1-2v_',num2str(vertPatches),'x',num2str(horzPatches),fileSpecifier],sizeGoodIm)
            close(H)
        end
    end
    
    
    
    
    
    % Eigenvectors 1-3
    if( bestPerf.maxDMM.Eig3(1)>0 | bestPerf.maxDSM.Eig3(1)>0 )
        [H, Pmtric, DivMarg] =  scatter_plot_DistPairs(Eig3DistsPW,SMDistsPW_L,dirEigMats,patchesDir,'BSDS_patch','51x51_ds1','Eig1-3',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,GtNumSegs,1);
        if(H~=0)
            saveGoodImg(H,[scatVizJpgsDir,'Scat2D_Eig1-3',fileSpecifier],sizeGoodIm)
            close(H)
        end
        %
        [H] =  plot64patches(Eig3DistsPW,SMDistsPW_L,SMDistsPW_L,dirEigMats,patchesDir,'BSDS_patch','51x51_ds1','Eig1-3',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,vertPatches,horzPatches);
        if(H~=0)
            saveGoodImg(H,[scatVizJpgsDir,'VizPatches_Eig1-3_',num2str(vertPatches),'x',num2str(horzPatches),fileSpecifier],sizeGoodIm)
            close(H)
        end
    end
    
    
    
    
    % Eigenvectors 1-3 with Vizualiztion Nonlinearity
    if( bestPerf.maxDMM.Eig3v(1)>0 | bestPerf.maxDSM.Eig3v(1)>0 )
        [H, Pmtric, DivMarg] =  scatter_plot_DistPairs(Eig3vDistsPW,SMDistsPW_L,dirEigMats,patchesDir,'BSDS_patch','51x51_ds1','Eig1-3v',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,GtNumSegs,1);
        if(H~=0)
            saveGoodImg(H,[scatVizJpgsDir,'Scat2D_Eig1-3v',fileSpecifier],sizeGoodIm)
            close(H)
        end
        %
        [H] =  plot64patches(Eig3vDistsPW,SMDistsPW_L,SMDistsPW_L,dirEigMats,patchesDir,'BSDS_patch','51x51_ds1','Eig1-3v',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,vertPatches,horzPatches);
        if(H~=0)
            saveGoodImg(H,[scatVizJpgsDir,'VizPatches_Eig1-3v_',num2str(vertPatches),'x',num2str(horzPatches),fileSpecifier],sizeGoodIm)
            close(H)
        end
    end
        
    
    
    
    
%         Eig2DistsPW,  
%         Eig2vDistsPW,  
%         Eig3DistsPW, 
%         Eig3vDistsPW, 
        
        
    
    
    
    disp('FINISHED SAVING PLOTS!')
    clock
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
%     
%     [Y,Itop] = sort(metric_Kur,'descend');                          % best performing patches
%     [Y,Ibot] = sort(metric_Kur,'ascend');                           % worst performing patches
%     [Y,Imn] = sort(abs(metric_Kur - mean(metric_Kur) ),'ascend');   % patches closest to mean performance
%     TBM = {'Itop'}; % 'Ibot', ,'Imn'
% 
%     disp('Mean Metric Value')
%     mean(metric_Kur)
%     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




    %% Plot 2x2D scatter.  Strawman vs. Kuramoto.  And Strawman vs. Eigenvector.
    if(0) % Just for Now.  TURN THIS BACK TO 1 FOR A NICE PLOT !!!
        
        disp('Plotting and Saving 2D Histogram of DivMarg Scatter Plots')
        tic

        Xaxis = linspace(0,1,20);
        Yaxis = linspace(0,1,20);
        
        % Plot
        hSc = figure;   
        %
        hKS = subplot(121); 
        scatter( SMDivMarg_C , KurDivMarg(end,:), 3 ), hold on 
        %Hist2 = hist2d([SMDivMarg_C ; KurDivMarg(end,:)],20,20,[0 1],[0 1]); close
        %imagesc(Xaxis, Yaxis, Hist2./sum(Hist2(:))), hold on, axis xy
        %cb=colorbar('NorthOutside');
        %set(cb,'FontSize',16,'FontWeight','Bold');
        plot([0 1], [0 1],'k--','LineWidth',1)
        axis([0 1 0 1])
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg of Coupled Oscillator Model','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg of Pixel Contrast Model','FontSize',18,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{+}','FontSize',24,'FontWeight','Bold','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{-}','FontSize',24,'FontWeight','Bold','HorizontalAlignment','center')
        text(0.1,0.1,'\color{black}{0}','FontSize',24,'FontWeight','Bold','HorizontalAlignment','center')
        text(0.5,0.95,['\color{red}{<metric>=',num2str(mean(metric_Kur),2),'}'],'FontSize',18,'FontWeight','Bold','HorizontalAlignment','center')
        %
        hES = subplot(122); 
        scatter( SMDivMarg_L , Eig2ovDivMarg_L , 3 ), hold on
        %Hist2 = hist2d([SMDivMarg_L ; Eig2ovDivMarg_L],20,20,[0 1],[0 1]); close
        %imagesc(Xaxis, Yaxis, Hist2./sum(Hist2(:))), hold on, axis xy
        %cb=colorbar('NorthOutside');
        %set(cb,'FontSize',16,'FontWeight','Bold');
        plot([0 1], [0 1],'k--','LineWidth',1)
        axis([0 1 0 1])
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg of Evec1 Viz Computation','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg of Pixel Contrast Model','FontSize',18,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{+}','FontSize',24,'FontWeight','Bold','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{-}','FontSize',24,'FontWeight','Bold','HorizontalAlignment','center')
        text(0.1,0.1,'\color{black}{0}','FontSize',24,'FontWeight','Bold','HorizontalAlignment','center')
        text(0.5,0.95,['\color{red}{<metric>=',num2str(mean(metric_Eig),2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        annotation('textbox', [0 0.9 1 0.1],'String', ...
            [fileTypeStr,' ',fileSizeStr,' : ',paramsStr,' : ','# Pts= ',num2str(numel(ImgPtchID))], ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')
        %
        disp(['Saving ',imgKurDir,'DivMarg_Hist2d_',fileSpecifier])
        saveGoodImg(hSc,[imgKurDir,'DivMarg_Hist2d_',fileSpecifier],sizeGoodIm)
        close(hSc)
        
        toc
        
    end


    %% Sort points according to metric on 2D Kuramoto Strawman Divisive Margin Space.
    
    % if mean(metric_Kur) > some Threshold Value then do this stuff below.
    if( 0 ) % mean(metric_Kur) > 0.08 ) 
        
        num2viz = 3;  % number of points to visualize.
    
        for KK = 1:numel(TBM) % types of points to visualize (top, bottom, mean)
        
            % Keep track of which image patches I visualize in the num2viz loop below. Then I can avoid visualizing the same 
            % one several times when it occurs for a slightly different ground truth.
            alreadyVizzed = {};  % cell array to hold names of image patches we already visualized.
            vizNum = 1;          % counter to index into alreadyVizzed array


            % index into top or bottom or mean scatter points / image patches
            eval(['I=',TBM{KK},';'])
            
            
            
            % File names for image and movie looking at particular examples on single image patches
            figKurName = [imgKurDir,'Kur_DivMarg_Scatter_2D_Ts',Ts,'_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_',TBM{KK},num2str(num2viz)];
            figEigName = [imgEigDir,'Eig_DivMarg_Scatter_2D_Ts',Ts,'_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_',TBM{KK},num2str(num2viz)];
            movKurName = [imgKurDir,'MovPhaseEvoL_Ts',Ts,'_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_K',Ks,'_',TBM{KK},num2str(num2viz)];

            
            
            
            
            % Plot the 2D Histogram of Comparative Divisive Margin Scatter Plot
            H = figure; 
            subplot(121), 
            scatter( SMDivMarg_C , KurDivMarg(end,:),1 ), hold on 
            %Hist2 = hist2d([SMDivMarg_C ; KurDivMarg(end,:)],20,20,[0 1],[0 1]); close
            %imagesc(Xaxis, Yaxis, Hist2./sum(Hist2(:))), hold on, axis xy
            %freezeColors
            %cb=colorbar('NorthOutside');
            %set(cb,'FontSize',16,'FontWeight','Bold');
            %cbfreeze
            plot([0 1], [0 1],'w--','LineWidth',1.5)
            axis([0 1 0 1])
            set(gca,'FontSize',16,'FontWeight','Bold')
            ylabel('DivMarg of Coupled Oscillator Model','FontSize',18,'FontWeight','Bold')
            xlabel('DivMarg of Pixel Contrast Model','FontSize',18,'FontWeight','Bold')
            text(0.9,0.1,'\color{green}{+}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
            text(0.1,0.9,'\color{red}{-}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
            text(0.1,0.1,'\color{black}{0}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
            text(0.5,0.95,['\color{red}{<metric>=',num2str(mean(metric_Kur),2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')

            
            


            for i = 1:num2viz

                

                disp('Loading Mats')

                tic

                % (1). Load in Image Patch to plot it and ground truth.
                ImPtch = load([patchesDir,ImgPtchID{I(i)},'.mat']);


                % (2). Load in KurMC file and plot a.) ending solution and 
                %      b.) Kur DivMarg time evolution along with Strawman DivMarg.
                KurMC = load([dirKurMats,'KurMC_',ImgPtchID{I(i)},'_rM',rM,'_sD',sD,'_sP',sP,'_NF_60_',sW,'_kscale',Ks,'_tscale',Ts,'_runs1.mat']);
                %
                % DivMarg for Kuramoto at all times of simulation.
                DivMarg_K = KurMC.MC{1}.DistAvgPW(:,1,ImgGtID(I(i))) ./ KurMC.MC{1}.DistAvgPW(:,2,ImgGtID(I(i)));
                
                
                
                
                % DivMarg_K(:,g) =  Kur.MC{1}.DistAvgPW(:,1,g) ./ Kur.MC{1}.DistAvgPW(:,2,g); 
                
                
                


                % (3). Load in Evecs file to plot Eigenvector solution and Eig DivMarg.
                % EigDivMarg(I(i)) % Note: Not dealing with Eig right now...
                EigMC = load([dirEigMats,'Evecs_',ImgPtchID{I(i)},'_rM',rM,'_sD',sD,'_sP',sP,'.mat']);
                %
                % DivMarg for Kuramoto at all times of simulation.
                DivMarg_E1v = EigMC.MC.ev1v_l{ImgGtID(I(i))}.DistAvgPW(1) ./ EigMC.MC.ev1v_l{ImgGtID(I(i))}.DistAvgPW(2);
                DivMarg_E1 = EigMC.MC.ev1_l{ImgGtID(I(i))}.DistAvgPW(1) ./ EigMC.MC.ev1_l{ImgGtID(I(i))}.DistAvgPW(2);



                toc
                        
                        
                   
                
                
                
                
                % Plot Original Image Patch and Divisive Margin (Y-axis of 2D Hist)
                figure(H), 
                subplot(num2viz,6,(i-1)*6+4), % plot img
                imagesc(ImPtch.im)  
                axis square
                colormap('bone'), freezeColors
                if(i==1)
                    ylabel('Img','FontSize',20,'FontWeight','Bold')
                end
                set(gca,'Xtick',[],'Ytick',[],'FontSize',16,'FontWeight','Bold')
                xlabel([ num2str(DivMarg_K(1),2),' ',num2str(SMDivMarg_C(I(i)),2) ])
                
                
                % Plot Ground Truth Patch and Divisive Margin (Y-axis of 2D Hist)
                subplot(num2viz,6,(i-1)*6+5), % plot img
                imagesc(ImPtch.gT{ImgGtID(I(i))})  
                axis square
                colormap('jet'), freezeColors
                if(i==1)
                    ylabel('GT','FontSize',20,'FontWeight','Bold')
                end
                set(gca,'Xtick',[],'Ytick',[],'FontSize',16,'FontWeight','Bold')
                


                % Plot Kuramoto Phase Result and Divisive Margin (Y-axis of 2D Hist)
                subplot(num2viz,6,(i-1)*6+6), % plot kur
                KurPhaseFinl = visKurPhase_inHSV(KurMC.netParams.im, reshape(KurMC.metaCluster.phaseAtClk(:,end),KurMC.netParams.Ndims(1),KurMC.netParams.Ndims(2)));
                imagesc(KurPhaseFinl)  
                axis square
                colormap('hsv'), freezeColors
                if(i==1)
                    ylabel('Kur','FontSize',20,'FontWeight','Bold')
                end
                set(gca,'Xtick',[],'Ytick',[],'FontSize',16,'FontWeight','Bold')
                xlabel([ num2str(DivMarg_K(end),2),' ',num2str(KurDivMarg(end,I(i)),2) ])
                
                
%                 % Plot a red x on 2x2 DivMarg Scatter Plot to indicate where this image falls.
%                 subplot(121), hold on
%                 scatter( SMDivMarg_C(I(i)) , KurDivMarg(end,I(i)) , 'rx')
                
                
                
            end
              
            % Save Image with example patches from 2D Hist of Div Marg Scatter.
            disp(['Saving ',figKurName])
            saveGoodImg(H,figKurName,sizeGoodIm)
            close(H)    
            
            
%             keyboard
            
            
            
            
            
            
%             
%             
%                         
% 
%                         % (4). Plot all that shit.
%                         if ~exist([figKurName,'.jpg'],'file')
%                             disp('Plotting')
% 
%                             tic
%                             
%                             % KurPhaseInit = visKurPhase_inBone(KurMC.netParams.im, reshape(KurMC.metaCluster.phaseAtClk(:,1),KurMC.netParams.Ndims(1),KurMC.netParams.Ndims(2)));
%                             
%                             
% 
%  
% 
% 
% 
% 
% 
% 
%                             subplot(242), imagesc(ImPtch.im)                                           % Image Patch
%                             axis square
%                             title('Image Patch / Strawman','FontSize',20,'FontWeight','Bold')
%                             set(gca,'Xtick',[],'Ytick',[],'FontSize',16,'FontWeight','Bold')
%                             %ylabel(['DM = ',num2str(DivMarg_S,3)],'FontSize',18,'FontWeight','Bold')
%                             colormap('jet'), freezeColors
%                             colorbar('SouthOutside'), cbfreeze
% 
% 
%                             subplot(243), imagesc(ImPtch.gT{ImgGtID(I(i))})                            % Ground Truth
%                             axis square
%                             title(['Ground Truth #',num2str(ImgGtID(I(i)))],'FontSize',20,'FontWeight','Bold')
%                             set(gca,'Xtick',[],'Ytick',[],'FontSize',16,'FontWeight','Bold')
%                             colormap('bone'), freezeColors
%                             cb=colorbar('SouthOutside'); set(cb,'Visible', 'off'), cbfreeze
% 
% 
%                             subplot(244), imagesc(ImPtch.imFull), hold on                              % Full Image
%                             %axis square
%                             plot([ImPtch.pach.xpbeg, ImPtch.pach.xpbeg],[ImPtch.pach.ypbeg, ImPtch.pach.ypfin],'r')
%                             plot([ImPtch.pach.xpfin, ImPtch.pach.xpfin],[ImPtch.pach.ypbeg, ImPtch.pach.ypfin],'r')
%                             plot([ImPtch.pach.xpbeg, ImPtch.pach.xpfin],[ImPtch.pach.ypbeg, ImPtch.pach.ypbeg],'r')
%                             plot([ImPtch.pach.xpbeg, ImPtch.pach.xpfin],[ImPtch.pach.ypfin, ImPtch.pach.ypfin],'r') 
%                             x = ImgPtchID{I(i)};
%                             x(x=='_')=' ';
%                             title({'Full Image:',x},'FontSize',20,'FontWeight','Bold')
%                             set(gca,'Xtick',[],'Ytick',[],'FontSize',16,'FontWeight','Bold')
%                             colormap('bone'), freezeColors
% 
% 
%                             subplot(234), hold on
%                             axis square
%                             scatter( SMDivMarg_C,KurDivMarg, 'k.' )                    % all data points
%                             % find other points in the 2x2D scatter plot belonging to same image patch (different ground truth or ...)
%                             xx = find(strcmp(ImgPtchID,ImgPtchID{I(i)}));
%                             % ImgPtchID{xx}
%                             % ImgGtID(xx)
%                             scatter( SMDivMarg_C(xx),KurDivMarg(xx), 200, 'g.' )       % other data points from this same img patch
%                             scatter( SMDivMarg_C(I(i)),KurDivMarg(I(i)), 300, 'r.' )   % this data point (pertaining to kuramoto phase image)
%                             plot([0 1.2], [0 1.2],'r--','LineWidth',1.5)
%                             axis ([0 1.2 0 1.2])
%                             plot([0 1], [1 1],'r--','LineWidth',1.5)
%                             plot([1 1], [0 1],'r--','LineWidth',1.5)
%                             set(gca,'FontSize',16,'FontWeight','Bold')
%                             ylabel('Kuramoto DM','FontSize',18,'FontWeight','Bold')
%                             xlabel('Strawman DM','FontSize',18,'FontWeight','Bold')
%                             %text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','HorizontalAlignment','Right')
%                             %text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold')
%                             %text(0.2,0.25,'\color{cyan}{Neutral}','FontSize',16,'FontWeight','Bold','HorizontalAlignment','Center')
%                             title(['metric: (this pt = ',num2str(metric_Kur(I(i)),2),') & (avg = ',num2str(mean(metric_Kur),2),')'],'FontSize',20,'FontWeight','Bold')
% 
% 
%                             subplot(2,3,[5:6]), hold on                                                % DivMarg Kur & Straw
%                             plot(KurMC.kurParams.tau*(1:numel(DivMarg_K)),DivMarg_K,'b','LineWidth',2)
%                             %plot(KurMC.kurParams.tau*(1:numel(DivMarg_K)),repmat(DivMarg_S,size(DivMarg_K)),'r','LineWidth',2)
%                             xlabel('Time (sec)','FontSize',18,'FontWeight','Bold')
%                             ylabel('Divisive Margin','FontSize',18,'FontWeight','Bold')
%                             legend('Kuramoto Osc','Pixel Strawman','FontSize',16,'FontWeight','Bold')
%                             title(['Parameters: Tscale=',Ts,' | \sigma_P=',sP,' | Rmax=',rM,' | \sigma_D=',sD,' | \sigma_W=',sW,' | Kscale=',Ks],'FontSize',20,'FontWeight','Bold')
%                             set(gca,'FontSize',16,'FontWeight','Bold')
%                             
%                             toc
% 
%                             disp('Saving')
%                             tic
% 
%                             disp(['Saving ',num2str(i),' / ',num2str(num2viz),' : ',figName])
%                             saveGoodImg(H,figName,sizeGoodIm)
%                             close(H)
%                             
%                             toc
%                             
%                         end
%                         
%                         
%                         % Make a video of phase evolution of oscillators
%                         if(phaseMovie & ~exist([movKurName,'.avi'],'file') )
%                             
%                             disp('Making Movie File:')
%                             tic
%                         
%                             makeMovie_phaseEvolution( movKurName, KurMC, ImPtch, ImgPtchID{I(i)}, ImgGtID(I(i)) )
%                             
%                             toc
%                         
%                         end
%                         
% 
%                         % put ImgPtchID into the alreadyVizzed array after we made and saved plot.
%                         alreadyVizzed{vizNum} = ImgPtchID{I(i)};
%                         vizNum = vizNum+1;
% 
%                         
%                     
%                     %end % if the jpg file has been (or hasnt) generated already.
% 
%                 %end % if we havent already visualized this patch

%             %end % loop through number of patches to visualize
        
        end % types of points to visualize (top, bottom, mean)

    end % if the 2D scatter metric exceeds a threshold

% end % loop over all mat files (all parameter value combinations)