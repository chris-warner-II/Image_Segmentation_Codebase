function visualize_single_imgPtch_compare_all_methods_optimized_KurNEig


    % RIGHT ONE !!!

    % NOTE: This is an old script.  It is being replaced or update now in
    % visualize_pb_png_single_image_performance_method_optimized.m
    %
    % This script / function will loop through 6 different directories, grab
    % the corresponding file from each one (some png, some txt) and will
    % display the images and plot the data.  It will output ~500 jpg image
    % files that one can flip through to compare performance of 2 methods (Iso
    % & SK) with optimized parameters and compare them both to ImPix (a
    % strawman). 
    %
    % This can be made more general to input into this function the methods to
    % compare and parameters for each method.  And maybe Ill do that at some
    % point...

    [dirPre,sizeGoodIm] = onCluster;
    addpath([dirPre,'images/BSDS_images/BSR/bench/benchmarks/'])


    using_preBlur = 1;
    %
    if(using_preBlur)
        blurTag = 'blur_sig1';
        blurTagB = 'blur_sz13_sig1';
    else
        blurTag = '';
        blurTagB = '';
    end



    method = {'Mod_SKHAdj','Mod_N&G','AAnrm','GLnrm','IsoDiff'};

    which_F_computation = 'maxGT'; % ,'meanGT'
    relative_to_what = {'relImBlur'};

    rM = {'1','3','5','10'};
    ks =  {'sml','mid','lrg'}; 
    eV =  {'ev1','ev2o','ev3o','ev2','ev3','ev2w','ev3w'}; 


    % Optimal Parameters: (recorded by hand) ordered with method.  (for Kur)              
                     % maxGT              % meanGT
    rM_optimal_K = [ [3, 4, 3, 3, 1] ; [3, 3, 3, 3, 1] ];
    ks_optimal_K = [ [3, 3, 3, 3, 1] ; [3, 3, 3, 3, 1] ];
    
    
    % Optimal Parameters: (recorded by hand) ordered with method. (for Eig)               
                     % maxGT              % meanGT
    rM_optimal_E = [ [4, 4, 4, 4, 4] ; [4, 4, 4, 4, 4] ] ;
    eV_optimal_E = [ [7, 1, 2, 4, 3] ; [7, 1, 2, 4, 3] ];
    


    switch which_F_computation
        case 'maxGT'
            k=1;
        case 'meanGT'
            k=2;
    end



    % directory structure to pb.png files and benchmark.txt files.
    Kuramoto_pre = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1_',blurTag,'/data/Kur_PIF_Fourier1/'];
    spectral_pre = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1_',blurTag,'/data/spectral/'];
    %
    % Method = Mod_SK
    j=1;
    benchDirK_SK = [Kuramoto_pre,method{j},'/benchmark_results/rM',rM{rM_optimal_K(k,j)},'/sDInf/sP0p2/NF_60_0/ks',ks{ks_optimal_K(k,j)},'/'];
    pbDirK_SK = [Kuramoto_pre,method{j},'/pb_png/rM',rM{rM_optimal_K(k,j)},'/sDInf/sP0p2/NF_60_0/ks',ks{ks_optimal_K(k,j)},'/'];
    benchDirE_SK = [spectral_pre,method{j},'/benchmark_results/rM',rM{rM_optimal_E(k,j)},'/sDInf/sP0p2/',eV{eV_optimal_E(k,j)},'/'];
    pbDirE_SK = [spectral_pre,method{j},'/pb_png/rM',rM{rM_optimal_E(k,j)},'/sDInf/sP0p2/',eV{eV_optimal_E(k,j)},'/'];
    %
    % Method = Mod_NG
    j=2; 
    benchDirK_NG = [Kuramoto_pre,method{j},'/benchmark_results/rM',rM{rM_optimal_K(k,j)},'/sDInf/sP0p2/NF_60_0/ks',ks{ks_optimal_K(k,j)},'/'];
    pbDirK_NG = [Kuramoto_pre,method{j},'/pb_png/rM',rM{rM_optimal_K(k,j)},'/sDInf/sP0p2/NF_60_0/ks',ks{ks_optimal_K(k,j)},'/'];
    benchDirE_NG = [spectral_pre,method{j},'/benchmark_results/rM',rM{rM_optimal_E(k,j)},'/sDInf/sP0p2/',eV{eV_optimal_E(k,j)},'/'];
    pbDirE_NG = [spectral_pre,method{j},'/pb_png/rM',rM{rM_optimal_E(k,j)},'/sDInf/sP0p2/',eV{eV_optimal_E(k,j)},'/'];
    %
    % Method = AA
    j=3; 
    benchDirK_AA = [Kuramoto_pre,method{j},'/benchmark_results/rM',rM{rM_optimal_K(k,j)},'/sDInf/sP0p2/NF_60_0/ks',ks{ks_optimal_K(k,j)},'/'];
    pbDirK_AA = [Kuramoto_pre,method{j},'/pb_png/rM',rM{rM_optimal_K(k,j)},'/sDInf/sP0p2/NF_60_0/ks',ks{ks_optimal_K(k,j)},'/'];
    benchDirE_AA = [spectral_pre,method{j},'/benchmark_results/rM',rM{rM_optimal_E(k,j)},'/sDInf/sP0p2/',eV{eV_optimal_E(k,j)},'/'];
    pbDirE_AA = [spectral_pre,method{j},'/pb_png/rM',rM{rM_optimal_E(k,j)},'/sDInf/sP0p2/',eV{eV_optimal_E(k,j)},'/'];
    % 
    % Method = GL
    j=4; 
    benchDirK_GL = [Kuramoto_pre,method{j},'/benchmark_results/rM',rM{rM_optimal_K(k,j)},'/sDInf/sP0p2/NF_60_0/ks',ks{ks_optimal_K(k,j)},'/'];
    pbDirK_GL = [Kuramoto_pre,method{j},'/pb_png/rM',rM{rM_optimal_K(k,j)},'/sDInf/sP0p2/NF_60_0/ks',ks{ks_optimal_K(k,j)},'/'];
    benchDirE_GL = [spectral_pre,method{j},'/benchmark_results/rM',rM{rM_optimal_E(k,j)},'/sDInf/sP0p2/',eV{eV_optimal_E(k,j)},'/'];
    pbDirE_GL = [spectral_pre,method{j},'/pb_png/rM',rM{rM_optimal_E(k,j)},'/sDInf/sP0p2/',eV{eV_optimal_E(k,j)},'/'];
    % 
    % Method = IsoDiff
    j=5; 
    benchDirK_Iso = [Kuramoto_pre,method{j},'/benchmark_results/rM',rM{rM_optimal_K(k,j)},'/NF_60_0/ks',ks{ks_optimal_K(k,j)},'/'];
    pbDirK_Iso = [Kuramoto_pre,method{j},'/pb_png/rM',rM{rM_optimal_K(k,j)},'/NF_60_0/ks',ks{ks_optimal_K(k,j)},'/'];
    benchDirE_Iso = [spectral_pre,method{j},'/benchmark_results/rM',rM{rM_optimal_E(k,j)},'/',eV{eV_optimal_E(k,j)},'/'];
    pbDirE_Iso = [spectral_pre,method{j},'/pb_png/rM',rM{rM_optimal_E(k,j)},'/',eV{eV_optimal_E(k,j)},'/'];


    


    % directory to ground truth and original image patch dir that includes imPix pb files.
    benchDirPix = [dirPre,'images/BSDS_patch/101x101_ds1/benchmark_results/'];
    pbDirPix = [dirPre,'images/BSDS_patch/101x101_ds1/pb_png/'];
    gtDir = [dirPre,'images/BSDS_patch/101x101_ds1/groundTruth/'];
    imDir = [dirPre,'images/BSDS_patch/101x101_ds1/'];
    %
    benchDirBlur = [dirPre,'images/BSDS_patch/101x101_ds1/',blurTagB,'/benchmark_results/'];
    pbDirBlur = [dirPre,'images/BSDS_patch/101x101_ds1/',blurTagB,'/pb_png/'];







    % Directory to put output images into
    outDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1_',blurTag,'/imgs/Kur_PIF_Fourier1/compare_single_imgPtchs_all_methods_optimized_KurNEig/'];
    if ~exist(outDir,'dir')
        mkdir(outDir)
    end











    % Loop thru all ev1.txt files in the first directory
    files = dir([benchDirPix,'*_d4_ev2.txt']);

    for i = 1:numel(files)

        i
        imPtchName = files(i).name(1:end-8)


        % Load in all Probabalistic Boundary Image Files.
        try
            pbK_SK = double(imread([pbDirK_SK,imPtchName(1:end-3),'.png']))/255;
        catch
            pbK_SK=0;
        end
        %
        try
            pbK_NG = double(imread([pbDirK_NG,imPtchName(1:end-3),'.png']))/255;
        catch
            pbK_NG=0;
        end
        %
        try
            pbK_AA = double(imread([pbDirK_AA,imPtchName(1:end-3),'.png']))/255;
        catch
            pbK_AA=0;
        end
        %
        try
            pbK_GL = double(imread([pbDirK_GL,imPtchName(1:end-3),'.png']))/255;
        catch
            pbK_GL=0;
        end
        %
        try
            pbK_Iso = double(imread([pbDirK_Iso,imPtchName(1:end-3),'.png']))/255;
        catch
            pbK_Iso=0;
        end
        %
        try
            pbE_SK = double(imread([pbDirE_SK,imPtchName(1:end-3),'.png']))/255;
        catch
            pbE_SK=0;
        end
        %
        try
            pbE_NG = double(imread([pbDirE_NG,imPtchName(1:end-3),'.png']))/255;
        catch
            pbE_NG=0;
        end
        %
        try
            pbE_AA = double(imread([pbDirE_AA,imPtchName(1:end-3),'.png']))/255;
        catch
            pbE_AA=0;
        end
        %
        try
            pbE_GL = double(imread([pbDirE_GL,imPtchName(1:end-3),'.png']))/255;
        catch
            pbE_GL=0;
        end
        %
        try
            pbE_Iso = double(imread([pbDirE_Iso,imPtchName(1:end-3),'.png']))/255;
        catch
            pbE_Iso=0;
        end
        %
        try
            pbPix = double(imread([pbDirPix,imPtchName(1:end-3),'.png']))/255;
        catch
            pbPix=0;
        end
        %
        try
            pbBlur = double(imread([pbDirBlur,imPtchName(1:end-3),'.png']))/255;
        catch 
            pbBlur=0;
        end






        % Load in all Benchmark F-measure Performance Files.
        try
            benK_SK = dlmread([benchDirK_SK,files(i).name]);
            bestK_SK = find_max_mean_stats(benK_SK);
        catch
            benK_SK=0;
        end
        %
        try
            benK_NG = dlmread([benchDirK_NG,files(i).name]);
            bestK_NG = find_max_mean_stats(benK_NG);
        catch
            benK_NG=0;
        end
        %
        try
            benK_AA = dlmread([benchDirK_AA,files(i).name]);
            bestK_AA = find_max_mean_stats(benK_AA);
        catch
            benK_AA=0;
        end
        %
        try
            benK_GL = dlmread([benchDirK_GL,files(i).name]);
            bestK_GL = find_max_mean_stats(benK_GL);
        catch
            benK_GL=0;
        end
        %
        try
            benK_Iso = dlmread([benchDirK_Iso,files(i).name]);
            bestK_Iso = find_max_mean_stats(benK_Iso);
        catch
            benK_Iso=0;
        end
        %
        try
            benE_SK = dlmread([benchDirE_SK,files(i).name]);
            bestE_SK = find_max_mean_stats(benE_SK);
        catch
            benE_SK=0;
        end
        %
        try
            benE_NG = dlmread([benchDirE_NG,files(i).name]);
            bestE_NG = find_max_mean_stats(benE_NG);
        catch
            benE_NG=0;
        end
        %
        try
            benE_AA = dlmread([benchDirE_AA,files(i).name]);
            bestE_AA = find_max_mean_stats(benE_AA);
        catch
            benE_AA=0;
        end
        %
        try
            benE_GL = dlmread([benchDirE_GL,files(i).name]);
            bestE_GL = find_max_mean_stats(benE_GL);
        catch
            benE_GL=0;
        end
        %
        try
            benE_Iso = dlmread([benchDirE_Iso,files(i).name]);
            bestE_Iso = find_max_mean_stats(benE_Iso);
        catch
            benE_Iso=0;
        end
        %
        try
            benPix = dlmread([benchDirPix,files(i).name]);
            bestPix = find_max_mean_stats(benPix);
        catch
            benPix=0;
        end
        %
        try
            benBlur = dlmread([benchDirBlur,files(i).name]);
            bestBlur = find_max_mean_stats(benBlur);
        catch 
            benBlur=0;
        end







        % Load in Image Patch groundtruth matfile because we can plot groundtruth boundaries from it.
        load([gtDir,imPtchName(1:end-3),'.mat']);
        %
        bD = uint8(zeros(size(pbPix)));
        for j = 1:numel(groundTruth)
            bD = bD + uint8(groundTruth{j}.Boundaries);
        end



        % Load in original image patch mat file.
        load([imDir,imPtchName(1:end-3),'.mat']);



        % Construct the plot comparing F-measure and pb.png images...
        H0 = figure;
        
        ha = tight_subplot(2,8,[0 0.003],[0 0],[0.035 0.005]);

        % % % %

        % (#1). Plot Precision-Recall 2D Space with Iso-F curves and each method's performance.
        if(0)
            s = subplot(ha(1));
            h1 = openfig('isoF.fig','reuse'); % open figure
            ax1 = gca; % get handle to axes of figure
            fig1 = get(ax1,'children'); %get handle to all the children in the figure
            copyobj(fig1,s); %copy children to new parent axes i.e. the subplot axes
            close(h1)
            figure(H0)
            subplot(s)
            hold on
            scatter(0.70,0.90,1300,'w.','Linewidth',2) % cover the green dot in Iso-F contours
            %
            %scatter(bestE_NG.Rmn,bestE_NG.Pmn,'bv','Linewidth',2,'filled')
            scatter(bestE_NG.Rmax,bestE_NG.Pmax,100,'b^','Linewidth',4)
            %
            %scatter(bestE_AA.Rmn,bestE_AA.Pmn,'yv','Linewidth',2,'filled')
            scatter(bestE_AA.Rmax,bestE_AA.Pmax,100,'y^','Linewidth',4)
            %
            %scatter(bestE_GL.Rmn,bestE_GL.Pmn,'gv','Linewidth',2,'filled')
            scatter(bestE_GL.Rmax,bestE_GL.Pmax,100,'g^','Linewidth',4)
            %
            %scatter(bestE_Iso.Rmn,bestE_Iso.Pmn,'kv','Linewidth',2,'filled')
            scatter(bestE_Iso.Rmax,bestE_Iso.Pmax,100,'k^','Linewidth',4)
            %
            %scatter(bestPix.Rmn,bestPix.Pmn,'ko','Linewidth',2)
            scatter(bestPix.Rmax,bestPix.Pmax,150,'mo','Linewidth',4)
            %
            %scatter(bestBlur.Rmn,bestBlur.Pmn,'cv','Linewidth',2,'filled')
            scatter(bestBlur.Rmax,bestBlur.Pmax,150,'co','Linewidth',4)
            %
            %scatter(bestE_SK.Rmn,bestE_SK.Pmn,'rv','Linewidth',2,'filled')
            scatter(bestE_SK.Rmax,bestE_SK.Pmax,100,'r^','Linewidth',4)
            %
            xlabel('R','FontSize',16,'FontWeight','Bold')
            ylabel('P','FontSize',16,'FontWeight','Bold')
            title('Eig','FontSize',18,'FontWeight','Bold')
            axis square
            hold off
        end
        
        
        
        % (#9). Plot Precision-Recall 2D Space with Iso-F curves and each method's performance.
        if(0)
            s = subplot(ha(8+1));
            
            h1 = openfig('isoF.fig','reuse'); % open figure
            ax1 = gca; % get handle to axes of figure
            fig1 = get(ax1,'children'); %get handle to all the children in the figure
            copyobj(fig1,s); %copy children to new parent axes i.e. the subplot axes
            close(h1)
            figure(H0)
            subplot(s)
            hold on
            scatter(0.70,0.90,1300,'w.','Linewidth',2) % cover the green dot in Iso-F contours
            %
            %scatter(bestK_NG.Rmn,bestK_NG.Pmn,'bv','Linewidth',2,'filled')
            scatter(bestK_NG.Rmax,bestK_NG.Pmax,100,'b^','Linewidth',4)
            %
            %scatter(bestK_AA.Rmn,bestK_AA.Pmn,'yv','Linewidth',2,'filled')
            scatter(bestK_AA.Rmax,bestK_AA.Pmax,100,'y^','Linewidth',4)
            %
            %scatter(bestK_GL.Rmn,bestK_GL.Pmn,'gv','Linewidth',2,'filled')
            scatter(bestK_GL.Rmax,bestK_GL.Pmax,100,'g^','Linewidth',4)
            %
            %scatter(bestK_Iso.Rmn,bestK_Iso.Pmn,'kv','Linewidth',2,'filled')
            scatter(bestK_Iso.Rmax,bestK_Iso.Pmax,100,'k^','Linewidth',4)
            %
            %scatter(bestPix.Rmn,bestPix.Pmn,'ko','Linewidth',2)
            scatter(bestPix.Rmax,bestPix.Pmax,150,'mo','Linewidth',4)
            %
            %scatter(bestBlur.Rmn,bestBlur.Pmn,'cv','Linewidth',2,'filled')
            scatter(bestBlur.Rmax,bestBlur.Pmax,150,'co','Linewidth',4)
            %
            %scatter(bestK_SK.Rmn,bestK_SK.Pmn,'rv','Linewidth',2,'filled')
            scatter(bestK_SK.Rmax,bestK_SK.Pmax,100,'r^','Linewidth',4)
            %
            xlabel('R','FontSize',16,'FontWeight','Bold')
            ylabel('P','FontSize',16,'FontWeight','Bold')
            title('Kur','FontSize',18,'FontWeight','Bold')
            axis square
            hold off
        end
        
        % % % %
        
        % (#2). Show just Pixels (spatial gradients)
        subplot(ha(2)), 
        imagesc(pbPix), colormap(bone), freezeColors
        axis square
        set(gca,'XTick',[],'YTick',[])
        title('RawPix','Color','magenta','FontSize',18,'FontWeight','Bold')
        if(0) %bestPix.Fmax > bestBlur.Fmax)
            xlabel(['\color{red}F=[',num2str(bestPix.Fmax,'%5.2f'),',',num2str(bestPix.Fmn,'%5.2f'),']'],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestPix.Fmn,'%+5.2f'),
                    
        else
            xlabel(['F=[',num2str(bestPix.Fmax,'%5.2f'),',',num2str(bestPix.Fmn,'%5.2f'),']'],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestPix.Fmn,'%+5.2f'),
                    
        end

        % (#10). Show Blurring (spatial gradients)
        subplot(ha(8+2)), 
        imagesc(pbBlur), colormap(bone), freezeColors
        axis square
        set(gca,'XTick',[],'YTick',[])
        title('GaussRF','Color','cyan','FontSize',18,'FontWeight','Bold')
        if(0) %bestBlur.Fmax > bestPix.Fmax)
            xlabel(['\color{red}F=[',num2str(bestBlur.Fmax,'%5.2f'),',',num2str(bestBlur.Fmn,'%5.2f'),']'],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestPix.Fmn,'%+5.2f'),
        else
            xlabel(['F=[',num2str(bestBlur.Fmax,'%5.2f'),',',num2str(bestBlur.Fmn,'%5.2f'),']'],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestPix.Fmn,'%+5.2f'),
        end
        
        % % % %
        
        % (#3).  Show image patch
        subplot(ha(3)), 
        imagesc(im), colormap(bone), freezeColors
        axis square
        set(gca,'XTick',[],'YTick',[])
        title('Img','FontSize',18,'FontWeight','Bold')
        
        % (#11). Show ground truth boundaries.
        subplot(ha(8+3)), 
        imagesc(bD), colormap(bone), freezeColors
        axis square
        set(gca,'XTick',[],'YTick',[])
        title('gT','FontSize',18,'FontWeight','Bold')
        % xlabel(['[',num2str(F_tot_stats2(1),'%5.2f'),']']) % Dont Know what this is.

        % % % %

        % (#4). Show Mod SK (phase spatial gradients)
        subplot(ha(4)), 
        imagesc(pbE_SK), colormap(bone), freezeColors
        axis square
        set(gca,'XTick',[],'YTick',[])
        title('TM','Color','red','FontSize',18,'FontWeight','Bold')
        if(0) %bestE_SK.Fmax > bestBlur.Fmax)
            xlabel(['\color{red}F=[',num2str(bestE_SK.Fmax,'%5.2f'),',',num2str(bestE_SK.Fmn,'%5.2f'),']'],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestSK.Fmn,'%+5.2f'),
        else
            xlabel(['F=[',num2str(bestE_SK.Fmax,'%5.2f'),',',num2str(bestE_SK.Fmn,'%5.2f'),']'],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestSK.Fmn,'%+5.2f'),
        end

        % (#12). Show Mod SK (phase spatial gradients)
        subplot(ha(8+4)), 
        imagesc(pbK_SK), colormap(bone), freezeColors
        axis square
        set(gca,'XTick',[],'YTick',[])
        %title('TM','Color','red','FontSize',18,'FontWeight','Bold')
        if(0) %bestK_SK.Fmax > bestBlur.Fmax)
            xlabel(['\color{red}F=[',num2str(bestK_SK.Fmax,'%5.2f'),',',num2str(bestK_SK.Fmn,'%5.2f'),']'],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestSK.Fmn,'%+5.2f'),
        else
            xlabel(['F=[',num2str(bestK_SK.Fmax,'%5.2f'),',',num2str(bestK_SK.Fmn,'%5.2f'),']'],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestSK.Fmn,'%+5.2f'),
        end
        
        % % % %

        % (#5). Show Mod N&G (phase spatial gradients)
        subplot(ha(5)), 
        imagesc(pbE_NG), colormap(bone), freezeColors
        axis square
        set(gca,'XTick',[],'YTick',[])
        title('M','Color','blue','FontSize',18,'FontWeight','Bold')
        if(0) %bestE_NG.Fmax > bestBlur.Fmax)
            xlabel(['\color{red}F=[',num2str(bestE_NG.Fmax,'%5.2f'),',',num2str(bestE_NG.Fmn,'%5.2f'),']'],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestNG.Fmn,'%+5.2f'),
        else
            xlabel(['F=[',num2str(bestE_NG.Fmax,'%5.2f'),',',num2str(bestE_NG.Fmn,'%5.2f'),']'],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestNG.Fmn,'%+5.2f'),
        end
        
        % (#13). Show Mod N&G (phase spatial gradients)
        subplot(ha(8+5)), 
        imagesc(pbK_NG), colormap(bone), freezeColors
        axis square
        set(gca,'XTick',[],'YTick',[])
        %title('M','Color','blue','FontSize',18,'FontWeight','Bold')
        if(0) %bestK_NG.Fmax > bestBlur.Fmax)
            xlabel(['\color{red}F=[',num2str(bestK_NG.Fmax,'%5.2f'),',',num2str(bestK_NG.Fmn,'%5.2f'),']'],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestNG.Fmn,'%+5.2f'),
        else
            xlabel(['F=[',num2str(bestK_NG.Fmax,'%5.2f'),',',num2str(bestK_NG.Fmn,'%5.2f'),']'],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestNG.Fmn,'%+5.2f'),
        end
        
        
        % % % %

        % (#6). Show Avg Assoc (phase spatial gradients)
        subplot(ha(6)), 
        imagesc(pbE_AA), colormap(bone), freezeColors
        axis square
        set(gca,'XTick',[],'YTick',[])
        title('AA','Color','yellow','FontSize',18,'FontWeight','Bold')
        set(get(gca,'Title'),'Background',[.5 .5 .5]);
        if(0) %bestE_AA.Fmax > bestBlur.Fmax)
            xlabel(['\color{red}F=[',num2str(bestE_AA.Fmax,'%5.2f'),',',num2str(bestE_AA.Fmn,'%5.2f'),']'],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestAA.Fmn,'%+5.2f'),
        else
            xlabel(['F=[',num2str(bestE_AA.Fmax,'%5.2f'),',',num2str(bestE_AA.Fmn,'%5.2f'),']'],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestAA.Fmn,'%+5.2f'),
        end

        % (#14). Show Avg Assoc (phase spatial gradients)
        subplot(ha(8+6)), 
        imagesc(pbK_AA), colormap(bone), freezeColors
        axis square
        set(gca,'XTick',[],'YTick',[])
        %title('AA','Color','yellow','FontSize',18,'FontWeight','Bold')
        set(get(gca,'Title'),'Background',[.5 .5 .5]);
        if(0) %bestK_AA.Fmax > bestBlur.Fmax)
            xlabel(['\color{red}F=[',num2str(bestK_AA.Fmax,'%5.2f'),',',num2str(bestK_AA.Fmn,'%5.2f'),']'],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestAA.Fmn,'%+5.2f'),
        else
            xlabel(['F=[',num2str(bestK_AA.Fmax,'%5.2f'),',',num2str(bestK_AA.Fmn,'%5.2f'),']'],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestAA.Fmn,'%+5.2f'),
        end
        
        % % % %
        
        % (#7). Show Graph Laplacian (phase spatial gradients)
        subplot(ha(7)), 
        imagesc(pbE_GL), colormap(bone), freezeColors
        axis square
        set(gca,'XTick',[],'YTick',[])
        title('GL','Color','green','FontSize',18,'FontWeight','Bold')
        if(0) %bestE_GL.Fmax > bestBlur.Fmax)
            xlabel(['\color{red}F=[',num2str(bestE_GL.Fmax,'%5.2f'),',',num2str(bestE_GL.Fmn,'%5.2f'),']'],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestGL.Fmn,'%+5.2f'),
        else
            xlabel(['F=[',num2str(bestE_GL.Fmax,'%5.2f'),',',num2str(bestE_GL.Fmn,'%5.2f'),']'],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestGL.Fmn,'%+5.2f'),
        end
        
        % (#15). Show Graph Laplacian (phase spatial gradients)
        subplot(ha(8+7)), 
        imagesc(pbK_GL), colormap(bone), freezeColors
        axis square
        set(gca,'XTick',[],'YTick',[])
        %title('GL','Color','green','FontSize',18,'FontWeight','Bold')
        if(0) %bestK_GL.Fmax > bestBlur.Fmax)
            xlabel(['\color{red}F=[',num2str(bestK_GL.Fmax,'%5.2f'),',',num2str(bestK_GL.Fmn,'%5.2f'),']'],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestGL.Fmn,'%+5.2f'),
        else
            xlabel(['F=[',num2str(bestK_GL.Fmax,'%5.2f'),',',num2str(bestK_GL.Fmn,'%5.2f'),']'],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestGL.Fmn,'%+5.2f'),
        end
        
        % % % %

        % (#8). Show Iso Diff (phase spatial gradients)
        subplot(ha(8)), 
        imagesc(pbE_Iso), colormap(bone), freezeColors
        axis square
        set(gca,'XTick',[],'YTick',[])
        title('ISO','Color','black','FontSize',18,'FontWeight','Bold')
        if(0) %bestE_Iso.Fmax > bestBlur.Fmax)
            xlabel(['\color{red}F=[',num2str(bestE_Iso.Fmax,'%5.2f'),',',num2str(bestE_Iso.Fmn,'%5.2f'),']'],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestIso.Fmn,'%+5.2f'),
        else
            xlabel(['F=[',num2str(bestE_Iso.Fmax,'%5.2f'),',',num2str(bestE_Iso.Fmn,'%5.2f'),']'],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestIso.Fmn,'%+5.2f'),
        end

        % (#16). Show Iso Diff (phase spatial gradients)
        subplot(ha(8+8)), 
        imagesc(pbK_Iso), colormap(bone), freezeColors
        axis square
        set(gca,'XTick',[],'YTick',[])
        %title('ISO','Color','black','FontSize',18,'FontWeight','Bold')
        if(0) %bestK_Iso.Fmax > bestBlur.Fmax)
            xlabel(['\color{red}F=[',num2str(bestK_Iso.Fmax,'%5.2f'),',',num2str(bestK_Iso.Fmn,'%5.2f'),']'],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestIso.Fmn,'%+5.2f'),
        else
            xlabel(['F=[',num2str(bestK_Iso.Fmax,'%5.2f'),',',num2str(bestK_Iso.Fmn,'%5.2f'),']'],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestIso.Fmn,'%+5.2f'),
        end
        
        % % % %
        
        % plot2svg([outDir,'delFmax_SK_blur_',num2str(bestK_SK.Fmax-bestBlur.Fmax,'%+5.2f'),'_',imPtchName,'.svg'],H0)

        saveGoodImg(H0,[outDir,'delFmax_SK_blur_',num2str(bestK_SK.Fmax-bestBlur.Fmax,'%+5.2f'),'_',imPtchName,'.jpg'],[0 0.25 1 0.55])
        close(H0)





        %%











%         H=figure;
%         subplot(2,3,1)
%         hold on,
%         plot(thr,F1,'b','Linewidth',2)
%         plot(thr,F2,'r','Linewidth',2)
%         plot(thr,F3,'k--','Linewidth',2)
%         legend({[meth1,' (rM',rM1,',ks',ks1,')'],[meth2,' (rM',rM2,',ks',ks2,')'],'imPix'})
%         xlabel('threshold')
%         ylabel('f-measure')
%         axis square
%         %
%         subplot(2,3,2), imagesc(pb3), colormap(bone), axis square off
%         title(['ImPix',' (F=',num2str(max(F3),2),')']),
%         %
%         subplot(2,3,3), imagesc(pb1), colormap(bone), axis square off
%         title(['\color{blue}',meth1,' (F=',num2str(max(F1),2),')',' \color{black}(\DeltaF=',num2str(max(F1-F3),'%+2.2f'),')']),
%         %
%         subplot(2,3,5), imagesc(pb2), colormap(bone), axis square off
%         title(['\color{red}',meth2,' (F=',num2str(max(F2),2),')',' \color{black}(\DeltaF=',num2str(max(F2-F3),'%+2.2f'),')']), 
%         %
%         subplot(2,3,6), imagesc(bD), colormap(bone), title('Ground Truth'), axis square off
%         %
%         s = subplot(2,3,4);
%         h1 = openfig('isoF.fig','reuse'); % open figure
%         ax1 = gca; % get handle to axes of figure
%         fig1 = get(ax1,'children'); %get handle to all the children in the figure
%         copyobj(fig1,s); %copy children to new parent axes i.e. the subplot axes
%         figure(H)
%         subplot(s)
%         hold on
%         plot(R1,P1,'b','Linewidth',2)
%         plot(R2,P2,'r','Linewidth',2)
%         plot(R3,P3,'k--','Linewidth',2)
%         xlabel('Recall')
%         ylabel('Precision')
%         axis square
%         hold off
%         close(h1);
% 
% 
% 
%         saveGoodImg(H,[outDir,imPtchName,'.jpg'],sizeGoodIm)
%         close(H)


    end





end

% compute f-measure fromm recall and precision
function [f] = fmeasure(r,p)

    f = 2*p.*r./(p+r+((p+r)==0));

end



function [best] = find_max_mean_stats(ben)

    % [thresh, dt, mnR, stdR, Rmax, meanP, stdP, Pmax, meanF, stdF, Fmax, #GT, BestGT]
    thr     = ben(:, 1);
    dt      = ben(:, 2);
    %
    Rmn     = ben(:, 3);
    Rstd    = ben(:, 4);
    Rmax    = ben(:, 5);
    Pmn     = ben(:, 6);
    Pstd    = ben(:, 7);
    Pmax    = ben(:, 8);
    Fmn     = ben(:, 9);
    Fstd    = ben(:, 10);
    Fmax    = ben(:, 11);
    gT_num  = ben(:, 12);
    gT_best = ben(:, 13);

     %
    ind_mean = find(Fmn == max(Fmn));
    ind_mean = ind_mean(1);
    best.Fmn = Fmn(ind_mean); 
    best.Rmn = Rmn(ind_mean);
    best.Pmn = Pmn(ind_mean);
    best.Fstd = Fstd(ind_mean); 
    best.Rstd = Rstd(ind_mean);
    best.Pstd = Pstd(ind_mean);
    %
    ind_max = find(Fmax == max(Fmax));
    ind_max = ind_max(1);
    best.Fmax = Fmax(ind_max); 
    best.Rmax = Rmax(ind_max);
    best.Pmax = Pmax(ind_max);

end




% % interpolate to find best F and coordinates thereof
% function [bestT,bestR,bestP,bestF] = maxF(thresh,R,P)
% bestT = thresh(1);
% bestR = R(1);
% bestP = P(1);
% bestF = fmeasure(R(1),P(1));
% for i = 2:numel(thresh),
%   for d = linspace(0,1),
%     t = thresh(i)*d + thresh(i-1)*(1-d);
%     r = R(i)*d + R(i-1)*(1-d);
%     p = P(i)*d + P(i-1)*(1-d);
%     f = fmeasure(r,p);
%     if f > bestF,
%       bestT = t;
%       bestR = r;
%       bestP = p;
%       bestF = f;
%     end
%   end
% end
% 
% end