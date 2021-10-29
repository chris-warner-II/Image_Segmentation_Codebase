function visualize_single_imgPtch_compare_all_methods_optimized_Kur


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


    % Optimal Parameters: (recorded by hand) ordered with method.                
                     % maxGT              % meanGT
    % OLD
    %rM_optimal = [ [3, 4, 3, 3, 1] ; [3, 3, 3, 3, 1] ];
    %ks_optimal = [ [3, 3, 3, 3, 1] ; [3, 3, 3, 3, 1] ];
    
    rM_optimal = [ [4, 4, 2, 2, 1] ; [4, 3, 2, 2, 1] ];
    ks_optimal = [ [2, 3, 3, 3, 1] ; [2, 3, 3, 3, 1] ];


    switch which_F_computation
        case 'maxGT'
            k=1;
        case 'meanGT'
            k=2;
    end



    % directory structure to pb.png files and benchmark.txt files.
    Kuramoto_pre = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1_',blurTag,'/data/Kur_PIF_Fourier1/'];
    % Method = Mod_SK with maxGT
    j=1;
    benchDirSK = [Kuramoto_pre,method{j},'/benchmark_results/rM',rM{rM_optimal(k,j)},'/sDInf/sP0p2/NF_60_0/ks',ks{ks_optimal(k,j)},'/'];
    pbDirSK = [Kuramoto_pre,method{j},'/pb_png/rM',rM{rM_optimal(k,j)},'/sDInf/sP0p2/NF_60_0/ks',ks{ks_optimal(k,j)},'/'];
    %
    % Method = Mod_NG with maxGT
    j=2; 
    benchDirNG = [Kuramoto_pre,method{j},'/benchmark_results/rM',rM{rM_optimal(k,j)},'/sDInf/sP0p2/NF_60_0/ks',ks{ks_optimal(k,j)},'/'];
    pbDirNG = [Kuramoto_pre,method{j},'/pb_png/rM',rM{rM_optimal(k,j)},'/sDInf/sP0p2/NF_60_0/ks',ks{ks_optimal(k,j)},'/'];
    %
    % Method = AA with maxGT
    j=3; 
    benchDirAA = [Kuramoto_pre,method{j},'/benchmark_results/rM',rM{rM_optimal(k,j)},'/sDInf/sP0p2/NF_60_0/ks',ks{ks_optimal(k,j)},'/'];
    pbDirAA = [Kuramoto_pre,method{j},'/pb_png/rM',rM{rM_optimal(k,j)},'/sDInf/sP0p2/NF_60_0/ks',ks{ks_optimal(k,j)},'/'];
    % 
    % Method = GL with maxGT
    j=4; 
    benchDirGL = [Kuramoto_pre,method{j},'/benchmark_results/rM',rM{rM_optimal(k,j)},'/sDInf/sP0p2/NF_60_0/ks',ks{ks_optimal(k,j)},'/'];
    pbDirGL = [Kuramoto_pre,method{j},'/pb_png/rM',rM{rM_optimal(k,j)},'/sDInf/sP0p2/NF_60_0/ks',ks{ks_optimal(k,j)},'/'];
    % 
    % Method = IsoDiff with maxGT
    j=5; 
    benchDirIso = [Kuramoto_pre,method{j},'/benchmark_results/rM',rM{rM_optimal(k,j)},'/NF_60_0/ks',ks{ks_optimal(k,j)},'/'];
    pbDirIso = [Kuramoto_pre,method{j},'/pb_png/rM',rM{rM_optimal(k,j)},'/NF_60_0/ks',ks{ks_optimal(k,j)},'/'];






    % directory to ground truth and original image patch dir that includes imPix pb files.
    benchDirPix = [dirPre,'images/BSDS_patch/101x101_ds1/benchmark_results/'];
    pbDirPix = [dirPre,'images/BSDS_patch/101x101_ds1/pb_png/'];
    gtDir = [dirPre,'images/BSDS_patch/101x101_ds1/groundTruth/'];
    imDir = [dirPre,'images/BSDS_patch/101x101_ds1/'];
    %
    benchDirBlur = [dirPre,'images/BSDS_patch/101x101_ds1/',blurTagB,'/benchmark_results/'];
    pbDirBlur = [dirPre,'images/BSDS_patch/101x101_ds1/',blurTagB,'/pb_png/'];







    % Directory to put output images into
    outDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1_',blurTag,'/imgs/Kur_PIF_Fourier1/compare_single_imgPtchs_all_methods_optimized/'];
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
            pbSK = double(imread([pbDirSK,imPtchName(1:end-3),'.png']))/255;
        catch
            pbSK=0;
        end
        %
        try
            pbNG = double(imread([pbDirNG,imPtchName(1:end-3),'.png']))/255;
        catch
            pbNG=0;
        end
        %
        try
            pbAA = double(imread([pbDirAA,imPtchName(1:end-3),'.png']))/255;
        catch
            pbAA=0;
        end
        %
        try
            pbGL = double(imread([pbDirGL,imPtchName(1:end-3),'.png']))/255;
        catch
            pbGL=0;
        end
        %
        try
            pbIso = double(imread([pbDirIso,imPtchName(1:end-3),'.png']))/255;
        catch
            pbIso=0;
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
            benSK = dlmread([benchDirSK,files(i).name]);
            bestSK = find_max_mean_stats(benSK);
            
            

        catch
            benSK=0;
        end
        %
        try
            benNG = dlmread([benchDirNG,files(i).name]);
            bestNG = find_max_mean_stats(benNG);
        catch
            benNG=0;
        end
        %
        try
            benAA = dlmread([benchDirAA,files(i).name]);
            bestAA = find_max_mean_stats(benAA);
        catch
            benAA=0;
        end
        %
        try
            benGL = dlmread([benchDirGL,files(i).name]);
            bestGL = find_max_mean_stats(benGL);
        catch
            benGL=0;
        end
        %
        try
            benIso = dlmread([benchDirIso,files(i).name]);
            bestIso = find_max_mean_stats(benIso);
        catch
            benIso=0;
        end
        %
        try
            benPix = dlmread([benchDirPix,files(i).name]);
            bestPix = find_max_mean_stats(benPix);
        catch
            benPix=0;
        end
        %clr
        
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
        
        ha = tight_subplot(1,10,[0 0.003],[0 0],[0.02 0.005]);


        % (1/10). Plot Precision-Recall 2D Space with Iso-F curves and each method's performance.
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
            %scatter(bestNG.Rmn,bestNG.Pmn,'bv','Linewidth',2,'filled')
            scatter(bestNG.Rmax,bestNG.Pmax,100,'b^','Linewidth',4)
            %
            %scatter(bestAA.Rmn,bestAA.Pmn,'yv','Linewidth',2,'filled')
            scatter(bestAA.Rmax,bestAA.Pmax,100,'y^','Linewidth',4)
            %
            %scatter(bestGL.Rmn,bestGL.Pmn,'gv','Linewidth',2,'filled')
            scatter(bestGL.Rmax,bestGL.Pmax,100,'g^','Linewidth',4)
            %
            %scatter(bestIso.Rmn,bestIso.Pmn,'kv','Linewidth',2,'filled')
            scatter(bestIso.Rmax,bestIso.Pmax,100,'k^','Linewidth',4)
            %
            %scatter(bestPix.Rmn,bestPix.Pmn,'ko','Linewidth',2)
            scatter(bestPix.Rmax,bestPix.Pmax,150,'mo','Linewidth',4)
            %
            %scatter(bestBlur.Rmn,bestBlur.Pmn,'cv','Linewidth',2,'filled')
            scatter(bestBlur.Rmax,bestBlur.Pmax,150,'co','Linewidth',4)
            %
            %scatter(bestSK.Rmn,bestSK.Pmn,'rv','Linewidth',2,'filled')
            scatter(bestSK.Rmax,bestSK.Pmax,100,'r^','Linewidth',4)
            %
            xlabel('R','FontSize',16,'FontWeight','Bold')
            ylabel('P','FontSize',16,'FontWeight','Bold')
            axis square
            hold off
        end
        
        % (2/10).  Show image patch
        subplot(ha(2)), 
        imagesc(im), colormap(bone), freezeColors
        axis square
        set(gca,'XTick',[],'YTick',[])
        title('Img','FontSize',18,'FontWeight','Bold')

        % (3/10). Show just Pixels (spatial gradients)
        subplot(ha(6)), 
        imagesc(pbPix), colormap(bone), freezeColors
        axis square
        set(gca,'XTick',[],'YTick',[])
        title('RawPix','Color','magenta','FontSize',18,'FontWeight','Bold')
        if(0)%bestPix.Fmax > bestBlur.Fmax)
            xlabel(['\color{red}F=',num2str(bestPix.Fmax,'%+5.2f')],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestPix.Fmn,'%+5.2f'),
        else
            xlabel(['F=',num2str(bestPix.Fmax,'%+5.2f')],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestPix.Fmn,'%+5.2f'),
        end

        % (4/10). Show Blurring (spatial gradients)
        subplot(ha(4)), 
        imagesc(pbBlur), colormap(bone), freezeColors
        axis square
        set(gca,'XTick',[],'YTick',[])
        title('GaussRF','Color','cyan','FontSize',18,'FontWeight','Bold')
        if(0)%bestBlur.Fmax > bestPix.Fmax)
            xlabel(['\color{red}F=',num2str(bestBlur.Fmax,'%+5.2f')],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestPix.Fmn,'%+5.2f'),
        else
            xlabel(['F=',num2str(bestBlur.Fmax,'%+5.2f')],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestPix.Fmn,'%+5.2f'),
        end
        
        % (5/10). Show ground truth boundaries.
        subplot(ha(3)), 
        imagesc(bD), colormap(bone), freezeColors
        axis square
        set(gca,'XTick',[],'YTick',[])
        title('gT','FontSize',18,'FontWeight','Bold')
        % xlabel(['[',num2str(F_tot_stats2(1),'%5.2f'),']']) % Dont Know what this is.

        % (6/10). Show Mod SK (phase spatial gradients)
        subplot(ha(5)), 
        imagesc(pbSK), colormap(bone), freezeColors
        axis square
        set(gca,'XTick',[],'YTick',[])
        title('TM','Color','red','FontSize',18,'FontWeight','Bold')
        if(0)%bestSK.Fmax > bestBlur.Fmax)
            xlabel(['\color{red}F=',num2str(bestSK.Fmax,'%+5.2f')],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestSK.Fmn,'%+5.2f'),
        else
            xlabel(['F=',num2str(bestSK.Fmax,'%+5.2f')],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestSK.Fmn,'%+5.2f'),
        end

        % (7/10). Show Mod N&G (phase spatial gradients)
        subplot(ha(7)), 
        imagesc(pbNG), colormap(bone), freezeColors
        axis square
        set(gca,'XTick',[],'YTick',[])
        title('M','Color','blue','FontSize',18,'FontWeight','Bold')
        if(0)%bestNG.Fmax > bestBlur.Fmax)
            xlabel(['\color{red}F=',num2str(bestNG.Fmax,'%+5.2f')],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestNG.Fmn,'%+5.2f'),
        else
            xlabel(['F=',num2str(bestNG.Fmax,'%+5.2f')],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestNG.Fmn,'%+5.2f'),
        end

        % (8/10). Show Avg Assoc (phase spatial gradients)
        subplot(ha(8)), 
        imagesc(pbAA), colormap(bone), freezeColors
        axis square
        set(gca,'XTick',[],'YTick',[])
        title('AA','Color','yellow','FontSize',18,'FontWeight','Bold')
        set(get(gca,'Title'),'Background',[.5 .5 .5]);
        if(0)%bestAA.Fmax > bestBlur.Fmax)
            xlabel(['\color{red}F=',num2str(bestAA.Fmax,'%+5.2f')],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestAA.Fmn,'%+5.2f'),
        else
            xlabel(['F=',num2str(bestAA.Fmax,'%+5.2f')],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestAA.Fmn,'%+5.2f'),
        end

        % (9/10). Show Graph Laplacian (phase spatial gradients)
        subplot(ha(9)), 
        imagesc(pbGL), colormap(bone), freezeColors
        axis square
        set(gca,'XTick',[],'YTick',[])
        title('GL','Color','green','FontSize',18,'FontWeight','Bold')
        if(0)%bestGL.Fmax > bestBlur.Fmax)
            xlabel(['\color{red}F=',num2str(bestGL.Fmax,'%+5.2f')],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestGL.Fmn,'%+5.2f'),
        else
            xlabel(['F=',num2str(bestGL.Fmax,'%+5.2f')],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestGL.Fmn,'%+5.2f'),
        end

        % (10/10). Show Iso Diff (phase spatial gradients)
        subplot(ha(10)), 
        imagesc(pbIso), colormap(bone), freezeColors
        axis square
        set(gca,'XTick',[],'YTick',[])
        title('ISO','Color','black','FontSize',18,'FontWeight','Bold')
        if(0)%bestIso.Fmax > bestBlur.Fmax)
            xlabel(['\color{red}F=',num2str(bestIso.Fmax,'%+5.2f')],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestIso.Fmn,'%+5.2f'),
        else
            xlabel(['F=',num2str(bestIso.Fmax,'%+5.2f')],'FontSize',16,'FontWeight','Bold') % ',',num2str(bestIso.Fmn,'%+5.2f'),
        end

        
        % plot2svg([outDir,'delFmax_SK_blur_',num2str(bestSK.Fmax-bestBlur.Fmax,'%+5.2f'),'_',imPtchName,'.svg'],H0)

        saveGoodImg(H0,[outDir,'delFmax_SK_blur_',num2str(bestSK.Fmax-bestBlur.Fmax,'%+5.2f'),'_',imPtchName,'.jpg'],[0 0.25 1 0.25])
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