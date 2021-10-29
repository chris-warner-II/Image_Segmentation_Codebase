

% WRONG ONE: use visualize_single_imgPtch_compare_all_methods_optimized!!!


function visualize_pb_png_simgle_image_performance_method_optimized

%
% This script / function will loop through 6 different directories, grab
% the corresponding file from each one (some png, some txt) and will
% display the images and plot the data.  
%
% The user inputs the method and its optimized parameters, whether to look 
% at {'best','worst', OR 'mean'} and how many pb_png images to generate/visualize.  
%
% The code will display (1) the original image patch, (2) the pb_png for the method
% with optimized parameters, (3) the pb_png for imPix, (4) the pb_png for Im Blur
% with optimized parameters, (5) the ground truth segment boundaries, and 
% (6) some plot of F-measure or Precision&Recall.





[dirPre,sizeGoodIm] = onCluster;
addpath([dirPre,'images/BSDS_images/BSR/bench/benchmarks/'])







% Two methods we want to compare.
meth1 = 'Mod_SKHAdj'; 
rM1 = '10';
ks1 = 'mid';
% meth2 = 'IsoDiff';
% rM2 = '3';
% ks2 = 'sml';


using_preBlur = 1;
%
if(using_preBlur)
    blurTag = '_blur_sig1';
    blurTagB = 'blur_sz13_sig1';
else
    blurTag = '';
    blurTagB = '';
end


whichF = 'max'; % either: 'mean' OR 'max'.



% directory structure to pb.png files and benchmark.txt files.
benchDirMethod1 = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blurTag,'/data/Kur_PIF_Fourier1/',meth1,'/benchmark_results/rM',rM1,'/sDInf/sP0p2/NF_60_0/ks',ks1,'/'];
pbDirMethod1 = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blurTag,'/data/Kur_PIF_Fourier1/',meth1,'/pb_png/rM',rM1,'/sDInf/sP0p2/NF_60_0/ks',ks1,'/'];



% directory to ground truth and original image patch dir that includes imBlur pb files.
benchDirImBlur = [dirPre,'images/BSDS_patch/101x101_ds1/',blurTagB,'/benchmark_results/'];
pbDirImBlur = [dirPre,'images/BSDS_patch/101x101_ds1/',blurTagB,'/pb_png/'];



% directory to ground truth and original image patch dir that includes imPix pb files.
benchDirImPix = [dirPre,'images/BSDS_patch/101x101_ds1/benchmark_results/'];
pbDirImPix = [dirPre,'images/BSDS_patch/101x101_ds1/pb_png/'];
gtDir = [dirPre,'images/BSDS_patch/101x101_ds1/groundTruth/'];
imDir = [dirPre,'images/BSDS_patch/101x101_ds1/'];


% CHANGE THIS: Directory to put output images into
outDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blurTag,'/imgs/Kur_PIF_Fourier1/visualize_single_image_patches/',meth1,'_rM',rM1,'_ks',ks1,'/F_',whichF,'/'];
if ~exist(outDir,'dir')
    mkdir(outDir)
end


% Loop thru all ev1.txt files in the first directory
files = dir([benchDirMethod1,'*_d4_ev2.txt']);



% TO DO: WANT TO SORT THRU FILES AND ORDER THEM BY F-MEASURE PERFORMANCE SO
% I CAN JUST LOOK AT THE K BEST, OR WORST, OR MEAN.



for i = 1:numel(files)

    i
    imPtchName = files(i).name(1:end-4) % take off .txt

    % Load in 3 Probabalistic Boundary Image Files.
    pb_Iso = double(imread([pbDirMethod1,imPtchName(1:end-7),'.png']))/255;
    
    
    
    
    pbImBlur = double(imread([pbDirImBlur,imPtchName(1:end-7),'.png']))/255;
    pbImPix = double(imread([pbDirImPix,imPtchName(1:end-7),'.png']))/255;
    
    
    % Load in 1st benchmark text file and compute from it Recall, Precision & F-measure.
    benMethod = dlmread([benchDirMethod1,files(i).name]);
    % [thresh, dt, mnR, stdR, Rmax, meanP, stdP, Pmax, meanF, stdF, Fmax, #GT, BestGT]
    %   
    thr            = benMethod(:, 1);
    dt   = benMethod(:, 2);
    %
    Rmn_Method   = benMethod(:, 3);
    Rstd_Method  = benMethod(:, 4);
    Rmax_Method   = benMethod(:, 5);
    Pmn_Method   = benMethod(:, 6);
    Pstd_Method  = benMethod(:, 7);
    Pmax_Method   = benMethod(:, 8);
    Fmn_Method   = benMethod(:, 9);
    Fstd_Method   = benMethod(:, 10);
    Fmax_Method   = benMethod(:, 11);
    gT_num_Method  = benMethod(:, 12);
    gT_best_Method  = benMethod(:, 13);
    %
%     Fmn_MethodB = fmeasure(Rmn_Method,Pmn_Method);
%     Fstd_MethodB = fmeasure(Rstd_Method,Pstd_Method);
%     Fmax_MethodB =  fmeasure(Rmax_Method,Pmax_Method);
    

    
    
    % Load in 2nd benchmark text file and compute from it Recall, Precision & F-measure.
    benImBlur = dlmread([benchDirImBlur,files(i).name]);
    %
    Rmn_ImBlur   = benImBlur(:, 3);
    Rstd_ImBlur  = benImBlur(:, 4);
    Rmax_ImBlur   = benImBlur(:, 5);
    Pmn_ImBlur  = benImBlur(:, 6);
    Pstd_ImBlur  = benImBlur(:, 7);
    Pmax_ImBlur   = benImBlur(:, 8);
    Fmn_ImBlur   = benImBlur(:, 9);
    Fstd_ImBlur   = benImBlur(:, 10);
    Fmax_ImBlur   = benImBlur(:, 11);
    gT_num_ImBlur  = benImBlur(:, 12);
    gT_best_ImBlur  = benImBlur(:, 13);
    %
%     Fmn_ImBlurB = fmeasure(Rmn_ImBlur,Pmn_ImBlur);
%     Fstd_ImBlurB = fmeasure(Rstd_ImBlur,Pstd_ImBlur);
%     Fmax_ImBlurB =  fmeasure(Rmax_ImBlur,Pmax_ImBlur);
    
    
    
    
    % Load in 3rd benchmark text file and compute from it Recall, Precision & F-measure.
    benImPix = dlmread([benchDirImPix,files(i).name]);
    %
    Rmn_ImPix   = benImPix(:, 3);
    Rstd_ImPix  = benImPix(:, 4);
    Rmax_ImPix   = benImPix(:, 5);
    Pmn_ImPix  = benImPix(:, 6);
    Pstd_ImPix  = benImPix(:, 7);
    Pmax_ImPix   = benImPix(:, 8);
    Fmn_ImPix   = benImPix(:, 9);
    Fstd_ImPix   = benImPix(:, 10);
    Fmax_ImPix   = benImPix(:, 11);
    gT_num_ImPix  = benImPix(:, 12);
    gT_best_ImPix  = benImPix(:, 13);
    %
%     Fmn_ImPixB = fmeasure(Rmn_ImPix,Pmn_ImPix);
%     Fstd_ImPixB = fmeasure(Rstd_ImPix,Pstd_ImPix);
%     Fmax_ImPixB =  fmeasure(Rmax_ImPix,Pmax_ImPix);
    
    
    
    
    
    % Compute Change in F-measure from Method to ImPix or ImBlur
    
    delF_mean_imPix  = max(Fmn_Method) - max(Fmn_ImPix);
    delF_mean_imBlur = max(Fmn_Method) - max(Fmn_ImBlur);
    
    delF_max_imPix  = max(Fmax_Method) - max(Fmax_ImPix);
    delF_max_imBlur = max(Fmax_Method) - max(Fmax_ImBlur);
    
    
    
    
    
    
    % Load in Image Patch groundtruth matfile because we can plot groundtruth boundaries from it.
    load([gtDir,imPtchName(1:end-7),'.mat']);
    %
    bD = uint8(zeros(size(pbMethod)));
    for j = 1:numel(groundTruth)
        bD = bD + uint8(groundTruth{j}.Boundaries);
    end
    
    
    % Load in Image Patch matfile so we can display original image patch
    load([imDir,imPtchName(1:end-7),'.mat']);
    
    
    
    
    
    
    % Choose which F-measure to use (old, mean, max)
    switch whichF
        case 'mean'
            F_pix = Fmn_ImPix;
            F_blur = Fmn_ImBlur;
            F_meth = Fmn_Method;
            %
            R_pix = Rmn_ImPix;
            R_blur = Rmn_ImBlur;
            R_meth = Rmn_Method;
            %
            P_pix = Pmn_ImPix;
            P_blur = Pmn_ImBlur;
            P_meth = Pmn_Method;
        case 'max'
            F_pix = Fmax_ImPix;
            F_blur = Fmax_ImBlur;
            F_meth = Fmax_Method;
            %
            R_pix = Rmax_ImPix;
            R_blur = Rmax_ImBlur;
            R_meth = Rmax_Method;
            %
            P_pix = Pmax_ImPix;
            P_blur = Pmax_ImBlur;
            P_meth = Pmax_Method;
    end
    %
    delF_pix = max(F_meth) - max(F_pix);
    delF_blur = max(F_meth) - max(F_blur);
    
    i_pix = find(F_pix==max(F_pix));
    i_blur = find(F_blur==max(F_blur));
    i_meth = find(F_meth==max(F_meth));
    
    

    

    H=figure;
    
    
    % Plot original image patch.
    subplot(2,3,1), imagesc(im)
    axis square off
    title(['Image Patch : ',imPtchName])
    
    
    
    
    
    % I dont want to plot both F-measure and Precision&Recall (do one or the other)
    if(0)
        subplot(2,3,4)
        hold on,
        plot(thr,F_meth,'k--','Linewidth',2)
        plot(thr,F_blur,'r','Linewidth',2)
        plot(thr,F_pix,'b','Linewidth',2)
        legend({meth1,'imBlur','imPix'})
        xlabel('threshold')
        ylabel('f-measure')
        axis square
    end
    %
    
    
    
    % I dont want to plot both F-measure and Precision&Recall (do one or the other)
    if(0)
        s = subplot(2,3,4);
        h1 = openfig('isoF.fig','reuse'); % open figure
        ax1 = gca; % get handle to axes of figure
        fig1 = get(ax1,'children'); %get handle to all the children in the figure
        copyobj(fig1,s); %copy children to new parent axes i.e. the subplot axes
        figure(H)
        subplot(s)
        hold on
        plot(R_meth,P_meth,'k--','Linewidth',2)
        plot(R_meth(i_meth),P_meth(i_meth),'kx','Linewidth',2)
        plot(R_blur,P_blur,'r','Linewidth',2)
        plot(R_blur(i_blur),P_blur(i_blur),'rx','Linewidth',2)
        plot(R_pix,P_pix,'b','Linewidth',2)
        plot(R_pix(i_pix),P_pix(i_pix),'bx','Linewidth',2)
        xlabel('Recall')
        ylabel('Precision')
        axis square
        hold off
        close(h1);
    end
    
    
    
    
    
    
    
    
    
    
    
    subplot(2,3,2), imagesc(pb_Iso), colormap(bone), axis square off
    title([meth1,' (F=',num2str(max(F_meth),2),')']),
    %
    subplot(2,3,3), imagesc(pbImPix), colormap(bone), axis square off
    title(['\color{black}ImPix (F=',num2str(max(F_pix),2),')',' \color{blue}(\DeltaF=',num2str(delF_pix,'%+2.2f'),')']),
    %
    subplot(2,3,5), imagesc(pbImBlur), colormap(bone), axis square off
    title(['\color{black}ImBlur (F=',num2str(max(F_blur),2),')',' \color{red}(\DeltaF=',num2str(delF_blur,'%+2.2f'),')']), 
    %
    subplot(2,3,6), imagesc(bD), colormap(bone), title('Ground Truth'), axis square off
    %
    
    
    
    
    
    

    
    saveGoodImg(H,[outDir,'delF_blur',num2str(delF_blur,'%+2.2f'),'_',imPtchName,'.jpg'],sizeGoodIm)
    close(H)
    
    
end




end




% compute f-measure fromm recall and precision
function [f] = fmeasure(r,p)
f = 2*p.*r./(p+r+((p+r)==0));

end


% interpolate to find best F and coordinates thereof
function [bestT,bestR,bestP,bestF] = maxF(thresh,R,P)
bestT = thresh(1);
bestR = R(1);
bestP = P(1);
bestF = fmeasure(R(1),P(1));
for i = 2:numel(thresh),
  for d = linspace(0,1),
    t = thresh(i)*d + thresh(i-1)*(1-d);
    r = R(i)*d + R(i-1)*(1-d);
    p = P(i)*d + P(i-1)*(1-d);
    f = fmeasure(r,p);
    if f > bestF,
      bestT = t;
      bestR = r;
      bestP = p;
      bestF = f;
    end
  end
end

end