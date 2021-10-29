function plot_single_imgPtch_compare_ImPix_Iso_SK_optimized

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


% Two methods we want to compare.
meth1 = 'IsoDiff'; % will have to make edits to benchDir1 & pbDir1 if not IsoDiff.
rM1 = '3';
ks1 = 'sml';
meth2 = 'Mod_SKHAdj';
rM2 = '3';
ks2 = 'lrg';




% directory structure to pb.png files and benchmark.txt files.
benchDir1 = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/',meth1,'/benchmark_results/rM',rM1,'/NF_60_0/ks',ks1,'/'];
pbDir1 = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/',meth1,'/pb_png/rM',rM1,'/NF_60_0/ks',ks1,'/'];
%
benchDir2 = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/',meth2,'/benchmark_results/rM',rM2,'/sDInf/sP0p2/NF_60_0/ks',ks1,'/'];
pbDir2 = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/',meth2,'/pb_png/rM',rM2,'/sDInf/sP0p2/NF_60_0/ks',ks1,'/'];



% directory to ground truth and original image patch dir that includes imPix pb files.
benchDir3 = [dirPre,'images/BSDS_patch/101x101_ds1/benchmark_results/'];
pbDir3 = [dirPre,'images/BSDS_patch/101x101_ds1/pb_png/'];
gtDir = [dirPre,'images/BSDS_patch/101x101_ds1/groundTruth/'];


% Directory to put output images into
outDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/imgs/Kur_PIF_Fourier1/compare_single_imgs_Iso_SK_Pix/'];
if ~exist(outDir,'dir')
    mkdir(outDir)
end


% Loop thru all ev1.txt files in the first directory
files = dir([benchDir1,'*_ev1.txt']);

for i = 1:numel(files)

    i
    imPtchName = files(i).name(1:end-8)


    % Load in 3 Probabalistic Boundary Image Files.
    pb1 = double(imread([pbDir1,imPtchName,'.png']))/255;
    pb2 = double(imread([pbDir2,imPtchName,'.png']))/255;
    pb3 = double(imread([pbDir3,imPtchName,'.png']))/255;
    
    
    % Load in 1st benchmark text file and compute from it Recall, Precision & F-measure.
    ben1 = dlmread([benchDir1,files(i).name]); % [thresh cntR sumR cntP sumP cntR0 cntP0]
    thr  = ben1(:, 1);
    sumR1 = ben1(:, 3);     % sumR
    sumP1 = ben1(:, 5);     % sumP
    cntR1 = ben1(:, 6);     % cntR0               
    cntP1 = ben1(:, 7);     % cntP0

    R1 = cntR1 ./ (sumR1 + (sumR1==0)); 
    P1 = cntP1 ./ (sumP1 + (sumP1==0)); 
    F1 = fmeasure(R1,P1);

    %[thr, R1, P1, F1]
    
    
    
    
    % Load in 2nd benchmark text file and compute from it Recall, Precision & F-measure.
    ben2 = dlmread([benchDir2,files(i).name]);
    sumR2 = ben2(:, 3);     % sumR
    sumP2 = ben2(:, 5);     % sumP
    cntR2 = ben2(:, 6);     % cntR0           
    cntP2 = ben2(:, 7);     % cntP0

    R2 = cntR2 ./ (sumR2 + (sumR2==0)); 
    P2 = cntP2 ./ (sumP2 + (sumP2==0)); 
    F2 = fmeasure(R2,P2);

    %[thr, R2, P2, F2]
    
    
    
    
    % Load in 3rd benchmark text file and compute from it Recall, Precision & F-measure.
    ben3 = dlmread([benchDir3,files(i).name]);
    sumR3 = ben3(:, 3);     % sumR
    sumP3 = ben3(:, 5);     % sumP
    cntR3 = ben3(:, 6);     % cntR0           
    cntP3 = ben3(:, 7);     % cntP0

    R3 = cntR3 ./ (sumR3 + (sumR3==0)); 
    P3 = cntP3 ./ (sumP3 + (sumP3==0)); 
    F3 = fmeasure(R3,P3);

    %[thr, R3, P3, F3]
    
    
    
    % Load in Image Patch groundtruth matfile because we can plot groundtruth boundaries from it.
    load([gtDir,imPtchName,'.mat']);
    
    
    
    bD = uint8(zeros(size(pb1)));
    for j = 1:numel(groundTruth)
    
        bD = bD + uint8(groundTruth{j}.Boundaries);
        
    end
    

    

    H=figure;
    subplot(2,3,1)
    hold on,
    plot(thr,F1,'b','Linewidth',2)
    plot(thr,F2,'r','Linewidth',2)
    plot(thr,F3,'k--','Linewidth',2)
    legend({[meth1,' (rM',rM1,',ks',ks1,')'],[meth2,' (rM',rM2,',ks',ks2,')'],'imPix'})
    xlabel('threshold')
    ylabel('f-measure')
    axis square
    %
    subplot(2,3,2), imagesc(pb3), colormap(bone), axis square off
    title(['ImPix',' (F=',num2str(max(F3),2),')']),
    %
    subplot(2,3,3), imagesc(pb1), colormap(bone), axis square off
    title(['\color{blue}',meth1,' (F=',num2str(max(F1),2),')',' \color{black}(\DeltaF=',num2str(max(F1-F3),'%+2.2f'),')']),
    %
    subplot(2,3,5), imagesc(pb2), colormap(bone), axis square off
    title(['\color{red}',meth2,' (F=',num2str(max(F2),2),')',' \color{black}(\DeltaF=',num2str(max(F2-F3),'%+2.2f'),')']), 
    %
    subplot(2,3,6), imagesc(bD), colormap(bone), title('Ground Truth'), axis square off
    %
    s = subplot(2,3,4);
    h1 = openfig('isoF.fig','reuse'); % open figure
    ax1 = gca; % get handle to axes of figure
    fig1 = get(ax1,'children'); %get handle to all the children in the figure
    copyobj(fig1,s); %copy children to new parent axes i.e. the subplot axes
    figure(H)
    subplot(s)
    hold on
    plot(R1,P1,'b','Linewidth',2)
    plot(R2,P2,'r','Linewidth',2)
    plot(R3,P3,'k--','Linewidth',2)
    xlabel('Recall')
    ylabel('Precision')
    axis square
    hold off
    close(h1);
        

    
    saveGoodImg(H,[outDir,imPtchName,'.jpg'],sizeGoodIm)
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