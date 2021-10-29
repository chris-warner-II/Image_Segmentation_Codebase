function plot_single_imgPtch_compare_ImPix_Iso_SK_optimized

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
%pbDir1 = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/',meth1,'/pb_png/rM',rM1,'/NF_60_0/ks',ks1,'/'];
%
benchDir2 = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/',meth2,'/benchmark_results/rM',rM2,'/sDInf/sP0p2/NF_60_0/ks',ks1,'/'];
%pbDir2 = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/',meth2,'/pb_png/rM',rM2,'/sDInf/sP0p2/NF_60_0/ks',ks1,'/'];



% directory to ground truth and original image patch dir that includes imPix pb files.
benchDir3 = [dirPre,'images/BSDS_patch/101x101_ds1/benchmark_results/'];
%pbDir3 = [dirPre,'images/BSDS_patch/101x101_ds1/pb_png/'];
%gtDir = [dirPre,'images/BSDS_patch/101x101_ds1/groundTruth/'];


% Directory to put output images into
outDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/imgs/Kur_PIF_Fourier1/compare_single_imgs_Iso_SK_Pix/'];
if ~exist(outDir,'dir')
    mkdir(outDir)
end




% Loop thru all ev1.txt files in the first directory
files = dir([benchDir1,'*_ev1.txt']);

% Preallocate memory for vector fields.
x = zeros(1,numel(files));
y = zeros(1,numel(files));
%
u1 = zeros(1,numel(files));
v1 = zeros(1,numel(files));
delF1 = zeros(1,numel(files));
%
u2 = zeros(1,numel(files));
v2 = zeros(1,numel(files));
delF2 = zeros(1,numel(files));





for i = 1:numel(files)

    i
    imPtchName = files(i).name(1:end-8)


%     % Load in 3 Probabalistic Boundary Image Files.
%     pb1 = double(imread([pbDir1,imPtchName,'.png']))/255;
%     pb2 = double(imread([pbDir2,imPtchName,'.png']))/255;
%     pb3 = double(imread([pbDir3,imPtchName,'.png']))/255;
    
    
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
%     load([gtDir,imPtchName,'.mat']);
%     
%     
%     
%     bD = uint8(zeros(size(pb1)));
%     for j = 1:numel(groundTruth)
%     
%         bD = bD + uint8(groundTruth{j}.Boundaries);
%         
%     end
    


    
    
    % set up (x,y,u,v) data arrays for the quiver function to plot a vector flow field
    
    [bestT1,bestR1,bestP1,bestF1] = maxF(thr,R1,P1);
    [bestT2,bestR2,bestP2,bestF2] = maxF(thr,R2,P2);
    [bestT3,bestR3,bestP3,bestF3] = maxF(thr,R3,P3);
    
    
    
    x(i) = bestR3;
    y(i) = bestP3;
    %
    u1(i) = bestR1-bestR3;
    v1(i) = bestP1-bestP3;
    delF1(i) = bestF1-bestF3;
    %
    u2(i) = bestR2-bestR3;
    v2(i) = bestP2-bestP3;
    delF2(i) = bestF2-bestF3;
    
    
    
    
    
end










h1 = openfig('isoF.fig','reuse'); % open figure
ax1 = gca; % get handle to axes of figure
fig1 = get(ax1,'children'); %get handle to all the children in the figure
H = figure; %create new figure
s1 = subplot(3,2,[1,3]); axis square %create and get handle to the subplot axes
s2 = subplot(2,2,2); axis square
s3 = subplot(2,2,4); axis square

s4 = subplot(6,2,9);
s5 = subplot(6,2,11);

copyobj(fig1,s1); colormap(bone); %copy children to new parent axes i.e. the subplot axes


close(h1);


% Plot vector flow field of how method changed P-R-F from ImPix.
subplot(s1)
hold on
quiver(x,y,u1,v1,2,'r','LineWidth',1.3)
quiver(x,y,u2,v2,2,'b','LineWidth',1.3)
hold off
legend({'human?',['\Delta w/ ',meth1],['\Delta w/ ',meth2]})
xlabel('Recall')
ylabel('Precision')
%
subplot(s2), 
hold on
quiver(zeros(1,numel(files)),zeros(1,numel(files)),u1,v1,'r')
quiver(0,0,mean(u1),mean(v1),'k','LineWidth',2)
scatter(0,0,'k.')
grid on
hold off
axis([min([u1,u2]) max([u1,u2]) min([v1,v2]) max([v1,v2])])
xlabel('\Delta Recall')
ylabel('\Delta Precision')
title(['Change From ImPix w/ ',meth1])
legend({'indiv. imgs','average'})
%
subplot(s3), 
hold on
quiver(zeros(1,numel(files)),zeros(1,numel(files)),u2,v2,'b')
quiver(0,0,mean(u2),mean(v2),'y','LineWidth',2)
scatter(0,0,'y.')
grid on
hold off
axis([min([u1,u2]) max([u1,u2]) min([v1,v2]) max([v1,v2])])
xlabel('\Delta Recall')
ylabel('\Delta Precision')
title(['Change From ImPix w/ ',meth2])
legend({'indiv. imgs','average'})

numbins = 25;
bins = linspace( min([delF1,delF2]), max([delF1,delF2]), numbins );

hitsF1 = hist(delF1,bins);
hitsF2 = hist(delF2,bins);
subplot(s4), bar(bins,hitsF1./sum(hitsF1),'r')
ylabel(meth1)
subplot(s5), bar(bins,hitsF2./sum(hitsF2),'b')
xlabel(['\Delta F-measure from ImPix'])
ylabel(meth2)




saveGoodImg(H,[outDir,'PRdelF_vectorField_',meth1,'_vs_',meth2,'.jpg'],sizeGoodIm)
close(H)



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
