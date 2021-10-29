% This script will loop through all image patch files in a directory and
% will compute the number of pixels in each.  This is important because the
% MLE estimation of parameters for mean and std of clusters performs better
% for larger clusters.  So I want to see what regime of size of clusters I
% am dealing with in the image patches.

[dirPre,sizeGoodIm] = onCluster;

dirSave = [dirPre,'output/truncatedNormalFittingPlots/'];


if ~exist('dirSave','dir')
    mkdir(dirSave)
end

patchDir = [dirPre,'images/BSDS_patch/51x51_ds1/'];

patchFiles = dir([patchDir,'*.mat']);

k=0;

for F = 1:numel(patchFiles)

    load([patchDir,patchFiles(F).name])

    for i = 1:numel(gT) % loop thru ground truths

        C = unique(gT{i});
        for j = 1:numel(C)
            k=k+1;
            s(k) = numel(find(gT{i}==C(j)));

        end
        
    end

    F
    
end



H=figure, hist(s,100)
xlabel('Size of Cluster','FontSize',18,'FontWeight','Bold')
ylabel('# in all im patches & gTs','FontSize',18,'FontWeight','Bold')
title(['Cluster Size Distribution in 51x51 BSDS Image Patches'],'FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',16,'FontWeight','Bold')

saveGoodImg(H,[dirSave,'gT_clusterSizeDist'],sizeGoodIm)
close(H);