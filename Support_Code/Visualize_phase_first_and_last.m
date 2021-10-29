


hmap(1:256,1) = linspace(0,1,256); 
hmap(:,[2 3]) = 0.7; %brightness 
huemap = hsv2rgb(hmap); 
%colormap(huemap)

dirPre ='./output/Kuramoto/NetsFromImgs/BSDS_full_481x321_ds1_blur_sig1/';
%
dataDir = [dirPre,'data/Kur_PIF_Fourier1/Mod_SKHAdj/'];
dataFiles = dir([dataDir,'*kslrg.mat']);
%
imgDir = [dirPre,'imgs/Kur_PIF_Fourier1/Mod_SKHAdj/'];
if ~exist(imgDir,'dir')
    mkdir(imgDir)
end


for i = 1:numel(dataFiles)
    
    
    load([dataDir,dataFiles(i).name])
    


    H = figure;
    subplot(121), imagesc(reshape(metaCluster.phaseAtClk(:,1),netParams.Ndims)), colormap(huemap), freezeColors, axis off, colorbar('Location','SouthOutside')
    subplot(122), imagesc(reshape(metaCluster.phaseAtClk(:,19),netParams.Ndims)), colormap(huemap), freezeColors, axis off, colorbar('Location','SouthOutside')

    saveGoodImg(H,[imgDir,dataFiles(i).name(1:end-4)],[0 0 1 1])
    close(H)



end


