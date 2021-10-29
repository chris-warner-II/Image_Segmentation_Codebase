% This function will take all the different patches extracted from a single
% image and will stitch them back together into one large image.

% 


[dirPre,sizeGoodIm] = onCluster;


% Directory to image patch files that tile the image.
dirTile = [dirPre,'images/BSDS_tile/51x51_ds1/'];


% Directory to KurMC files output from ImgSeg jobs.  They contain Kuramoto oscillator phases as function of time.
dirKur = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_tile_51x51_ds1/data/Kur_PIF_Fourier1/Mod_SKHAdj/'];

BSDSimage = '202012';

% Loop through all files in dirKur.  Grab them.  Use info in tile files to place them in larger picture.
filesKur = dir([dirKur,'KurMC_',BSDSimage,'*sDInf*.mat']);


% matrix of zeros that I will fill with the phase info from running simulation on patches
imKur = zeros(321,481);


for i  = 1: numel(filesKur)
    
    % load in KurMC file
    load([dirKur,filesKur(i).name])
    
    % Get patch # from file name so I can load in proper tile file
    x = strfind(filesKur(i).name,'ptch')+4;
    y = strfind(filesKur(i).name,'_rM')-1;
    patchTileStr = filesKur(i).name(x:y);
    patchTileNum = str2double(patchTileStr);
    
    % Load tile file
    load([dirTile,BSDSimage,'_ptch',patchTileStr,'.mat'])
    
    x1 = pach.xpbeg(patchTileNum);
    x2 = pach.xpfin(patchTileNum);
    y1 = pach.ypbeg(patchTileNum);
    y2 = pach.ypfin(patchTileNum);
    imKur(y1:y2,x1:x2) = reshape(metaCluster.phaseAtClk(:,end),ximg,yimg);
    
    
end





% Plot Image and Phases of Oscillators side by side.
H=figure; 
subplot(121), imagesc(imFull), title('Original Image','FontSize',20,'FontWeight','Bold')
colormap('bone'), freezeColors
colorbar('SouthOutside'), %cbfreeze
set(gca,'xtick',[],'ytick',[])
subplot(122), imagesc(imKur), title('Oscillator Phases','FontSize',20,'FontWeight','Bold')
colormap('hsv'), freezeColors
colorbar('SouthOutside'), %cbfreeze
set(gca,'xtick',[],'ytick',[])

% Save the plot.
saveGoodImg(H,'stitch_patches_info_image',sizeGoodIm)
close(H)