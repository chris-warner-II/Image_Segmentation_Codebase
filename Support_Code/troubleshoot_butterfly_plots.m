% This scripts is to look at KurMC mat files in both GL & AA methods and
% compare them to the mat files in image pixels. I want to plot d' for
% these things to compare them.  I expect them to be the same and that the 
% butterfly plot should hug the diagonal line.  But in practice I find some
% spread to the butterfly plots and I want to investigate what causes that
% or if it is real...



[dirPre,sizeGoodIm] = onCluster;

% dirImgSave = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/imgs/methodsComparePlots/'];
% if ~exist(dirImgSave,'dir')
%     mkdir(dirImgSave);
% end


method = 'AAnrm';




%% Load KurMC & Evecs mat files with Simulation Results and Phase Distributions.

files = dir([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/',method,'/KurMC*rM10_*kssml*.mat']);


hc = figure;
subplot(131), axis square, hold on
subplot(132), axis square, hold on
subplot(133), axis square, hold on


for F = 1:numel(files)
    
    
    disp(['File # ',num2str(F),' / ',num2str(numel(files))])

    mthSml = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/',method,'/',files(F).name]);
    mthMid = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/',method,'/',files(F).name(1:end-7),'mid.mat']);
    mthLrg = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/',method,'/',files(F).name(1:end-7),'lrg.mat']);
    %
    st = 7;
    fn = strfind(files(F).name,'_rM')-1;
    imgPix = load([dirPre,'images/BSDS_patch/101x101_ds1/',files(F).name(st:fn),'.mat']);

    

    % [imgPix.MC.D,mthSml.MC.D,mthMid.MC.D,mthLrg.MC.D,]
    
    
    
    figure(hc)
    subplot(131), scatter(imgPix.MC.D(1), mthSml.MC.D(1), 'bx')
    subplot(132), scatter(imgPix.MC.D(1), mthMid.MC.D(1), 'gx')
    subplot(133), scatter(imgPix.MC.D(1), mthLrg.MC.D(1), 'rx')
    
    
end