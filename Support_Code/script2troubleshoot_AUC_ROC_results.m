function script2troubleshoot_AUC_ROC_results(methodType,rM)

% This script will input a mat file from AUC_ROC_results directory and
% analyze how often and under what conditions the 2 problems with Kur GS
% tend to happen.  Maybe at some point I will combine this with the
% compare_methods_AUC_ROC_res.m file.



fileType = 'BSDS_patch';
fileSize = '51x51_ds1';
%methodType = 'GLnrm';
%
%rM = 3;
sW = 0; 
%
sP = 0.2; sPx = num2str(sP); sPx(sPx=='.')='p';
sD = inf;


[dirPre,sizeGoodIm] = onCluster;

% directory to AUC_ROC_results mat file
dirAUCMats = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/AUC_ROC_results/',methodType,'/'];

if strcmp(methodType,'IsoDiff')
    load([dirAUCMats,'AUCdata_',methodType,'_rM',num2str(rM),'_NF_60_',num2str(sW),'.mat']);
else
    load([dirAUCMats,'AUCdata_',methodType,'_rM',num2str(rM),'_sD',num2str(sD),'_sP',sPx,'_NF_60_',num2str(sW),'.mat'])
end


% directory to KurMC mat files
dirKurMats = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/Kur_PIF_Fourier1_C/',methodType,'/'];


% % directory to the Evecs mat files
% dirEigMats = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/spectral_C/',methodType,'/'];





%% Problem #1: when the AUC for kur1D search is better than AUC for kur Gridsearch.
ind1 = find(kur_ROC_AUC_1D_accum - kur_ROC_AUC_GS_accum > 1e-6);

if isempty(ind1)
    disp('Problem #1 (where AUC 1D > AUC GS) solved.')
else
    disp('Problem #1 : AUC 1D > AUC GS.')
    methodType
    disp('ImgPatchID')
    ImgPatchID(ind1)
    disp('GndTruthID')
    GndTruthID(ind1)
    disp('ClusterPairID')
    ClusterPairID(:,ind1)
    disp('ClusterSizes')
    ClusterSizes(:,ind1)
    disp('GS & 1D Kur AUC values')
    [kur_ROC_AUC_1D_accum(ind1);kur_ROC_AUC_GS_accum(ind1)]
end


% Loop through ind1 and open KurMC & Evec files to see what is happening...
for i = 1:-1 %numel(ind1)

    i
    ImgPatchID{ind1(i)}
    GndTruthID(ind1(i))
    
    if strcmp(methodType,'IsoDiff')
        load([dirKurMats,'KurMC_',ImgPatchID{ind1(i)},'_rM',num2str(rM),'_NF_60_',num2str(sW),'.mat'])
    else
        load([dirKurMats,'KurMC_',ImgPatchID{ind1(i)},'_rM',num2str(rM),'_sD',num2str(sD),'_sP',sPx,'_NF_60_',num2str(sW),'.mat'])
    end

    
    AUC_ROC_1D.kur{GndTruthID(ind1(i))}
    AUC_ROC_GS.kur{GndTruthID(ind1(i))}

    
    (AUC_ROC_1D.kur{GndTruthID(ind1(i))} > AUC_ROC_GS.kur{GndTruthID(ind1(i))})
    disp('1D Search AUC Kuramoto')
    AUC_ROC_1D.kur{GndTruthID(ind1(i))}(AUC_ROC_1D.kur{GndTruthID(ind1(i))} > AUC_ROC_GS.kur{GndTruthID(ind1(i))})
    disp('Grid Search AUC Kuramoto')
    AUC_ROC_GS.kur{GndTruthID(ind1(i))}(AUC_ROC_1D.kur{GndTruthID(ind1(i))} > AUC_ROC_GS.kur{GndTruthID(ind1(i))})
    
    
    ClusterPairID(:,ind1(i))
    ClusterSizes(:,ind1(i))
    
    
    
    
    
    if(1)
        figure, 
        subplot(141),imagesc(netParams.gT{GndTruthID(ind1(i))}), colormap('jet'), colorbar('SouthOutside'), 
        freezeColors, cbfreeze,
        axis off square, title(['GT#',num2str(GndTruthID(ind1(i))),' - -',num2str(ClusterPairID(1,ind1(i))),' & ',num2str(ClusterPairID(2,ind1(i)))])
        %
        subplot(142),imagesc(netParams.im), colormap('bone'), colorbar('SouthOutside'), 
        freezeColors, cbfreeze
        axis off square, title(['image'])
        %
        phaseInit = visKurPhase_inHSV(netParams.im , reshape(metaCluster.phaseAtClk(:,1),51,51) );
        subplot(143),imagesc(phaseInit), caxis([0 2*pi]), colormap('hsv'), colorbar('SouthOutside'),  
        freezeColors, cbfreeze,
        axis off square, title([methodType,' Kuramoto Init Phase'])
        %
        phaseFin = visKurPhase_inHSV(netParams.im , reshape(metaCluster.phaseAtClk(:,end),51,51) );
        subplot(144),imagesc(phaseFin), caxis([0 2*pi]), colormap('hsv'), colorbar('SouthOutside'),  
        freezeColors, cbfreeze,
        axis off square, title([methodType,' Kuramoto Final Phase'])
    end
    
    

    keyboard






end


%% Problem #2: when Convex Hull computation for kur Gridsearch failed.
ind2 = find(isnan(kur_ROC_AUC_GS_accum));



if isempty(ind2)
    disp('Problem #2 (where AUC GS is NAN because Convex Hull failed) solved.')
else
    disp('Problem #2 : AUC GS is NAN because Convex Hull failed.')
    methodType
    disp('ImgPatchID')
    ImgPatchID(ind2)
    disp('GndTruthID')
    GndTruthID(ind2)
    disp('ClusterPairID')
    ClusterPairID(:,ind2)
    disp('ClusterSizes')
    ClusterSizes(:,ind2)
    disp('GS & 1D Kur AUC values')
    [kur_ROC_AUC_1D_accum(ind2);kur_ROC_AUC_GS_accum(ind2);...
        ev1_ROC_AUC_1D_accum(ind2); ev2o_ROC_AUC_1D_accum(ind2); ev3o_ROC_AUC_1D_accum(ind2); im_ROC_AUC_1D_accum(ind2)]
end





%% So many small cluster pairs. (Doesnt seem to matter in results frankly)
if(1)
    figure, hist2d(ClusterSizes)
    title('Sizes of Clusters in Pair')
end



ind10pix = 1:numel(ImgPatchID);                        % all data points regardless of cluster size
% find(ClusterSizes(1,:)<10 & ClusterSizes(2,:)<10);   % When both clusters are less than 10 pixels.

%[kur_ROC_AUC_1D_accum(ind10pix);kur_ROC_AUC_GS_accum(ind10pix);...
%        ev1_ROC_AUC_1D_accum(ind10pix); ev2o_ROC_AUC_1D_accum(ind10pix); ev3o_ROC_AUC_1D_accum(ind10pix); im_ROC_AUC_1D_accum(ind10pix)]
    
    
    
figure,
subplot(321), hist(kur_ROC_AUC_1D_accum(ind10pix)), title('Kur 1D'), axis([0.5 1 0 numel(ind10pix)]), grid on
subplot(323), hist(kur_ROC_AUC_GS_accum(ind10pix)), title('Kur GS'), axis([0.5 1 0 numel(ind10pix)]), grid on
subplot(322), hist(ev1_ROC_AUC_1D_accum(ind10pix)), title('Ev 1'), axis([0.5 1 0 numel(ind10pix)]), grid on
subplot(324), hist(ev2o_ROC_AUC_1D_accum(ind10pix)), title('Ev 2'), axis([0.5 1 0 numel(ind10pix)]), grid on
subplot(326), hist(ev3o_ROC_AUC_1D_accum(ind10pix)), title('Ev 3'), axis([0.5 1 0 numel(ind10pix)]), grid on
subplot(325), hist(im_ROC_AUC_1D_accum(ind10pix)), title('Im Pix'), axis([0.5 1 0 numel(ind10pix)]), grid on
xlabel('ROC AUC'), ylabel('# counts')
    
    

%keyboard



% 3D Scatter plot size of clusters vs AUC performance of different segMethods.
figure, hold on,
scatter3(ClusterSizes(1,:),ClusterSizes(2,:),kur_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum,30,'b.')
scatter3(ClusterSizes(1,:),ClusterSizes(2,:),im_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum,30,'g.')
scatter3(ClusterSizes(1,:),ClusterSizes(2,:),ev2o_ROC_AUC_1D_accum - im_ROC_AUC_1D_accum,30,'r.')

legend({'kur','img','ev2'})


%% Split data into 3 cases based on AUC_imgPix.

ind_AUCi_sml = find(im_ROC_AUC_1D_accum<0.7);
ind_AUCi_mid = find( (im_ROC_AUC_1D_accum>=0.7) & (im_ROC_AUC_1D_accum<=0.9) );
ind_AUCi_lrg = find(im_ROC_AUC_1D_accum>0.9);



keyboard

