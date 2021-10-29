function plot_phaseAtClk(netMethod)

% Script to look at Kuramoto Phases when using Isotropic Diffusion - uninformed by image. 

[dirPre,sizeGoodIm] = onCluster;

fileType = 'BSDS_patch';
fileSize = '101x101_ds1';

%netMethod = 'IsoDiff'; % ,'AAnrm','Mod_N&G','Mod_SKHAdj'}; {,'Mod_SKHAdj'

matKurFilesDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/Kur_PIF_Fourier1/',netMethod,'/'];
matEigFilesDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/spectral/',netMethod,'/'];
    

phase_time_vector = [1,2,3,4,5,6,7,8,9];



files = dir([matKurFilesDir,'*100075_ptch2_rM10_*kslrg.mat']); % different Kscale values

LogPhaseChange = zeros(numel(files),18);



h=figure;





for k = 1:numel(files)
    
    
    nf = strfind(files(k).name,'_NF');
    params_deets = files(k).name(7:nf-1);
    

    kur = load([matKurFilesDir,files(k).name]);
%    eig = load([matEigFilesDir,'Evecs_',params_deets,'.mat']);
  
    
    % How much is the phase of oscillators changing with simulation time?
    x = diff(kur.metaCluster.phaseAtClk,1,2);
    LogPhaseChange(k,:) = log(mean(abs(x)));
    
    
    
       
%    % Compute mean pairwise Area under ROC Curve to put on figure.
%    for j = 1:numel(kur.netParams.gT)
%         kur1D(k,j) = mean( kur.AUC_ROC_1D.kur{j}(kur.AUC_ROC_1D.kur{j}>0) );
%         kurUB(k,j) = mean( kur.AUC_ROC_GS.kur{j}(kur.AUC_ROC_GS.kur{j}>0) );
%    end



   for j = 1:numel(phase_time_vector)
       figure(h), subplot(numel(files)+2,numel(phase_time_vector),numel(phase_time_vector)*(k-1)+j), 
       imagesc( visKurPhase_inHSV( kur.netParams.im, reshape(kur.metaCluster.phaseAtClk(:,j),kur.netParams.Ndims) ) ),
       caxis([0 2*pi])
       set(gca,'XTick',[],'YTick',[]), axis square, colormap('hsv'), freezeColors
       if(k==1); title(['(t=',num2str(phase_time_vector(j)),')']); end
       if(j==1); ylabel(['(Ks=',num2str(kur.kurParams.Kscale,2),')']); end
   end
       
    
    
%     for j = 1:numel(eig.netParams.gT)
%         im(j) = mean( eig.AUC_ROC_1D.im{j}(eig.AUC_ROC_1D.im{j}>0) );
%     end
    
end    




% Plot Phase Change vs. Kuramoto Simulation Time-Step for Different Kscales.
subplot(numel(files)+2,1,numel(files)+1),
plot(LogPhaseChange')
xlabel('Simulation Time Step')
ylabel('Log \Delta in Phase')
title(['(Ks=',num2str(kur.kurParams.Kscale,2),')'])
legend({'1','2','3'})



subplot(numel(files)+2,6,6*(numel(files)+1)+1), 
imagesc(eig.netParams.im), set(gca,'XTick',[],'YTick',[]), 
axis square, colormap('bone'),freezeColors
ylabel(['Img'])
xlabel({ 'im pix', 'ks=1', 'ks=2', 'ks=3'})
%
for j = 1:numel(eig.netParams.gT)
    subplot(numel(files)+2,6,6*(numel(files)+1)+j+1), 
    imagesc(eig.netParams.gT{j}), set(gca,'XTick',[],'YTick',[]), 
    axis square, colormap('jet'), freezeColors
    ylabel(['gT',num2str(j)])
    xlabel({num2str(im(j),2),num2str(kur1D(:,j),2)})
end


keyboard



