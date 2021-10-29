% Script to look at Kuramoto Phases when using Isotropic Diffusion - uninformed by image. 

[dirPre,sizeGoodIm] = onCluster;

fileType = 'BSDS_patch';
fileSize = '51x51_ds1';

netMethod = 'Mod_SKHAdj'; % ,'AAnrm','Mod_N&G','Mod_SKHAdj'}; {'IsoDiff',

matKurFilesDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/Kur_PIF_Fourier1/',netMethod,'/'];
matEigFilesDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/spectral/',netMethod,'/'];
    


files = dir([matKurFilesDir,'*.mat']);




    
% im_patch = {'10081_ptch1'}; %{'100007_ptch1','100007_ptch2','100007_ptch3','100039_ptch1','100039_ptch2', ...
%     %'100075_ptch1','100080_ptch1','100098_ptch1','100099_ptch1','101027_ptch1','101084_ptch1',...
%     %'101085_ptch1','101087_ptch1','102061_ptch1','102062_ptch1'};
% 
% rM = {'1'}; %}; % ,'Inf', ,'3','5','10'


%%

Q_isotrope = zeros(2601,2601,numel(files)); 


phase_time_vector = [1:19];


for k = 1:numel(files)
    
    
    nf = strfind(files(k).name,'_NF');
    params_deets = files(k).name(7:nf-1);
    

    h=figure;

    for i = 1%:numel(rM)

%         try   
            kur = load([matKurFilesDir,files(k).name]);
            eig = load([matEigFilesDir,'Evecs_',params_deets,'.mat']);
%         catch
%             disp(['meh, file not there? : ',im_patch{k}])
%             continue
%         end
       
       % for rM=X, do Q matrices match?  What is Eigenvalue Spectrum? Evectors do not match.  There is a degeneracy.
       if(i==2)
           Q_isotrope(:,:,k) = kur.netParams.Q;
       end
           
       
       
       
       % Compute mean pairwise Area under ROC Curve to put on figure.
       for j = 1:numel(kur.netParams.gT)
            kur1D(j) = mean( kur.AUC_ROC_1D.kur{j}(kur.AUC_ROC_1D.kur{j}>0) );
            kurUB(j) = mean( kur.AUC_ROC_GS.kur{j}(kur.AUC_ROC_GS.kur{j}>0) );
       end

       
       % Plot Eigenvectors 1-3 for each rMax value
       
       % Use visKurPhase_inHSV ...
       
       for j = 1:numel(phase_time_vector)
           subplot(1+1,numel(phase_time_vector),numel(phase_time_vector)*(i-1)+j), 
           imagesc( visKurPhase_inHSV( kur.netParams.im, reshape(kur.metaCluster.phaseAtClk(:,j),kur.netParams.Ndims) ) ), 
           set(gca,'XTick',[],'YTick',[]), axis square, colormap('hsv'), freezeColors
           %ylabel(['rM=',rM{i}]) 
           if(i==1); title(['(t=',num2str(phase_time_vector(j)),')']); end
       end
       
%        %
%        subplot(1+1,numel(phase_time_vector),5*(i-1)+2), imagesc(reshape(kur.metaCluster.phaseAtClk(:,2),51,51)), set(gca,'XTick',[],'YTick',[]), axis square, colormap('hsv'), freezeColors
%        if(i==1); title('(t=2)'); end
%        %
%        subplot(1+1,numel(phase_time_vector),5*(i-1)+3), imagesc(reshape(kur.metaCluster.phaseAtClk(:,3),51,51)), set(gca,'XTick',[],'YTick',[]), axis square, colormap('hsv'), freezeColors
%        if(i==1); title('(t=3)'); end
%        %
%        subplot(1+1,numel(phase_time_vector),5*(i-1)+4), imagesc(reshape(kur.metaCluster.phaseAtClk(:,9),51,51)), set(gca,'XTick',[],'YTick',[]), axis square, colormap('hsv'), freezeColors
%        if(i==1); title('(t=9)'); end
%        %
%        subplot(1+1,numel(phase_time_vector),5*(i-1)+5), imagesc(reshape(kur.metaCluster.phaseAtClk(:,end),51,51)), set(gca,'XTick',[],'YTick',[]), axis square, colormap('hsv'), freezeColors
%        xlabel(['(',num2str(mean(kur1D),2),',',num2str(mean(kurUB),2),')'])
%        if(i==1); title('(t=19)'); end

    end
    
    
    for j = 1:numel(eig.netParams.gT)
        im(j) = mean( eig.AUC_ROC_1D.im{j}(eig.AUC_ROC_1D.im{j}>0) );
    end
    
    

    subplot(1+1,5,5*(i)+1), imagesc(eig.netParams.im), set(gca,'XTick',[],'YTick',[]), axis square, colormap('bone'),freezeColors
    ylabel(['Img'])
    xlabel([num2str(mean(im),2)])
    %
    subplot(1+1,5,5*(i)+2), imagesc(eig.netParams.gT{1}), set(gca,'XTick',[],'YTick',[]), axis square, colormap('jet'), freezeColors
    ylabel(['gT1'])
    %
    subplot(1+1,5,5*(i)+3), imagesc(eig.netParams.gT{2}), set(gca,'XTick',[],'YTick',[]), axis square, colormap('jet'), freezeColors
    ylabel(['gT2'])
    %
    subplot(1+1,5,5*(i)+4), imagesc(eig.netParams.gT{3}), set(gca,'XTick',[],'YTick',[]), axis square, colormap('jet'), freezeColors
    ylabel(['gT3'])
    %
    subplot(1+1,5,5*(i)+5), imagesc(eig.netParams.gT{4}), set(gca,'XTick',[],'YTick',[]), axis square, colormap('jet'), freezeColors
    ylabel(['gT4'])
    
    keyboard

end


keyboard
    