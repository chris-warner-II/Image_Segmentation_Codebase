% Script to look at Eigenvectors of Matrices that are doing Isotropic
% Diffusion.  For a given image, the first 3 eigenvectors are plotted for
% different rMax values.  For a given rMax value, an eigenvector should
% look identical regardless of this image.  

[dirPre,sizeGoodIm] = onCluster;

fileType = 'BSDS_patch';
fileSize = '51x51_ds1';

matFilesDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/spectral_C/IsoDiff/'];
    
    
im_patch = {'100007_ptch1','100007_ptch2','100007_ptch3','100039_ptch1','100039_ptch2', ...
    '100075_ptch1','100080_ptch1','100098_ptch1','100099_ptch1','101027_ptch1','101084_ptch1',...
    '101085_ptch1','101087_ptch1','102061_ptch1','102062_ptch1'};

rM = {'1','3','5','10'}; % ,'Inf'


%%

Q_isotrope = zeros(2601,2601,numel(im_patch)); 

for k = 1:numel(im_patch)

    h=figure;

    for i = 1:numel(rM)

        try
            load([matFilesDir,'Evecs_',im_patch{k},'_rM',rM{i},'_sDInf_sPInf.mat']);
        catch
            disp('meh, file not there?')
            continue
        end
       
       % for rM=X, do Q matrices match?  What is Eigenvalue Spectrum? Evectors do not match.  There is a degeneracy.
       if(i==2)
           Q_isotrope(:,:,k) = netParams.Q;
           EValsML(1:3)
           diff(EValsML(1:3))
       end
           
       
%        % Compute Eigenvalues and Eigenvectors again to see if they look different.  Does the degeneracy result in different answers?
%        [xxx,yyy] = eigs(netParams.Q);
%        [xxx2,yyy2] = eigs(netParams.Q);
%        figure, imagesc(abs(xxx-xxx2))
       
       
       
       
       % Compute mean pairwise Area under ROC Curve to put on figure.
       for j = 1:numel(netParams.gT)
            ev1(j) = mean( AUC_ROC_1D.ev1{j}(AUC_ROC_1D.ev1{j}>0) );
            ev2(j) = mean( AUC_ROC_1D.ev2o{j}(AUC_ROC_1D.ev2o{j}>0) );
            ev3(j) = mean( AUC_ROC_1D.ev3o{j}(AUC_ROC_1D.ev3o{j}>0) );
       end

       
       % Plot Eigenvectors 1-3 for each rMax value
       subplot(numel(rM)+1,3,3*(i-1)+1), imagesc(reshape(EVecsML(:,1),51,51)), set(gca,'XTick',[],'YTick',[]), axis square, colormap('jet'), freezeColors
       ylabel(['rM=',rM{i}]) 
       xlabel([num2str(mean(ev1),2)])
       if(i==1); title('Ev1'); end
       %subplot(numel(rM)+1,6,6*(i-1)+2), imagesc(reshape(xxx(:,1),51,51)), set(gca,'XTick',[],'YTick',[]), axis square
       %
       subplot(numel(rM)+1,3,3*(i-1)+2), imagesc(reshape(EVecsML(:,2),51,51)), set(gca,'XTick',[],'YTick',[]), axis square, colormap('jet'), freezeColors
       xlabel([num2str(mean(ev2),2)])
       if(i==1); title('Ev2'); end
       %subplot(numel(rM)+1,6,6*(i-1)+4), imagesc(reshape(xxx(:,2),51,51)), set(gca,'XTick',[],'YTick',[]), axis square
       %
       subplot(numel(rM)+1,3,3*(i-1)+3), imagesc(reshape(EVecsML(:,3),51,51)), set(gca,'XTick',[],'YTick',[]), axis square, colormap('jet'), freezeColors
       xlabel([num2str(mean(ev3),2)])
       if(i==1); title('Ev3'); end
       %subplot(numel(rM)+1,6,6*(i-1)+6), imagesc(reshape(xxx(:,3),51,51)), set(gca,'XTick',[],'YTick',[]), axis square

    end
    
    
    for j = 1:numel(netParams.gT)
        im(j) = mean( AUC_ROC_1D.im{j}(AUC_ROC_1D.im{j}>0) );
    end
    
    

    subplot(numel(rM)+1,3,3*(i)+1), imagesc(netParams.im), set(gca,'XTick',[],'YTick',[]), axis square, colormap('bone'),freezeColors
    ylabel(['Img'])
    xlabel([num2str(mean(im),2)])

    subplot(numel(rM)+1,3,3*(i)+2), imagesc(netParams.gT{1}), set(gca,'XTick',[],'YTick',[]), axis square, colormap('jet'), freezeColors
    ylabel(['gT1'])

    subplot(numel(rM)+1,3,3*(i)+3), imagesc(netParams.gT{2}), set(gca,'XTick',[],'YTick',[]), axis square, colormap('jet'), freezeColors
    ylabel(['gT2'])
    
    %keyboard

end


keyboard
    