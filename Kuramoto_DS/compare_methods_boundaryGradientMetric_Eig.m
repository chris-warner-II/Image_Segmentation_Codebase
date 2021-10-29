function [] = compare_methods_boundaryGradientMetric_Eig(method, ptch, rM)

% weird.



% save flags.
sav_plt1 = 1;


[dirPre,sizeGoodIm] = onCluster;

dirImgSave = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/imgs/spectral/',method];
if ~exist(dirImgSave,'dir')
    mkdir(dirImgSave);
end


% dirMatSave = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/'];
% if ~exist(dirMatSave,'dir')
%     mkdir(dirMatSave);
% end


%ptch = '118015_ptch2'; % Mod_SKH performs much better than IsoDiff
%ptch = '117054_ptch2'; % Iso Diff performs well here
%ptch = '100080_ptch2'; % Iso Diff performs well here too



%% Load KurMC & Evecs mat files with Simulation Results and Phase Distributions.
if strmatch(method,'IsoDiff')
    load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/spectral/',method,'/Evecs_',ptch,'_rM',rM,'.mat']);
else
    load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/spectral/',method,'/Evecs_',ptch,'_rM',rM,'_sDInf_sP0p2.mat']);
end

GTFile = load([dirPre,'images/BSDS_patch/101x101_ds1/',ptch,'.mat']); % Ground truth file containing Boundary information







% Image Patch and Various Eigenvector Combinations for visualization and plotting.
im = netParams.im;
ev1 = reshape(EVecsML(:,1), netParams.Ndims);
ev2o = reshape(EVecsML(:,2), netParams.Ndims);
ev3o = reshape(EVecsML(:,3), netParams.Ndims);

ev2 = reshape(sum(EVecsML(:,1:2),2), netParams.Ndims);
ev3 = reshape(sum(EVecsML(:,1:3),2), netParams.Ndims);
ev2w = reshape(EVecsML(:,1).*EValsML(1) + EVecsML(:,2).*EValsML(2), netParams.Ndims);
ev3w = reshape(EVecsML(:,1).*EValsML(1) + EVecsML(:,2).*EValsML(2) + EVecsML(:,3).*EValsML(3), netParams.Ndims);






%% PLOT #1: Plot Combinations of Eigenvector & Gradients & d'-metric.
if(sav_plt1)

    
    Hc=figure;
    %ha = tight_subplot(2, 8, [0.05 0.01], 0.05, 0.01);
    
    % Image Pixels
    %axes(ha([1:2])), 
    subplot(2,3,1)
    imagesc(im), colormap(bone), freezeColors, axis square, 
    title(['\color{black}Img Pix'],'FontSize',20,'FontWeight','Bold') 
    set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
    %
    % Ground Truth Boundaries
    subplot(2,3,4)
    imagesc(max(GTFile.bD(:))-GTFile.bD), colormap(bone), freezeColors, axis square, 
    title(['GT Boundaries'],'FontSize',20,'FontWeight','Bold') 
    set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
    %
    % Text Information about Method & Parameters.
    subplot(4,6,3) 
    text(0.5,0.5,{[ptch],[method],['rM=',num2str(rM)]},'HorizontalAlignment','Center','VerticalAlignment','Middle','FontSize',16,'FontWeight','Bold')  
    axis off
    
    %
    % Image Pixel Gradients
    %axes(ha(9)), 
    subplot(4,6,6+3)   
    imagesc(MC.im.F), colormap(jet), freezeColors, axis square, 
    title(['Image'],'FontSize',20,'FontWeight','Bold')
    ylabel(['Gradients'],'FontSize',20,'FontWeight','Bold')
    text(netParams.Ndims(1),netParams.Ndims(2),['\color{white}d''= ',num2str(MC.im.D(1),2)],'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',16,'FontWeight','Bold')  
    set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
    xlabel(['\color{blue}',num2str(min(MC.im.F(:)),'%10.0e'),'    \color{red}',num2str(max(MC.im.F(:)),2)],'FontSize',16,'FontWeight','Bold')
    
    % Eigenvector 1
    %axes(ha(2)), 
    subplot(4,6,4)
    imagesc(ev1), colormap(jet), freezeColors, axis square, 
    ylabel(['\color{black}Evec 1'],'FontSize',20,'FontWeight','Bold') 
    set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
    xlabel(['\color{blue}',num2str(min(ev1(:)),'%10.0e'),'    \color{red}',num2str(max(ev1(:)),'%10.0e')],'FontSize',16,'FontWeight','Bold')
    
    %
    % Eigenvector 1 Gradients
    %axes(ha(10)), 
    subplot(4,6,6+4)
    imagesc(MC.ev1.F), colormap(jet), freezeColors, axis square, 
    text(netParams.Ndims(1),netParams.Ndims(2),['\color{white}d''= ',num2str(MC.ev1.D(1),2)],'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',16,'FontWeight','Bold') 
    ylabel(['\Delta\Phi'],'FontSize',20,'FontWeight','Bold')
    set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
    xlabel(['\color{blue}',num2str(min(MC.ev1.F(:)),'%10.0e'),'    \color{red}',num2str(max(MC.ev1.F(:)),2)],'FontSize',16,'FontWeight','Bold')
    
    % Eigenvector 2
    %axes(ha(3)), 
    subplot(4,6,5)
    imagesc(ev2o), colormap(jet), freezeColors, axis square, 
    ylabel(['\color{black}Evec 2'],'FontSize',20,'FontWeight','Bold') 
    set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
    xlabel(['\color{blue}',num2str(min(ev2o(:)),'%10.0e'),'    \color{red}',num2str(max(ev2o(:)),'%10.0e')],'FontSize',16,'FontWeight','Bold')
    %
    % Eigenvector 2 Gradients
    %axes(ha(11)), 
    subplot(4,6,6+5)
    imagesc(MC.ev2o.F), colormap(jet), freezeColors, axis square, 
    text(netParams.Ndims(1),netParams.Ndims(2),['\color{white}d''= ',num2str(MC.ev2o.D(1),2)],'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',16,'FontWeight','Bold') 
    ylabel(['\Delta\Phi'],'FontSize',20,'FontWeight','Bold')
    set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
    xlabel(['\color{blue}',num2str(min(MC.ev2o.F(:)),'%10.0e'),'    \color{red}',num2str(max(MC.ev2o.F(:)),2)],'FontSize',16,'FontWeight','Bold')
    
    % Eigenvector 3
    %axes(ha(4)), 
    subplot(4,6,6)
    imagesc(ev3o), colormap(jet), freezeColors, axis square, 
    ylabel(['\color{black}Evec 3'],'FontSize',20,'FontWeight','Bold') 
    set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
    xlabel(['\color{blue}',num2str(min(ev3o(:)),'%10.0e'),'    \color{red}',num2str(max(ev3o(:)),'%10.0e')],'FontSize',16,'FontWeight','Bold')
    %
    % Eigenvector 3 Gradients
    %axes(ha(12)), 
    subplot(4,6,6+6)
    imagesc(MC.ev3o.F), colormap(jet), freezeColors, axis square, 
    text(netParams.Ndims(1),netParams.Ndims(2),['\color{white}d''= ',num2str(MC.ev3o.D(1),2)],'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',16,'FontWeight','Bold') 
    ylabel(['\Delta\Phi'],'FontSize',20,'FontWeight','Bold')
    set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
    xlabel(['\color{blue}',num2str(min(MC.ev3o.F(:)),'%10.0e'),'    \color{red}',num2str(max(MC.ev3o.F(:)),2)],'FontSize',16,'FontWeight','Bold')
    
    % Eigenvector 1&2
    %axes(ha(5)), 
    subplot(4,6,2*6+3)
    imagesc(ev2), colormap(jet), freezeColors, axis square, 
    ylabel(['\color{black}Evec 1-2'],'FontSize',20,'FontWeight','Bold') 
    set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
    xlabel(['\color{blue}',num2str(min(ev2(:)),'%10.0e'),'    \color{red}',num2str(max(ev2(:)),'%10.0e')],'FontSize',16,'FontWeight','Bold')
    %
    % Eigenvector 1&2 Gradients
    %axes(ha(13)), 
    subplot(4,6,3*6+3)
    imagesc(MC.ev2.F), colormap(jet), freezeColors, axis square, 
    text(netParams.Ndims(1),netParams.Ndims(2),['\color{white}d''= ',num2str(MC.ev2.D(1),2)],'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',16,'FontWeight','Bold')
    ylabel(['\Delta\Phi'],'FontSize',20,'FontWeight','Bold')
    set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
    xlabel(['\color{blue}',num2str(min(MC.ev2.F(:)),'%10.0e'),'    \color{red}',num2str(max(MC.ev2.F(:)),2)],'FontSize',16,'FontWeight','Bold')
    
    % Eigenvector 1-3
    %axes(ha(6)), 
    subplot(4,6,2*6+4)
    imagesc(ev3), colormap(jet), freezeColors, axis square, 
    ylabel(['\color{black}Evec 1-3'],'FontSize',20,'FontWeight','Bold') 
    set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
    xlabel(['\color{blue}',num2str(min(ev3(:)),'%10.0e'),'    \color{red}',num2str(max(ev3(:)),'%10.0e')],'FontSize',16,'FontWeight','Bold')
    %
    % Eigenvector 1-3 Gradients
    %axes(ha(14)), 
    subplot(4,6,3*6+4)
    imagesc(MC.ev3.F), colormap(jet), freezeColors, axis square, 
    text(netParams.Ndims(1),netParams.Ndims(2),['\color{white}d''= ',num2str(MC.ev3.D(1),2)],'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',16,'FontWeight','Bold')
    ylabel(['\Delta\Phi'],'FontSize',20,'FontWeight','Bold')
    set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
    xlabel(['\color{blue}',num2str(min(MC.ev3.F(:)),'%10.0e'),'    \color{red}',num2str(max(MC.ev3.F(:)),2)],'FontSize',16,'FontWeight','Bold')
    
    % Eigenvector 1&2 Weighted
    %axes(ha(7)), 
    subplot(4,6,2*6+5)
    imagesc(ev2w), colormap(jet), freezeColors, axis square, 
    ylabel(['\color{black}Evec 1-2 W'],'FontSize',20,'FontWeight','Bold') 
    set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
    xlabel(['\color{blue}',num2str(min(ev2w(:)),'%10.0e'),'    \color{red}',num2str(max(ev2w(:)),'%10.0e')],'FontSize',16,'FontWeight','Bold')
    %
    % Eigenvector 1&2 Weighted Gradients
    %axes(ha(15)), 
    subplot(4,6,3*6+5)
    imagesc(MC.ev2w.F), colormap(jet), freezeColors, axis square, 
    text(netParams.Ndims(1),netParams.Ndims(2),['\color{white}d''= ',num2str(MC.ev2w.D(1),2)],'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',16,'FontWeight','Bold')
    ylabel(['\Delta\Phi'],'FontSize',20,'FontWeight','Bold')
    set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
    xlabel(['\color{blue}',num2str(min(MC.ev2w.F(:)),'%10.0e'),'    \color{red}',num2str(max(MC.ev2w.F(:)),2)],'FontSize',16,'FontWeight','Bold')
    
    % Eigenvector 1-3 Weighted
    %axes(ha(8)), 
    subplot(4,6,2*6+6)
    imagesc(ev3w), colormap(jet), freezeColors, axis square, 
    ylabel(['\color{black}Evec 1-3 W'],'FontSize',20,'FontWeight','Bold') 
    set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
    xlabel(['\color{blue}',num2str(min(ev3w(:)),'%10.0e'),'    \color{red}',num2str(max(ev3w(:)),'%10.0e')],'FontSize',16,'FontWeight','Bold')
    %
    % Eigenvector 1-3 Weighted Gradients
    %axes(ha(16)), 
    subplot(4,6,3*6+6)
    imagesc(MC.ev3w.F), colormap(jet), freezeColors, axis square, 
    text(netParams.Ndims(1),netParams.Ndims(2),['\color{white}d''= ',num2str(MC.ev3w.D(1),2)],'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',16,'FontWeight','Bold')
    ylabel(['\Delta\Phi'],'FontSize',20,'FontWeight','Bold')
    set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3),
    xlabel(['\color{blue}',num2str(min(MC.ev3w.F(:)),'%10.0e'),'    \color{red}',num2str(max(MC.ev3w.F(:)),2)],'FontSize',16,'FontWeight','Bold')
    


    saveGoodImg(Hc,[dirImgSave,'/Evecs_',ptch,'_rM',rM],sizeGoodIm)
    close(Hc)

end



