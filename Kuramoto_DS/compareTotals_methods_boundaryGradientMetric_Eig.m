

% save flags.
sav_plt_butterfly = 1; % Make butterfly scatter plots of d' Pix vs. d' Method for all 5 methods and different blurs.
sav_plt6 = 0;          % 

sav_mat = 1;


[dirPre,sizeGoodIm] = onCluster;

dirImgSave = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/imgs/methodsComparePlots/'];
if ~exist(dirImgSave,'dir')
    mkdir(dirImgSave);
end





%% Load KurMC & Evecs mat files with Simulation Results and Phase Distributions.
IsoE = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/spectral/IsoDiff/BoundaryGradientMetricEig_rM3_files500.mat']);
IsoE.rM = 3;
%
SK_E = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/spectral/Mod_SKHAdj/BoundaryGradientMetricEig_rM3_files500.mat']);
SK_E.rM = 3;
%
% NG_E = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/spectral/Mod_N&G/BoundaryGradientMetricEig_rM3_files500.mat']);
% NG_E.rM = 3;
%
AA_E = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/spectral/AAnrm/BoundaryGradientMetricEig_rM3_files500.mat']);
AA_E.rM = 3;
%
GL_E = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/spectral/GLnrm/BoundaryGradientMetricEig_rM3_files500.mat']);
GL_E.rM = 3;






%% Compute Discriminability (d') of boundaries for all 6 different methods

% %PixK = computeDiscriminability(PixK); % SHOULD BE THE SAME AS D_PRIME ABOVE.  CHECK AND BE SURE.  THEN USE THIS ONE.
% IsoK = computeDiscriminability(IsoK);
% AA_K = computeDiscriminability(AA_K);
% GL_K = computeDiscriminability(GL_K);
% %NG_K = computeDiscriminability(NG_K);
% SK_K = computeDiscriminability(SK_K);

% I THINK IT IS NOW BEING DONE RIGHT SO ONLY LOOK INTO IF IT SEEMS WRONG.






%% PLOT # 5: Make a butterfly plot or a scatter plot that compares directly d' of method vs pixels.
if(sav_plt_butterfly)
    
    which_single_patch = '';
    indicate_scatter_single_patch=0;
    show_single_gradient_inset=0;
    
    
    % find which data point pertains to the single patch of interest.
    for i = 1:numel(AA_E.ImgPtchID)
        x=strmatch(which_single_patch,AA_E.ImgPtchID{i});
        if(x)
           patch_ind = i; 
        end
    end
    
    
    
    
    
    
    
    if(show_single_gradient_inset)
        
        
    end
    
    
    
    
    

    for i=[1]  % 1:numel(PixK.bDc_blurs_info)
        

        Hc=figure;
        ha = tight_subplot(2, 3, [0.05 0], 0.05, 0.01);
        maxx = 1.05.*max([AA_E.boundaryDiscriminability.im(i,:),AA_E.boundaryDiscriminability.ev1(i,:),AA_E.boundaryDiscriminability.ev2o(i,:),AA_E.boundaryDiscriminability.ev3o(i,:),...
                          IsoE.boundaryDiscriminability.im(i,:),IsoE.boundaryDiscriminability.ev1(i,:),IsoE.boundaryDiscriminability.ev2o(i,:),IsoE.boundaryDiscriminability.ev3o(i,:),...
                          GL_E.boundaryDiscriminability.im(i,:),GL_E.boundaryDiscriminability.ev1(i,:),GL_E.boundaryDiscriminability.ev2o(i,:),GL_E.boundaryDiscriminability.ev3o(i,:),...
                          SK_E.boundaryDiscriminability.im(i,:),SK_E.boundaryDiscriminability.ev1(i,:),SK_E.boundaryDiscriminability.ev2o(i,:),SK_E.boundaryDiscriminability.ev3o(i,:),...
                          ]);
        set(ha,'FontSize',16,'FontWeight','Bold','XTick',[0:0.5:round(maxx)],'YTick',[0:0.5:round(maxx)],'XtickLabel',[0:0.5:round(maxx)],'YtickLabel',[0:0.5:round(maxx)])

        
        
        
        % #1. Average Association (AA)
        axes(ha(1)), hold on
        scatter(AA_E.boundaryDiscriminability.im(i,:),AA_E.boundaryDiscriminability.ev1(i,:),50,'r.')
        scatter(AA_E.boundaryDiscriminability.im(i,:),AA_E.boundaryDiscriminability.ev2o(i,:),50,'g.')
        scatter(AA_E.boundaryDiscriminability.im(i,:),AA_E.boundaryDiscriminability.ev3o(i,:),50,'b.')
        %scatter(AA_E.boundaryDiscriminability.im(i,:),AA_E.boundaryDiscriminability.ev2(i,:),50,'k.')
        %scatter(AA_E.boundaryDiscriminability.im(i,:),AA_E.boundaryDiscriminability.ev2w(i,:),50,'c.')
        plot([0 maxx], [0 maxx],'k--','LineWidth',2)
        axis([0 maxx 0 maxx])
        axis square
        %xlabel(['d'' Pix'],'FontSize',18,'FontWeight','Bold')
        ylabel(['d'' \color{magenta}AA'],'FontSize',18,'FontWeight','Bold')
        text( maxx, 0, {['rM=',num2str(AA_E.rM)]},'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',16,'FontWeight','Bold','EdgeColor','k')
        grid on
        
        
        
        
%         % #2. Non-Topographic Modularity (NG)
        axes(ha(2)), hold on
%         scatter(NG_E.boundaryDiscriminability.im(i,:),NG_E.boundaryDiscriminability.ev1(i,:),100,'r.')
%         scatter(NG_E.boundaryDiscriminability.im(i,:),NG_E.boundaryDiscriminability.ev2o(i,:),100,'g.')
%         scatter(NG_E.boundaryDiscriminability.im(i,:),NG_E.boundaryDiscriminability.ev3o(i,:),100,'b.')
        plot([0 maxx], [0 maxx],'k--','LineWidth',2)
        axis([0 maxx 0 maxx])
        axis square
%         %xlabel(['d'' Pix'],'FontSize',18,'FontWeight','Bold')
        ylabel(['d'' \color{cyan}NG'],'FontSize',18,'FontWeight','Bold')
%         text( maxx, 0, {['rM=',num2str(NG_E.rM)]},'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',16,'FontWeight','Bold','EdgeColor','k')
        grid on
        


        % # 3. Isotropic Diffusion (Iso)
        axes(ha(3)), hold on
        scatter(IsoE.boundaryDiscriminability.im(i,:),IsoE.boundaryDiscriminability.ev1(i,:),100,'r.')
        scatter(IsoE.boundaryDiscriminability.im(i,:),IsoE.boundaryDiscriminability.ev2o(i,:),100,'g.')
        scatter(IsoE.boundaryDiscriminability.im(i,:),IsoE.boundaryDiscriminability.ev3o(i,:),100,'b.')
        plot([0 maxx], [0 maxx],'k--','LineWidth',2)
        axis([0 maxx 0 maxx])
        axis square
        %xlabel(['d'' Pix'],'FontSize',18,'FontWeight','Bold')
        ylabel(['d'' \color{blue}Iso'],'FontSize',18,'FontWeight','Bold')
        text( maxx, 0, {['rM=',num2str(IsoE.rM)]},'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',16,'FontWeight','Bold','EdgeColor','k')
        grid on


        % #4. Graph Laplacian (GL)
        axes(ha(4)), hold on
        scatter(GL_E.boundaryDiscriminability.im(i,:),GL_E.boundaryDiscriminability.ev1(i,:),100,'r.')
        scatter(GL_E.boundaryDiscriminability.im(i,:),GL_E.boundaryDiscriminability.ev2o(i,:),100,'g.')
        scatter(GL_E.boundaryDiscriminability.im(i,:),GL_E.boundaryDiscriminability.ev3o(i,:),100,'b.')
        plot([0 maxx], [0 maxx],'k--','LineWidth',2)
        axis([0 maxx 0 maxx])
        axis square
        %xlabel(['d'' Pix'],'FontSize',18,'FontWeight','Bold')
        ylabel(['d'' \color{green}GL'],'FontSize',18,'FontWeight','Bold')
        text( maxx, 0, {['rM=',num2str(GL_E.rM)]},'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',16,'FontWeight','Bold','EdgeColor','k')
        grid on
        
        
        % #5. Topographic Modularity (SK)
        axes(ha(5)), hold on
        scatter(SK_E.boundaryDiscriminability.im(i,:),SK_E.boundaryDiscriminability.ev1(i,:),100,'r.')
        scatter(SK_E.boundaryDiscriminability.im(i,:),SK_E.boundaryDiscriminability.ev2o(i,:),100,'g.')
        scatter(SK_E.boundaryDiscriminability.im(i,:),SK_E.boundaryDiscriminability.ev3o(i,:),100,'b.')
        plot([0 maxx], [0 maxx],'k--','LineWidth',2)
        axis([0 maxx 0 maxx])
        axis square
        %xlabel(['d'' Pix'],'FontSize',18,'FontWeight','Bold')
        ylabel(['d'' \color{red}SK'],'FontSize',18,'FontWeight','Bold')
        text( maxx, 0, {['rM=',num2str(SK_E.rM)]},'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',16,'FontWeight','Bold','EdgeColor','k')
        grid on


        
        
        % #6. Clear the 6th axis
        axes(ha(6)),
        axis off
        
        
    end
    

    disp('Note: Important Remaining Problem.  What to do with negative d'' values?')
    

end



%% PLOT # 6: Plot Full Pair of 2D Histograms (mean & std of gradients on & off boundary)
if(sav_plt6)
    hist2d([SK_K.meanGradientOnBoundary(1,:)./pi;SK_K.stdGradientOnBoundary(1,:)./pi],20,20,[0 0.5],[0 0.5])
    xlabel('mean \Delta')
    ylabel('std \Delta')
    title('SK On'),view(0,90)
    %
    hist2d([SK_K.meanGradientOffBoundary(1,:)./pi;SK_K.stdGradientOffBoundary(1,:)./pi],20,20,[0 0.5],[0 0.5])
    xlabel('mean \Delta')
    ylabel('std \Delta')
    title('SK Off'),view(0,90)

    hist2d([PixK.meanGradientOnBoundary(1,:);PixK.stdGradientOnBoundary(1,:)],20,20,[0 0.5],[0 0.5])
    xlabel('mean \Delta')
    ylabel('std \Delta')
    title('Pix On'),view(0,90)
    %
    hist2d([PixK.meanGradientOffBoundary(1,:);PixK.stdGradientOffBoundary(1,:)],20,20,[0 0.5],[0 0.5])
    xlabel('mean \Delta')
    ylabel('std \Delta')
end
title('Pix Off'),view(0,90)

keyboard