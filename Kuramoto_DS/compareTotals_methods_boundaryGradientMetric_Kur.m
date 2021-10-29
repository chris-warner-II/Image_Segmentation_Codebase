

% save flags.
sav_plt1 = 0;          % Error bar plot (mean & std) of gradient on & off boundary for 6 different methods {pix,iso,AA,GL,NG,SK}
                       % Maybe better would be some violin plots or something. [Note: for 6 different boundary blurrs.]
sav_plt2 = 0;          % Error bar plot for boundary gradient (on & off) for 6 methods with single (ie. no) blur.
sav_plt3 = 0;          % Gradient on vs. off boundary.  Not a useful plot any more.
sav_plt4 = 0;          % Full Histograms and mean&std of Boundary Discriminability (d') 
                       % for each method and Discriminability Improvement over Image Pixels
sav_plt_butterfly = 1; % Make butterfly scatter plots of d' Pix vs. d' Method for all 5 methods and different blurs.
sav_plt6 = 0;          % 

sav_mat = 1;


which_bD = 1;


[dirPre,sizeGoodIm] = onCluster;

dirImgSave = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/imgs/Kur_PIF_Fourier1/compareMethods/'];
if ~exist(dirImgSave,'dir')
    mkdir(dirImgSave);
end


% dirMatSave = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/'];
% if ~exist(dirMatSave,'dir')
%     mkdir(dirMatSave);
% end




%% Load KurMC & Evecs mat files with Simulation Results and Phase Distributions.
IsoK = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/IsoDiff/BoundaryGradientMetricKur_rM1_KSmid_files500.mat']);
IsoK.rM = 1;
IsoK.KS = 'mid';
%
AA_K = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/AAnrm/BoundaryGradientMetricKur_rM1_KSlrg_files500.mat']);
AA_K.rM = 1;
AA_K.KS = 'lrg';
%
GL_K = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/GLnrm/BoundaryGradientMetricKur_rM1_KSlrg_files500.mat']);
GL_K.rM = 1;
GL_K.KS = 'lrg';
%
NG_K = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/Mod_N&G/BoundaryGradientMetricKur_rM3_KSmid_files500.mat']);
NG_K.rM = 3;
NG_K.KS = 'mid';
%
SK_K = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/Mod_SKHAdj/BoundaryGradientMetricKur_rM3_KSmid_files500.mat']);
SK_K.rM = 3;
SK_K.KS = 'mid';
%
PixK = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/ImPix/BoundaryGradientMetricKur_rM_KS_files500.mat']);







% have to find intersection between PixK and various methods for butterfly plots.
[x,indIso,indPixIso] = intersect(IsoK.ImgPtchID,PixK.ImgPtchID);
[x,indAA,indPixAA] = intersect(AA_K.ImgPtchID,PixK.ImgPtchID);
[x,indGL,indPixGL] = intersect(GL_K.ImgPtchID,PixK.ImgPtchID);
[x,indNG,indPixNG] = intersect(NG_K.ImgPtchID,PixK.ImgPtchID);
[x,indSK,indPixSK] = intersect(SK_K.ImgPtchID,PixK.ImgPtchID);




%% PLOT #1: Visualize Gradients (phase,contrast or otherwise)
% for different competitor methods and plot avg gradient/pixel
% on boundaries vs. off for different boundary blurs.
if(sav_plt1)

    
    on(1,1,:) = mean(PixK.meanGradientOnBoundary,2);
    on(1,2,:) = std(PixK.meanGradientOnBoundary,[],2);
    %on(1,3,:) = mean(PixK.stdGradientOnBoundary,2);
    %
    off(1,1,:) = mean(PixK.meanGradientOffBoundary,2);
    off(1,2,:) = std(PixK.meanGradientOffBoundary,[],2);
    %off(1,3,:) = mean(PixK.stdGradientOffBoundary,2);

    
    on(2,1,:) = mean(IsoK.meanGradientOnBoundary,2)./pi;
    on(2,2,:) = std(IsoK.meanGradientOnBoundary,[],2)./pi;
    %on(2,3,:) = mean(IsoK.stdGradientOnBoundary,2)./pi;
    %
    off(2,1,:) = mean(IsoK.meanGradientOffBoundary,2)./pi;
    off(2,2,:) = std(IsoK.meanGradientOffBoundary,[],2)./pi;
    %off(2,3,:) = mean(IsoK.stdGradientOffBoundary,2)./pi;
    
    
    on(3,1,:) = mean(SK_K.meanGradientOnBoundary,2)./pi;
    on(3,2,:) = std(SK_K.meanGradientOnBoundary,[],2)./pi;
    %on(3,3,:) = mean(SK_K.stdGradientOnBoundary,2)./pi;
    %
    off(3,1,:) = mean(SK_K.meanGradientOffBoundary,2)./pi;
    off(3,2,:) = std(SK_K.meanGradientOffBoundary,[],2)./pi;
    %off(3,3,:) = mean(SK_K.stdGradientOffBoundary,2)./pi;
    
    
    on(4,1,:) = mean(NG_K.meanGradientOnBoundary,2)./pi;
    on(4,2,:) = std(NG_K.meanGradientOnBoundary,[],2)./pi;
    %on(4,3,:) = mean(NG_K.stdGradientOnBoundary,2)./pi;
    %
    off(4,1,:) = mean(NG_K.meanGradientOffBoundary,2)./pi;
    off(4,2,:) = std(NG_K.meanGradientOffBoundary,[],2)./pi;
    %off(4,3,:) = mean(NG_K.stdGradientOffBoundary,2)./pi;
    
    
    on(5,1,:) = mean(AA_K.meanGradientOnBoundary,2)./pi;
    on(5,2,:) = std(AA_K.meanGradientOnBoundary,[],2)./pi;
    %on(5,3,:) = mean(AA_K.stdGradientOnBoundary,2)./pi;
    %
    off(5,1,:) = mean(AA_K.meanGradientOffBoundary,2)./pi;
    off(5,2,:) = std(AA_K.meanGradientOffBoundary,[],2)./pi;
    %off(5,3,:) = mean(AA_K.stdGradientOffBoundary,2)./pi;
    
    
    on(6,1,:) = mean(GL_K.meanGradientOnBoundary,2)./pi;
    on(6,2,:) = std(GL_K.meanGradientOnBoundary,[],2)./pi;
    %on(6,3,:) = mean(GL_K.stdGradientOnBoundary,2)./pi;
    %
    off(6,1,:) = mean(GL_K.meanGradientOffBoundary,2)./pi;
    off(6,2,:) = std(GL_K.meanGradientOffBoundary,[],2)./pi;
    %off(6,3,:) = mean(GL_K.stdGradientOffBoundary,2)./pi;
    
    
    


    % Find max value for y-axis of error bar plots (M_mn+M_std or S_mn+S_std)
    del_mn_std = [max(max(sum(on,2))), max(max(sum(off,2)))];


    % Figure with Subplots !
    Hc=figure;
    %
    % Model 1: Image Pixels only (straw man model)
    subplot(231), hold on, errorbar(on(1,1,:),on(1,2,:),'b.-'), errorbar(off(1,1,:),off(1,2,:),'r.-')
    set(gca,'FontSize',16,'FontWeight','Bold'), axis([0.5 5+0.5 0 1.05*max(del_mn_std)])
    title('Pix','FontSize',20,'FontWeight','Bold')
    %
    % Model 2: Isotropic Diffusion (Oscillator Relaxation)
    subplot(234), hold on, errorbar(on(2,1,:),on(2,2,:),'b.-'), errorbar(off(2,1,:),off(2,2,:),'r.-')
    set(gca,'FontSize',16,'FontWeight','Bold'), axis([0.5 5+0.5 0 1.05*max(del_mn_std)])
    legend({'on boundary','off boundary'})
    title('Iso','FontSize',20,'FontWeight','Bold')
    xlabel('gT boundary blur')
    ylabel({'mean \Delta','(% of D.R.)'})
    %
    % Model 3: SKH Topographic Modularity (Oscillator Relaxation)
    subplot(232), hold on, errorbar(on(3,1,:),on(3,2,:),'b.-'), errorbar(off(3,1,:),off(3,2,:),'r.-')
    set(gca,'FontSize',16,'FontWeight','Bold'), axis([0.5 5+0.5 0 1.05*max(del_mn_std)])
    title('SK','FontSize',20,'FontWeight','Bold')
    %
    % Model 4: Newman & Girvan Modularity (Oscillator Relaxation)
    subplot(235), hold on, errorbar(on(4,1,:),on(4,2,:),'b.-'), errorbar(off(4,1,:),off(4,2,:),'r.-')
    set(gca,'FontSize',16,'FontWeight','Bold'), axis([0.5 5+0.5 0 1.05*max(del_mn_std)])
    title('NG','FontSize',20,'FontWeight','Bold')
    %
    % Model 5: Average Association (Oscillator Relaxation)
    subplot(233), hold on, errorbar(on(5,1,:),on(5,2,:),'b.-'), errorbar(off(5,1,:),off(5,2,:),'r.-')
    set(gca,'FontSize',16,'FontWeight','Bold'), axis([0.5 5+0.5 0 1.05*max(del_mn_std)])
    title('AA','FontSize',20,'FontWeight','Bold')
    %
    % Model 6: Graph Laplacian (Oscillator Relaxation)
    subplot(236), hold on, errorbar(on(6,1,:),on(6,2,:),'b.-'), errorbar(off(6,1,:),off(6,2,:),'r.-')
    set(gca,'FontSize',16,'FontWeight','Bold'), axis([0.5 5+0.5 0 1.05*max(del_mn_std)])
    title('GL','FontSize',20,'FontWeight','Bold')
    

    saveGoodImg(Hc,[dirImgSave,'Plot1'],sizeGoodIm)
    close(Hc)

end










%% PLOT #2: Plot 6 Different methods with errorbars with no boundary blurring.
if(sav_plt2)
    % Figure with Subplots !
    H2=figure; hold on
    %
    % Model 1: Image Pixels only (straw man model)
    errorbar(1,mean(PixK.meanGradientOnBoundary(1,:)),mean(PixK.stdGradientOnBoundary(1,:)),'b.-')
    errorbar(1,mean(PixK.meanGradientOffBoundary(1,:)),mean(PixK.stdGradientOffBoundary(1,:)),'r.-')
    errorbar(1.2,mean(PixK.meanGradientOnBoundary(1,:)),std(PixK.meanGradientOnBoundary(1,:)),'b.-') 
    errorbar(1.2,mean(PixK.meanGradientOffBoundary(1,:)),std(PixK.meanGradientOffBoundary(1,:)),'r.-')
    %
    % Model 2: Isotropic Diffusion (Oscillator Relaxation)
    errorbar(2,mean(IsoK.meanGradientOnBoundary(1,:)),mean(IsoK.stdGradientOnBoundary(1,:)),'b.-')
    errorbar(2.2,mean(IsoK.meanGradientOnBoundary(1,:)),std(IsoK.meanGradientOnBoundary(1,:)),'b.-') 
    errorbar(2,mean(IsoK.meanGradientOffBoundary(1,:)),mean(IsoK.stdGradientOffBoundary(1,:)),'r.-')
    errorbar(2.2,mean(IsoK.meanGradientOffBoundary(1,:)),std(IsoK.meanGradientOffBoundary(1,:)),'r.-') 
    %
    % Model 3: SKH Topographic Modularity (Oscillator Relaxation)
    errorbar(3,mean(SK_K.meanGradientOnBoundary(1,:)),mean(SK_K.stdGradientOnBoundary(1,:)),'b.-')
    errorbar(3.2,mean(SK_K.meanGradientOnBoundary(1,:)),std(SK_K.meanGradientOnBoundary(1,:)),'b.-') 
    errorbar(3,mean(SK_K.meanGradientOffBoundary(1,:)),mean(SK_K.stdGradientOffBoundary(1,:)),'r.-')
    errorbar(3.2,mean(SK_K.meanGradientOffBoundary(1,:)),std(SK_K.meanGradientOffBoundary(1,:)),'r.-')
    %
    % Model 4: Newman & Girvan Modularity (Oscillator Relaxation)
    errorbar(4,mean(NG_K.meanGradientOnBoundary(1,:)),mean(NG_K.stdGradientOnBoundary(1,:)),'b.-')
    errorbar(4.2,mean(NG_K.meanGradientOnBoundary(1,:)),std(NG_K.meanGradientOnBoundary(1,:)),'b.-') 
    errorbar(4,mean(NG_K.meanGradientOffBoundary(1,:)),mean(NG_K.stdGradientOffBoundary(1,:)),'r.-')
    errorbar(4.2,mean(NG_K.meanGradientOffBoundary(1,:)),std(NG_K.meanGradientOffBoundary(1,:)),'r.-')
    %
    % Model 5: Average Association (Oscillator Relaxation)
    errorbar(5,mean(AA_K.meanGradientOnBoundary(1,:)),mean(AA_K.stdGradientOnBoundary(1,:)),'b.-')
    errorbar(5.2,mean(AA_K.meanGradientOnBoundary(1,:)),std(AA_K.meanGradientOnBoundary(1,:)),'b.-') 
    errorbar(5,mean(AA_K.meanGradientOffBoundary(1,:)),mean(AA_K.stdGradientOffBoundary(1,:)),'r.-')
    errorbar(5.2,mean(AA_K.meanGradientOffBoundary(1,:)),std(AA_K.meanGradientOffBoundary(1,:)),'r.-')
    %
    % Model 6: Graph Laplacian (Oscillator Relaxation)
    errorbar(6,mean(GL_K.meanGradientOnBoundary(1,:)),mean(GL_K.stdGradientOnBoundary(1,:)),'b.-')
    errorbar(6.2,mean(GL_K.meanGradientOnBoundary(1,:)),std(GL_K.meanGradientOnBoundary(1,:)),'b.-') 
    errorbar(6,mean(GL_K.meanGradientOffBoundary(1,:)),mean(GL_K.stdGradientOffBoundary(1,:)),'r.-')
    errorbar(6.2,mean(GL_K.meanGradientOffBoundary(1,:)),std(GL_K.meanGradientOffBoundary(1,:)),'r.-')
    %
    set(gca,'FontSize',16,'FontWeight','Bold','XTick',[1:6],'XTickLabel',{'Pix','Iso','SK','NG','AA','GL'}) 
    axis([0.5 6+0.5 0 1])
   
    title('Average Gradient (150 Image Patches)','FontSize',20,'FontWeight','Bold')
    legend({'on boundary','off boundary'})
    xlabel('Method','FontSize',16,'FontWeight','Bold')
    ylabel('% Dynamic Range','FontSize',16,'FontWeight','Bold')
    grid on
    

    saveGoodImg(H2,[dirImgSave,'Plot2'],sizeGoodIm)
    close(H2)
end



%% PLOT #3: Plot (a). Average Gradient at boundaries in ground truth w/ error bars and 
%                (b). Ratio of Avg Gradient Step at Boundaries vs within Segments.
if(sav_plt3)

    % Plot ratio of mean gradients for different methods (on boundary vs off)
    gg1 = [M.im(:,1)./S.im(:,1), ...
          M.iso(:,1)./S.iso(:,1), ...
          M.sk(:,1)./S.sk(:,1), ...
          M.ng(:,1)./S.ng(:,1), ...
          M.gl(:,1)./S.gl(:,1), ...
          M.aa(:,1)./S.aa(:,1)];

    Hl = figure; 
    subplot(121), plot(gg1,'LineWidth',2)
    title('Ratio of Mean Gradient Value On vs Off Boundaries','Fontweight','Bold','FontSize',20)
    xlabel('blur','Fontweight','Bold','FontSize',18)
    ylabel('M_{B}:/:M_{~B}','Fontweight','Bold','FontSize',18)
    set(gca,'FontSize',16,'FontWeight','Bold')
    %legend({'pix','iso','sk','ng','gl','aa'}) 
    text(size(gg1,1)-1, max(gg1(:)), [ptch_str],'Fontweight','Bold','FontSize',16,'HorizontalAlignment','Center')



    gg2 = [M.im(:,1), ...
          M.iso(:,1)./pi, ...
          M.sk(:,1)./pi, ...
          M.ng(:,1)./pi, ...
          M.gl(:,1)./pi, ...
          M.aa(:,1)./pi];

    subplot(122), plot(gg2,'LineWidth',2)
    title('Mean Gradient Value On Boundaries','Fontweight','Bold','FontSize',20)
    xlabel('blur','Fontweight','Bold','FontSize',18)
    ylabel('M_{B} (% of Dynamic Range)','Fontweight','Bold','FontSize',18)
    set(gca,'FontSize',16,'FontWeight','Bold')
    legend({'pix','iso','sk','ng','gl','aa'}) 

    saveGoodImg(Hl,[dirImgSave,ptch,'_grads_OnNoffBoundary'],sizeGoodIm)
    close(Hl)

end
       






%% PLOT #4: Lots of plots:  Full Histograms of Boundary Discriminability (d') 
%  for each method and Discriminability Improvement over Image Pixels and
%  also mean & std (over # image patches) of d' for different methods.

if(sav_plt4)
    
    % Plot full Histograms for 6 methods for each blur.
    bins = [-0.5:0.05:1.5];
    
    for i = 1:5

        figure,
        subplot(321), hist( PixK.boundaryDiscriminability(i,:), bins ), title('Pix d''','FontSize',20,'FontWeight','Bold'), xlim([min(bins), max(bins)]), set(gca,'FontSize',16,'FontWeight','Bold','YTick',[])
        subplot(322), hist( IsoK.boundaryDiscriminability(i,:), bins ), title('Iso d''','FontSize',20,'FontWeight','Bold'), xlim([min(bins), max(bins)]), set(gca,'FontSize',16,'FontWeight','Bold','YTick',[])
        subplot(323), hist( SK_K.boundaryDiscriminability(i,:), bins ), title('SK d''','FontSize',20,'FontWeight','Bold'), xlim([min(bins), max(bins)]), set(gca,'FontSize',16,'FontWeight','Bold','YTick',[])
        subplot(324), hist( NG_K.boundaryDiscriminability(i,:), bins ), title('NG d''','FontSize',20,'FontWeight','Bold'), xlim([min(bins), max(bins)]), set(gca,'FontSize',16,'FontWeight','Bold','YTick',[])
        subplot(325), hist( GL_K.boundaryDiscriminability(i,:), bins ), title('GL d''','FontSize',20,'FontWeight','Bold'), xlim([min(bins), max(bins)]), set(gca,'FontSize',16,'FontWeight','Bold','YTick',[]), ylabel('150 Im Patches','FontSize',18,'FontWeight','Bold')
        subplot(326), hist( AA_K.boundaryDiscriminability(i,:), bins ), title('AA d''','FontSize',20,'FontWeight','Bold'), xlim([min(bins), max(bins)]), set(gca,'FontSize',16,'FontWeight','Bold','YTick',[])

    end
    
    
    % Plot blur vs mean/std d' with errorbars
    figure, hold on,
    errorbar(mean(PixK.boundaryDiscriminability'), std(PixK.boundaryDiscriminability'),'b.-')
    errorbar(mean(IsoK.boundaryDiscriminability'), std(IsoK.boundaryDiscriminability'),'g.-')
    errorbar(mean(SK_K.boundaryDiscriminability'), std(SK_K.boundaryDiscriminability'),'r.-')
    errorbar(mean(NG_K.boundaryDiscriminability'), std(NG_K.boundaryDiscriminability'),'c.-')
    errorbar(mean(GL_K.boundaryDiscriminability'), std(GL_K.boundaryDiscriminability'),'m.-')
    errorbar(mean(AA_K.boundaryDiscriminability'), std(AA_K.boundaryDiscriminability'),'y.-')
    legend({'Pix','Iso','SK','NG','GL','AA'})
    
    % Plot blur vs d' mean & std separately
    figure, 
    subplot(211), hold on,
    plot(mean(PixK.boundaryDiscriminability'),'b.-','LineWidth',2)
    plot(mean(IsoK.boundaryDiscriminability'),'g.-','LineWidth',2)
    plot(mean(SK_K.boundaryDiscriminability'),'r.-','LineWidth',2)
    plot(mean(NG_K.boundaryDiscriminability'),'c.-','LineWidth',2)
    plot(mean(GL_K.boundaryDiscriminability'),'m.-','LineWidth',2)
    plot(mean(AA_K.boundaryDiscriminability'),'y.-','LineWidth',2)
    ylabel('mean(d'')','FontSize',18,'FontWeight','Bold')
    set(gca,'FontSize',16,'FontWeight','Bold','XTick',[1:size(PixK.boundaryDiscriminability,1)])
    grid on
    title(['Boundary Discriminability (d'') Metric - Across ',num2str(size(AA_K.boundaryDiscriminability,2)),' Image Patches'],'FontSize',20,'FontWeight','Bold')
    %
    subplot(212), hold on,
    plot(std(PixK.boundaryDiscriminability'),'b.-','LineWidth',2)
    plot(std(IsoK.boundaryDiscriminability'),'g.-','LineWidth',2)
    plot(std(SK_K.boundaryDiscriminability'),'r.-','LineWidth',2)
    plot(std(NG_K.boundaryDiscriminability'),'c.-','LineWidth',2)
    plot(std(GL_K.boundaryDiscriminability'),'m.-','LineWidth',2)
    plot(std(AA_K.boundaryDiscriminability'),'y.-','LineWidth',2)
    ylabel('std(d'')','FontSize',18,'FontWeight','Bold')
    xlabel(' Progressive Blur Applied to Ground Truth Boundary','FontSize',18,'FontWeight','Bold')
    set(gca,'FontSize',16,'FontWeight','Bold','XTick',[1:size(PixK.boundaryDiscriminability,1)],'XTickLabel',PixK.bDc_blurs_info)
    grid on
    legend({'Pix','Iso','SK','NG','GL','AA'})

    
    % plot SK d' full histograms for different blurs
    figure,
    subplot(511), hist( SK_K.boundaryDiscriminability(1,:), bins ), ylim([0 20]), ylabel('blur1'),title('SK d'' Distribution (150 img patches)')
    subplot(512), hist( SK_K.boundaryDiscriminability(2,:), bins ), ylim([0 20]), ylabel('blur2')
    subplot(513), hist( SK_K.boundaryDiscriminability(3,:), bins ), ylim([0 20]), ylabel('blur3')
    subplot(514), hist( SK_K.boundaryDiscriminability(4,:), bins ), ylim([0 20]), ylabel('blur4')
    subplot(515), hist( SK_K.boundaryDiscriminability(5,:), bins ), ylim([0 20]), ylabel('blur5')
    
    
    
    % plot SK d' improvement over img pixels full histograms for different blurs
    figure,
    subplot(511), hist( SK_K.boundaryDiscriminability(1,:) - PixK.boundaryDiscriminability(1,:), bins ), ylim([0 20]), ylabel('blur1'),title('SK d'' Improvement over Pixels Distribution (150 img patches)')
    subplot(512), hist( SK_K.boundaryDiscriminability(2,:) - PixK.boundaryDiscriminability(2,:), bins ), ylim([0 20]), ylabel('blur2')
    subplot(513), hist( SK_K.boundaryDiscriminability(3,:) - PixK.boundaryDiscriminability(3,:), bins ), ylim([0 20]), ylabel('blur3')
    subplot(514), hist( SK_K.boundaryDiscriminability(4,:) - PixK.boundaryDiscriminability(4,:), bins ), ylim([0 20]), ylabel('blur4')
    subplot(515), hist( SK_K.boundaryDiscriminability(5,:) - PixK.boundaryDiscriminability(5,:), bins ), ylim([0 20]), ylabel('blur5')
    
    keyboard
    
end



%% PLOT # 5: Make a butterfly plot or a scatter plot that compares directly d' of method vs pixels.
if(sav_plt_butterfly)
    
    which_single_patch = '118015_ptch2';
    indicate_scatter_single_patch=1;
    show_single_gradient_inset=0;
    
    
    % find which data point pertains to the single patch of interest.
    for i = 1:numel(AA_K.ImgPtchID)
        x=strmatch(which_single_patch,AA_K.ImgPtchID{i});
        if(x)
           patch_ind = i; 
        end
    end
    
    
    
    
    
    
    
    if(show_single_gradient_inset)
        
        
    end
    
    
    
    
    

    for i=[which_bD]  % 1:numel(PixK.bDc_blurs_info)
        

        
        % Find average improvement of d' after method from image pixels.
        ind = find( (PixK.boundaryDiscriminability(i,indPixAA)>=0) & (AA_K.boundaryDiscriminability(i,indAA)>=0) );
        disp('Number of AA & Pix positive')
        numel(ind)
        mn_del_dp.AA = mean( (AA_K.boundaryDiscriminability(i,indAA(ind)) - PixK.boundaryDiscriminability(i,indPixAA(ind))) );
        %
        ind = find( (PixK.boundaryDiscriminability(i,indPixNG)>=0) & (NG_K.boundaryDiscriminability(i,indNG)>=0) );
        disp('Number of NG & Pix positive')
        numel(ind)
        mn_del_dp.NG = mean( (NG_K.boundaryDiscriminability(i,indNG(ind)) - PixK.boundaryDiscriminability(i,indPixNG(ind))) );
        %
        ind = find( (PixK.boundaryDiscriminability(i,indPixGL)>=0) & (GL_K.boundaryDiscriminability(i,indGL)>=0) );
        numel(ind)
        mn_del_dp.GL = mean( (GL_K.boundaryDiscriminability(i,indGL(ind)) - PixK.boundaryDiscriminability(i,indPixGL(ind))) );
        %
        ind = find( (PixK.boundaryDiscriminability(i,indPixSK)>=0) & (SK_K.boundaryDiscriminability(i,indSK)>=0) );
        numel(ind)
        mn_del_dp.SK = mean( (SK_K.boundaryDiscriminability(1,indSK(ind)) - PixK.boundaryDiscriminability(1,indPixSK(ind))) );
        %
        ind = find( (PixK.boundaryDiscriminability(i,indPixIso)>=0) & (IsoK.boundaryDiscriminability(i,indIso)>=0) );
        numel(ind)
        mn_del_dp.Iso = mean( (IsoK.boundaryDiscriminability(i,indIso(ind)) - PixK.boundaryDiscriminability(i,indPixIso(ind))) );
        
        
        
        
        Hc=figure;
        ha = tight_subplot(2, 3, [0.05 0], 0.05, 0.01);
        maxx = 1.05.*max([PixK.boundaryDiscriminability(i,:),SK_K.boundaryDiscriminability(i,:),NG_K.boundaryDiscriminability(i,:),AA_K.boundaryDiscriminability(i,:),GL_K.boundaryDiscriminability(i,:),IsoK.boundaryDiscriminability(i,:),]);
        minn = 1.00.*min([PixK.boundaryDiscriminability(i,:),SK_K.boundaryDiscriminability(i,:),NG_K.boundaryDiscriminability(i,:),AA_K.boundaryDiscriminability(i,:),GL_K.boundaryDiscriminability(i,:),IsoK.boundaryDiscriminability(i,:),]);
        set(ha,'FontSize',16,'FontWeight','Bold','XTick',[floor(minn):1:ceil(maxx)],'YTick',[floor(minn):1:ceil(maxx)],'XtickLabel',[floor(minn):1:ceil(maxx)],'YtickLabel',[floor(minn):1:ceil(maxx)])

        
        
        
        % #1. Average Association (AA)
        axes(ha(1)), hold on
        scatter(PixK.boundaryDiscriminability(i,indPixIso),IsoK.boundaryDiscriminability(i,indIso),50,'ko')
        scatter(PixK.boundaryDiscriminability(i,indPixAA),AA_K.boundaryDiscriminability(i,indAA),100,'m.')
        plot([minn maxx], [minn maxx],'k--','LineWidth',2)

        if(indicate_scatter_single_patch)
            scatter(PixK.boundaryDiscriminability(i,patch_ind),AA_K.boundaryDiscriminability(i,patch_ind),300,'kd','filled')
            scatter(PixK.boundaryDiscriminability(i,patch_ind),AA_K.boundaryDiscriminability(i,patch_ind),100,'m.')
        end

        axis([minn maxx minn maxx])
        axis square
        %xlabel(['d'' Pix'],'FontSize',18,'FontWeight','Bold')
        ylabel(['d'' \color{magenta}AA'],'FontSize',18,'FontWeight','Bold')
        text( maxx, minn, {['rM=',num2str(AA_K.rM)],['Ks=',AA_K.KS]},'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',16,'FontWeight','Bold','EdgeColor','k')
        text( maxx, (1/3)*maxx, ['< \Deltad'' > =',num2str(mn_del_dp.AA,2)],'HorizontalAlignment','Right','VerticalAlignment','Middle','FontSize',16,'FontWeight','Bold','EdgeColor','m','LineWidth',2)
        grid on
        
        if(show_single_gradient_inset)
            
            
            
            
        end
        
        
        
        % #2. Non-Topographic Modularity (NG)
        axes(ha(2)), hold on
        scatter(PixK.boundaryDiscriminability(i,indPixIso),IsoK.boundaryDiscriminability(i,indIso),50,'ko')
        scatter(PixK.boundaryDiscriminability(i,indPixNG),NG_K.boundaryDiscriminability(i,indNG),100,'c.')
        plot([minn maxx], [minn maxx],'k--','LineWidth',2)
        
        if(indicate_scatter_single_patch)
            scatter(PixK.boundaryDiscriminability(i,patch_ind),NG_K.boundaryDiscriminability(i,patch_ind),300,'kd','filled')
            scatter(PixK.boundaryDiscriminability(i,patch_ind),NG_K.boundaryDiscriminability(i,patch_ind),100,'c.')
        end
        
        axis([minn maxx minn maxx])
        axis square
        title([PixK.bDc_blurs_info(i)],'FontSize',20,'FontWeight','Bold')
        %xlabel(['d'' Pix'],'FontSize',18,'FontWeight','Bold')
        ylabel(['d'' \color{cyan}NG'],'FontSize',18,'FontWeight','Bold')
        text( maxx, minn, {['rM=',num2str(NG_K.rM)],['Ks=',NG_K.KS]},'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',16,'FontWeight','Bold','EdgeColor','k')
        text( maxx, (1/3)*maxx, ['< \Deltad'' > =',num2str(mn_del_dp.NG,2)],'HorizontalAlignment','Right','VerticalAlignment','Middle','FontSize',16,'FontWeight','Bold','EdgeColor','c','LineWidth',2)
        grid on
        


        % # 3. Isotropic Diffusion (Iso)
        axes(ha(3)), hold on
        scatter(PixK.boundaryDiscriminability(i,indPixIso),IsoK.boundaryDiscriminability(i,indIso),50,'ko')
        scatter(PixK.boundaryDiscriminability(i,indPixIso),IsoK.boundaryDiscriminability(i,indIso),100,'b.')
        plot([minn maxx], [minn maxx],'k--','LineWidth',2)
        
        if(indicate_scatter_single_patch)
            scatter(PixK.boundaryDiscriminability(i,patch_ind),IsoK.boundaryDiscriminability(i,patch_ind),300,'kd','filled')
            scatter(PixK.boundaryDiscriminability(i,patch_ind),IsoK.boundaryDiscriminability(i,patch_ind),100,'b.')
        end
        
        axis([minn maxx minn maxx])
        axis square
        xlabel(['d'' Pix'],'FontSize',18,'FontWeight','Bold')
        ylabel(['d'' \color{blue}Iso'],'FontSize',18,'FontWeight','Bold')
        text( maxx, minn, {['rM=',num2str(IsoK.rM)],['Ks=',IsoK.KS]},'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',16,'FontWeight','Bold','EdgeColor','k')
        text( maxx, (1/3)*maxx, ['< \Deltad'' > =',num2str(mn_del_dp.Iso,2)],'HorizontalAlignment','Right','VerticalAlignment','Middle','FontSize',16,'FontWeight','Bold','EdgeColor','b','LineWidth',2)
        grid on


        % #4. Graph Laplacian (GL)
        axes(ha(4)), hold on
        scatter(PixK.boundaryDiscriminability(i,indPixIso),IsoK.boundaryDiscriminability(i,indIso),50,'ko')
        scatter(PixK.boundaryDiscriminability(i,indPixGL),GL_K.boundaryDiscriminability(i,indGL),100,'g.')
        plot([minn maxx], [minn maxx],'k--','LineWidth',2)
        
        if(indicate_scatter_single_patch)
            scatter(PixK.boundaryDiscriminability(i,patch_ind),GL_K.boundaryDiscriminability(i,patch_ind),300,'kd','filled')
            scatter(PixK.boundaryDiscriminability(i,patch_ind),GL_K.boundaryDiscriminability(i,patch_ind),100,'g.')
        end
        
        axis([minn maxx minn maxx])
        axis square
        %xlabel(['d'' Pix'],'FontSize',18,'FontWeight','Bold')
        ylabel(['d'' \color{green}GL'],'FontSize',18,'FontWeight','Bold')
        text( maxx, minn, {['rM=',num2str(GL_K.rM)],['Ks=',GL_K.KS]},'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',16,'FontWeight','Bold','EdgeColor','k')
        text( maxx, (1/3)*maxx, ['< \Deltad'' > =',num2str(mn_del_dp.GL,2)],'HorizontalAlignment','Right','VerticalAlignment','Middle','FontSize',16,'FontWeight','Bold','EdgeColor','g','LineWidth',2)
        grid on
        
        
        % #5. Topographic Modularity (SK)
        axes(ha(5)), hold on
        scatter(PixK.boundaryDiscriminability(i,indPixIso),IsoK.boundaryDiscriminability(i,indIso),50,'ko')
        scatter(PixK.boundaryDiscriminability(i,indPixSK),SK_K.boundaryDiscriminability(i,indSK),100,'r.')
        plot([minn maxx], [minn maxx],'k--','LineWidth',2)
        
        if(indicate_scatter_single_patch)
            scatter(PixK.boundaryDiscriminability(i,patch_ind),SK_K.boundaryDiscriminability(i,patch_ind),300,'kd','filled')
            scatter(PixK.boundaryDiscriminability(i,patch_ind),SK_K.boundaryDiscriminability(i,patch_ind),100,'r.')
        end
        
        axis([minn maxx minn maxx])
        axis square
        %xlabel(['d'' Pix'],'FontSize',18,'FontWeight','Bold')
        ylabel(['d'' \color{red}SK'],'FontSize',18,'FontWeight','Bold')
        text( maxx, minn, {['rM=',num2str(SK_K.rM)],['Ks=',SK_K.KS]},'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',16,'FontWeight','Bold','EdgeColor','k')
        text( maxx, (1/3)*maxx, ['< \Deltad'' > =',num2str(mn_del_dp.SK,2)],'HorizontalAlignment','Right','VerticalAlignment','Middle','FontSize',16,'FontWeight','Bold','EdgeColor','r','LineWidth',2)
        grid on

        
        
        
        % #6. Clear the 6th axis
        axes(ha(6)),
        axis off
        
        
        saveGoodImg(Hc,[dirImgSave,'Butterfly_optimized_params_bD',num2str(which_bD)],sizeGoodIm)
        close(Hc)
        
        
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