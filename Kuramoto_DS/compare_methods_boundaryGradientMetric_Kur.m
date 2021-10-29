function [] = compare_methods_boundaryGradientMetric_Kur(ptch)

% weird.



% save flags.
sav_degDists = 0;
sav_mov0 = 0;
sav_plt1 = 1;
sav_plt1p5 = 0;
sav_plt2 = 0;
sav_mat = 0;


[dirPre,sizeGoodIm] = onCluster;

dirImgSave = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/imgs/Kur_PIF_Fourier1/compareMethods/'];
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
IsoK = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/IsoDiff/KurMC_',ptch,'_rM1_NF_60_0_ksmid.mat']);
%IsoE = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/spectral/IsoDiff/Evecs_',ptch,'_rM1.mat']);
%
SK_K = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/Mod_SKHAdj/KurMC_',ptch,'_rM3_sDInf_sP0p2_NF_60_0_ksmid.mat']);
%
NG_K = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/Mod_N&G/KurMC_',ptch,'_rM3_sDInf_sP0p2_NF_60_0_ksmid.mat']);
%
AA_K = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/AAnrm/KurMC_',ptch,'_rM3_sDInf_sP0p2_NF_60_0_ksmid.mat']);
%
GL_K = load([dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/GLnrm/KurMC_',ptch,'_rM3_sDInf_sP0p2_NF_60_0_kslrg.mat']);


% PUT PARAMS RM & KS ON PLOTS !!!

GTFile = load([dirPre,'images/BSDS_patch/101x101_ds1/',ptch,'.mat']); % Ground truth file containing Boundary informatuon



% Investigate modifying ground truth by requiring consensus & blurring boundaries.
bD = GTFile.bD;
MC.bDc_blurs_info = {'(bD>0)','(bD>0)blur_binom2','(bD>0)blur_binom3','(bD>1)','(bD>1)blur_binom2','(bD>1)blur_binom3'};
bDc_blurs(:,:,1) = (bD>0);
bDc_blurs(:,:,2) = logical( blur( single(bD>0),1,'binom2' ) );
bDc_blurs(:,:,3) = logical( blur( single(bD>0),1,'binom3' ) );
bDc_blurs(:,:,4) = (bD>1);
bDc_blurs(:,:,5) = logical( blur( single(bD>1),1,'binom2' ) );
bDc_blurs(:,:,6) = logical( blur( single(bD>1),1,'binom3' ) ); % note: 'power' increases with increased blurring





% As a sanity check, create a segmentation that consists of a
% blurred version of the ground truth.  What do M & S look like?
if(0)
    disp('Sanity Check: Blur 1')
    X = single(bD_consensus);
    [F,M.b1,S.b1] = compute_BoundaryGradientMetric(X,bDc_blurs);
    %
    disp('Sanity Check: Blur 2')
    X = single(bD2c);
    [F,M.b2,S.b2] = compute_BoundaryGradientMetric(X,bDc_blurs);
    %
    disp('Sanity Check: Blur 3')
    X = single(bD3c);
    [F,M.b3,S.b3] = compute_BoundaryGradientMetric(X,bDc_blurs);
    %
    disp('Sanity Check: Blur 4')
    X = single(bD4c);
    [F,M.b4,S.b4] = compute_BoundaryGradientMetric(X,bDc_blurs);
    %
    disp('Sanity Check: Blur 5clr')
    X = single(bD5c);
    [F,M.b5,S.b5] = compute_BoundaryGradientMetric(X,bDc_blurs);

    % Plot ratio of mean gradients for different methods (on boundary vs off)
    gg = [M.b1(:,1)./S.b1(:,1), ...
          M.b2(:,1)./S.b2(:,1), ...
          M.b3(:,1)./S.b3(:,1), ...
          M.b4(:,1)./S.b4(:,1), ...
          M.b5(:,1)./S.b5(:,1)];

    figure, 
    subplot(121), plot(gg,'LineWidth',2)
    title('Ratio of Mean Gradient Value On vs Off Boundaries')
    xlabel('blur')
    %legend({'pix','iso','sk','ng','gl','aa'}) 



    gg = [M.b1(:,1), ...
          M.b2(:,1), ...
          M.b3(:,1), ...
          M.b4(:,1), ...
          M.b5(:,1)];

    subplot(122), plot(gg,'LineWidth',2)
    title('Mean Gradient Value On Boundaries')
    xlabel('blur')
    legend({'b1','2','3','4','5'}) 

end



%% Plot -1: For a given Image patch, look at the degree distribution
if(sav_degDists)
    
    figure

    deg = sum(IsoK.netParams.Q);
    subplot(231), imagesc(reshape(deg,101,101)), 
    colormap('jet'), colorbar, freezeColors, cbfreeze
    set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','Bold')
    axis square
    title([],'FontSize',20,'FontWeight','Bold')
    %
    deg = sum(SK_K.netParams.Q);
    subplot(232), imagesc(reshape(deg,101,101)), 
    colormap('jet'), colorbar, freezeColors, cbfreeze
    set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','Bold'), 
    axis square
    title('SK','FontSize',20,'FontWeight','Bold')
    %
    deg = sum(NG_K.netParams.Q);
    subplot(233), imagesc(reshape(deg,101,101)), 
    colormap('jet'), colorbar, freezeColors, cbfreeze
    set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','Bold'), 
    axis square
    title('NG','FontSize',20,'FontWeight','Bold')
    %
    deg = sum(AA_K.netParams.Q);
    subplot(234), imagesc(reshape(deg,101,101)), 
    colormap('jet'), colorbar, freezeColors, cbfreeze
    set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','Bold'), 
    axis square
    title('AA','FontSize',20,'FontWeight','Bold')
    %
    deg = sum(GL_K.netParams.Q);
    subplot(235), imagesc(reshape(deg,101,101)), 
    colormap('jet'), colorbar, freezeColors, cbfreeze
    set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','Bold'), 
    axis square
    title('GL','FontSize',20,'FontWeight','Bold')
    %
    subplot(4,3,9), imagesc(GTFile.im), colormap('bone'), freezeColors,
    set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','Bold')
    axis square
    title('Image','FontSize',18,'FontWeight','Bold')
    %
    subplot(4,3,12), imagesc(max(GTFile.bD(:)) - GTFile.bD), colormap('bone'), freezeColors,
    set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','Bold')
    axis square
    title('GT Boundaries','FontSize',18,'FontWeight','Bold')
    
    
end



%% Plot #0:  Make a movie of phase relaxation
if(sav_mov0)
    
    makeMovie_phaseEvolution([dirImgSave,ptch,'_SK_v_Iso_phase'],GTFile,SK_K,IsoK)
    
    makeMovie_phaseEvolution([dirImgSave,ptch,'_NG_v_Iso_phase'],GTFile,NG_K,IsoK)
    
    makeMovie_phaseEvolution([dirImgSave,ptch,'_GL_v_Iso_phase'],GTFile,GL_K,IsoK)
    
    makeMovie_phaseEvolution([dirImgSave,ptch,'_AA_v_Iso_phase'],GTFile,AA_K,IsoK)
    
    % makeMovie_phaseEvolution([dirImgSave,ptch,'_Iso_phase'],GTFile,IsoK)
    
    
    %keyboard
    
end


%% PLOT #1: Visualize Gradients (phase,contrast or otherwise)
% for different competitor methods and plot avg gradient/pixel
% on boundaries vs. off for different boundary blurs.
if(sav_plt1)

    numTHs=20;


    % Compute Mean & Std of Gradient Fields On & Off Ground Truth Boundaries
    X = GTFile.im;
    %[F.im,M.im,S.im] = compute_BoundaryGradientMetric(X,bDc_blurs,0,'Img Pix');
    Ph.im = X;
    %[precision.im,recall.im,f_measure.im] = compute_BoundaryGradientPrecRec(X,bDc_blurs,numTHs,'Img Pix');
    %
    X = visKurPhase_inHSV(IsoK.netParams.im, reshape(IsoK.metaCluster.phaseAtClk(:,end),IsoK.netParams.Ndims));
    %[F.iso,M.iso,S.iso] = compute_BoundaryGradientMetric(X,bDc_blurs,1,'Iso Diff');
    Ph.iso = X;
    %[precision.iso,recall.iso,f_measure.iso] = compute_BoundaryGradientPrecRec(X,bDc_blurs,numTHs,'Iso Diff');
    %
    X = visKurPhase_inHSV(NG_K.netParams.im, reshape(NG_K.metaCluster.phaseAtClk(:,end),NG_K.netParams.Ndims));
    %[F.ng,M.ng,S.ng] = compute_BoundaryGradientMetric(X,bDc_blurs,1,'Modularity NG');
    Ph.ng =X;
    %[precision.ng,recall.ng,f_measure.ng] = compute_BoundaryGradientPrecRec(X,bDc_blurs,numTHs,'Modularity NG');
    %
    X = visKurPhase_inHSV(SK_K.netParams.im, reshape(SK_K.metaCluster.phaseAtClk(:,end),SK_K.netParams.Ndims));
    %[F.sk,M.sk,S.sk] = compute_BoundaryGradientMetric(X,bDc_blurs,1,'Modularity SK');
    Ph.sk = X;
    %[precision.sk,recall.sk,f_measure.sk] = compute_BoundaryGradientPrecRec(X,bDc_blurs,numTHs,'Modularity SK');
    %
    X = visKurPhase_inHSV(GL_K.netParams.im, reshape(GL_K.metaCluster.phaseAtClk(:,end),GL_K.netParams.Ndims));
    %[F.gl,M.gl,S.gl] = compute_BoundaryGradientMetric(X,bDc_blurs,1,'GL');
    Ph.gl = X;
    %[precision.gl,recall.gl,f_measure.gl] = compute_BoundaryGradientPrecRec(X,bDc_blurs,numTHs,'GL');
    %
    X = visKurPhase_inHSV(AA_K.netParams.im, reshape(AA_K.metaCluster.phaseAtClk(:,end),AA_K.netParams.Ndims));
    %[F.aa,M.aa,S.aa] = compute_BoundaryGradientMetric(X,bDc_blurs,1,'AA');
    Ph.aa = X;
    %[precision.aa,recall.aa,f_measure.aa] = compute_BoundaryGradientPrecRec(X,bDc_blurs,numTHs,'AA');
    
    
    
    
% % NOTE: THese are incorrect right now because STD on bottom are not squared.    
%     % From M (mean & std on boundary) and S (mean & std off boundary), compute d' for each blur.
%     PixK.d_prime = ( M.im(:,1) - S.im(:,1) ) ./ sqrt( 0.5*( M.im(:,2) + S.im(:,2) ) );
%     IsoK.d_prime = ( M.iso(:,1) - S.iso(:,1) ) ./ sqrt( 0.5*( M.iso(:,2) + S.iso(:,2) ) );
%     AA_K.d_prime = ( M.aa(:,1) - S.aa(:,1) ) ./ sqrt( 0.5*( M.aa(:,2) + S.aa(:,2) ) );
%     GL_K.d_prime = ( M.gl(:,1) - S.gl(:,1) ) ./ sqrt( 0.5*( M.gl(:,2) + S.gl(:,2) ) );
%     NG_K.d_prime = ( M.ng(:,1) - S.ng(:,1) ) ./ sqrt( 0.5*( M.ng(:,2) + S.ng(:,2) ) );
%     SK_K.d_prime = ( M.sk(:,1) - S.sk(:,1) ) ./ sqrt( 0.5*( M.sk(:,2) + S.sk(:,2) ) );
    
    
    


%     % Find max value for y-axis of error bar plots (M_mn+M_std or S_mn+S_std)
%     del_mn_std = [sum(M.im,2); sum(M.iso,2)./pi; sum(M.ng,2)./pi; sum(M.sk,2)./pi; sum(M.gl,2)./pi; sum(M.aa,2)./pi;...
%         sum(S.im,2); sum(S.iso,2)./pi; sum(S.ng,2)./pi; sum(S.sk,2)./pi; sum(S.gl,2)./pi; sum(S.aa,2)./pi;...
%         PixK.d_prime; IsoK.d_prime; AA_K.d_prime; GL_K.d_prime; NG_K.d_prime; SK_K.d_prime];


    % Figure with Subplots !
    
    
    showImgNgT = 1; % flag to include Image Patch and Ground Truth Boundaries in this plot.
    
    if(showImgNgT)
        subplotsW = 4;
        gap_h = 0.02;
    else
        subplotsW = 3;
        gap_h = 0;
    end
    
%     % for plotting phases        
%     F = Ph;                         
%     cmap = 'hsv';      
%     ftag = '_visPhases';  
    
    
    % for plotting gradients
    cmap = 'jet';
    ftag = '_visGradients';
    
    
    Hc=figure;
    ha = tight_subplot(2, subplotsW, [0.05 gap_h], 0.05, 0.02);
    
    % Model 1: Average Association (Oscillator Relaxation)
    axes(ha(1)), imagesc(AA_K.MC.F), colormap(cmap), axis square, title(['\color{magenta}Average Association (AA)'],'FontSize',20,'FontWeight','Bold') 
    set(gca,'XTick',[],'YTick',[],'XColor','m','YColor','m','LineWidth',3), freezeColors, % colorbar, cbfreeze
    text(AA_K.netParams.Ndims(1),AA_K.netParams.Ndims(2),['\color{white}d''= ',num2str(AA_K.MC.D(1),2)],'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',18,'FontWeight','Bold')
    ylabel(['( \color{blue}{min = ',num2str( min(AA_K.MC.F(:)) ,2),'} \color{black}{,} \color{red}{max = ',num2str( max(AA_K.MC.F(:)) ,2),'}\color{black}{)}'],'FontSize',18,'FontWeight','Bold')

    % Model 3: Newman & Girvan Modularity (Oscillator Relaxation)
    axes(ha(2)), imagesc(NG_K.MC.F), colormap(cmap), axis square, title(['\color{cyan}Non-Topo Modularity (NG)'],'FontSize',20,'FontWeight','Bold')
    set(gca,'XTick',[],'YTick',[],'XColor','c','YColor','c','LineWidth',3), freezeColors, % colorbar, cbfreeze 
    text(NG_K.netParams.Ndims(1),NG_K.netParams.Ndims(2),['\color{white}d''= ',num2str(NG_K.MC.D(1),2)],'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',18,'FontWeight','Bold')
    ylabel(['( \color{blue}{min = ',num2str( min(NG_K.MC.F(:)) ,2),'} \color{black}{,} \color{red}{max = ',num2str( max(NG_K.MC.F(:)) ,2),'}\color{black}{)}'],'FontSize',18,'FontWeight','Bold')
    
    % Model 5: Isotropic Diffusion (Oscillator Relaxation)
    axes(ha(3)), imagesc(IsoK.MC.F), colormap(cmap),  axis square, title(['\color{blue}Isotropic Diffusion (Iso)'],'FontSize',20,'FontWeight','Bold') 
    set(gca,'XTick',[],'YTick',[],'XColor','b','YColor','b','LineWidth',3), freezeColors, % colorbar, cbfreeze
    text(IsoK.netParams.Ndims(1),IsoK.netParams.Ndims(2),['\color{white}d''= ',num2str(IsoK.MC.D(1),2)],'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',18,'FontWeight','Bold')
    ylabel(['( \color{blue}{min = ',num2str( min(IsoK.MC.F(:)) ,2),'} \color{black}{,} \color{red}{max = ',num2str( max(IsoK.MC.F(:)) ,2),'}\color{black}{)}'],'FontSize',18,'FontWeight','Bold')
    
    % Model 2: Graph Laplacian (Oscillator Relaxation)
    axes(ha(subplotsW+1)), imagesc(GL_K.MC.F), colormap(cmap), axis square, title(['\color{green}Graph Laplacian (GL)'],'FontSize',20,'FontWeight','Bold') 
    set(gca,'XTick',[],'YTick',[],'XColor','g','YColor','g','LineWidth',3), freezeColors, % colorbar, cbfreeze
    text(GL_K.netParams.Ndims(1),GL_K.netParams.Ndims(2),['\color{white}d''= ',num2str(GL_K.MC.D(1),2)],'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',18,'FontWeight','Bold')
    ylabel(['( \color{blue}{min = ',num2str( min(GL_K.MC.F(:)) ,2),'} \color{black}{,} \color{red}{max = ',num2str( max(GL_K.MC.F(:)) ,2),'}\color{black}{)}'],'FontSize',18,'FontWeight','Bold')
    
    % Model 4: SKH Topographic Modularity (Oscillator Relaxation)
    axes(ha(subplotsW+2)), imagesc(SK_K.MC.F), colormap(cmap), axis square, title(['\color{red}Topographic Modularity (SK)'],'FontSize',20,'FontWeight','Bold') 
    set(gca,'XTick',[],'YTick',[],'XColor','r','YColor','r','LineWidth',3), freezeColors, % colorbar, cbfreeze
    text(SK_K.netParams.Ndims(1),SK_K.netParams.Ndims(2),['\color{white}d''= ',num2str(SK_K.MC.D(1),2)],'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',18,'FontWeight','Bold')
    ylabel(['( \color{blue}{min = ',num2str( min(SK_K.MC.F(:)) ,2),'} \color{black}{,} \color{red}{max = ',num2str( max(SK_K.MC.F(:)) ,2),'}\color{black}{)}'],'FontSize',18,'FontWeight','Bold')
    
    % Model 6: Just Image Pixels
    axes(ha(subplotsW+3)), imagesc(GTFile.MC.F), colormap(cmap), axis square, title(['\color{black}Image Pixels (Pix)'],'FontSize',20,'FontWeight','Bold'),
    set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3), freezeColors, % colorbar, cbfreeze
    text(AA_K.netParams.Ndims(1),AA_K.netParams.Ndims(2),['\color{white}d''= ',num2str(GTFile.MC.D(1),2)],'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',18,'FontWeight','Bold')
    ylabel(['( \color{blue}{min = ',num2str( min(GTFile.MC.F(:)) ,2),'} \color{black}{,} \color{red}{max = ',num2str( max(GTFile.MC.F(:)) ,2),'}\color{black}{)}'],'FontSize',18,'FontWeight','Bold')
    
    
    if(showImgNgT)
        
        ptch_str = ptch;
        ptch_str(ptch_str=='_')=' ';
        
        % Plot Image Patch just for Observer Orientation.
        axes(ha(4)), imagesc(GTFile.im), colormap('bone'), axis square, title(['Im Patch : ',ptch_str],'FontSize',20,'FontWeight','Bold')
        set(gca,'XTick',[],'YTick',[]), freezeColors,

        % Show Ground Truth Boundaries
        axes(ha(subplotsW+4)), imagesc(1-(GTFile.bD>0)), colormap('bone'), axis square, title('gT boundaries','FontSize',20,'FontWeight','Bold')
        set(gca,'XTick',[],'YTick',[]), freezeColors,
    end


    saveGoodImg(Hc,[dirImgSave,ptch,ftag],sizeGoodIm)
    close(Hc)

end




%% Plot # 1.5:  Plot Random stuff I am wanting to include in poster & presentation.
if (sav_plt1p5)
    
    % Plot Final Phase Result for a given method (here AA)
    figure, imagesc( visKurPhase_inHSV( GTFile.im, reshape( AA_K.metaCluster.phaseAtClk(:,end), AA_K.netParams.Ndims ) ) ), 
    colormap('hsv'), caxis([0 2*pi]), colorbar, freezeColors, cbfreeze
    set(gca,'XTick',[],'YTick',[],'FontSize',20,'FontWeight','Bold'), axis square
    title('AA Phase','FontSize',20,'FontWeight','Bold')
    
    
    
    % Plot Boundary ground truth & consensus and blurring processing we do on it.
    figure, colormap('bone')
    subplot(231), imagesc(max(GTFile.bD(:)) - GTFile.bD), title('GT Boundary','FontSize',20,'FontWeight','Bold'), 
    set(gca,'XTick',[],'YTick',[]), axis square
    for i = 1:5
        subplot(2,3,i+1), imagesc(1 - bDc_blurs(:,:,i) ), title(bDc_blurs_info{i},'FontSize',20,'FontWeight','Bold'), 
        set(gca,'XTick',[],'YTick',[]), axis square
    end
    
    
    
end





%% PLOT #2: Plot (a). Average Gradient at boundaries in ground truth w/ error bars and 
%                (b). Ratio of Avg Gradient Step at Boundaries vs within Segments.
if(sav_plt2)

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
         

%% Save KurMC matfile in new directory...

if(0)
    figure, imagesc(F.im), colormap('jet'), colorbar, set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','Bold'),title('Pixel Gradient','FontSize',20,'FontWeight','Bold')
    figure, imagesc(mod(F.sk,pi)), colormap('jet'), colorbar, set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','Bold'), title('SK Phase Gradient','FontSize',20,'FontWeight','Bold')
    figure, imagesc(mod(F.ng,pi)), colormap('jet'), colorbar, set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','Bold'), title('NG Phase Gradient','FontSize',20,'FontWeight','Bold')

    figure, imagesc(mod(F.aa,pi)), colormap('jet'), colorbar, set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','Bold'), title('AA Phase Gradient','FontSize',20,'FontWeight','Bold')
    figure, imagesc(mod(F.gl,pi)), colormap('jet'), colorbar, set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','Bold'), title('GL Phase Gradient','FontSize',20,'FontWeight','Bold')
    figure, imagesc(mod(F.iso,pi)), colormap('jet'), colorbar, set(gca,'XTick',[],'YTick',[],'FontSize',16,'FontWeight','Bold'), title('Isotropic Diffusion Phase Gradient','FontSize',20,'FontWeight','Bold')
end
