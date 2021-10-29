function [] = compare_allParams_boundaryGradientMetric_Eig(method,ptch)

% This function will make a single diagnostic plot for a method and image
% patch.  It will show 



    % save flags.
    sav_mov0 = 0;
    sav_plt1 = 1;

    which_bD = 7;


    [dirPre,sizeGoodIm] = onCluster;

    dirImgSave = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/imgs/spectral/allParams/',method,'/'];
    if ~exist(dirImgSave,'dir')
        mkdir(dirImgSave);
    end

    ftag = '_visGradients';

    if exist([dirImgSave,ptch,ftag,'_bD',num2str(which_bD),'.jpg'],'file')
        disp('This jpg image already exists: Not gonna replace it.')
        disp([dirImgSave,ptch,ftag,'_bD',num2str(which_bD),'.jpg'])
        return
    end


    %% Load KurMC & Evecs mat files with Simulation Results and Phase Distributions.

    eig_data_dir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/spectral/',method,'/'];
    img_data_dir = [dirPre,'images/BSDS_patch/101x101_ds1/'];

    if ~isempty(strmatch(method,'IsoDiff'))
        fname_fill = '';
    else
        fname_fill = '_sDInf_sP0p2';
    end

    % disp('rM1')
    rM1 = load([eig_data_dir,'Evecs_',ptch,'_rM1',fname_fill,'.mat']);
    %rM1_KSsml = rmfield(rM1_KSsml,'netParams'); % doing this so I dont overload memory with Mod_N&G ( can get rid of this)

    % disp('rM3')
    rM3 = load([eig_data_dir,'Evecs_',ptch,'_rM3',fname_fill,'.mat']);
    %rM3_KSsml = rmfield(rM3_KSsml,'netParams'); % doing this so I dont overload memory with Mod_N&G ( can get rid of this

    % disp('rM5')
    rM5 = load([eig_data_dir,'Evecs_',ptch,'_rM5',fname_fill,'.mat']);
    %rM5_KSsml = rmfield(rM5_KSsml,'netParams'); % doing this so I dont overload memory with Mod_N&G ( can get rid of this)

    % disp('rM10')
    rM10 = load([eig_data_dir,'Evecs_',ptch,'_rM10',fname_fill,'.mat']);
    %rM10_KSsml = rmfield(rM10_KSsml,'netParams'); % doing this so I dont overload memory with Mod_N&G ( can get rid of this)


    GTFile = load([img_data_dir,ptch,'.mat']);


    %% Plot #0:  Make a movie of phase relaxation
    if(sav_mov0)

        makeMovie_phaseEvolution([dirImgSave,ptch,'_SK_v_Iso_phase'],GTFile,SK_K,IsoK)
        % MAYBE STILL DO THIS.  USEFUL TO WATCH ALL PARAM VALUES OVER TIME...?
    end


    %% PLOT #1: Visualize Gradients (phase,contrast or otherwise)
    % for different competitor methods and plot avg gradient/pixel
    % on boundaries vs. off for different boundary blurs.
    if(sav_plt1)

    %     numTHs=20;
    % 
    %
    %     % Compute Mean & Std of Gradient Fields On & Off Ground Truth Boundaries
    %     X = GTFile.im;
    %     %[F.im,M.im,S.im] = compute_BoundaryGradientMetric(X,bDc_blurs,0,'Img Pix');
    %     Ph.im = X;
    %     %[precision.im,recall.im,f_measure.im] = compute_BoundaryGradientPrecRec(X,bDc_blurs,numTHs,'Img Pix');
    %     %
    %     X = visKurPhase_inHSV(IsoK.netParams.im, reshape(IsoK.metaCluster.phaseAtClk(:,end),IsoK.netParams.Ndims));
    %     %[F.iso,M.iso,S.iso] = compute_BoundaryGradientMetric(X,bDc_blurs,1,'Iso Diff');
    %     Ph.iso = X;
    %     %[precision.iso,recall.iso,f_measure.iso] = compute_BoundaryGradientPrecRec(X,bDc_blurs,numTHs,'Iso Diff');
    %     %
    %     X = visKurPhase_inHSV(NG_K.netParams.im, reshape(NG_K.metaCluster.phaseAtClk(:,end),NG_K.netParams.Ndims));
    %     %[F.ng,M.ng,S.ng] = compute_BoundaryGradientMetric(X,bDc_blurs,1,'Modularity NG');
    %     Ph.ng =X;
    %     %[precision.ng,recall.ng,f_measure.ng] = compute_BoundaryGradientPrecRec(X,bDc_blurs,numTHs,'Modularity NG');
    %     %
    %     X = visKurPhase_inHSV(SK_K.netParams.im, reshape(SK_K.metaCluster.phaseAtClk(:,end),SK_K.netParams.Ndims));
    %     %[F.sk,M.sk,S.sk] = compute_BoundaryGradientMetric(X,bDc_blurs,1,'Modularity SK');
    %     Ph.sk = X;
    %     %[precision.sk,recall.sk,f_measure.sk] = compute_BoundaryGradientPrecRec(X,bDc_blurs,numTHs,'Modularity SK');
    %     %
    %     X = visKurPhase_inHSV(GL_K.netParams.im, reshape(GL_K.metaCluster.phaseAtClk(:,end),GL_K.netParams.Ndims));
    %     %[F.gl,M.gl,S.gl] = compute_BoundaryGradientMetric(X,bDc_blurs,1,'GL');
    %     Ph.gl = X;
    %     %[precision.gl,recall.gl,f_measure.gl] = compute_BoundaryGradientPrecRec(X,bDc_blurs,numTHs,'GL');
    %     %
    %     X = visKurPhase_inHSV(AA_K.netParams.im, reshape(AA_K.metaCluster.phaseAtClk(:,end),AA_K.netParams.Ndims));
    %     %[F.aa,M.aa,S.aa] = compute_BoundaryGradientMetric(X,bDc_blurs,1,'AA');
    %     Ph.aa = X;
    %     %[precision.aa,recall.aa,f_measure.aa] = compute_BoundaryGradientPrecRec(X,bDc_blurs,numTHs,'AA');



        % Figure with Subplots !


    %     % for plotting phases        
    %     F = Ph;                         
    %     cmap = 'hsv';      
    %     ftag = '_visPhases';  


        % for plotting gradients
        cmap = 'jet';



        Hc=figure;
        ha = tight_subplot(6, 5, [0.05 0.02], 0.05, 0.02);

        % rM1 & evec1
        axes(ha(1)), imagesc_gradients(rM1.MC.ev1,'magenta',cmap,which_bD)
        xlabel(['\color{black} \{rM=1 ; \color{magenta} ev1 \}'],'FontSize',16,'FontWeight','Bold')
        %
        % rM3 & evec1
        axes(ha(2)), imagesc_gradients(rM3.MC.ev1,'magenta',cmap,which_bD)
        xlabel(['\color{black} \{rM=3\}'],'FontSize',16,'FontWeight','Bold')
        %
        % rM5 & evec1
        axes(ha(3)), imagesc_gradients(rM5.MC.ev1,'magenta',cmap,which_bD)
        xlabel(['\color{black} \{rM=5\}'],'FontSize',16,'FontWeight','Bold')
        %
        % rM10 & evec1
        axes(ha(4)), imagesc_gradients(rM10.MC.ev1,'magenta',cmap,which_bD)
        xlabel(['\color{black} \{rM=10\}'],'FontSize',16,'FontWeight','Bold')
        %
        % Image
        axes(ha(5)), imagesc_gradients(GTFile.MC,'black',cmap,which_bD)
        xlabel(['\{Img\}'],'FontSize',16,'FontWeight','Bold')
        
        
        
        
        % mean & std of phase gradient distributions on & off boundaries.
        axes(ha(5+1)), plot_B_notB_dists(rM1.MC.ev1,which_bD)
        %
        axes(ha(5+2)), plot_B_notB_dists(rM3.MC.ev1,which_bD)
        %
        axes(ha(5+3)), plot_B_notB_dists(rM5.MC.ev1,which_bD)
        %
        axes(ha(5+4)), plot_B_notB_dists(rM10.MC.ev1,which_bD)
        %
        axes(ha(5+5)), plot_B_notB_dists(GTFile.MC,which_bD)
        
        
        
        
        


        % rM1 & evec2
        axes(ha(2*5+1)), imagesc_gradients(rM1.MC.ev2,'green',cmap,which_bD)
        xlabel(['\color{green} \{ev2\}'],'FontSize',16,'FontWeight','Bold')
        %
        % rM3 & evec2
        axes(ha(2*5+2)), imagesc_gradients(rM3.MC.ev2,'green',cmap,which_bD)
        %
        % rM5 & evec2
        axes(ha(2*5+3)), imagesc_gradients(rM5.MC.ev2,'green',cmap,which_bD)
        %
        % rM10 & evec2
        axes(ha(2*5+4)), imagesc_gradients(rM10.MC.ev2,'green',cmap,which_bD)
        %
        % Ground Truth Boundaries
        axes(ha(2*5+5)), imagesc(max(GTFile.bD(:))-GTFile.bD), colormap(bone), axis square, 
        set(gca,'XTick',[],'YTick',[],'XColor','k','YColor','k','LineWidth',3), freezeColors, % colorbar, cbfreeze
        xlabel(['\{GT BD ',num2str(which_bD),'\}'],'FontSize',16,'FontWeight','Bold')

        
        
        
        
        
        % mean & std of phase gradient distributions on & off boundaries.
        axes(ha(3*5+1)), plot_B_notB_dists(rM1.MC.ev2,which_bD)
        %
        axes(ha(3*5+2)), plot_B_notB_dists(rM3.MC.ev2,which_bD)
        %
        axes(ha(3*5+3)), plot_B_notB_dists(rM5.MC.ev2,which_bD)
        %
        axes(ha(3*5+4)), plot_B_notB_dists(rM10.MC.ev2,which_bD)
        %
        axes(ha(3*5+5)), axis off
        
        
        
        
        
        

        % rM1 & evec3
        axes(ha(4*5+1)), imagesc_gradients(rM1.MC.ev3,'cyan',cmap,which_bD)
        xlabel(['\color{cyan} \{ev3\}'],'FontSize',16,'FontWeight','Bold')

        % rM3 & evec3
        axes(ha(4*5+2)), imagesc_gradients(rM3.MC.ev3,'cyan',cmap,which_bD)
        
        % rM5 & evec3
        axes(ha(4*5+3)), imagesc_gradients(rM5.MC.ev3,'cyan',cmap,which_bD)
        
        % rM10 & evec3
        axes(ha(4*5+4)), imagesc_gradients(rM10.MC.ev3,'cyan',cmap,which_bD)
        % Ranking of d'
        axes(ha(4*5+5)), hold on
        plot([1,4],[GTFile.MC.D(which_bD),GTFile.MC.D(which_bD)],'k--')
        scatter(1,rM1.MC.ev1.D(which_bD),50,'m','filled')
        scatter(1,rM1.MC.ev2.D(which_bD),50,'g','filled')
        scatter(1,rM1.MC.ev3.D(which_bD),50,'c','filled')
        scatter(2,rM3.MC.ev1.D(which_bD),50,'m','filled')
        scatter(2,rM3.MC.ev2.D(which_bD),50,'g','filled')
        scatter(2,rM3.MC.ev3.D(which_bD),50,'c','filled')
        scatter(3,rM5.MC.ev1.D(which_bD),50,'m','filled')
        scatter(3,rM5.MC.ev2.D(which_bD),50,'g','filled')
        scatter(3,rM5.MC.ev3.D(which_bD),50,'c','filled')
        scatter(4,rM10.MC.ev1.D(which_bD),50,'m','filled')
        scatter(4,rM10.MC.ev2.D(which_bD),50,'g','filled')
        scatter(4,rM10.MC.ev3.D(which_bD),50,'c','filled')
        %
        axis square, 
        set(gca,'XTick',[1:4],'XTickLabel',{'1','3','5','10'},'FontSize',16,'FontWeight','Bold') 
        ylabel('d''','Fontsize',16,'FontWeight','Bold')
        xlabel('rM','Fontsize',16,'FontWeight','Bold')
        
        
        
        
        
        
        % mean & std of phase gradient distributions on & off boundaries.
        axes(ha(5*5+1)), plot_B_notB_dists(rM1.MC.ev3,which_bD)
        %
        axes(ha(5*5+2)), plot_B_notB_dists(rM3.MC.ev3,which_bD)
        %
        axes(ha(5*5+3)), plot_B_notB_dists(rM5.MC.ev3,which_bD)
        %
        axes(ha(5*5+4)), plot_B_notB_dists(rM10.MC.ev3,which_bD)
        %
        axes(ha(5*5+5)), axis off
        
        
        
        
        
        

        method_tag = method;
        method_tag(method_tag=='_')=' ';
        ptch_tag = ptch;
        ptch_tag(ptch_tag=='_')=' ';
        annotation('textbox', [0 0.9 1 0.1],'String', [method_tag,' : ',ptch_tag], 'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',18,'FontWeight','Bold')




        saveGoodImg(Hc,[dirImgSave,ptch,ftag,'_bD',num2str(which_bD)],sizeGoodIm)
        close(Hc)

    end


    disp('finished')

end


%% Function to make single pane in Imagesc Phase Gradients (F)
function imagesc_gradients(RmEv,boxColor,cmap,which_bD)
    
    imagesc(RmEv.F), colormap(cmap), axis square, 
    set(gca,'XTick',[],'YTick',[],'XColor',['',boxColor,''],'YColor',['',boxColor,''],'LineWidth',3), freezeColors, % colorbar, cbfreeze
    text(size(RmEv.F,1),size(RmEv.F,2),['\color{white}d''= ',num2str(RmEv.D(which_bD),2)],'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',16,'FontWeight','Bold')
    ylabel(['( \color{blue}{',num2str( min(RmEv.F(:)) ,2),'} \color{black}{,} \color{red}{',num2str( max(RmEv.F(:)) ,2),'}\color{black}{)}'],'FontSize',16,'FontWeight','Bold')
    
end


%% Function to make single pane in Mean & Std Errorbar Plots
function plot_B_notB_dists(RmEv,which_bD)

    hold on,
    herrorbar( RmEv.M(which_bD,1), 0.95, RmEv.M(which_bD,2) ,'ro') % on boundary
    herrorbar( RmEv.S(which_bD,1), 1.05, RmEv.S(which_bD,2) ,'bo') % off boundary
    del_mu = RmEv.M(which_bD,1) - RmEv.S(which_bD,1);
    text(double(RmEv.M(which_bD,1)), 0.93,{['\color{red}\sigma=',num2str(RmEv.M(which_bD,2)./del_mu,2)]},'VerticalAlignment','top','HorizontalAlignment','center')
    text( double((RmEv.S(which_bD,1)+RmEv.M(which_bD,1) )./2) , 1.00,{['\color{black}\Delta\mu=',num2str(del_mu,2)]},'VerticalAlignment','middle','HorizontalAlignment','center')
    text(double(RmEv.S(which_bD,1)), 1.07,{['\color{blue}%\sigma=',num2str(RmEv.S(which_bD,2)./del_mu,2)]},'VerticalAlignment','bottom','HorizontalAlignment','center')
    text(double(RmEv.S(which_bD,1)-RmEv.S(which_bD,2)), 1.02,{['\color{blue}\\B:']},'VerticalAlignment','bottom','HorizontalAlignment','right')
    text(double(RmEv.S(which_bD,1)-RmEv.S(which_bD,2)), 0.99,{['\color{red}B:']},'VerticalAlignment','top','HorizontalAlignment','right')
    ylim([0.8 1.2])
    axis off

end
