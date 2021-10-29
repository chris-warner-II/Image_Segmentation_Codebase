% This script is largely hardcoded to produce a figure for the Image Seg
% paper. It loads in .mat data files produced in compute_null_modelB which
% is called from Loop_ImgSegMethodsD(6,5,0.2,{[11,11]},'BSDS_patch',[2,2],0);
%
% Now, the goal is to plot the results for the 3 different modularity null
% models (NG, 1D and 2D TMs) on one figure and make it nice for the paper.



%% Load in data
consis_dir = './output/consistency_plots_mat_files/';

M       = load([consis_dir,'Consistency_AvsNM_Mod_NG.mat']);              % N&G
TM1   = load([consis_dir,'Consistency_AvsNM_Mod_topo_1D.mat']);     % TM 1D
TM2   = load([consis_dir,'Consistency_AvsNM_Mod_topo_2D.mat']);     % TM 2D








%% Plot figure.
    H=figure; 
    
        
    minn = min( [ full(M.W(:)); full(M.NM(:)); full(TM1.NM(:)); full(TM2.NM(:)) ]  );
    maxx = max( [ full(M.W(:)); full(M.NM(:)); full(TM1.NM(:)); full(TM2.NM(:)) ]  );

    
    % ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
    ha = tight_subplot(3, 5, [0.02, 0], [0.04, 0.04], [0.02,0.02]);
    hb = tight_subplot(3, 1, [0.10, 0.08], [0.10, 0.05], [0.04,0.04]);
    
	subplot(ha(1)), imagesc( M.imagin ), colormap(bone), freezeColors 
    title(['Image'],'FontSize',20,'FontWeight','Bold'), axis square
    set(gca,'XTick',[],'YTick',[1, size(M.imagin,1) ],'FontSize',16,'FontWeight','Bold'), 
    
    subplot(ha(2)), imagesc(M.W, [minn, maxx] ), colormap(jet) %colorbar, 
    title(['Adjacency'],'FontSize',20,'FontWeight','Bold'), axis square
    set(gca,'XTick',[],'YTick',[1, size(M.W,1)],'FontSize',16,'FontWeight','Bold'), 
    %ylabel(['$\sum_{ij} A_{ij} = $',num2str(sum(M.W(:)),5)],'FontSize',14,'FontWeight','Bold','Interpreter','LaTex')

    subplot(ha(3)), imagesc(M.NM, [minn, maxx]), colormap(jet) %colorbar, 
    title(['\color{blue}Mod N&G NM'],'FontSize',20,'FontWeight','Bold'), axis square
    set(gca,'XTick',[],'YTick',[1, size(M.NM,1)],'FontSize',16,'FontWeight','Bold'), 
   %ylabel(['$\sum_{ij} N_{ij} = $',num2str(sum(M.NM(:)),5)],'FontSize',14,'FontWeight','Bold','Interpreter','LaTex')
    
    subplot(ha(4)), imagesc(TM1.NM, [minn, maxx]), colormap(jet) %colorbar, 
    title(['\color{green}1D tMod NM'],'FontSize',20,'FontWeight','Bold'), axis square
    set(gca,'XTick',[],'YTick',[1, size(TM1.NM,1)],'FontSize',16,'FontWeight','Bold'), 
    %ylabel(['$\sum_{ij} N_{ij} = $',num2str(sum(TM1.NM(:)),5)],'FontSize',14,'FontWeight','Bold','Interpreter','LaTex')
    
    subplot(ha(5)), imagesc(TM2.NM, [minn, maxx]), colormap(jet), colorbar, 
    title(['\color{red}2D tMod NM'],'FontSize',20,'FontWeight','Bold'), axis square
    set(gca,'XTick',[],'YTick',[1, size(TM2.NM,1)],'FontSize',16,'FontWeight','Bold'), 
    %ylabel(['$\sum_{ij} N_{ij} = $',num2str(sum(TM2.NM(:)),5)],'FontSize',14,'FontWeight','Bold','Interpreter','LaTex')
    
    
    
    %set(ha(1), 'Visible', 'off')
    set(hb(1), 'Visible', 'off')
    for i = 6:15
        set(ha(i), 'Visible', 'off')
    end

    x = mean(M.W); % sort nodes based on degree and plot sorted nodes (easier to see)
    yM = mean(M.NM);
    y1 = mean(TM1.NM);
    y2 = mean(TM2.NM);
    YM = yM; % not sorted.
    Y1 = y1; % not sorted.
    Y2 = y2; % not sorted.
    X = x;      % not sorted.
    [X,I] = sort(x);  % or sorted.
    YM = yM(I);
    Y1 = y1(I);
    Y2 = y2(I);
    %
    subplot(hb(2)),hold on,  
    
    plot(YM,'b','LineWidth',3),
    plot(Y1,'g','LineWidth',3),
    plot(Y2,'r','LineWidth',3),
    plot(X,'k--','LineWidth',2), 
    xlabel(['Node #'],'FontSize',18,'FontWeight','Bold') % 'Sorted Node #'
    ylabel(['node degree'],'FontSize',18,'FontWeight','Bold')
    axis('tight')
    grid('on')
    %legend('Adjacency', 'N&G Mod', '1D top Mod', '2D top Mod','Location','SouthEast')
    set(gca,'FontSize',16,'FontWeight','Bold')

    subplot(hb(3)),hold on, 
    plot([M.MaskFull.distance],M.distNM,'b','LineWidth',3), 
    plot([TM1.MaskFull.distance],TM1.distNM,'g','LineWidth',3), 
    plot([TM2.MaskFull.distance],TM2.distNM,'r','LineWidth',3), 
    plot([M.MaskFull.distance],M.distW,'k--','LineWidth',2),
    xlabel(['Distance (in pixels) between nodes'],'FontSize',18,'FontWeight','Bold')
    ylabel(['avg. edge-weight'],'FontSize',18,'FontWeight','Bold')
    axis('tight')
    grid('on')
    set(gca,'FontSize',16,'FontWeight','Bold')

%     tagB = tag;
%     tagB(tagB==' ')='_';
  
    saveGoodImg(H,['./Consistency_Adj_vs_3ModNMs.eps'],[0 0 1 1])