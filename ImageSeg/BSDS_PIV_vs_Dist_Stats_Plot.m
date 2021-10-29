function BSDS_PIV_vs_Dist_Stats_Plot(ds_fctr, ptch_sz, sigpix, sigdist, r_range_range, plt_flg)

% Words of Wizdumb
%
% syntax: BSDS_PIV_vs_Dist_Stats_Plot(ds_fctr, ptch_sz, sigpix, sigdist, r_range_range, plt_flg);
%
% INPUTS:
% ------
%
% ds_fctr
%
% ptch_sz
%
% sigpix
%
% sigdist
%
% r_range_range = [rmin, rmax]
%
% plt_flg



% Make a directories to put output images and data in
[dirPre,sizeGoodIm] = onCluster;

bins = [0:0.1:1];                 % bin centers for histogram (Hr & Hr_tot)
cmap = jet(numel(bins));          % colormap for line histogram plots.






% Construct outer image directory name
sP = num2str(sigpix);  % string of sigpix, for output file name
sP(sP=='.')='p';
%
sD = num2str(sigdist); % string of sigdist, for output file name
sD(sD=='.')='p';
% 



% Loop through different r_ranges and build up r_tot vector which will
% include all r values across all ranges in r_range_range
iid = [];
r_tot = []; % string together the r vector from each distance range.
%
for R = 1:size(r_range_range,1)
    
    % convert r_range_range into rmin & rmax iteratively.
    rmin = r_range_range(R,1);
    rmax = r_range_range(R,2);
    %
    rMin = num2str(rmin);  % string of rmin, for output file name
    rMin(rMin=='.')='p';
    %
    rMax = num2str(rmax);  % string of rmax, for output file name
    rMax(rMax=='.')='p';
    
    im_dir = [dirPre,'output/ImgSeg/simpleExamples/BSDS/PIV_v_Dist/',...
        num2str(ptch_sz),'x',num2str(ptch_sz),'_ds',num2str(ds_fctr),'/'...
        'sP',sP,'_sD',sD,'_r',rMin,'-',rMax,'/'];
    
    iids = dir([im_dir,'*.mat']);

    if isempty(iid)
        iid = iids(1).name;
    end
    
    load([im_dir,iid]);
    
    r_tot = [r_tot, r];
    
end

% if isempty( find(r_tot==1) ) % add distance r=1 so we can plot fall off relative to it in 3rd plot.
%     r_tot = [1, r_tot];
% end


% Prealloate memory for summary arrays that will include Adjacency Weight
% value as a function of distance r across all image patches.
Pr_tot = zeros( numel(iids) , numel(r_tot) );
Sr_tot = zeros( numel(iids) , numel(r_tot) );
edges_tot = zeros( 1        , numel(r_tot) );
Hr_tot = zeros( numel(bins) , numel(r_tot) );


% Loop through different r ranges & 
for R = 1:size(r_range_range,1)

    im_st = 1;
    im_fin = numel(iids);
    
    num_ims = 1 + im_fin - im_st;
    
    % convert r_range_range into rmin & rmax iteratively.
    rmin = r_range_range(R,1);
    rmax = r_range_range(R,2);
    %
    rMin = num2str(rmin);  % string of rmin, for output file name
    rMin(rMin=='.')='p';
    %
    rMax = num2str(rmax);  % string of rmax, for output file name
    rMax(rMax=='.')='p';
    
    im_dir = [dirPre,'output/ImgSeg/simpleExamples/BSDS/PIV_v_Dist/',...
        num2str(ptch_sz),'x',num2str(ptch_sz),'_ds',num2str(ds_fctr),'/'...
        'sP',sP,'_sD',sD,'_r',rMin,'-',rMax,'/'];
    
    iids = dir([im_dir,'*.mat']);

    for i = im_st:im_fin
        
        [R,i]

        iid = iids(i).name;
        xx = load([im_dir,iid]);
        load([im_dir,iid])
        
        
        rmin_ind = find(r_tot==r(1));
        rmax_ind = find(r_tot==r(end));

        % Accumulate Pr, Sr, Hr & edges across all image patches in reasonable ways.
        Pr_tot(i,rmin_ind:rmax_ind) = Pr;
        Sr_tot(i,rmin_ind:rmax_ind) = Sr;
        Hr_tot(:,rmin_ind:rmax_ind) = Hr_tot(:,rmin_ind:rmax_ind) + Hr;
        if(i==1)
            edges_tot(i,rmin_ind:rmax_ind) = xx.edges;
        end


%         % make plots I would have made above for single image patches.
%         if(plt_flg)
% 
%             ints = find( mod(r,1)==0 );
% 
%             bins_leg = num2cell(bins); % for plot legend
%             for k = 1:numel(bins_leg)
%                 bins_leg{k} = num2str(bins_leg{k});
%             end
% 
%             % (1). Errorbar (mean,std) plot of weights (W) for pixel pairs separated by distance (r).
%             figure('units','normalized','outerposition',[0 0 1 1])
%             subplot(311), hold on
%             errorbar(1:numel(r),Pr,Sr)
%             plot(1:numel(r),Pr,'LineWidth',2,'Color','r')
%             plot(1:numel(r),ones(size(r)),'k--')
%             set(gca,'FontSize',16,'FontWeight','Bold','Xtick',ints,'XtickLabel',r(ints))
%             ylabel('< W >','FontSize',18,'FontWeight','Bold')
%             xlabel('r','FontSize',18,'FontWeight','Bold')
%             h = legend({'\sigma','\mu'},'Location','EastOutside'); %'Orientation','Horizontal');
%             set(get(h,'title'),'string','W','FontSize',18,'FontWeight','Bold');
%             title('Expectation of W as a function of r','FontSize',18,'FontWeight','Bold')
%             axis tight  
% 
%             % (2). Histogram w stacked of weight distribution at each distance (not just mean & std)
%             subplot(312)
%             bar(Hr' ./ repmat(edges,numel(bins),1)', 'stacked' ), colormap(jet)
%             set(gca,'FontSize',16,'FontWeight','Bold','Xtick',ints,'XtickLabel',r(ints),'Ytick',[0 1])
%             ylabel({'Distribution of available edges','at distance (r) with W value'},'FontSize',18,'FontWeight','Bold')
%             xlabel('r','FontSize',18,'FontWeight','Bold')
%             h = legend(bins_leg,'Location','EastOutside'); %'Orientation','Horizontal');
%             set(get(h,'title'),'string','W','FontSize',18,'FontWeight','Bold');
%             title(['Histogram of W as a function of r'],'FontSize',18,'FontWeight','Bold')
%             axis([0 numel(r) 0 1]) 
% 
%             % (3). 
%             subplot(313), hold on
%             for i = 1:numel(bins)
%                 plot(Hr_tot(i,:) ./ (num_ims.*edges),'LineWidth',2,'Color',cmap(i,:))
%             end
%             set(gca,'FontSize',16,'FontWeight','Bold','Xtick',ints,'XtickLabel',r(ints),'Ytick',[0:0.2:1])
%             ylabel({'Distribution of available edges','at distance (r) with W value'},'FontSize',18,'FontWeight','Bold')
%             xlabel('r','FontSize',18,'FontWeight','Bold')
%             h = legend(bins_leg,'Location','EastOutside'); %'Orientation','Horizontal');
%             set(get(h,'title'),'string','W bins','FontSize',18,'FontWeight','Bold');
%             title(['Histogram of W as a function of r'],'FontSize',18,'FontWeight','Bold')
%             axis([0 numel(r) 0 1])
%             grid on
%             freezeColors
% 
%             keyboard
% 
%         end


    end % looping over images in a given r range.

end % looping r range values.




%% Make plots of activity (mean, std & histogram) across all image patches from im_st to im_fin.
if(1)
        
    ints = find( mod(r_tot,1)==0 );

    % (1). Mean of weights (W) for pixel pairs separated by distance (r).
    H = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(211), hold on
    plot(1:numel(r_tot),Pr_tot,'LineWidth',1)
    %plot(1:numel(r_tot),ones(size(r_tot)),'k--')
    set(gca,'FontSize',16,'FontWeight','Bold','Xtick',ints,'XtickLabel',r_tot(ints))
    ylabel('\mu W (within single img patch)','FontSize',18,'FontWeight','Bold')
    xlabel('r','FontSize',18,'FontWeight','Bold')
    title(['W as a function of r (across ',num2str(num_ims),' img patches)'],'FontSize',18,'FontWeight','Bold')
    axis tight
    grid on
    %
    shadedErrorBar(1:numel(r_tot), mean(Pr_tot,1), std(Pr_tot))
    plot(1:numel(r_tot), mean(Pr_tot,1),'k', 'LineWidth',2)
    %
    % (2). Std of weights (W) for pixel pairs separated by distance (r).
    subplot(212), hold on
    plot(1:numel(r_tot),Sr_tot,'LineWidth',1)
    set(gca,'FontSize',16,'FontWeight','Bold','Xtick',ints,'XtickLabel',r_tot(ints))
    ylabel('\sigma W (within single img patch)','FontSize',18,'FontWeight','Bold')
    xlabel('r','FontSize',18,'FontWeight','Bold')
    axis tight 
    grid on
    %
    shadedErrorBar(1:numel(r_tot), mean(Sr_tot,1), std(Sr_tot))
    plot(1:numel(r_tot), mean(Sr_tot,1),'k', 'LineWidth',2)
    %
    saveGoodImg(H,[im_dir,'Plot1.jpg'],[0 0 1 1])
    close(H)
    
    
    
    
    
    
    
    
    % (3). Histogram w stacked bars of weight distribution at each distance (r)
    bins_leg = num2cell(bins); % for plot legend
    for k = 1:numel(bins_leg)
        bins_leg{k} = num2str(bins_leg{k});
    end
    %
    H = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(211),
    bar(Hr_tot' ./ repmat(num_ims.*edges_tot,numel(bins),1)', 'stacked' ), colormap(jet)
    set(gca,'FontSize',16,'FontWeight','Bold','Xtick',ints,'XtickLabel',r_tot(ints),'Ytick',[0 1])
    ylabel({'Distribution of available edges','at distance (r) with W value'},'FontSize',18,'FontWeight','Bold')
    xlabel('r','FontSize',18,'FontWeight','Bold')
    h = legend(bins_leg,'Location','EastOutside'); %'Orientation','Horizontal');
    set(get(h,'title'),'string','W bins','FontSize',18,'FontWeight','Bold');
    title(['Histogram of W as a function of r'],'FontSize',18,'FontWeight','Bold')
    axis([0 numel(r_tot) 0 1])
    freezeColors
    %
    % (4). Histogram w lines of weight distribution at each distance (r)
    subplot(212), hold on
    for i = 1:numel(bins)
        plot(Hr_tot(i,:) ./ (num_ims.*edges_tot),'LineWidth',2,'Color',cmap(i,:))
    end
    set(gca,'FontSize',16,'FontWeight','Bold','Xtick',ints,'XtickLabel',r_tot(ints),'Ytick',[0:0.2:1])
    ylabel({'Distribution of available edges','at distance (r) with W value'},'FontSize',18,'FontWeight','Bold')
    xlabel('r','FontSize',18,'FontWeight','Bold')
    h = legend(bins_leg,'Location','EastOutside'); %'Orientation','Horizontal');
    set(get(h,'title'),'string','W bins','FontSize',18,'FontWeight','Bold');
    title(['Histogram of W as a function of r'],'FontSize',18,'FontWeight','Bold')
    axis([0 numel(r_tot) 0 1])
    grid on
    freezeColors
    saveGoodImg(H,[im_dir,'Plot2.jpg'],[0 0 1 1])
    close(H)


    
    
    
    
    
    
    
    
    
    % Plot falloff of weight relative to value at 1 pixel distance.
    ind = find(r_tot==1);
    if ~isempty(ind)
        Pr_tot_rel = Pr_tot ./ repmat( Pr_tot(:,ind), 1, numel(r_tot) );

        H = figure('units','normalized','outerposition',[0 0 1 1]); hold on,
        shadedErrorBar(1:numel(r_tot), mean(Pr_tot_rel), std(Pr_tot_rel))
        plot(1:numel(r_tot), mean(Pr_tot_rel),'k', 'LineWidth',2)
        errorbar( 1:numel(r_tot), mean(Pr_tot_rel), std(Pr_tot_rel), 'bx', 'LineWidth',2) % round(numel(r_tot)./2)
        %errorbar( r_tot(ind), mean(Pr_tot(:,ind)), std(Pr_tot(:,ind)), 'gx', 'LineWidth',2) % round(numel(r_tot)./2)
        legend({'\mu \sigma abs','\mu \sigma rel'},'Location','Best')
        set(gca,'FontSize',16,'FontWeight','Bold','Xtick',ints,'XtickLabel',r_tot(ints))
        ylabel({'Weight at r, relative to Weight at r=1.'},'FontSize',18,'FontWeight','Bold')
        xlabel('r','FontSize',18,'FontWeight','Bold')
        title(['W(r) relative to W(r=1) (\mu & \sigma across ',num2str(num_ims),' img patches)'],'FontSize',18,'FontWeight','Bold')
        axis([0 numel(r_tot) 0 1])
        grid on
        saveGoodImg(H,[im_dir,'Plot3.jpg'],[0 0 1 1])
        close(H)
    end
    
end





