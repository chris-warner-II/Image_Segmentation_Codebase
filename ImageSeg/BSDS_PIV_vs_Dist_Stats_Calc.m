function BSDS_PIV_vs_Dist_Stats_Calc(im_st, im_fin, ds_fctr, ptch_sz, sigpix, sigdist, r_range, plt_flg)

% Words of Wizdumb
%
% syntax: BSDS_PIV_vs_Dist_Stats_Calc(im_st, im_fin, ds_fctr, ptch_sz, sigpix, sigdist, r_range, plt_flg);
%
% INPUTS:
% ------
%
% im_st
%
% im_fin
%
% ds_fctr
%
% ptch_sz
%
% sigpix
%
% sigdist
%
% r_range = [rmin, rmax] or just rmax (then rmin defaults to 1)
%
% plt_flg


%% Make a directories to put output images and data in
    % which onCluster.m
    [dirPre,sizeGoodIm] = onCluster;



%% convert r_range into rmin & rmax, whether it is 1 number of 2 numbers.
if( numel(r_range)==2 )
    rmin = r_range(1);
    rmax = r_range(2);
elseif( numel(r_range)==1 )
    rmin = 0;
    rmax = r_range;
else
    disp('r_range variable does not make sense.')
    keyboard
end


%% Read through image patch matlab files saved in Projects/images/BSDS_patch/101x101_ds1 directory
maskFlg = 1;
topo = 2;

bins = [0:0.1:1];                 % bin centers for histogram (Hr)
cmap = jet(numel(bins));          % colormap for line histogram plots.

im_dir = [dirPre,'images/BSDS_patch/',num2str(ptch_sz),'x',num2str(ptch_sz),'_ds',num2str(ds_fctr),'/'];
iids = dir([im_dir,'*.mat']);


sP = num2str(sigpix);  % string of sigpix, for output file name
sP(sP=='.')='p';
%
sD = num2str(sigdist); % string of sigdist, for output file name
sD(sD=='.')='p';
%
rMin = num2str(rmin);  % string of rmin, for output file name
rMin(rMin=='.')='p';
%
rMax = num2str(rmax);  % string of rmax, for output file name
rMax(rMax=='.')='p';
%
dir_out = [dirPre,'output/ImgSeg/simpleExamples/BSDS/PIV_v_Dist/',...
    num2str(ptch_sz),'x',num2str(ptch_sz),'_ds',num2str(ds_fctr),'/'...
    'sP',sP,'_sD',sD,'_r',rMin,'-',rMax,'/'];
%
if ~exist(dir_out,'dir')
    mkdir(dir_out);
end






%% (1). Loop through BSDS input images, compute Weight as a function of 
%        distance between pixel pairs, & save that in a mat file.
for i = im_st:im_fin

    i
    iid = iids(i).name
    
    % Check if output save file is already there. If so, dont repeat calculation.
    if exist([dir_out,iid],'file')
        disp([num2str(i),' - file already exists.'])
        continue
    end
    
    load([im_dir,iid]);
    
    disp('Calc Weights')
    
%     im,
%     sigpix,
%     sigdist,
%     rmin,
%     rmax,
%     maskFlg,
%     topo
%     
%     which calc_weightsB.m
    tic
    [W, Wdist, Mask] = calc_weightsB(im,sigpix,sigdist,[rmin,rmax],maskFlg,topo);
    
    disp('Calc PIV Distance Dependencies')
    [B] = PIV_vs_MaskDistB(W,Wdist,Mask);
    tic


    % Compute Weight statistics (mean, std, histogram)
    r = [Mask.distance];              % distances in image plane between 0 & rmax
    %
    r_ind = find(r>=rmin);            % indecies to r between rmin & rmax
    r = r(r_ind);                     % now r goes between rmin & rmax
    Mask = Mask(r_ind);               % now Mask goes between rmin & rmax
    %
    edges = zeros(1,numel(r));        % # of edges possible at each distance 
    Pr = zeros(1,numel(r));           % mean edge weight at each distance
    Sr = zeros(1,numel(r));           % std of edge weight at each distance
    Hr = zeros(numel(bins),numel(r)); % histogram of how many weight values 
                                      % fell into each bin at each distance
                                      
    bins_leg = num2cell(bins); % for plot legend
    for k = 1:numel(bins_leg)
        bins_leg{k} = num2str(bins_leg{k});
    end
    
    
    
    for j = 1:numel(Mask)
        edges(j) = numel(find(Mask(j).mask));
        Pr(j) = mean(W(find(Mask(j).mask)));
        Sr(j) = std(W(find(Mask(j).mask)));
        Hr(:,j) = hist(W(find(Mask(j).mask)),bins);
    end

    
    
    
    
    % Make some plots of Weight statistics
    if(plt_flg)
        ints = find( mod(r,1)==0 );
        
        % (1). Errorbar plot of weights (W) for pixel pairs separated by distance (r).
        figure('units','normalized','outerposition',[0 0 1 1])
        subplot(221), hold on
        errorbar(1:numel(r),Pr,Sr)
        plot(1:numel(r),Pr,'LineWidth',2,'Color','r')
        plot(1:numel(r),ones(size(r)),'k--')
        set(gca,'FontSize',16,'FontWeight','Bold','Xtick',ints,'XtickLabel',r(ints))
        ylabel('< W >','FontSize',18,'FontWeight','Bold')
        xlabel('r','FontSize',18,'FontWeight','Bold')
        h = legend({'\sigma','\mu'},'Location','EastOutside'); %'Orientation','Horizontal');
        set(get(h,'title'),'string','W','FontSize',18,'FontWeight','Bold');
        title('Expectation of W as a function of r','FontSize',18,'FontWeight','Bold')
        axis tight
        freezeColors
        %
        subplot(122)
        imagesc(im), colormap('bone') %,colorbar
        set(gca,'FontSize',16,'FontWeight','Bold','Xtick',[1:10:ptch_sz],'Ytick',[1:10:ptch_sz],'YtickLabel',[ptch_sz:-10:1])
        title(['Image ',iid],'FontSize',20,'FontWeight','Bold')
        axis square
        freezeColors, %cbfreeze
   
    
        % (2). Display histogram of weight distribution at each distance (not just mean & std)
        subplot(223)
        bar(Hr' ./ repmat(edges,numel(bins),1)', 'stacked' ), colormap(jet)
        set(gca,'FontSize',16,'FontWeight','Bold','Xtick',ints,'XtickLabel',r(ints),'Ytick',[0 1])
        ylabel({'Distribution of available edges','at distance (r) with W value'},'FontSize',18,'FontWeight','Bold')
        xlabel('r','FontSize',18,'FontWeight','Bold')
        h = legend(bins_leg,'Location','EastOutside'); %'Orientation','Horizontal');
        set(get(h,'title'),'string','W','FontSize',18,'FontWeight','Bold');
        title(['Histogram of W as a function of r'],'FontSize',18,'FontWeight','Bold')
        axis([0 numel(r) 0 1])
        freezeColors
        
        pause
        close
        
    end
    
    % Save output mat file with important statistics for single image.
    save([dir_out,iid],'r','Pr','Sr','Hr','edges','bins','iid')
    
    clear Mask B W Wdist bDfull bD gT gTfull imFull 

    
end



