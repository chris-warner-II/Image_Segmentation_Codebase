% This script will take in a pair of mat files (1 for Kur & 1 for Eig) with
% the outcome of each's segmentation.  It will then plot a single figure
% for a single image patch.  It will display different segmentations in the
% first column {im, ev1, ev2o, ev3o, kur, ...} and will display the
% different ground truths across the top row.  At the intersection of a
% segmentation & ground truth, it will display a matrix of the pairwise 
% hit-rate metric.


%% Set up Directories to mat files.
netMethod = 'GLnrm'; % ,,'Mod_SKHAdj','Mod_N&G'

[dirPre,sizeGoodIm] = onCluster;

dirBase = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/'];

eigDir = [dirBase,'spectral_fixed_B/',netMethod,'/'];

kurDir = [dirBase,'Kur_PIF_Fourier1_fixed_B/',netMethod,'/'];

imgOutDir = [dirBase,'../imgs/visHitsMetric/',netMethod,'/'];
%
if ~exist(imgOutDir,'dir')
    mkdir(imgOutDir)
end






%% Loop thru Eig directory to input mat file and corresponding mat file from Kur dir.
eigFiles = dir([eigDir,'*.mat']);
% kurFiles = dir([kurDir,'*.mat']);


for i = 1:numel(eigFiles)
    
    i
    eigFname = eigFiles(i).name;
    eig = load([eigDir,eigFname]);
    
    try
        st = 6; %'Evecs' is start of filenames
        nd = length(eigFname)-4;
        kur = load([kurDir,'KurMC',eigFname(st:nd),'_NF_60_0_kscale300_tscale1_runs1.mat']);
    catch
        disp('Looks like KurMC file doesnt exist.  Skipping...')
        ['KurMC',eigFname(st:nd),'_NF_60_0_kscale300_tscale1_runs1.mat']
        continue
    end
    
    
    
    
    
    % In one figure, plot segmentations, ground truths, & hits results.
    H=figure;
    
    numCols = numel(kur.netParams.gT)+1;
    numRows = 6; % {im, ev1,ev2o,ev3o, kur} + 1
    
    % ground truth
    for j = 1:numel(kur.netParams.gT)
        subplot(numRows,numCols,j+1), imagesc(kur.netParams.gT{j})
        colormap('jet'), freezeColors,
        title(['gT#',num2str(j)],'FontSize',16,'FontWeight','Bold')
        set(gca,'Xtick',[],'YTick',[])
    end
    text(1.1*eig.netParams.Ndims(1),0.5*eig.netParams.Ndims(2),{'mean','across','gTs&pairs'},'HorizontalAlignment','Left','FontSize',16,'FontWeight','Bold')
    
    
    
    % im
    subplot(numRows,numCols,numCols+1), imagesc(kur.netParams.im)
    colormap('jet'), freezeColors,
    ylabel(['im'],'FontSize',16,'FontWeight','Bold')
    set(gca,'Xtick',[],'YTick',[])
    %
    sum_hits = [];
    for j = 1:numel(kur.netParams.gT)
        subplot(numRows,numCols,numCols+j+1),
        imagesc(eig.hits.im{j}), caxis([0.5,1])
        colormap('bone'), freezeColors,
        set(gca,'Xtick',[],'YTick',[])
        %
        sum_hits = [sum_hits; eig.hits.im{j}(find(eig.hits.im{j}))]; % find average pairwise hit-rate (weight it by cluster size?)
    end
    text(size(eig.hits.im{j},1)+0.7,size(eig.hits.im{j},2),num2str(mean(sum_hits),2),'HorizontalAlignment','Left','FontSize',16,'FontWeight','Bold')
    
    % kur
    subplot(numRows,numCols,2*numCols+1), imagesc(reshape(kur.metaCluster.phaseAtClk(:,end),kur.netParams.Ndims))
    colormap('hsv'), freezeColors,
    ylabel(['kur'],'FontSize',16,'FontWeight','Bold')
    set(gca,'Xtick',[],'YTick',[])
    %
    sum_hits = [];
    for j = 1:numel(kur.netParams.gT)
        subplot(numRows,numCols,2*numCols+j+1),
        imagesc(kur.hits.kur{j}), caxis([0.5,1])
        colormap('bone'), freezeColors,
        set(gca,'Xtick',[],'YTick',[])
        %
        sum_hits = [sum_hits; kur.hits.kur{j}(find(kur.hits.kur{j}))]; % find average pairwise hit-rate (weight it by cluster size?)
    end
    text(size(kur.hits.kur{j},1)+0.7,size(kur.hits.kur{j},2),num2str(mean(sum_hits),2),'HorizontalAlignment','Left','FontSize',16,'FontWeight','Bold')
    
    % ev1
    subplot(numRows,numCols,3*numCols+1), imagesc(reshape(eig.EVecsML(:,1),kur.netParams.Ndims))
    colormap('jet'), freezeColors,
    ylabel(['ev1'],'FontSize',16,'FontWeight','Bold')
    set(gca,'Xtick',[],'YTick',[])
    %
    sum_hits = [];
    for j = 1:numel(kur.netParams.gT)
        subplot(numRows,numCols,3*numCols+j+1),
        imagesc(eig.hits.ev1{j}), caxis([0.5,1])
        colormap('bone'), freezeColors,
        set(gca,'Xtick',[],'YTick',[])
        %
        sum_hits = [sum_hits; eig.hits.ev1{j}(find(eig.hits.ev1{j}))]; % find average pairwise hit-rate (weight it by cluster size?)
    end
    text(size(eig.hits.ev1{j},1)+0.7,size(eig.hits.ev1{j},2),num2str(mean(sum_hits),2),'HorizontalAlignment','Left','FontSize',16,'FontWeight','Bold')
    
    
    % ev2
    subplot(numRows,numCols,4*numCols+1), imagesc(reshape(eig.EVecsML(:,2),kur.netParams.Ndims))
    colormap('jet'), freezeColors,
    ylabel(['ev2'],'FontSize',16,'FontWeight','Bold')
    set(gca,'Xtick',[],'YTick',[])
    %
    sum_hits = [];
    for j = 1:numel(kur.netParams.gT)
        subplot(numRows,numCols,4*numCols+j+1),
        imagesc(eig.hits.ev2o{j}), caxis([0.5,1])
        colormap('bone'), freezeColors,
        set(gca,'Xtick',[],'YTick',[])
        %
        sum_hits = [sum_hits; eig.hits.ev2o{j}(find(eig.hits.ev2o{j}))]; % find average pairwise hit-rate (weight it by cluster size?)
    end
    text(size(eig.hits.ev2o{j},1)+0.7,size(eig.hits.ev2o{j},2),num2str(mean(sum_hits),2),'HorizontalAlignment','Left','FontSize',16,'FontWeight','Bold')
    
    
    % ev3
    subplot(numRows,numCols,5*numCols+1), imagesc(reshape(eig.EVecsML(:,3),kur.netParams.Ndims))
    colormap('jet'), freezeColors,
    ylabel(['ev3'],'FontSize',16,'FontWeight','Bold')
    set(gca,'Xtick',[],'YTick',[])
    %
    sum_hits = [];
    for j = 1:numel(kur.netParams.gT)
        subplot(numRows,numCols,5*numCols+j+1),
        imagesc(eig.hits.ev3o{j}), caxis([0.5,1])
        colormap('bone'), freezeColors,
        set(gca,'Xtick',[],'YTick',[])
        %
        sum_hits = [sum_hits; eig.hits.ev3o{j}(find(eig.hits.ev3o{j}))]; % find average pairwise hit-rate (weight it by cluster size?)
    end
    text(size(eig.hits.ev3o{j},1)+0.7,size(eig.hits.ev3o{j},2),num2str(mean(sum_hits),2),'HorizontalAlignment','Left','FontSize',16,'FontWeight','Bold')
    
    
    % colorbar up in {1,1}
    subplot(numRows,numCols,1)
    colormap('bone'), caxis([0.5,1])
    cb=colorbar('Location','South'); axis off
    title({'pairwise','hit-rate'},'FontSize',16,'FontWeight','Bold')
    set(cb,'XTick',[0.5,1],'FontSize',16,'FontWeight','Bold')
    
    
    % Include image patch information and segmentation parameters in an annotation.
    annotation('textbox', [0 0 1 0.1],'String', [kur.netflags.method,'   ',kur.kurflags.fname,'   ',kur.kurflags.KurTitleTag,'   '], ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')
    

    % Save the figure so I can flip thru them.
    saveGoodImg(H,[imgOutDir,kur.kurflags.fname,kur.kurflags.KurParamsTag],sizeGoodIm)
    close(H);
    
    
end