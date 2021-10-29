
% THIS IS OLD.  SEE PLOT_SEPARATION_VS_PARAMETERS

[dirPre] = onCluster;

imageIn = 'grad_21x21_Box_0_1_0_1';

methods = {'AA','AAnrm','Mod N&G', 'Mod SKHAdj'};


cmapRWB = rd_plotColorbar('redwhiteblue',256);
cmapRW = rd_plotColorbar('whitered',128);



%  Looping over Methods:  For each method, Look at all parameter settings.
for i = 1 %1:numel(methods)
    
    
    method = methods{i};
    dataDir = [dirPre,'output/Kuramoto/NetsFromImgs/',imageIn,'/data/spectral/',method,'/'];
    files = dir([dataDir,'Evecs_*']);
    
    
    % Loop over different parameter settings for one method.
    for j = 1:numel(files)
       
        load([dataDir,files(j).name]);
        
        Sep(j)  = MC.phase.meanCSep;
        Rmax(j) = netParams.Rmax;
        SigD(j) = netParams.sigD;
        SigP(j) = netParams.sigP;
        
        if(0) % Plot each eigenvector
            figure, imagesc(EvecsML{1}), colorbar
            title(['Sep ',num2str(Sep(j)),' -- Rmax ',num2str(Rmax(j)),' -- SigD ',num2str(SigD(j)),' -- SigP ',num2str(SigP(j))])
        end
        
    end
    

    
    
    
    % Order Separation values into a 2D array that we can image to see best parameter settings
    RmU = unique(Rmax);
    SpU = unique(SigP);
    % This assumes that SdU is a single value (inf).  It doesnt have to be.
    
    for j = 1:numel(RmU)
        for k = 1:numel(SpU)
            ind = find( Rmax == RmU(j) & SigP == SpU(k));
            SepOrdr(j,k) = Sep(ind);
        end
    end
    
    % find parameter values with maximum Separation.
    maxInd = find(Sep==max(Sep));
    [maxSp,maxRm] = find(SepOrdr==max(Sep));

    
    
    
    
    h=figure; 
    
    % Plot Separation value vs. Parameters for an image / method.
    subplot(121), imagesc(SepOrdr), colorbar, 
    hold on, scatter(maxRm,maxSp,200,'kx')
    xlabel('\sigma Pix','FontSize',18,'FontWeight','Bold')
    ylabel('Rmax','FontSize',18,'FontWeight','Bold')
    title(['Separation for ',imageIn,' with ',method,' Evecs'],'FontSize',20,'FontWeight','Bold')
    set(gca,'XTick',1:numel(SpU),'YTick',1:numel(RmU),'XTickLabel',SpU,'YTickLabel',RmU,'FontSize',16,'FontWeight','Bold')
    caxis([-1 1])
    colormap(cmapRWB)
    
    % Plot Eigenvector with maximum separation value
    load([dataDir,files(maxInd).name]);
    subplot(222),imagesc(EvecsML{1}),
    title(['Sep ',num2str(Sep(maxInd)),' -- Rmax ',num2str(Rmax(maxInd)),' -- SigD ',num2str(SigD(maxInd)),' -- SigP ',num2str(SigP(maxInd))])
    
    % Plot Ground Truth
    subplot(224),imagesc(netParams.gndTruth{1}),
    title(['Ground Truth'])
    
    % Note: For any square, I can plot how Separation changes with
    % nonLinearity on visualization using plot_Separation_vs_EvecVizNonLin
    
    
    
    
    
    
    

end