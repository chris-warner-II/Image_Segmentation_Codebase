function [] = explore_Separation_vs_Parameters(fileGeneral,fileSize,fileSubset,tscale)


% This function will take output mat file from main_Kuramoto function that is titled
% metaClusterKur_...  This will plot the phase of each oscillator at a 60Hz
% clock and will colorcode them by ground truth segment.
%
%
% Inputs to function must be of this form
%
% fileGeneral = 'GradientBox'; % 'BSDS_patch'; % 
% fileSize = '21x21_ds1';
% fileSubset = '_'; % '_100075_ptch7_';  %          % certain BSDS image & patch.  Can leave blank too to get all files in dir..
% tscale = '1';                                     % scale for phase initialization of oscillators based on input image



plot2x2Dscatter = 1; % flag to make plot of 2x2D scatter Kur vs Straw & Eig vs Straw for each individual patch. 
                     % If this flag is zero, then this function is just creating the Kur_metaSummary files if they do not exist.

                     
                     
% Directories to Input Images or Data :: Change these paths below to get a different images or patches
[dirPre,sizeGoodIm] = onCluster;
dirKur = [dirPre,'output/Kuramoto/NetsFromImgs/',fileGeneral,'_',fileSize,'/data/Kur_PIF_Fourier1/Mod_SKHAdj/']; %Theta Scale ',tscale,'/'];
dirEig = [dirPre,'output/Kuramoto/NetsFromImgs/',fileGeneral,'_',fileSize,'/data/spectral/Mod_SKHAdj/'];
dirImg = [dirPre,'images/',fileGeneral,'/',fileSize,'/'];

filesKur = dir([dirKur,'KurMC',fileSubset,'*tscale',tscale,'*']); % can loop through these later
filesEig = dir([dirEig,'Evecs',fileSubset,'*tscale',tscale,'*']); % will be in a different directory later.

% Directory to save images into (histograms, DivMarg Distributions for top & bottom parameter settings, etc...)
imgKur = [dirPre,'output/Kuramoto/NetsFromImgs/',fileGeneral,'_',fileSize,'/imgs/Kur_PIF_Fourier1/Mod_SKHAdj/',fileSubset,'/'];
%
if ~exist(imgKur,'dir')
    mkdir(imgKur)
end

% A Picture of a colorwheel I will use as a colorbar for oscillator phase
colorwheel = imread([dirPre,'images/HSV_colorwheel.jpeg']);

load([dirKur,filesKur(1).name])

% Note:  This will work ok up to 7 clusters (cludge)
colors = ['bgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmy',...
          'bgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmybgrkcmy'];
stylee = cell(1,numel(netParams.gT{1}));
for i = 1:numel(netParams.gT{1})
    stylee{i} = colors(netParams.gT{1}(i)); % This only uses 1st groundtruth segmentation.
end





% If Patch from BSDS images then this shows full image, patch and segmentation at pixel mean
if(0)
    
    fKur = filesKur(1).name;
    st = 7;
    nd = strfind(fKur,'_rM')-1;
    fImg = fKur(st:nd); 
    load([dirImg,fImg]);
    %
    try
        explore_BSDS_patches % figure handle in script is hIm
    catch
        hIm=figure; imagesc(im), colormap('bone'), title('Gradient Box Input Image')
    end
    
    % save this image
    saveGoodImg(hIm,[imgKur,'Input_Image'],sizeGoodIm)
    close(hIm)

end


%% Look Into Statistics of Separability across Many Runs of Kuramoto Simulation.  
if(1)
    
    % Check if mat file looping thru single parameter metaCluster analysis
    % mat files exists.  If not, loop thru all those single parameter files
    % and generate this accumulated statistics file.
    if exist([dirKur,'Kur_metaSummary_',fileGeneral,fileSubset,'_tscale',tscale,'_',num2str(numel(filesKur)),'files.mat'],'file') 
        
        % This is kinda hacky and may throw errors if I start plotting distribution of results over multiple runs later.  But for now, 
        % I only want to load the metaSummary file if I am also plotting. Otherwise, if I am not plotting, leave metaSummary file alone
        % if it already exists and make it if it doesnt. This should make it run quicker since I dont really care about these individual
        % patch scatter plots anyway.
        if(plot2x2Dscatter)
            disp(['Loading file from previous Summary run using ',num2str(numel(filesKur)),' files.'])
            load([dirKur,'Kur_metaSummary_',fileGeneral,fileSubset,'_tscale',tscale,'_',num2str(numel(filesKur)),'files.mat'])  
        end
        
    else
    
        disp('Time to loop through all metaCluster mat files:')
        tic

        % Histogram (across simulation runs) of DivM / Relative Margin / Pairwise distance measure
        numBins = 10;
        HistBins = linspace(0,1,numBins);

        % Preallocate Arrays to hold statistics across runs for all parameter settings.
        sigW = zeros(1,numel(filesKur));
        Kscale = zeros(1,numel(filesKur));
        sigP = zeros(1,numel(filesKur));
        sigD = zeros(1,numel(filesKur));
        Rmax = zeros(1,numel(filesKur));
        %
        DivM_End_Hist = zeros( numel(filesKur), numBins, numel(netParams.gT) );
        DivM_Mean_Hist = zeros( numel(filesKur), numBins, numel(netParams.gT) );
        %
        runs=1;
        DivM_End = ones( numel(filesKur), numel(netParams.gT), runs ); %numel(MC) );
        DivM_Mean = ones( numel(filesKur), numel(netParams.gT), runs ); %numel(MC) ); % replace runs with numel(MC) or runs.
        DivM_EV = ones( numel(filesKur), numel(netParams.gT));
        

        for i = 1:numel(filesKur)

            % Load file containing metaClustering Analysis results for Kuramoto & StrawMan for parameter combination
            try
                load([dirKur,filesKur(i).name])
            catch
                disp('Kuramoto file is corrupted:  Deleting it and Skipping to next...')
                [dirKur,filesKur(i).name]
                delete([dirKur,filesKur(i).name])
                continue
                % note: files can get corrupted on cluster if I kill jobs
                % right in the middle of a run (after file has been created
                % but before it is finished or closed properly)
            end

            % Load in Eigenvector Segmentation results too so we can compare all 3 (Kuramoto, Eigenvector, Strawman)
            try
                Eig = load([dirEig,'Evecs_',kurflags.fname,'_rM',netflags.rM,' sD',netflags.sD,' sP',netflags.sP,'.mat']);
            catch
                disp('Eigen file is corrupted:  Deleting it and Skipping to next...')
                [dirEig,'Evecs_',kurflags.fname,'_rM',netflags.rM,' sD',netflags.sD,' sP',netflags.sP,'.mat']
                delete([dirEig,'Evecs_',kurflags.fname,'_rM',netflags.rM,' sD',netflags.sD,' sP',netflags.sP,'.mat'])
                continue
            end
            
            
%             MC{2} = MC{1};
%             MC{3} = MC{1}; % to test with multiple runs.
%             MC{4} = MC{1}; % delete this eventually when
%             MC{5} = MC{1}; % I am satisfied it works consistently.

            
            fKur = filesKur(i).name;
            st = 7; % strfind(fKur,'_rM')+1;
            fImg = fKur(st:end); 

            runs = numel(MC);

            DivMarg = zeros( size(MC{1}.DistAvgPW,1), numel(netParams.gT), runs );     % Divisive Margin for Coupled Osc. Model
            DivMarg_SM = zeros( 1, numel(netParams.gT));                               % Divisive Margin for Straw Man Model
            % DivMarg_EV = zeros( 1, numel(netParams.gT));                               % Divisive Margin for Eigenvector Computation
            
            for g = 1:numel(netParams.gT)
            	DivMarg_SM(1,g) = MCsm{g}.DistAvgPW(1,1) ./ MCsm{g}.DistAvgPW(1,2);   
                DivM_EV(i,g) = Eig.MC{g}.DistAvgPW(1,1) ./ Eig.MC{g}.DistAvgPW(1,2); 
                
                for r = 1:runs
                    DivMarg(:,g,r) =  MC{r}.DistAvgPW(:,1,g) ./ MC{r}.DistAvgPW(:,2,g);
                end
            end
          
            DivMarg(isnan(DivMarg)) = 1;        % set to 1 if DivMarg = 0/0
            DivMarg_SM(isnan(DivMarg_SM)) = 1;
            DivM_EV(isnan(DivM_EV)) = 1;

            DivM_mean_of_run = reshape( mean(DivMarg,1), numel(netParams.gT), runs ); % value averaged across entire simulation time
            DivM_end_of_run = reshape( DivMarg(end,:,:), numel(netParams.gT), runs ); % final value at end of simulation

                
            % distribution across all r runs
            if runs>1
                [cntMn] = hist(DivM_mean_of_run',HistBins);
                [cntNd] = hist(DivM_end_of_run',HistBins);
            else
                for g = 1:numel(netParams.gT)
                    cntMn(:,g) = hist(DivM_mean_of_run(g),HistBins);
                    cntNd(:,g) = hist(DivM_end_of_run(g),HistBins);
                    % makes no sense to calculate histogram with single run but I do it for continuity
                    % sake because I am saving results of histogram calculation
                end
            end
            
            

            % find a run indicative of each bar in histogram (for Mean & End of Simulation)
            for k = 1:numel(HistBins)
                
                for g = 1:numel(netParams.gT)
                
                    diffFrmPeakMn = abs(DivM_mean_of_run(g,:)-HistBins(k));
                    tmp = find(diffFrmPeakMn == min(diffFrmPeakMn));
                    runInd_histMn(k) = tmp(1); % take only 1st if multiple

                    diffFrmPeakNd = abs(DivM_end_of_run(g,:)-HistBins(k));
                    tmp = find(diffFrmPeakNd == min(diffFrmPeakNd));
                    runInd_histNd(k) = tmp(1); % take only 1st if multiple
                
                end
                
            end
            
            DivM_End_Hist(i,:,:) = cntNd;
            DivM_Mean_Hist(i,:,:) = cntMn;
            

            if(0)
                % Plot time evolution of DivMarg for all 100 simulations.
                % Also, plot histogram of final value reached.  And plot
                % exemplar time evolution & final phase distributions for each
                % bin.

                plot_DivMarg_Distrib

            end

            disp(['File #',num2str(i),' / ',num2str(numel(filesKur))])
            
            % Can analyze these below like we do CS & AUC (only works for single run - if deterministic).
            DivM_End(i,:,:) = DivM_end_of_run; % want to save a vector with Divisive Margin for each combination of parameters.
            DivM_Mean(i,:,:) = DivM_mean_of_run;
            %
            % DivM_SM(i,:) = DivMarg_SM; % have to do a different one for each ground truth (or average over).
            
            % Save parameter values into array so we can visualize performance vs parameters.
            sigW(i) = kurParams.sigW;
            Kscale(i) = kurParams.Kscale;
            sigP(i) = netParams.sigP;
            sigD(i) = netParams.sigD;  % note sigD is related to Rmax (either its inf or 1/4*Rmax)
            Rmax(i) = netParams.Rmax;
            
      
        end % loop over files

        toc


        CS = cumsum(DivM_End_Hist,2);
        AUC = squeeze( trapz(CS,2) ./ trapz(repmat(runs,1,numBins)) );
        
        
        uSW = unique(sigW);
        uKS = unique(Kscale);
        uSP = unique(sigP);
        uSD = unique(sigD); % note sigD is related to Rmax (either its inf or 1/4*Rmax)
        uRM = unique(Rmax);
    
    

        % Running this loop over metaCluster mat files takes a long time (~20 mins for ~216 files with 100 runs in each).  
        % So I should save a mat file containing what I am generating inside there.
        save([dirKur,'Kur_metaSummary_',fileGeneral,fileSubset,'_tscale',tscale,'_',num2str(numel(filesKur)),'files.mat'],...
            'HistBins','numBins','sigW','Kscale','sigP','sigD','Rmax','uSW','uKS','uSP','uSD','uRM','runs',...
            'DivM_End_Hist','DivM_Mean_Hist','CS','AUC','DivM_End','DivM_Mean','DivMarg_SM','DivM_EV','filesKur') % 'CS_cube',
        
        
        % If some other version of Kur_metaSummary file exists with fewer files (parameter combinations), delete it.
        %
        % NOT SURE HOW TO DO THIS, BUT WORK ON IT.
        %
        
        
    end % if Kur_metaSummary exists
    

    
%     if ~exist('numBins','var')
%         numBins = numel(HistBins); 
%     end
    if ~exist('runs','var')
        runs = 1; 
    end
    

    
    % Plot Statistics (mean Std) across 100 runs of max CSep value reached
    % during a simulation for different parameter settings.
    if(plot2x2Dscatter)
        CS_mngT = mean(CS,3);                    % Average over different Ground Truths for plotting.
        AUC_mngT = mean(AUC,2);                  %
        DivM_End_mnRn = mean(DivM_End,3);        % average DivM_End across runs (better to look at CS & AUC for multiple runs)
        DivM_End_mnRnGt = mean(DivM_End_mnRn,2); % average DivM End across runs, then across ground truths (for plots only)
    end
    
    if(0)
        
        % (0). Plot Cumulative Sum (averaged over Ground Truths)
        h1=figure;
        subplot(3,4,[1:3,5:7]), imagesc( CS_mngT ), colorbar
        set(gca,'XTickLabel',sprintf('%1.1f|',HistBins))
        xlabel('Relative Margin - Pairwise Distance Within vs. Across Cluster')
        ylabel('Single Instantiation of Parameter Value Combination')
        title({'Probability that Simulation does better than given RM Value',...
            'Brightness for Lower RM values indicates better Segmentation'},...
            'FontSize',14,'FontWeight','Bold')
        %
        subplot(3,4,[4,8]), hold on,
        plot( 1-AUC_mngT ,1:numel(filesKur) ,'b'), 
        plot( DivM_End_mnRnGt ,1:numel(filesKur) ,'r'), 
        axis ij tight
        title('\color{blue}{(1-AUC)} \color{black}{&} \color{red}{DivMarg}')
        xlabel('smaller is better \in  (0,1)')
        %
        subplot(3,1,3), hist( AUC_mngT )
        title('Histogram of Area Under CumSum Values')
        xlabel('smaller is better \in  (0,1)')
        %
        saveGoodImg(h1,[imgKur,'cumSumStats_all_tscale',tscale],sizeGoodIm)
        close(h1)

        
       
        
        
        
        % (1). Histogram when we Vary SigW
        x = zeros(size(CS_mngT));
        y = zeros(size(AUC_mngT));
        z = ones(size(DivM_End_mnRnGt));
        cur=0; 
        brk(1)=cur;
        for k = 1:numel(uSW)
            ind = find(sigW==uSW(k));
            x(brk(k)+1:brk(k)+numel(ind),:) = CS_mngT(ind,:);
            y(brk(k)+1:brk(k)+numel(ind)) = AUC_mngT(ind);
            z(brk(k)+1:brk(k)+numel(ind)) = DivM_End_mnRnGt(ind);
            cur = numel(ind); 
            brk(k+1) = brk(k)+cur;
        end
        %
        for k = 1:numel(brk)-1
            lbl(k) = mean([brk(k),brk(k+1)]);
        end
        %
        h2=figure;
        subplot(1,4,1:3), imagesc(x), colorbar
        hold on, 
        for k = 1:numel(brk)
            plot([0 numBins], [brk(k) brk(k)],'w--','LineWidth',2)
        end
        %
        set(gca,'XTickLabel',sprintf('%1.1f|',HistBins),'YTick',lbl,'YTickLabel',sprintf('%1.1f|',uSW))
        xlabel('Relative Margin - Pairwise Distance Within vs. Across Cluster')
        ylabel('Separating by \sigma_w value')
        title({'Probability that Simulation does better than given RM Value',...
            'Brightness for Lower RM values indicates better Segmentation'},...
            'FontSize',14,'FontWeight','Bold')
        %
        subplot(1,4,4), hold on,
        plot( 1-y, 1:numel(filesKur), 'b'), 
        plot( z, 1:numel(filesKur), 'r'), 
        axis ij tight
        for k = 1:numel(brk)
            plot([0 max([(1-y);z])], [brk(k) brk(k)],'k--','LineWidth',2)
        end
        set(gca,'YTick',[])
        title('\color{blue}{(1-AUC)} \color{black}{&} \color{red}{DivMarg}')
        xlabel('smaller is better \in  (0,1)')
        %
        clear brk lbl
        saveGoodImg(h2,[imgKur,'cumSumStats_sigW_tscale',tscale],sizeGoodIm)
        close(h2)
        
        
        % (2). Histogram when we Vary KScale
        x = zeros(size(CS_mngT));
        y = zeros(size(AUC_mngT));
        z = ones(size(DivM_End_mnRnGt));
        cur=0; 
        brk(1)=cur;
        for k = 1:numel(uKS)
            ind = find(Kscale==uKS(k));
            x(brk(k)+1:brk(k)+numel(ind),:) = CS_mngT(ind,:);
            y(brk(k)+1:brk(k)+numel(ind)) = AUC_mngT(ind);
            z(brk(k)+1:brk(k)+numel(ind)) = DivM_End_mnRnGt(ind);
            cur = numel(ind); 
            brk(k+1) = brk(k)+cur;
        end
        %
        for k = 1:numel(brk)-1
            lbl(k) = mean([brk(k),brk(k+1)]);
        end
        %
        h3=figure;
        subplot(1,4,1:3), imagesc(x), colorbar
        hold on, 
        for k = 1:numel(brk)
            plot([0 numBins], [brk(k) brk(k)],'w--','LineWidth',2)
        end
        %
        set(gca,'XTickLabel',sprintf('%1.1f|',HistBins),'YTick',lbl,'YTickLabel',sprintf('%1.1f|',uKS))
        xlabel('Relative Margin - Pairwise Distance Within vs. Across Cluster')
        ylabel('Separating by Kscale value')
        title({'Probability that Simulation does better than given RM Value',...
            'Brightness for Lower RM values indicates better Segmentation'},...
            'FontSize',14,'FontWeight','Bold')
        %
        subplot(1,4,4), hold on,
        plot( 1-y, 1:numel(filesKur), 'b'), 
        plot( z, 1:numel(filesKur), 'r'), 
        axis ij tight
        for k = 1:numel(brk)
            plot([0 max([(1-y);z])], [brk(k) brk(k)],'k--','LineWidth',2)
        end
        set(gca,'YTick',[])
        title('\color{blue}{(1-AUC)} \color{black}{&} \color{red}{DivMarg}')
        xlabel('smaller is better \in  (0,1)')
        %
        clear brk lbl
        saveGoodImg(h3,[imgKur,'cumSumStats_Kscale_tscale',tscale],sizeGoodIm)
        close(h3)
        
        
        
        % (3). Histogram when we Vary sigP
        x = zeros(size(CS_mngT));
        y = zeros(size(AUC_mngT));
        z = ones(size(DivM_End_mnRnGt));
        cur=0; 
        brk(1)=cur;
        for k = 1:numel(uSP)
            ind = find(sigP==uSP(k));
            x(brk(k)+1:brk(k)+numel(ind),:) = CS_mngT(ind,:);
            y(brk(k)+1:brk(k)+numel(ind)) = AUC_mngT(ind);
            z(brk(k)+1:brk(k)+numel(ind)) = DivM_End_mnRnGt(ind);
            cur = numel(ind); 
            brk(k+1) = brk(k)+cur;
        end
        %
        for k = 1:numel(brk)-1
            lbl(k) = mean([brk(k),brk(k+1)]);
        end
        %
        h4=figure;
        subplot(1,4,1:3), imagesc(x), colorbar
        hold on, 
        for k = 1:numel(brk)
            plot([0 numBins], [brk(k) brk(k)],'w--','LineWidth',2)
        end
        %
        set(gca,'XTickLabel',sprintf('%1.1f|',HistBins),'YTick',lbl,'YTickLabel',sprintf('%1.1f|',uSP))
        xlabel('Relative Margin - Pairwise Distance Within vs. Across Cluster')
        ylabel('Separating by \sigma_P value')
        title({'Probability that Simulation does better than given RM Value',...
            'Brightness for Lower RM values indicates better Segmentation'},...
            'FontSize',14,'FontWeight','Bold')
        %
        subplot(1,4,4), hold on,
        plot( 1-y, 1:numel(filesKur), 'b'), 
        plot( z, 1:numel(filesKur), 'r'), 
        axis ij tight
        for k = 1:numel(brk)
            plot([0 max([(1-y);z])], [brk(k) brk(k)],'k--','LineWidth',2)
        end
        set(gca,'YTick',[])
        title('\color{blue}{(1-AUC)} \color{black}{&} \color{red}{DivMarg}')
        xlabel('smaller is better \in  (0,1)')
        %
        clear brk lbl
        saveGoodImg(h4,[imgKur,'cumSumStats_sigP_tscale',tscale],sizeGoodIm)
        close(h4)
        
        
        % (4). Histogram when we Vary Rmax with binary distance dependence (sigD=inf)
        x = zeros(ceil(size(CS_mngT,1)/2), size(CS_mngT,2));
        y = zeros(ceil(size(AUC_mngT,1)/2), 1);
        z = ones(ceil(numel(DivM_End_mnRnGt)/2), 1);
        cur=0; 
        brk(1)=cur;
        for k = 1:numel(uRM)
            ind = find(Rmax==uRM(k) & sigD==inf);
            x(brk(k)+1:brk(k)+numel(ind),:) = CS_mngT(ind,:);
            y(brk(k)+1:brk(k)+numel(ind)) = AUC_mngT(ind);
            z(brk(k)+1:brk(k)+numel(ind)) = DivM_End_mnRnGt(ind);
            cur = numel(ind); 
            brk(k+1) = brk(k)+cur;
        end
        %
        for k = 1:numel(brk)-1
            lbl(k) = mean([brk(k),brk(k+1)]);
        end
        %
        h5=figure;
        subplot(1,4,1:3), imagesc(x), colorbar
        hold on, 
        for k = 1:numel(brk)
            plot([0 numBins], [brk(k) brk(k)],'w--','LineWidth',2)
        end
        %
        set(gca,'XTickLabel',sprintf('%1.1f|',HistBins),'YTick',lbl,'YTickLabel',sprintf('%1.1f|',uRM))
        xlabel('Relative Margin - Pairwise Distance Within vs. Across Cluster')
        ylabel('Separating by Rmax value with sD = \infty')
        title({'Probability that Simulation does better than given RM Value',...
            'Brightness for Lower RM values indicates better Segmentation'},...
            'FontSize',14,'FontWeight','Bold')
        %
        subplot(1,4,4), hold on,
        plot( (1-y)', 1:numel(y), 'b'), 
        plot( z, 1:numel(z), 'r'), 
        axis ij tight
        for k = 1:numel(brk)
            plot([0 max([(1-y);z])], [brk(k) brk(k)],'k--','LineWidth',2)
        end
        set(gca,'YTick',[])
        title('\color{blue}{(1-AUC)} \color{black}{&} \color{red}{DivMarg}')
        xlabel('smaller is better \in  (0,1)')
        %
        clear brk lbl
        saveGoodImg(h5,[imgKur,'cumSumStats_Rmax_sDInf_tscale',tscale],sizeGoodIm)
        close(h5)
        
        
        
        
        
        
        % (4b). Histogram when we Vary Rmax with continuous falloff in distance dependence (sigD=(1/4)*Rmax)
        x = zeros(ceil(size(CS_mngT,1)/2), size(CS_mngT,2));
        y = zeros(ceil(numel(AUC_mngT)/2), 1);
        z = ones(ceil(numel(DivM_End_mnRnGt)/2), 1);
        cur=0; 
        brk(1)=cur;
        for k = 1:numel(uRM)
            ind = find(Rmax==uRM(k) & sigD~=inf);
            x(brk(k)+1:brk(k)+numel(ind),:) = CS_mngT(ind,:);
            y(brk(k)+1:brk(k)+numel(ind)) = AUC_mngT(ind);
            z(brk(k)+1:brk(k)+numel(ind)) = DivM_End_mnRnGt(ind);
            cur = numel(ind); 
            brk(k+1) = brk(k)+cur;
        end
        %
        for k = 1:numel(brk)-1
            lbl(k) = mean([brk(k),brk(k+1)]);
        end
        %
        h6=figure;
        subplot(1,4,1:3), imagesc(x), colorbar
        hold on, 
        for k = 1:numel(brk)
            plot([0 numBins], [brk(k) brk(k)],'w--','LineWidth',2)
        end
        %
        set(gca,'XTickLabel',sprintf('%1.1f|',HistBins),'YTick',lbl,'YTickLabel',sprintf('%1.1f|',uRM))
        xlabel('Relative Margin - Pairwise Distance Within vs. Across Cluster')
        ylabel('Separating by Rmax value with sD = Rmax/4')
        title({'Probability that Simulation does better than given RM Value',...
            'Brightness for Lower RM values indicates better Segmentation'},...
            'FontSize',14,'FontWeight','Bold')
        %
        subplot(1,4,4), hold on,
        plot( (1-y)', 1:numel(y), 'b'), 
        plot( z, 1:numel(z), 'r'), 
        axis ij tight
        for k = 1:numel(brk)
            plot([0 max([(1-y);z])], [brk(k) brk(k)],'k--','LineWidth',2)
        end
        set(gca,'YTick',[])
        title('\color{blue}{(1-AUC)} \color{black}{&} \color{red}{DivMarg}')
        xlabel('smaller is better \in  (0,1)')
        %
        clear brk lbl
        saveGoodImg(h6,[imgKur,'cumSumStats_Rmax_sDQtr_tscale',tscale],sizeGoodIm)
        close(h6)
        
    end % if(1) to plot these histograms for different dimensions in param Space.
    
    
    
    
    % Scatter Plot Straw-Man Model Performance vs. Coupled Oscillator Model Performance
    if(plot2x2Dscatter)
        
        sigP_shapes = 'o^sdxxxx'; % adding x's to avoid out of bounds errors
        Rmax_colors = 'rbgkcmmm'; % adding m's to avoid out of bounds errors
        %
        StrawMan = repmat(DivMarg_SM,numel(filesKur),1);
        
        
        % Scatter Plot Kuramoto vs Strawman & Kuramoto vs. Eigenvector
        hSc=figure; 
        subplot(121), hold on, subplot(122), hold on
        %
        %
        for a = 1:numel(uSP)
            for b = 1:numel(uRM)
                ind = find(sigP==uSP(a) & Rmax==uRM(b));
                x = StrawMan(ind,:);
                y = DivM_End_mnRn(ind,:);
                z = DivM_EV(ind,:);
                subplot(121), scatter(x(:), y(:), 100, [Rmax_colors(b),sigP_shapes(a)],'LineWidth',2)
                subplot(122), scatter(z(:), y(:), 100, [Rmax_colors(b),sigP_shapes(a)],'LineWidth',2)
            end
        end
        %
        %
        subplot(121)
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg of Coupled Oscillator Model','FontSize',18,'FontWeight','Bold')
    	xlabel('DivMarg of Pixel Contrast Model','FontSize',18,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold')
        text(0.2,0.25,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold')
        %
        text(0.1,1.0,['*Note: Averaging across ',num2str(runs),' run.'])
        text(0.1,0.98,['*Legend.'])
        text(0.1,0.96,['*SigP Shapes = ',sigP_shapes])
        text(0.1,0.94,['*Rmax Colors = ',Rmax_colors])
        %
        subplot(122)
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg of Coupled Oscillator Model','FontSize',18,'FontWeight','Bold')
    	xlabel('DivMarg of Eigenvector Computation','FontSize',18,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold')
        text(0.2,0.25,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold')
        %

        title([fileGeneral,' ',fileSubset,' ',fileSize],'FontSize',20,'FontWeight','Bold')
        
        %
        saveGoodImg(hSc,[imgKur,'DivMarg_Scatter_2D_tscale',tscale],sizeGoodIm)
        close(hSc)
        
        
    end
    
    

    
    % Look DivM distribution plots for top (and bottom) couple parameter settings.
    if(0)
        
        [p,d] = sort(AUC_mngT,'descend');

        numtop = 0;
        numbot = 0; % Change these from zero to plot some example best and worst runs.

        top = d(1:numtop); 
        bot = d(end-numbot+1:end);

        % save this AUC value distribution image
        h=figure; stem(p), title(['Distribution of Performance for Tscale = ',tscale])
        xlabel(['Different Parameter Combinations (sigW, sigD, sigP, Rmax, Kscale)'])
        ylabel(['AUC of CDF of Histogram of Ending DivMarg'])
        %
        saveGoodImg(h,[imgKur,'AUC_Dist_sorted_tscale',tscale],sizeGoodIm)
        close(h)





        % Plot DivMarg_Distrib code for the K best parameter settings.
        for i = 1:numel(top)

            load([dirKur,filesKur(top(i)).name])

            fKur = filesKur(top(i)).name;
            st = 7; % strfind(fKur,'_rM')+1;
            fImg = fKur(st:end); 

            runs = numel(MC);

            DivMarg = zeros( runs, size(MC{1}.DistAvgPW,1) );     % Divisive Margin
            StdAcrossGt = zeros( size(MC{1}.DistAvgPW,1), size(MC{1}.DistAvgPW,2), runs ); % should be low if things are work.
            AvgAcrossGt = zeros( size(MC{1}.DistAvgPW,1), size(MC{1}.DistAvgPW,2), runs );

            % mean and std across groundtruth segmentations (good if std << mean)
            for r = 1:runs
                MC{r}.DistAvgPW_StdGT = std(MC{r}.DistAvgPW,[],3);
                MC{r}.DistAvgPW_AvgGT = mean(MC{r}.DistAvgPW,3);
                %
                DivMarg(r,:) =  MC{r}.DistAvgPW_AvgGT(:,1) ./ MC{r}.DistAvgPW_AvgGT(:,2);
                %
                StdAcrossGt(:,:,r) = MC{r}.DistAvgPW_StdGT;
                AvgAcrossGt(:,:,r) = MC{r}.DistAvgPW_AvgGT;
                
                % NOTE: I AM TRYING TO GET AWAY FROM CALCULATING & STD OF DIVMARG ACROSS GROUND TRUTHS

            end

            DivMarg(isnan(DivMarg)) = 1; % set to 1 if DivMarg = 0/0


%             % Variation of performance across different Human Segmentations.
%             % Think of a better way to do this later... Maybe include this
%             % single number in the plot_DivMarg_Distrib figure.
%             if(0)
%                 disp('Distance between Oscillator Pairs')
%                 disp('Std Across different Human Ground Truth Segmentations')
%                 disp('Averaging across Runs, Time & whether same or different cluster.')
%                 disp('Characterizing similarity of results using different Ground Truths.')
%                 disp('Good if Small')
%                 mean(StdAcrossGt(:))
%             end


            DivM_mean_of_run = mean(DivMarg,2); % value averaged across entire simulation time
            DivM_end_of_run = DivMarg(:,end);   % final value at end of simulation

            % distribution across all 100 runs
            [cntMn] = hist(DivM_mean_of_run,HistBins);
            [cntNd] = hist(DivM_end_of_run,HistBins);

            % find a run indicative of each bar in histogram (for Mean & End of Simulation)
            for g = 1:numel(HistBins)
                diffFrmPeakMn = abs(DivM_mean_of_run-HistBins(g));
                tmp = find(diffFrmPeakMn == min(diffFrmPeakMn));
                runInd_histMn(g) = tmp(1);

                diffFrmPeakNd = abs(DivM_end_of_run-HistBins(g));
                tmp = find(diffFrmPeakNd == min(diffFrmPeakNd));
                runInd_histNd(g) = tmp(1);
            end

            if(1)
                % Plot time evolution of DivMarg for all 100 simulations.
                % Also, plot histogram of final value reached.  And plot
                % exemplar time evolution & final phase distributions for each
                % bin.

                plot_DivMarg_Distrib

            end

        end % loop over files with best performance (top)







        % Plot DivMarg_Distrib code for the K worst parameter settings.
        for i = 1:numel(bot)

            load([dirKur,filesKur(bot(i)).name])

            fKur = filesKur(bot(i)).name;
            st = 7; % strfind(fKur,'_rM')+1;
            fImg = fKur(st:end); 

            runs = numel(MC);

            DivMarg = zeros( runs, size(MC{1}.DistAvgPW,1) );     % Divisive Margin
            StdAcrossGt = zeros( size(MC{1}.DistAvgPW,1), size(MC{1}.DistAvgPW,2), runs ); % should be low if things are work.
            AvgAcrossGt = zeros( size(MC{1}.DistAvgPW,1), size(MC{1}.DistAvgPW,2), runs );

            % mean and std across groundtruth segmentations (good if std << mean)
            for r = 1:runs
                MC{r}.DistAvgPW_StdGT = std(MC{r}.DistAvgPW,[],3);
                MC{r}.DistAvgPW_AvgGT = mean(MC{r}.DistAvgPW,3);
                %
                DivMarg(r,:) =  MC{r}.DistAvgPW_AvgGT(:,1) ./ MC{r}.DistAvgPW_AvgGT(:,2);
                %
                StdAcrossGt(:,:,r) = MC{r}.DistAvgPW_StdGT;
                AvgAcrossGt(:,:,r) = MC{r}.DistAvgPW_AvgGT;
                
                % NOTE: I AM TRYING TO GET AWAY FROM CALCULATING & STD OF DIVMARG ACROSS GROUND TRUTHS
                
            end

            DivMarg(isnan(DivMarg)) = 1; % set to 1 if DivMarg = 0/0


%             % Variation of performance across different Human Segmentations.
%             % Think of a better way to do this later... Maybe include this
%             % single number in the plot_DivMarg_Distrib figure.
%             if(0)
%                 disp('Distance between Oscillator Pairs')
%                 disp('Std Across different Human Ground Truth Segmentations')
%                 disp('Averaging across Runs, Time & whether same or different cluster.')
%                 disp('Characterizing similarity of results using different Ground Truths.')
%                 disp('Good if Small')
%                 mean(StdAcrossGt(:))
%             end


            DivM_mean_of_run = mean(DivMarg,2); % value averaged across entire simulation time
            DivM_end_of_run = DivMarg(:,end);   % final value at end of simulation

            % distribution across all 100 runs
            [cntMn] = hist(DivM_mean_of_run,HistBins);
            [cntNd] = hist(DivM_end_of_run,HistBins);

            % find a run indicative of each bar in histogram (for Mean & End of Simulation)
            for g = 1:numel(HistBins)
                diffFrmPeakMn = abs(DivM_mean_of_run-HistBins(g));
                tmp = find(diffFrmPeakMn == min(diffFrmPeakMn));
                runInd_histMn(g) = tmp(1);

                diffFrmPeakNd = abs(DivM_end_of_run-HistBins(g));
                tmp = find(diffFrmPeakNd == min(diffFrmPeakNd));
                runInd_histNd(g) = tmp(1);
            end

            if(1)
                % Plot time evolution of DivMarg for all 100 simulations.
                % Also, plot histogram of final value reached.  And plot
                % exemplar time evolution & final phase distributions for each
                % bin.

                plot_DivMarg_Distrib

            end

        end % loop over files with worst performance (bot)
    end
    
    
    % Write text file: List Parameters, AUC value & of Performance Equivalence Classes.
    if(0)
        
        xx = [AUC_mngT(d)'; DivM_End_mnRnGt   sigD(d); sigP(d); sigW(d); Rmax(d); Kscale(d)];
        fid = fopen([imgKur,'AUC_Dist_sorted_tscale',tscale,'.txt'],'w');
        %
        fprintf(fid,['Parameter Combinations for Tscale = ',tscale,' sorted from best to worst\n']);
        fprintf(fid,' \n');
        fprintf(fid,['   AUC    sigD    sigP    sigW    Rmax  Kscale\n']);
        fprintf(fid,' \n');
        fprintf(fid,'%6.2f  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f\n',xx);
        fclose(fid);

        % maybe put DivMarg here in addition to AUC or instead of.
        keyboard

        DivM_End_mnRnGt
        DivM_End_mnRn
        DivMarg_SM

        % What to save to txt file or mat file for next analysis step?

    end


end % if(1) Look at Kuramoto Stats big if statement...














% 
% 
% 
% 
% %% Plots for a single run.  Choose Run Number Here.
% if(0)
%     
%     load([dirPre,dirRest,'metaClusterKur_rM3_sDInf_sP0p4_NF_60_0_kscale300_runs100.mat'])
%     
%     
%     
%     % Plot statistics of max Csep & mean Csep.
%     for r = 1:runs
%         CSep_mean_of_run(r) = max(MC{r}.phase.meanCSep);
%         CSep_max_of_run(r) = mean(MC{r}.phase.meanCSep);
%     end
%     
%     
%     
%     r=50; % run number
% 
% 
% 
% 
%     %% (1). Plot Center at cluster distance and shaded errorbar showing Cluster Extent.  Good if shaded region is positive.
%     if(0)
%         figure, shadedErrorBar([1:kurParams.T]./(kurParams.spp.*kurParams.muW), MC{r}.phase.meanCDist, MC{r}.meanCExt)
%         hold on, plot([1:kurParams.T]./(kurParams.spp.*kurParams.muW), zeros(1,kurParams.T), 'k--')
%     end
% 
% 
% 
% 
% 
% 
% 
%     %% (2). Plot Phase of Oscillators at 60Hz Clock with Separability Metric below.
%     h=figure; 
%     subplot(2,1,1), hold on  %  
%     
%     
%     % netParams, metaCluster, kurParams, MC, stylee
% 
%     for i = 1:netParams.N
% 
%         % Dont want to have plots wrap around because its distracting
%         dude = 2*pi*mod(metaCluster(r).t_2pi{i},1/kurParams.muW)/(1/kurParams.muW);
%         tt = 0:numel(dude);
%         %
%         dl = diff(dude);
%         jumpind = [abs(dl)>pi]; % now if jumpind(i) = true, we know that the
%         ind = find(jumpind);
%         ind = [0, ind, numel(dude)];
%         %
%         for j = 2:numel(ind)
%             dudesub = dude(ind(j-1)+1:ind(j));
%             tsub = (1./MC{r}.w(i)).*tt(ind(j-1)+1:ind(j));
%             plot(tsub,dudesub,stylee{i},'LineWidth',1)
%         end
% 
%         % Annotate oscillator number and natural frequency on plot
%         %text(t(end), dude(end),['#',num2str(i),' - \omega=',num2str(w(i),4)],'Color',stylee{i},'Fontsize',14,'FontWeight','Bold')
%     end
%     %
%     plot([0 (1./kurParams.muW).*numel(dude)], [0 0],'k--')
%     plot([0 (1./kurParams.muW).*numel(dude)], [2*pi 2*pi],'k--')
% 
% 
%     title(['Phase of 2\pi Crossing w.r.t. ',num2str(kurParams.muW),' Hz Oscillation'],'Fontsize',20,'FontWeight','Bold')
%     % xlabel(['Period # of ',num2str(muW),'Hz oscillation'],'Fontsize',18,'FontWeight','Bold')
%     ylabel(['Phase of Oscillation'],'Fontsize',18,'FontWeight','Bold')
%     set(gca,'Fontsize',16,'FontWeight','Bold','Ytick',[0 pi/2 pi 3*pi/2 2*pi],'YTickLabel',{'0','90','180','270','360'}) % ,'Interpreter','Latex'
%     axis([0 (1./kurParams.muW).*numel(dude) 0 2*pi])
% 
% 
% 
% 
%     % Plot "Cluster Separation" - Clusterability Metric.
%     subplot(2,1,2), hold on
%     shadedErrorBar([1:kurParams.T]./(kurParams.spp.*kurParams.muW), MC{r}.phase.meanCSep, MC{r}.phase.stdCSep )
% 
%     %         plot([1:kurParams.T]./(kurParams.spp.*kurParams.muW), MCout.phase.minCSep-0.2,'k', 'LineWidth', 5)
%     %         for i = 1:kurParams.T
%     %             scatter(i./(kurParams.spp.*kurParams.muW), MCout.phase.minCSep(i), 30, 'Filled', 'MarkerFaceColor',colors(MCout.phase.minCSepID(i,1)))
%     %             scatter(i./(kurParams.spp.*kurParams.muW), MCout.phase.minCSep(i) + 0.05*max(MCout.phase.minCSep), 30, 'Filled', 'MarkerFaceColor',colors(MCout.phase.minCSepID(i,2)))
%     %         end
% 
%     %         for i = 1:C
%     %             for j = i+1:C
%     %                 plot([1:T]./(kurParams.spp.*kurParams.muW), squeeze(ClusterSep(i,j,:)), colors(i), 'LineWidth', 1.2)
%     %                 plot([1:T]./(kurParams.spp.*kurParams.muW), squeeze(ClusterSep(i,j,:)), colors(j), 'LineStyle', '--', 'LineWidth', 1.0)
%     %             end
%     %         end
% 
%     plot([1:kurParams.T]./(kurParams.spp.*kurParams.muW), MC{r}.phase.meanCSep,'k--', 'LineWidth', 5)
%     title(['Cluster Separation'],'FontSize',20,'FontWeight','Bold')
%     xlabel(['Simulation Time (sec)'],'FontSize',18,'FontWeight','Bold')
%     ylabel({['Distance Between Cluster Centers /'],['Avg Extent across Clusters']},'FontSize',18,'FontWeight','Bold')
%     set(gca,'FontSize',16,'FontWeight','Bold')
% 
%     % saveGoodImg(h,[imgsKurDir,'PhasePrecessionRelToClk_',KurParamsTag,'_run',num2str(r)],sizeGoodIm)
%     % close(h);
% 
% end
% 
% 
% 
% 
% 




