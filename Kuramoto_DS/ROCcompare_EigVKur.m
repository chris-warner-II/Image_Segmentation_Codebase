function ROCcompare_EigVKur(vecOfDirsK)


% This function takes in the output from main_Kuramoto function which
% includes results from segmentation at each time step of the Kuramoto
% simulation, for an array of cosine distance threshold values and a number
% of runs.  This uses results (true pos, false pos, true neg, false neg) to
% compute Sensitivity and Selectivity metrics and uses them to compute a
% Receiver Operating Characteristic (ROC) Curve at each time point of the
% simulation.  This yields a curve of AUC vs. t for each parameter ensemble 
% telling how accurate a segmentation from the current phase clustering of 
% oscillators would be.  From this curve, we extract two bits of
% information: (1). The maximum AUC and (2). The time it takes for the
% network to reach 95% of that maximum value.  These two values tell us
% about the quality and speed of segmentation.  We save these results along
% with the parameters for each run that produced them. On top of this we
% save the time course of AUC measure and chng measure (Hamming distance
% between segmentation at previous time using threshold and current
% segmentation).

dirPre = onCluster;

ROCmovFlg = 0;        % Flag to plot a movie of the ROC curve as it changes with time.
AUCvT_singParams = 0; % Flag to plot images of AUC vs T for single parameters settings.


grey = [0.6 0.6 0.6]; % 3x1 color vector - grey line for plotting.
colors = 'rgbkcmy';

% Directory to save output plots & data into...
imgsDir = [dirPre,'output/Kuramoto/HandCookedNetwork/imgs/SegStatsParamSearch/'];
if ~exist(imgsDir,'dir') 
    mkdir(imgsDir)
end
%
dataDistDir = [dirPre,'output/Kuramoto/HandCookedNetwork/data/ROC_KurEig_Distil/'];
if ~exist(dataDistDir,'dir') 
    mkdir(dataDistDir)
end

% Directory to save Single Parameter AUC vs Time Plots
if ( AUCvT_singParams && ~exist([imgsDir,'AUCvT_singParams/'],'dir') )
    mkdir([imgsDir,'AUCvT_singParams/'])
end

% Directory to look for data output from Kuramoto main to analyze here
dataKurDir = [dirPre,'output/Kuramoto/HandCookedNetwork/data/SegResKur/'];
%
dataEigDir = [dirPre,'output/Kuramoto/HandCookedNetwork/data/SegResEig/'];


BB=0; % a counter to keep track of number of files.

NsubsetSpec = 'N1x48'; % to look at only certain files or directories
numRunsSubsetSpec = 'runs100.'; % to look at only files with certain number of Kuramoto simulation runs

dirsK = dir([dataKurDir,'*',NsubsetSpec,'*']);

% This vecOfDirsK is input into function - a subset of directories to run thru.  If not, run through them all.
if ~exist('vecOfDirsK','var')
    vecOfDirsK = 1:numel(dirsK);
end


% If the data file already exists.  Dont rerun.
if exist([dataDistDir,'KurEig_ROC_Distil_',num2str(vecOfDirsK(1)),'_',num2str(vecOfDirsK(end)),'.mat'],'file')
    disp(['Data File - KurEig_ROC_Distil_',num2str(vecOfDirsK(1)),'_',num2str(vecOfDirsK(end)),' - already exists.  Next...'])
    return
end


% (1). Loop through directories containing Kuramoto segmentation results for different (N,C,Rmax,Pfar) network parameters
for hh = vecOfDirsK
    
    
    % (2). find directory containing Eigenvector segmentation results with matching network parameter settings (N,C,Rmax,Pfar)
    EigSubsetSpec = dirsK(hh).name (strfind(dirsK(hh).name,NsubsetSpec) : end );
    dirE = [dataEigDir,EigSubsetSpec];
    filesE = dir(dirE);   
    
    % (3).  Find (Win,Wout) connectivity strengths for given network parameter settings (N,C,Rmax,Pfar)
    for ii = 1:numel(filesE)
        
        if ~isempty(strfind(filesE(ii).name,'SegResEig_'))

            
            % (4).  For given (Win,Wout) connectivity strengths, find different settings of sigW in Kuramoto files that match
            filesK = dir([dataKurDir,dirsK(hh).name,'/SegResKur*',filesE(ii).name( numel('SegResEig_') : end-4 ),'*',numRunsSubsetSpec,'*']); % num runs here!
                        
        else
            filesK = [];
            continue
        end
    
        
        
        % Load file saved from main_Kuramoto run with given parameter settings
        SRE = load([dirE,'/',filesE(ii).name]);
        
        if(0)
            params = SRE.runParams;
            flags = SRE.runflags;
            params.thresh_seg = 100;
            flags.illROC = 1; 
            flags.saveSegRes = 0;
            [Vdom_pair_EB, K_EB, segres_EB] = main_EvecSeg(params, flags);  % not using any output.  only for ill_ROC plots.
        end
        
        
        eWin = SRE.runParams.Strng;
        eWout = SRE.runParams.Weak;
        eRmax = SRE.runParams.rmax;
        ePfar = SRE.runParams.PconnFar;
        eN = SRE.runParams.N;
        eC = SRE.runParams.C; 
        networkParams = {['N=',num2str(eN)],['C=',num2str(eC)],['Win=',num2str(eWin)],['Wout=',num2str(eWout)],['Rmax=',num2str(eRmax)],['Pfar=',num2str(ePfar)]};
        
        % Do I need to regenerate K?  It is just for plotting.  It will not
        % be exact same K as what generated results in main_Kuramoto &
        % main_EigSeg if pfar is nonzero.  Best to generate K & K_true and
        % to use same K for both main_Kuramoto & main_EigSeg.
        [K] = build_K2(SRE.runParams,SRE.runflags);
        
        %
        tp = SRE.segres.tp; % true positive: 
        fp = SRE.segres.fp; % false positive:
        tn = SRE.segres.tn; % true negative:
        fn = SRE.segres.fn; % false negative:
        %
        sens = tp ./ (tp+fn); % 
        spec = tn ./ (fp+tn); % 
        
        thresh_seg = SRE.runParams.thresh_seg;
        for b = 2:thresh_seg
            AUCe(b) = abs( trapz(1-spec(1:b),sens(1:b)) );
        end
        
        % "Robustness" is characterized as the area underneath the accumulated 
        % ROC AUC curve plotted vs threshold.  Small value is more robust.
        robustoE = abs( trapz(1:thresh_seg,AUCe./max(AUCe) ) )./thresh_seg;
        
        % "Accuracy" measure: Not used currently because is heavily biased
        % by skewed distributions (unequal values of positives & negatives)
        accE = (tp + tn) ./ (tp + fp + tn + fn);
        
        
        % NOTE:  Not sure I will need to calculate AUCe & Robusto etc here
        % because I am now doing it inside main_EvecSeg code.  I should
        % redo it here and confirm that the results are the same or
        % comparable.
        
        
        % (4).  Find (Win,Wout,SigW) for given network parameter settings (N,C,Rmax,Pfar)
        for jj = 1:numel(filesK)
            
            
            
            % Load file saved from main_Kuramoto run with given parameter settings
            SRK(jj) = load([dataKurDir,dirsK(hh).name,'/',filesK(jj).name]);
            
            if (SRK(jj).runParams.runs~=100)
               disp(['Skip. Only ',num2str(SRK(jj).runParams.runs),' runs.'])
               continue
            end
            
            BB=BB+1; % increment file number
            
            kWin = SRK(jj).runParams.Strng;
            kWout = SRK(jj).runParams.Weak;
            kRmax = SRK(jj).runParams.rmax;
            kPfar = SRK(jj).runParams.PconnFar;
            kN = SRK(jj).runParams.N;
            kC = SRK(jj).runParams.C;
            kMuW = SRK(jj).runParams.muW;
            kSigW = SRK(jj).runParams.sigW;
            
            [kMuW, kSigW];
            
            % an error check.
            if ~( kWin==eWin & kWout==eWout & kRmax==eRmax & kPfar==ePfar & kN==eN & kC==eC )
            
                disp('Something is wrong with Kuramoto & Eigenvector file matching!')
                %
                disp('Eigen')
                [eWin, eWout, eRmax, ePfar, eN, eC]
                %
                disp(['Kuramoto ',num2str(jj)])
                [kWin, kWout, kRmax, kPfar, kN, kC]
                
                keyboard
            end
            
            kurParams = ['N=',num2str(kN),' C=',num2str(kC),'  --  Win=',num2str(kWin),' Wout=',num2str(kWout),...
                ' --  Rmax=',num2str(kRmax),' Pfar=',num2str(kPfar),' --  MuW=',num2str(kMuW),' SigW=',num2str(kSigW)];
        
            
            
            %% Rerun Kuramoto for only single run so I can have pairwise cosine phase distance!
            if(0)
                params = SRK(jj).runParams;
                params.runs = 1;
                params.thresh_seg=100;
                flags = SRK(jj).runflags;
                flags.illROC = 0; 
                flags.saveSegRes = 0;
                [Rpair3_acc_KB,K_KB,segres_KB] = main_Kuramoto(params,flags); % not using any output.  only for ill_ROC plots.
                
                NoSelfLoopsMask = ones(kN) - eye(kN);
                
                tp = segres_KB.tp; % true positive: 
                fp = segres_KB.fp; % false positive:
                tn = segres_KB.tn; % true negative:
                fn = segres_KB.fn; % false negative:
                chng = segres_KB.chng;
                chng(:,1,:) = 0;
                
                sens = tp ./ (tp+fn); % 
                spec = tn ./ (fp+tn); % 
                
                for t = 1:size(sens,2)
                    spec_t = spec(:,t);
                    sens_t = sens(:,t);
                    %
                    % AUCk(t) = abs( trapz(1-spec_t,sens_t) );
                    %
                    for b = 2:params.thresh_seg
                        AUCk(t,b) = abs( trapz(1-spec_t(1:b),sens_t(1:b)) );
                    end
                    
                    robustoK(t) = abs( trapz(1:params.thresh_seg,AUCk(t,:)./max(AUCk(t,:)) ) )./params.thresh_seg;
                    
                end
        
        
                % "Accuracy" measure
                accK = (tp + tn) ./ (tp + fp + tn + fn);
                
                
                
                
                % Plot Eigenvector similarity & Kuramoto Pairwise Cosine Distance with ROC, AUC, Accuracy, etc.
                if(0)
                    
                    % Want to find peaks in AUC curve and plot pairwise cosine phase distance and ROC curve for those times
                    x = AUCk(:,end);
                    y = find(diff(x)<=0);
                    z = [1; find(diff(y)>1)+1];
                    t_vec = [1; y(z)];
                    
                    % How to make sure we have grabbed the right peaks and
                    % that there are the right amount of them?
                    if (numel(t_vec)>6)
                        
                    else
                        
                    end
                    
                    keyboard
                    
                    figure,
                    subplot(3,6,1), imagesc(K), set(gca,'XTick',[],'YTick',[])
                    title('Coupling Matrix (K)')
                    subplot(3,6,2), imagesc(SRE.K_true), set(gca,'XTick',[],'YTick',[])
                    title('Ground Truth (K true)')
                    subplot(3,6,3), imagesc(max(SRE.Vdom_pair(:)) - SRE.Vdom_pair), set(gca,'XTick',[],'YTick',[])
                    title('Eigenvalue Pairwise Similarity')
                    subplot(6,6,4:6),
                    hold on
                    plot(AUCk(:,end),'r','LineWidth',2)
                    scatter(t_vec,AUCk(t_vec,end),'bo','LineWidth',2)
                    title('AUC vs. T')
                    %
                    subplot(3,6,7), imagesc(squeeze(Rpair3_acc_KB(t_vec(1),:,:)).*NoSelfLoopsMask), caxis([-1 1]), set(gca,'XTick',[],'YTick',[])
                    ylabel('Pairwise Cosine Dist')
                    subplot(3,6,8), imagesc(squeeze(Rpair3_acc_KB(t_vec(2),:,:)).*NoSelfLoopsMask), caxis([-1 1]), set(gca,'XTick',[],'YTick',[])
                    subplot(3,6,9), imagesc(squeeze(Rpair3_acc_KB(t_vec(3),:,:)).*NoSelfLoopsMask), caxis([-1 1]), set(gca,'XTick',[],'YTick',[])
                    subplot(3,6,10), imagesc(squeeze(Rpair3_acc_KB(t_vec(4),:,:)).*NoSelfLoopsMask), caxis([-1 1]), set(gca,'XTick',[],'YTick',[])
                    subplot(3,6,11), imagesc(squeeze(Rpair3_acc_KB(t_vec(5),:,:)).*NoSelfLoopsMask), caxis([-1 1]), set(gca,'XTick',[],'YTick',[])
                    subplot(3,6,12), imagesc(squeeze(Rpair3_acc_KB(t_vec(6),:,:)).*NoSelfLoopsMask), caxis([-1 1]),  set(gca,'XTick',[],'YTick',[])
                    %
                    subplot(3,6,13), plot(1-spec(:,t_vec(1)), sens(:,t_vec(1)), 'bx'), set(gca,'XTick',[],'YTick',[])
                    ylabel('ROC Curve')
                    xlabel(['t=',num2str(t_vec(1))])
                    subplot(3,6,14), plot(1-spec(:,t_vec(2)), sens(:,t_vec(2)), 'bx'), set(gca,'XTick',[],'YTick',[])
                    xlabel(['t=',num2str(t_vec(2))])
                    subplot(3,6,15), plot(1-spec(:,t_vec(3)), sens(:,t_vec(3)), 'bx'), set(gca,'XTick',[],'YTick',[])
                    xlabel(['t=',num2str(t_vec(3))])
                    subplot(3,6,16), plot(1-spec(:,t_vec(4)), sens(:,t_vec(4)), 'bx'), set(gca,'XTick',[],'YTick',[])
                    xlabel(['t=',num2str(t_vec(4))])
                    subplot(3,6,17), plot(1-spec(:,t_vec(5)), sens(:,t_vec(5)), 'bx'), set(gca,'XTick',[],'YTick',[])
                    xlabel(['t=',num2str(t_vec(5))])
                    subplot(3,6,18), plot(1-spec(:,t_vec(6)), sens(:,t_vec(6)), 'bx'), set(gca,'XTick',[],'YTick',[])
                    xlabel(['t=',num2str(t_vec(6))])
                end
                
            
                
            end

           
            
            %
            tp = SRK(jj).segres.tp; % true positive: 
            fp = SRK(jj).segres.fp; % false positive:
            tn = SRK(jj).segres.tn; % true negative:
            fn = SRK(jj).segres.fn; % false negative:
            %
            sens = tp ./ (tp+fn); % 
            spec = tn ./ (fp+tn); % 
            
            % size of all these is = (thresholds) x (time points) x (runs)
            
            runs = SRK(jj).runParams.runs;
            tpts = SRK(jj).runParams.T ./ SRK(jj).runParams.spp;
            thrs = SRK(jj).runParams.thresh_seg;
            
            
            
            


            % calculating AUC of ROC individually on each run to get error bars
            for r = 1:runs
                
                for j = 1:tpts

                    AUCr{jj}(1,j,r) = 0;
                    
                    for b = 2:thrs
                        
                        % Accumulated Area under ROC Curve used to calculate performance
                        % and robustness of threshold results.
                        AUCr{jj}(b,j,r) = abs( trapz( 1-spec(1:b,j,r), sens(1:b,j,r) ) );
                 
                    end % loop over thresholds used for ROC analysis
                    
                    % "Robustness" is characterized as the area underneath the accumulated 
                    % ROC AUC curve plotted vs threshold.  Small value is more robust.
                    robustoK_r(j,r,jj) = abs( trapz(1:thrs, AUCr{jj}(:,j,r)./max(AUCr{jj}(:,j,r)) ) )./thrs;
                    % jj = # of Kur files (different sigma NF) 
                    
                end % loop over timepoints in Kuramoto simulation
                
               
                % Characterize the volatility of individual runs.  How
                % fluctuating is AUC with time? This varies with sigW.
                AUCrun = AUCr{jj}(end,1:end-1,r);
                
                
                % Fritz says we dont care other attracting basins - too complicated.  
                % We just care about when it is near perfect. Look when derivative / change is small to find peaks.
                maxAUC = max(AUCrun);
                nearPerf  = find( AUCrun >= 0.99.*maxAUC );
                %
                breaks = find(diff(nearPerf)>1);
                breaks = [0, breaks, numel(nearPerf)];
                %
                for a = 1:numel(breaks)-1
                    epoch = nearPerf(breaks(a)+1:breaks(a+1));
                    
                    OnDur(a) = numel(epoch);
                    
                    if(a==1)
                        OffDur(a) = epoch(1)-1;
                    else
                        OffDur(a) = nearPerf(breaks(a)+1) - nearPerf(breaks(a));
                    end
                    
                    FirstOnT(a) = epoch(1);
                    robusto(a) =  mean( robustoK_r(epoch,r,jj) );
                    
                end
                
                AUCrunStats.maxAUC(r) = maxAUC;
                AUCrunStats.numVisits(r) = numel(breaks)-1;
                AUCrunStats.meanOnDuration(r) = mean(OnDur);
                AUCrunStats.meanOffDuration(r) = mean(OffDur);
                AUCrunStats.pctTimeOn(r) = numel(nearPerf) ./ numel(AUCrun);
                AUCrunStats.FirstTimeOn(r) = FirstOnT(1);
                AUCrunStats.robustoDuringVisit(r) = mean(robusto);
                
                
                % Plot AUC value vs. derivative of AUC value
                if(0)
                    
                    figure, subplot(121), hold on,
                    plot(AUCrun(1:end-1), diff(AUCrun))
                    scatter(AUCrun(1:end-1), diff(AUCrun),'rx')
                    xlabel('AUC val')
                    ylabel('change')

                    subplot(122), hold on,
                    plot(AUCrun,'k','LineWidth',2)
                    scatter(nearPerf, AUCrun(nearPerf), 'gs', 'LineWidth',2)

                    keyboard
                
                end
                

                
            end % loop over runs of Kuramoto simulation
            
            
            % Plot Histograms of different statistics gathered over the many runs
            if(0)
                figure,
                subplot(511), hist(AUCrunStats.maxAUC)
                title([kurParams,' : Statistics across 100 runs'])
                ylabel('max AUC')
                subplot(512), hist(AUCrunStats.pctTimeOn)
                ylabel('% time on')
                subplot(513), hist(AUCrunStats.FirstTimeOn)
                ylabel('1st onset time')
                subplot(514), hist(AUCrunStats.meanOnDuration)
                ylabel('Duration On')
                subplot(515), hist(AUCrunStats.meanOffDuration)
                ylabel('Duration Off')
            end
            
            % Plot average AUC & Robusto measures along with examples from 7 random trials. 
            % (Illustrates that the average AUC value does not reflect single run AUC activity)
            if(0)

                pickRands = ceil(runs*rand(1,6));
                AUC = squeeze( AUCr{jj}(end,1:end-1,:) );
                
                % (1). Plot AUC & Robusto vs time both averaged across runs and for randomly chosen individual runs.
                figure, 
                subplot(2,5,1:4),hold on 
                shadedErrorBar( 1:tpts-1, mean(AUC,2), std(AUC,[],2) ); % mean and errorbar of AUC
                for i=1:6
                    plot(AUCr{jj}(end,1:end-1,pickRands(i)),colors(i),'LineWidth',2) % AUC of randomly chosen runs
                end
                plot(AUCk(:,end),colors(7),'LineWidth',2)   % AUC from our just generated sample run with main_Kuramoto
                scatter(tpts,AUCe(end),300,'k*','LineWidth',2)
                text(tpts, AUCe(end)-0.1, {'Eigen','Soln'},'FontSize',16,'FontWeight','Bold')
                axis([1 tpts 0.49 1.01])
                set(gca,'FontSize',16,'FontWeight','Bold')
                title(['6 Random Runs : ',kurParams],'FontSize',18,'FontWeight','Bold')
                ylabel(['Area Under ROC Curve (AUC)'],'FontSize',18,'FontWeight','Bold')
                %


                subplot(2,5,6:9),hold on 
                shadedErrorBar( 1:tpts, mean(robustoK_r(:,:,jj),2), std(robustoK_r(:,:,jj),[],2) ); % mean and errorbar of Robusto
                for i=1:6
                    plot(robustoK_r(1:end-1,pickRands(i),jj),colors(i),'LineWidth',2) % Robusto of randomly chosen runs
                end
                plot(robustoK,colors(7),'LineWidth',2)   % Robusto from our just generated sample run with main_Kuramoto
                scatter(tpts,robustoE,300,'k*','LineWidth',2)
                text(tpts, robustoE-0.1, {'Eigen','Soln'},'FontSize',16,'FontWeight','Bold')
                axis([1 tpts 0 1])
                set(gca,'FontSize',16,'FontWeight','Bold')
                xlabel(['Time (Period of 60Hz cycle)'],'FontSize',18,'FontWeight','Bold')
                ylabel({'Robustness','(Area under accumulated AUC)'},'FontSize',18,'FontWeight','Bold')

                %
                subplot(3,5,5), imagesc(SRE.K_true), colorbar
                set(gca,'XTick',[],'YTick',[])
                axis square
                title('Ground Truth','FontSize',16,'FontWeight','Bold')

                subplot(3,5,10), imagesc(K), colorbar
                set(gca,'XTick',[],'YTick',[])
                axis square
                title('Coupling Matrix (K)','FontSize',16,'FontWeight','Bold')

                subplot(3,5,15), imagesc( max(Vdom_pair_EB(:)) - Vdom_pair_EB ), colorbar
                set(gca,'XTick',[],'YTick',[])
                axis square
                title('Eigenvector Similarity','FontSize',16,'FontWeight','Bold')


                % (2). Plot in 2D space AUC vs Robusto trajectory for some randomly chosen runs and for mean run.
                figure, hold on,
                for i=1:6
                    plot(AUCr{jj}(end,1:end-1,pickRands(i)), robustoK_r(1:end-1,pickRands(i),jj),colors(i),'LineWidth',2) % AUC &Robusto of randomly chosen runs
                    scatter(AUCr{jj}(end,end-1,pickRands(i)), robustoK_r(end-1,pickRands(i),jj),100,colors(i),'LineWidth',2)
                end
                plot(AUCk(:,end), robustoK ,colors(7),'LineWidth',2) % AUC & Robusto from recently generated data
                scatter(AUCk(end,end), robustoK(end),100,colors(i),'LineWidth',2) % end point of recently generated data
                %
                
                plot(mean(AUC,2), mean(robustoK_r(1:end-1,:,jj),2) ,'k--','LineWidth',1.5)
                scatter(mean(AUC(end,:)), mean(robustoK_r(end-1,:,jj),2) ,'ks','LineWidth',1.5)
                %
                scatter(AUCe(end),robustoE,300,'rx','LineWidth',2)
                scatter(AUCe(end),robustoE,300,'rs','LineWidth',2)
                axis([0.49 1.01 0 1])
                xlabel('ROC AUC','FontSize',18,'FontWeight','Bold')
                ylabel('Robusto','FontSize',18,'FontWeight','Bold')
                title(['6 Random Runs : ',kurParams],'FontSize',20,'FontWeight','Bold')
                
                keyboard
                % Should also plot mean AUC / robusto trajectory
                
            end


            % Load distilled performance statistics into a vector that can be indexed into later.
            %
            Distil.params.Win(BB) = kWin;
            Distil.params.Wout(BB) = kWout;
            Distil.params.Rmax(BB) = kRmax;
            Distil.params.Pfar(BB) = kPfar;
            Distil.params.N(BB) = kN;
            Distil.params.C(BB) = kC;
            Distil.params.MuW(BB) = kMuW;
            Distil.params.SigW(BB) = kSigW;
            %
            Distil.params.runParams{BB} = SRK(jj).runParams;
            Distil.params.runflags{BB} = SRK(jj).runflags;
            %
            %
            Distil.Eigen.AUCe(BB) = AUCe(end);
            Distil.acrossRuns.AUCk_mn(:,BB) = mean(squeeze(AUCr{jj}(end,:,:)),2);
            Distil.acrossRuns.AUCk_std(:,BB) = std(squeeze(AUCr{jj}(end,:,:)),[],2);
            %
            % 
            Distil.Eigen.robustoE(BB) = robustoE;
            Distil.acrossRuns.robustoK_mn(:,BB) = mean(robustoK_r(:,:,jj),2);
            Distil.acrossRuns.robustoK_std(:,BB) = std(robustoK_r(:,:,jj),[],2);
            %
            %
            Distil.singleRuns.maxAUC_mn(BB) = mean(AUCrunStats.maxAUC);
            Distil.singleRuns.pctTimeOn_mn(BB) = mean(AUCrunStats.pctTimeOn);
            Distil.singleRuns.FirstTimeOn_mn(BB) = mean(AUCrunStats.FirstTimeOn);
            Distil.singleRuns.meanOnDuration(BB) = mean(AUCrunStats.meanOnDuration);
            Distil.singleRuns.meanOffDuration(BB) = mean(AUCrunStats.meanOffDuration);
                

            % Display Progresss through for loops.
            disp([num2str(hh),'/',num2str(numel(dirsK)),' Directories containing N,C,Rmax,Pfar'])
            disp([num2str(ii),'/',num2str(numel(filesE)),' Eigen-Files additionally containing Win, Wout'])
            disp([num2str(jj),'/',num2str(numel(filesK)),' Kuramoto Files additionally containing SigW'])
            disp(['------------------------------------------------------------------------------------'])

        end % loop over saved Kuramoto data files (jj) 
        
        
        
        
        
        if(0)
            figure, 
            subplot(121), hold on,
            scatter(Distil.params.SigW,Distil.singleRuns.meanOnDuration,'bo')
            scatter(Distil.params.SigW,Distil.singleRuns.meanOffDuration,'rx')
            title('Duration on (o) Depends on sigW. Duration off (x) Doesnt.')
            xlabel('SigW')
            ylabel('Duration (# periods of 60Hz)')
            subplot(122), hold on,
            scatter(Distil.params.SigW,Distil.singleRuns.pctTimeOn_mn,'bo')
            scatter(Distil.params.SigW,Distil.singleRuns.maxAUC_mn,'rx')
            title('% Time on (o) Depends on sigW. Sometimes max AUC (x) does too.')
            xlabel('SigW')
            ylabel('Percent of Time or Quality')
        end
        

    
    end % Loop over Eigen-files (ii)
    
    
    
    
end % look over Kuramoto Directories (hh)
    
    
    
    
    
    
    






%% save an output mat file.
save([dataDistDir,'KurEig_ROC_Distil_',num2str(vecOfDirsK(1)),'_',num2str(vecOfDirsK(end))],'Distil'); 







    
    
%     
%     
%     
% 
%     % 2nd Loop through output files saved from main_Kuramoto code and plot statistics
%     for i = 1:numel(files)
% 
% 
%         % Load data file saved from main_Kuramoto_HandCookedNetwork
%         if strfind(files(i).name,'.mat')
%             load([dataKurDir,'/',dirsK(hh).name,'/',files(i).name]);
%             k=k+1
%         else
%             continue
%         end
% 
% 
%         %% 
%         tp = segres.tp; % true positive: 
%         fp = segres.fp; % false positive:
%         tn = segres.tn; % true negative:
%         fn = segres.fn; % false negative:
%         chng = segres.chng;
%         chng(:,1,:) = 0;
% 
% 
% %         % Get rid of last time point (something weird.  no big deal)
% %         tp = tp(:,(1:end-1),:);
% %         fp = fp(:,(1:end-1),:);
% %         tn = tn(:,(1:end-1),:);
% %         fn = fn(:,(1:end-1),:);
% %         chng = chng(:,(1:end-1),:);
% 
% 
% 
% %         if ~isfield(runParams,'PconnFar')
% %             runParams.PconnFar = 0; % if not other wise stated probability of a connection out beyond Rmax is zero.
% %         end
% 
%         %
%         thresh_cosDist_seg = linspace(-1,1,runParams.thresh_seg);
%         for j = 1:runParams.thresh_seg
%             threshCell{j} = num2str(thresh_cosDist_seg(j),2);
%         end
%         
%         
% 
% 
%         %
%         runParamsTag = {['PIF ',runParams.PIFlg,' N',num2str(runParams.Ndims(1)),'x',num2str(runParams.Ndims(2)),' - C',num2str(runParams.C)],...
%           [' - rmax',num2str(runParams.rmax),' - NF ',num2str(runParams.muW),' ',num2str(runParams.sigW)],[' - Win ',num2str(runParams.Strng),' ',...
%                     num2str(runParams.sigSg),' - Wout ',num2str(runParams.Weak),' ',num2str(runParams.sigWk)]};
% 
% 
%         %% ROC Curve using Sensitivity & Specificity !!
% 
%         sens = tp ./ (tp+fn); % 
%         spec = tn ./ (fp+tn); % 
% 
% 
%         % calculating AUC of ROC individually on each run to get error bars
%         for r = 1:runParams.runs
% 
%             xr = sens(:,:,r); 
%             yr = spec(:,:,r);
% 
%             if(ROCmovFlg)
%                 ROCwriterObj = VideoWriter([imgsDir,'ROC_time_lapse_',num2str(r),'.avi']);
%                 ROCwriterObj.FrameRate = 2;
%                 open(ROCwriterObj);
%                 h=figure;
%             end
%             %
%             for j = 1:size(xr,2)
%                 xpr = [1, xr(:,j)', 0];
%                 ypr = 1 - [0, yr(:,j)', 1];
%                 AUCr(j,r) = abs( trapz(ypr,xpr) );
% 
%                 % to plot ROC curve that changes with time.
%                 if(ROCmovFlg) % can make a movie out of it.
%                     
%                     figure(h), subplot(3,1,1:2), hold on
%                     plot(ypr,xpr,'b','LineWidth',2) 
%                     plot([0 1],[0 1],'r--','LineWidth',2)
%                     text(0.02, 1, 'Best','Color','red','VerticalAlignment','Top','FontSize',16,'FontWeight','Bold')
%                     text(0.5, 0.5, 'Uninformative','Color','red','VerticalAlignment','Top','FontSize',16,'FontWeight','Bold')
%                     text(0.6, 0.2, runParamsTag,'VerticalAlignment','Top','FontSize',16,'FontWeight','Bold')
%                     title(['Run ',num2str(r),' Time ',num2str(j),' ROC Curve'],'FontSize',18,'FontWeight','Bold')
%                     ylabel('Sensitivity','FontSize',18,'FontWeight','Bold'), 
%                     xlabel('1 - Specificity','FontSize',18,'FontWeight','Bold')
%                     set(gca,'FontSize',18,'FontWeight','Bold')
%                     set(h,'Units','Normalized','Position',[0.1 0.1 0.5 0.9])
%                     %
%                     subplot(3,1,3), hold on, cla
%                     plot([0.5 size(xr,2)], [0.5 0.5], 'r--')
%                     text(size(xr,2)./2, 0.5, 'Uninformative','Color','red','VerticalAlignment','Bottom','FontSize',16,'FontWeight','Bold')
%                     text(size(xr,2)./2, 1, 'Best','Color','red','VerticalAlignment','Top','FontSize',16,'FontWeight','Bold')
%                     plot(AUCr(1:j,r),'b','LineWidth',2), 
%                     xlabel('Time - Period of 60Hz osc.','FontSize',18,'FontWeight','Bold')
%                     ylabel('Area Under ROC Curve','FontSize',18,'FontWeight','Bold')
%                     axis([0 size(xr,2) 0.4 1])
%                     set(gca,'FontSize',16,'FontWeight','Bold')
%                     %
%                     mov = getframe(h);
%                     writeVideo(ROCwriterObj,mov);
%                     %
%                     subplot(3,1,1:2), plot(ypr,xpr,'Color',grey,'LineWidth',2)
%                     
%                 end
%             end
%             
%             
%             if(ROCmovFlg)
%                 close(ROCwriterObj);
%             end
%             
% %             keyboard
%             
%             
%         end
%         
%         if(1)
% %             h=figure; hold on, shadedErrorBar(1:size(AUCr,1), mean(AUCr'), std(AUCr'))
%         end
%         
%         
% %         % HOW IT WAS.  THIS IS FINE IF I DONT WANT ERRORBARS        
% %         x = mean(sens,3); 
% %         y = mean(spec,3);
% % 
% %         %figure,
% %         for j = 1:size(x,2)
% %             xp = [1, x(:,j)', 0];
% %             yp = 1 - [0, y(:,j)', 1];
% %             AUC(j) = abs( trapz(yp,xp) );
% %             
% %             % to plot ROC curve that changes with time.
% %             if(0) % can make a movie out of it.
% %                 plot(yp,xp), title(num2str(j))
% %                 pause
% %             end
% %             
% %         end
% 
% 
%         % Parameters attached to performance (may need to grab more).
%         %
%         ROCparams.Strng(k) = runParams.Strng;
%         ROCparams.Weak(k) = runParams.Weak;
%         ROCparams.C(k) = runParams.C;
%         ROCparams.sigW(k) = runParams.sigW;
%         ROCparams.rmax(k) = runParams.rmax;
%         ROCparams.pfar(k) = runParams.PconnFar;
%         
%         AUC = mean(AUCr');
% 
%         AUCm2{k} = mean(AUCr'); % Curve of AUC vs T averaged over simulations for each parameter setting
%         AUCs2{k} = std(AUCr'); % Curve of AUC vs T averaged over simulations for each parameter setting
%         chng2{k} = mean(mean(chng,3)); % avg first across runs, then across thresholds
%         %
%         best(k) = max(AUC);
%         bestestT = find(AUC == max(AUC) );
%         bestT(k) = bestestT(1);
%         %
%         fastest = find(AUC > 0.98*max(AUC));
%         try
%             fast(k) = fastest(1);
%         catch
%             fast(k) = find(AUC == max(AUC));
%         end
%         %
%         stab(k) = chng2{k}(fast(k)-1);
%         %
%         %
%          dif = diff(AUCm2{k});
%          cut = find(dif < 0);
%          try
%             tFlat(k) = cut(1);
%          catch
%             tFlat(k) = numel(dif);
%          end
%          qFlat(k) = AUCm2{k}(tFlat(k));
% 
% 
% 
%         
%         
%         
%         
% 
%         % plot for single parameter setting, time course of AUC
%         if(AUCvT_singParams)
%             h=figure; hold on, 
%             shadedErrorBar(1:size(AUCr,1), mean(AUCr'), std(AUCr'))
%             plot(AUC,'k','LineWidth',2),
%             scatter(fast(k),AUC(fast(k)),100,'gx','LineWidth',2)
%             text(fast(k),0.98*AUC(fast(k)),'t_{98%}','Color','g','FontSize',16,'FontWeight','Bold','VerticalAlignment','Top')
%             scatter(bestT(k),best(k),100,'ro','LineWidth',2)
%             text(bestT(k),1.02*best(k),'Best','Color','r','FontSize',16,'FontWeight','Bold','VerticalAlignment','Bottom')
%             title([runParamsTag],'FontSize',20,'FontWeight','Bold') % ['File# ',num2str(i)],
%             xlabel('Time - period of 60Hz Osc (60=1sec)','FontSize',18,'FontWeight','Bold')
%             ylabel('Area under ROC Curve (100 runs)','FontSize',18,'FontWeight','Bold')
%             axis([1 numel(AUC) 0.5 1])
%             set(gca,'FontSize',16,'FontWeight','Bold')
%             grid on
%             %
%             % Save image
%             rpt = ['PIF_',runParams.PIFlg,'_N',num2str(runParams.Ndims(1)),'x',num2str(runParams.Ndims(2)),'_C',...
%                 num2str(runParams.C),'_rmax',num2str(runParams.rmax),'_NF_',num2str(runParams.muW),'_',...
%                 num2str(runParams.sigW),'_Win ',num2str(runParams.Strng),'_',num2str(runParams.sigSg),...
%                 '_Wout ',num2str(runParams.Weak),'_',num2str(runParams.sigWk)];
%             saveGoodImg(h,[imgsDir,'AUCvT_singParams/',rpt],[0 0 1 1])
%             close(h) 
%         end
%         
% 
% %         keyboard
%         clear AUC % get rid of this too.
%         
%         disp([num2str(hh),'/',num2str(numel(dirsK))])
%         disp([num2str(i),'/',num2str(numel(files))])
%         
%         
% 
%     end % Loop over files in a given directory
%     
% 
% %end % Loop over directories in data dir
% 
% 
% 
% 
% %% save an output mat file.
% save([dataKurDir,'AUC_ROC_stats'],'ROCparams','runParams','best','bestT','tFlat','qFlat','fast','stab','AUCm2','AUCs2','chng2','subsetSpec'); 
