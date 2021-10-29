% This script will loop through the data .mat files that were saved from
% single parameter settings (multiple simulation runs) using
% loop_main_Kuramoto_HandCookedNetwork.m.  
%
% The variables saved from that code are 'MC', 'metaCluster', 'runParams',
% 'runflags'.  From MC, I compute Separation clustering metric for each run
% and look at statistics across runs.
%
% There are options to make plot
% of dynamics on single runs.  Average values of 



% NOTE: I have update variable names in loop_main_Kuramoto_HandCookedNetwork 
% so things in here will have to be changed if I ever rerun that code to 
% produce new matfiles.
%
% 'netParams','netflags','kurParams','kurflags','metaCluster','MC'





%% Make directory to save output data files and images produced.
dirPre = onCluster;
saveDir = [dirPre,'output/Kuramoto/HandCookedNetwork/'];

SeparationInFlg = 'phase'; % 'location' or 'phase'

h5r_flg = 0; % flag to produce plots for 5 runs of oscillator phase along side Separation measure plot.

SepStructFlg = 1; % flag to create and save structure containing run parameters and Separation stats to be analyzed in another function


% Images Directory
imgsKurDir = [saveDir,'imgs/MetaCluster/'];
%
if ~exist(imgsKurDir,'dir')
    mkdir(imgsKurDir)
end
%
% Data Directory
dataKurDir = [saveDir,'data/MetaCluster/'];       
%
if ~exist(dataKurDir,'dir')
    mkdir(dataKurDir)
end


BB=0; % Counter of how many times we go through the innermost loop (to store calculated things)

for MM = [0 0.3 0.6 0.9 1.2] % Loop through SigW values

    for LL = [1, 4, inf] % Loop through Rmax value

        for KK = [2,3,4,6,8] % Loop through Cluster Sizes

            for JJ = 10 %[1:10] % Loop through Win

                for II = 1:10 %[-1:-1:-10] % Loop through Wout


                    [MM,LL,KK,JJ,II]
                    

                    % Load Data file
                    sigW_str = num2str(MM);
                    sigW_str(sigW_str=='.')='p';
                    
                    try
                        load([saveDir,'data/MetaCluster/PIF_Fourier1_N1x48_C',num2str(KK),'_rmax',num2str(LL),'_pfar0/metaClusterKur_NF_60_',sigW_str,'_Win_',num2str(JJ),'_0_Wout_',num2str(II),'_0_runs100.mat'])
                    catch
                    
                        disp(['File : PIF_Fourier1_N1x48_C',num2str(KK),'_rmax',num2str(LL),'_pfar0/metaClusterKur_NF_60_',sigW_str,'_Win_',num2str(JJ),'_0_Wout_',num2str(II),'_0_runs100.mat does not exist. Skipping.'])
                        continue
                    end

                    % Extract parameters used to generate the data file from runParams
                    muW        = runParams.muW;
                    sigSg      = runParams.sigSg;
                    sigWk      = runParams.sigWk;
                    runs       = runParams.runs;
                    Tsec       = runParams.Tsec;
                    spp        = runParams.spp;
                    tau        = runParams.tau;
                    T          = runParams.T;
                    PIFlg      = runParams.PIFlg;
                    PIFparams  = runParams.PIFparams;
                    numInt     = runParams.numInt;
                    thresh_seg = runParams.thresh_seg;
                    PconnFar   = runParams.PconnFar;
                    rmax       = runParams.rmax;
                    gndTruth   = runParams.gndTruth;
                    C          = runParams.C;
                    Ndims      = runParams.Ndims;
                    N          = runParams.N;
                    Weak       = runParams.Weak;
                    Strng      = runParams.Strng;
                    sigW       = runParams.sigW;

                    % Tag plots with parameters used just for easier analysis
                    paramsTit = ['PIF = ',PIFlg,' - N = ',num2str(Ndims(1)),'x',num2str(Ndims(2)),' - C = ',num2str(C),...
                        ' - rmax = ',num2str(rmax),' - pfar = ',num2str(PconnFar),' - NF = ',num2str(muW),' , ',num2str(sigW),' - Win = ',num2str(Strng),...
                        ' , ', num2str(sigSg),' - Wout = ',num2str(Weak),' , ',num2str(sigWk)];


                    % Note:  This will work ok up to 7 clusters without repeating colors
                    colors = 'kbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmy';
                    stylee = cell(1,numel(gndTruth));
                    for i = 1:numel(gndTruth)
                        stylee{i} = colors(gndTruth(i));
                    end



                    % % Is there a general trend across runs of meanCSep and stdCSep?
                    % for r = 1:runs
                    %     stdCSep(:,r)  = metaCluster(r).stdCSep;
                    %     meanCSep(:,r) = metaCluster(r).meanCSep;
                    % end


                    % Extract data from MC (metaCluster) data structure.
                    w         = MC.w;
                    meanCVar  = MC.meanCVar;
                    stdCVar   = MC.stdCVar;
                    meanCExt  = MC.meanCExt;
                    stdCExt   = MC.stdCExt;
                    %
                    switch SeparationInFlg

                        case 'location'
                            meanCSep  = MC.location.meanCSep;
                            stdCSep   = MC.location.stdCSep;
                            minCSep   = MC.location.minCSep;
                            minCSepID = MC.location.minCSepID;
                            meanCDist = MC.location.meanCDist;
                            stdCDist  = MC.location.stdCDist;

                        case 'phase'
                            meanCSep  = MC.phase.meanCSep;
                            stdCSep   = MC.phase.stdCSep;
                            minCSep   = MC.phase.minCSep;
                            minCSepID = MC.phase.minCSepID;
                            meanCDist = MC.phase.meanCDist;
                            stdCDist  = MC.phase.stdCDist;

                    end
                    
                    
                    
                    
                    % Collect mean and median (over runs) meanCSep value into a data structure along with parameter values.
                    if(SepStructFlg)
                        BB = BB +1;
                        SeparationStruct.sigW(BB) = sigW;
                        SeparationStruct.Win(BB) = Strng;
                        SeparationStruct.Wout(BB) = Weak;
                        SeparationStruct.C(BB) = C;
                        SeparationStruct.rmax(BB) = rmax;
                        %
                        SeparationStruct.meanOverRuns_CSep(BB) = mean(max(meanCSep'));
                        SeparationStruct.stdOverRuns_CSep(BB) = std(max(meanCSep'));
                        SeparationStruct.medianOverRuns_CSep(BB) = median(max(meanCSep'));
                        SeparationStruct.maxOverRuns_CSep(BB) = max(max(meanCSep'));
                    end

                    


%                     % find time when clusters have coalesced enough (to 5% of the 2pi/C space alotted for them)
%                     % NOTE: This is averaged across runs ... (maybe not bad) 
%                     try
%                         t_coalesce = find(mean(meanCVar)<0.05*2*pi/C);
%                         tc = t_coalesce(1);
%                     catch
%                         tc = T;
%                     end
% 
% 




                    % Plot mean (across runs) of CSep (mean/std across clusters of "Separation" value for cluster and its closest neighbor).
                    if(0)
                        figure
                        subplot(121), hold on, plot(meanCSep'), plot(mean(meanCSep),'k','LineWidth',2), title('mean Cluster Sep Across clusters')
                        subplot(122), hold on, plot(stdCSep'), plot(mean(stdCSep),'k','LineWidth',2), title('std Cluster Sep Across clusters')
                        % Put a Title over all subplots on the figure
                        annotation('textbox', [0 0.9 1 0.1],'String', ['Variability Across Runs: ',paramsTit], ...
                                'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')
                    end





                    % Plot mean (across runs) of CDist  & CVar.
                    if(0)
                        h=figure('units','normalized','outerposition',[0 0 1 1]);
                        subplot(221), hold on, plot(meanCDist'), plot(mean(meanCDist),'k','LineWidth',2), title('mean Cluster Dist Across clusters')
                        axis([0 T 0 2*pi/C])
                        plot([tc, tc], [0, 2*pi/C],'r--')
                        text(tc, 2*pi/C, 'tc', 'Color', 'r','HorizontalAlignment','Left','VerticalAlignment','Top')
                        subplot(222), hold on, plot(stdCDist'), plot(mean(stdCDist),'k','LineWidth',2), title('std Cluster Dist Across clusters')
                        axis([0 T 0 2*pi/C])
                        plot([tc, tc], [0, 2*pi/C],'r--')
                        text(tc, 2*pi/C, 'tc', 'Color', 'r','HorizontalAlignment','Left','VerticalAlignment','Top')
                        %
                        subplot(223), hold on, plot(meanCVar'), plot(mean(meanCVar),'k','LineWidth',2), title('mean Cluster Var Across clusters')
                        axis([0 T 0 2*pi/C])
                        plot([tc, tc], [0, 2*pi/C],'r--')
                        text(tc, 2*pi/C, {'coalesce','time'}, 'Color', 'r','HorizontalAlignment','Left','VerticalAlignment','Top')
                        xlabel('time')
                        ylabel('radians')


                        subplot(224), hold on, plot(stdCVar'), plot(mean(stdCVar),'k','LineWidth',2), title('std Cluster Var Across clusters')
                        axis([0 T 0 2*pi/C])
                        plot([tc, tc], [0, 2*pi/C],'r--')
                        text(tc, 2*pi/C, 'tc', 'Color', 'r','HorizontalAlignment','Left','VerticalAlignment','Top')
                        % Put a Title over all subplots on the figure
                        annotation('textbox', [0 0.9 1 0.1],'String', ['Variability Across Runs: ',paramsTit], ...
                                'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')
                    end


                    if(h5r_flg)
                        h5r = figure('units','normalized','outerposition',[0 0 1 1]);
                    end
                    
                    % Loop Through Runs
                    for r = 1:runs

                        % Extract "Clusterability" single-run data from metaCluster data structure
                        t_2pi      = metaCluster(r).t_2pi;
                        phaseAtClk = metaCluster(r).phaseAtClk;
                        w          = MC.w(r,:);
                        meanCVar  = MC.meanCVar(r,:);
                        stdCVar   = MC.stdCVar(r,:);
                        meanCExt  = MC.meanCExt(r,:);
                        stdCExt   = MC.stdCExt(r,:);
                        %
                        switch SeparationInFlg

                            case 'location'
                                meanCSep  = MC.location.meanCSep;
                                stdCSep   = MC.location.stdCSep;
                                minCSep   = MC.location.minCSep;
                                minCSepID = MC.location.minCSepID;
                                meanCDist = MC.location.meanCDist;
                                stdCDist  = MC.location.stdCDist;

                            case 'phase'
                                meanCSep  = MC.phase.meanCSep;
                                stdCSep   = MC.phase.stdCSep;
                                minCSep   = MC.phase.minCSep;
                                minCSepID = MC.phase.minCSepID;
                                meanCDist = MC.phase.meanCDist;
                                stdCDist  = MC.phase.stdCDist;

                        end


                        % find time when mean "Separation" measure of nearest in phase cluster pairs (subtractive one) goes positive.
                        % NOTE: Do not want to average this across runs...
                        try
                            t_pos_mean = find(meanCSep > 0);
                            tpo_mean = t_pos_mean(1);
                        catch
                            tpo_mean = T;
                        end
                        %
        %                 try
        %                     t_pos_min = find(meanCSep > 0);
        %                     tpo_min = t_pos_min(1);
        %                 catch
        %                     tpo_min = T;
        %                 end


                        if(h5r_flg)
                            figure(h5r)
                            subplot(runs,2,2*r-1), hold on
%                             subplot(211), hold on
                            for i = 1:N

                                % Dont want to have plots wrap around because its distracting
                                dude = 2*pi*mod(t_2pi{i},1/muW)/(1/muW);
                                tt = 0:numel(dude);
                                %
                                dl = diff(dude);
                                jumpind = [abs(dl)>pi]; % now if jumpind(i) = true, we know that the
                                ind = find(jumpind);
                                ind = [0, ind, numel(dude)];
                                %
                                for j = 2:numel(ind)
                                    dudesub = dude(ind(j-1)+1:ind(j));
                                    tsub = (1./w(i)).*tt(ind(j-1)+1:ind(j));
                                    plot(tsub,dudesub,stylee{i},'LineWidth',2)
                                end

                                % Annotate oscillator number and natural frequency on plot
                                text(Tsec, dude(end-2),['#',num2str(i),' - \omega=',num2str(w(i),4)],'Color',stylee{i},'Fontsize',14,'FontWeight','Bold')
                            end
                            %
                            % Mark time when mean Cluster Pair "Separation" turns positive
                            if(tpo_mean~=T)
                                plot([tpo_mean/T, tpo_mean/T], [0, 2*pi],'b--','LineWidth',2)
                                scatter(tpo_mean/T, 0, 100, 'bs','LineWidth',2)
                                text(tpo_mean/T, 0, {'posSep','time'}, 'Color', 'r','HorizontalAlignment','Left','VerticalAlignment','Bottom')
                            end
                            %
                        %     plot([tau*spp*LP2, tau*spp*LP2], [0, 2*pi],'r--')\
                            %
                            title(['Phase with Separation = ',num2str(max(meanCSep),3),' & tp = ',num2str(tpo_mean/T,3)],'Fontsize',20,'FontWeight','Bold')
                            % xlabel(['Period # of ',num2str(muW),'Hz oscillation'],'Fontsize',18,'FontWeight','Bold')
                            % ylabel(['Phase of Oscillation'],'Fontsize',18,'FontWeight','Bold')
                            set(gca,'Fontsize',16,'FontWeight','Bold','Ytick',[0 pi/2 pi 3*pi/2 2*pi],'YTickLabel',{'0','90','180','270','360'}) % ,'Interpreter','Latex'
                            axis([0 (1./muW).*numel(dude) 0 2*pi])
                            ylabel('Phase (degrees)','Fontsize',18,'FontWeight','Bold')



                            subplot(runs,2,2*r), hold on
%                             subplot(212), hold on
                            plot([1:T]/T, meanCSep,'k','LineWidth',2)
                            plot([1:T]/T, minCSep,'b','LineWidth',2)
                            if(r==1)
                                title(paramsTit,'Fontsize',18,'FontWeight','Bold')
                            end
                            %
                            if(tpo_mean~=T)
                                % plot([tpo/T, tpo/T], [0, 2*pi],'b--')
                                scatter(tpo_mean/T, 0, 100, 'bs')
                                text(tpo_mean/T, 0, {'tp'}, 'Color', 'r','HorizontalAlignment','Left','VerticalAlignment','Bottom')
                            end
                            
                            set(gca,'Fontsize',16,'FontWeight','Bold') % ,'Interpreter','Latex'
                            xlabel('Simulation Time (seconds)','Fontsize',18,'FontWeight','Bold')
                            ylabel('Separation Metric','Fontsize',18,'FontWeight','Bold')
                            axis tight
                            
                            keyboard
                            
                        end



                        if(0)
                            h=figure; 
                            subplot(4,1,1), hold on  %  

                            for i = 1:N

                                % Dont want to have plots wrap around because its distracting
                                dude = 2*pi*mod(t_2pi{i},1/muW)/(1/muW);
                                tt = 0:numel(dude);
                                %
                                dl = diff(dude);
                                jumpind = [abs(dl)>pi]; % now if jumpind(i) = true, we know that the
                                ind = find(jumpind);
                                ind = [0, ind, numel(dude)];
                                %
                                for j = 2:numel(ind)
                                    dudesub = dude(ind(j-1)+1:ind(j));
                                    tsub = (1./w(i)).*tt(ind(j-1)+1:ind(j));
                                    plot(tsub,dudesub,stylee{i},'LineWidth',2)
                                end

                                % Annotate oscillator number and natural frequency on plot
                                text(Tsec, dude(end),['#',num2str(i),' - \omega=',num2str(w(i),4)],'Color',stylee{i},'Fontsize',14,'FontWeight','Bold')
                            end
                            %
                        %     plot([tau*spp*LP2, tau*spp*LP2], [0, 2*pi],'r--')\
                            %
                            title(['Phase of 2\pi Crossing w.r.t. ',num2str(muW),' Hz Oscillation'],'Fontsize',20,'FontWeight','Bold')
                            % xlabel(['Period # of ',num2str(muW),'Hz oscillation'],'Fontsize',18,'FontWeight','Bold')
                            ylabel(['Phase of Oscillation'],'Fontsize',18,'FontWeight','Bold')
                            set(gca,'Fontsize',16,'FontWeight','Bold','Ytick',[0 pi/2 pi 3*pi/2 2*pi],'YTickLabel',{'0','90','180','270','360'}) % ,'Interpreter','Latex'
                            axis([0 (1./muW).*numel(dude) 0 2*pi])




                            % Plot "Cluster Separation" - Clusterability Metric.
                            subplot(4,1,2), hold on
                            shadedErrorBar([1:T]./(runParams.spp.*runParams.muW), meanCSep, stdCSep )
                            %plot([1:T]./(runParams.spp.*runParams.muW), minCSep-0.2,'k', 'LineWidth', 5)


                            for i = 1:T
                                scatter(i./(runParams.spp.*runParams.muW), minCSep(i), 30, 'Filled', 'MarkerFaceColor',colors(minCSepID(i,1)))
                                scatter(i./(runParams.spp.*runParams.muW), minCSep(i) + 0.05*max(minCSep), 30, 'Filled', 'MarkerFaceColor',colors(minCSepID(i,2)))
                            end


                            plot([1:T]./(runParams.spp.*runParams.muW), meanCSep,'k--', 'LineWidth', 5)

                        %     % Line at Launch Point (where "Separation" measure > threshold)
                        %     plot([LP/T, LP/T], [0, max(meanCSep)+max(stdCSep)],'r--')
                        %     text(0.99*(LP/T), 0.99*(max(meanCSep)+max(stdCSep)), 'LP', 'Color', 'r','HorizontalAlignment','Right','VerticalAlignment','Top')

                            title(['Cluster Separation'],'FontSize',20,'FontWeight','Bold')
                            xlabel(['Simulation Time (sec)'],'FontSize',18,'FontWeight','Bold')
                            ylabel({['Distance Between Cluster Centers /'],['sqrt Variance of Clusters']},'FontSize',18,'FontWeight','Bold')
                            set(gca,'FontSize',16,'FontWeight','Bold')
                            axis tight



                            subplot(4,1,3), hold on
                            shadedErrorBar([1:T]./(runParams.spp.*runParams.muW), meanCDist, stdCDist )
                            axis([0 Tsec 0 2*pi/C])
                            title('Distance between each Cluster and its Nearest neighbor')
                            set(gca,'FontSize',16,'FontWeight','Bold')

                            subplot(4,1,4), hold on
                            shadedErrorBar([1:T]./(runParams.spp.*runParams.muW), meanCVar, stdCVar )
                            axis([0 Tsec 0 2*pi/C])
                            title('Variance of Single Clusters')
                            set(gca,'FontSize',16,'FontWeight','Bold')
                        end





                        %     saveGoodImg(h,[imgsKurDir,'metaClusterAnalysis_run',num2str(r)],[0 0 1 1])
                        %     close(h);








                        % Plot Trajectory in 2D Cluster Variance vs. Distance between Cluster space
                        if(0)
                            figure('units','normalized','outerposition',[0 0 1 1]) 
                            subplot(221),hold on,

                            if(tc~=T)
                                shadedErrorBar(meanCVar(tc:end),meanCDist(tc:end),stdCDist(tc:end)) 
                                plot(meanCVar(tc:end),meanCDist(tc:end),'k','LineWidth',2)
                                scatter(meanCVar(tc),meanCDist(tc),100,'bd','LineWidth',2); % mark start of trajectory

            %                         scatter(meanCVar(1),meanCDist(1),200,'k*','LineWidth',2); % mark start of trajectory
                                text(meanCVar(tc),meanCDist(tc),{'coalesce','time'},'HorizontalAlignment','Left','VerticalAlignment','Top')

            %                         plot(meanCVar+stdCVar,meanCDist,'r--')
            %                         plot(meanCVar-stdCVar,meanCDist,'r--')
                                %
                                plot(meanCVar(tc:end)+stdCVar(tc:end),meanCDist(tc:end),'r--','LineWidth',1.5)
                                plot(meanCVar(tc:end)-stdCVar(tc:end),meanCDist(tc:end),'r--','LineWidth',1.5)

            %                     plot(meanCVar,meanCDist+stdCDist,'g--')
            %                     plot(meanCVar,meanCDist-stdCDist,'g--')
                                %
                                plot(meanCVar(tc:end),meanCDist(tc:end)+stdCDist(tc:end),'g--','LineWidth',1.5)
                                plot(meanCVar(tc:end),meanCDist(tc:end)-stdCDist(tc:end),'g--','LineWidth',1.5)
                            end

                            axis([0 2*pi/C 0 2*pi/C])
                            xlabel('Average Within Cluster Variance (Radians)')
                            ylabel('Average Distance between Nearest Pairs of Clusters (Radians)')
                            title(['Trajectory of Simulation after Cluster Coalesce Time (tc = ',num2str(tc/T),')'])



                            % Plot Trajectory in 2D Cluster Extent vs. Distance between Cluster space
                            subplot(222),hold on,

                            if(tc~=T)
                                shadedErrorBar(meanCExt(tc:end),meanCDist(tc:end),stdCDist(tc:end)) 
                                plot(meanCExt(tc:end),meanCDist(tc:end),'k','LineWidth',2)
                                scatter(meanCExt(tc),meanCDist(tc),100,'bd','LineWidth',2); % mark start of trajectory

            %                         scatter(meanCExt(1),meanCDist(1),200,'k*','LineWidth',2); % mark start of trajectory
                                text(meanCExt(tc),meanCDist(tc),{'coalesce','time'},'HorizontalAlignment','Left','VerticalAlignment','Top')

            %                         plot(meanCExt+stdCExt,meanCDist,'r--')
            %                         plot(meanCExt-stdCExt,meanCDist,'r--')
                                %
                                plot(meanCExt(tc:end)+stdCExt(tc:end),meanCDist(tc:end),'r--','LineWidth',1.5)
                                plot(meanCExt(tc:end)-stdCExt(tc:end),meanCDist(tc:end),'r--','LineWidth',1.5)

            %                     plot(meanCExt,meanCDist+stdCDist,'g--')
            %                     plot(meanCExt,meanCDist-stdCDist,'g--')
                                %
                                plot(meanCExt(tc:end),meanCDist(tc:end)+stdCDist(tc:end),'g--','LineWidth',1.5)
                                plot(meanCExt(tc:end),meanCDist(tc:end)-stdCDist(tc:end),'g--','LineWidth',1.5)
                            end

                            plot([0 2*pi/C],[0 2*pi/C],'k--') % plot diagonal line

                            axis([0 2*pi/C 0 2*pi/C])
                            xlabel('Average Within Cluster End-to-End Extent (Radians)')
                            ylabel('Average Distance between Nearest Pairs of Clusters (Radians)')
                            title(['Trajectory of Simulation after Cluster Coalesce Time (tc = ',num2str(tc/T),')'])






                            % Plot demodulated oscillator phase vs. T to watch clustering happen
                            subplot(413), hold on  %  
                            for i = 1:N

                                % Dont want to have plots wrap around because its distracting
                                dude = 2*pi*mod(t_2pi{i},1/muW)/(1/muW);
                                tt = 0:numel(dude);
                                %
                                dl = diff(dude);
                                jumpind = [abs(dl)>pi]; % now if jumpind(i) = true, we know that the
                                ind = find(jumpind);
                                ind = [0, ind, numel(dude)];
                                %
                                for j = 2:numel(ind)
                                    dudesub = dude(ind(j-1)+1:ind(j));
                                    tsub = (1./w(i)).*tt(ind(j-1)+1:ind(j));
                                    plot(tsub,dudesub,stylee{i},'LineWidth',2)
                                end

                                % Annotate oscillator number and natural frequency on plot
                                text(Tsec, dude(end),['#',num2str(i),' - \omega=',num2str(w(i),4)],'Color',stylee{i},'Fontsize',14,'FontWeight','Bold')
                            end
                            %
                        %     plot([tau*spp*LP2, tau*spp*LP2], [0, 2*pi],'r--')\
                            %
                            title(['Phase of 2\pi Crossing w.r.t. ',num2str(muW),' Hz Oscillation'],'Fontsize',20,'FontWeight','Bold')
                            % xlabel(['Period # of ',num2str(muW),'Hz oscillation'],'Fontsize',18,'FontWeight','Bold')
                            ylabel(['Phase of Oscillation'],'Fontsize',18,'FontWeight','Bold')
                            set(gca,'Fontsize',16,'FontWeight','Bold','Ytick',[0 pi/2 pi 3*pi/2 2*pi],'YTickLabel',{'0','$\pi$/2','$\pi$','3$\pi$/2','2$\pi$'}) % ,'Interpreter','Latex'
                            axis([0 (1./muW).*numel(dude) 0 2*pi])
                            if(tc~=T)
                                plot([tc/T, tc/T], [0, 2*pi],'b--')
                                scatter(tc/T, 0, 100, 'bd')
                                text(tc/T, 0, {'coalesce','time'}, 'Color', 'r','HorizontalAlignment','Left','VerticalAlignment','Bottom')
                            end
                            %
                            if(tpo~=T)
                                plot([tpo/T, tpo/T], [0, 2*pi],'b--')
                                scatter(tpo/T, 0, 100, 'bs')
                                text(tpo/T, 0, {'posSep','time'}, 'Color', 'r','HorizontalAlignment','Left','VerticalAlignment','Bottom')
                            end


                            % Plot the Cluster "Separation" metric vs. time.
                            subplot(414), hold on  %
                            shadedErrorBar([1:T]./(runParams.spp.*runParams.muW), meanCSep, stdCSep )
                                %plot([1:T]./(runParams.spp.*runParams.muW), minCSep-0.2,'k', 'LineWidth', 5)


                                for i = 1:T
                                    scatter(i./(runParams.spp.*runParams.muW), minCSep(i), 30, 'Filled', 'MarkerFaceColor',colors(minCSepID(i,1)))
                                    scatter(i./(runParams.spp.*runParams.muW), minCSep(i) + 0.05*max(minCSep), 30, 'Filled', 'MarkerFaceColor',colors(minCSepID(i,2)))
                                end


                                plot([1:T]./(runParams.spp.*runParams.muW), meanCSep,'k--', 'LineWidth', 5)

                            %     % Line at Launch Point (where "Separation" measure >
                            %     threshold)dbc
                            %     plot([LP/T, LP/T], [0, max(meanCSep)+max(stdCSep)],'r--')
                            %     text(0.99*(LP/T), 0.99*(max(meanCSep)+max(stdCSep)), 'LP', 'Color', 'r','HorizontalAlignment','Right','VerticalAlignment','Top')

                                title(['Cluster Separation'],'FontSize',20,'FontWeight','Bold')
                                xlabel(['Simulation Time (sec)'],'FontSize',18,'FontWeight','Bold')
                                ylabel({['Distance Between Cluster Centers /'],['sqrt Variance of Clusters']},'FontSize',18,'FontWeight','Bold')
                                set(gca,'FontSize',16,'FontWeight','Bold')
                                axis tight



                        end

        %                 keyboard
        
        
                        SeparationStruct.tp(BB,r) = tpo_mean; 
        

                    end % Loop Over Runs

                    if(h5r_flg)
                        % save and close figure showing sample 5 runs

                        if ~exist([imgsKurDir,'PIF_',PIFlg,'_N',num2str(Ndims(1)),'x',num2str(Ndims(2)),'_C',num2str(C),...
                            '_rmax',num2str(rmax),'_pfar',num2str(PconnFar),'/NF_',num2str(muW),'_',num2str(sigW),'/'],'dir')

                            mkdir([imgsKurDir,'PIF_',PIFlg,'_N',num2str(Ndims(1)),'x',num2str(Ndims(2)),'_C',num2str(C),...
                            '_rmax',num2str(rmax),'_pfar',num2str(PconnFar),'/NF_',num2str(muW),'_',num2str(sigW),'/'])

                        end

                        saveGoodImg(h5r,[imgsKurDir,'PIF_',PIFlg,'_N',num2str(Ndims(1)),'x',num2str(Ndims(2)),'_C',num2str(C),...
                            '_rmax',num2str(rmax),'_pfar',num2str(PconnFar),'/NF_',num2str(muW),'_',num2str(sigW),'/Win_',num2str(Strng),'_',...
                            num2str(sigSg),'_Wout_',num2str(Weak),'_',num2str(sigWk),'_5runs'],[0 0 1 1])
                        close(h5r);
                    end

%                     keyboard

                end % Loop Over Wout 

            end % Loop Over Win

        end % Loop Over C

    end % Loop Over Rmax 
   
end % Loop Over SigW


if(SepStructFlg)
    save([dataKurDir,'SeparationStruct_',SeparationInFlg],'SeparationStruct')
end

