
dirPre = onCluster;

% Directory to look for data output from Kuramoto main to analyze here
dataDir = [dirPre,'output/Kuramoto/HandCookedNetwork/data/ROC_KurEig_Distil/'];

matfile = 'Stitch_KurEig_ROC_Distil_1_112';
load([dataDir,matfile,'.mat'])




% Directory to save output plots into...
imgsDir = [dirPre,'output/Kuramoto/HandCookedNetwork/imgs/SegStatsParamSearch/',matfile,'/'];
if ~exist(imgsDir,'dir') 
    mkdir(imgsDir)
end

colorsJ = colormap('jet');
colorsB = flipud(colormap('bone'));

runParams = DTot.params.runParams{1};
tstep = runParams.tau*runParams.spp; 
% just indexing into first run since tau & spp didnt change.


mrkrScat = 'oxs^+';

colorScat(1,:) = [1 0 0];       % red 
colorScat(2,:) = [0 1 0];       % green
colorScat(3,:) = [0 0 1];       % blue
colorScat(4,:) = [0 0 0];       % black
colorScat(5,:) = [0 1 1];       % cyan
colorScat(6,:) = [1 0 1];       % magenta
colorScat(7,:) = [1 1 0];       % yellow
colorScat(8,:) = [0.6 0.6 0.6]; % grey



colorScat2 = {'red','green','blue','black','cyan','magenta','yellow','black'};



%% Make Plots of Slices through parameter Space that show interesting Trends in best and fast


% (0). extract parameters from ROCparams structure
MuW = DTot.params.MuW;
SigW = DTot.params.SigW;
C = DTot.params.C;
N = DTot.params.N;
Win = DTot.params.Win;
Wout = DTot.params.Wout;
Rmax = DTot.params.Rmax;
Pfar = DTot.params.Pfar;


uSW = unique(SigW); % find numbers of unique parameter settings to loop over below.
uC = unique(C);
uST = unique(Win);
uWK = unique(Wout);
uRM = unique(Rmax);
uPF = unique(Pfar);


% to Label axis of plot
for i = 1:numel(uWK)
    uWKc{i} = num2str(uWK(i));
end
%
for i = 1:numel(uST)
    uSTc{i} = num2str(uST(i));
end


% % fixing errors below.  Time to flatten out can only be as long as total time.
% tFlat(tFlat > runParams.T/runParams.spp) = runParams.T/runParams.spp;
% 
% 
% % What is the standard deviation (of AUC over the 100 runs) at time when AUC curve flattens out?
% vFlat = zeros(size(tFlat));
% for i = 1:numel(tFlat)
%     vFlat(i) = AUCs2{i}(tFlat(i));
% end

%% Plot AUC vs T. curves  (labeling within cluster weight strength by line width and across cluster weight by line color). 
%  KEEP THIS.  COULD BE NICE PLOTS TO SEE LATER.
if(0)

    for i = 1:numel(uSW)            % sigW - spread on natural frequency distribution
        for j = 1:numel(uRM)        % rmax - distance between oscillators overwhich influence can be felt
            for k = 1:numel(uC)     % number of clusters in the network
                
                x = find(SigW == uSW(i) & rmax == uRM(j) & C == uC(k) );

                h=figure; hold on
                for pp = 1:numel(x)
                    %
                    LC = round( size(colorsJ,1) * find(uWK == Weak(x(pp)) ) ./ numel(uWK) );
                    %
                    switch Strng(x(pp)) % labeling within cluster weight strength by line width
                        case 10
                            LW = 1.5; MK = '-';
                        case 5
                            LW = 1.5; MK = '--';
                        case 1
                            LW = 1.5; MK = '-.';
                    end
                    %
                    plot(AUCm2{x(pp)}(1:59), 'LineStyle',MK, 'LineWidth',LW, 'Color',colorsJ(LC,:))
                end
                %
                title(['Oscillator Phase Separation Quality -- \sigma_{\omega} = ',num2str(uSW(i)),...
                    ' -- rmax = ',num2str(uRM(j)),' -- #C = ',num2str(uC(k)) ],'FontSize',20,'FontWeight','Bold')
                xlabel(['Time (# period of 60Hz oscillation)'],'FontSize',18,'FontWeight','Bold')
                ylabel(['Area Under ROC Curve'],'FontSize',18,'FontWeight','Bold')
                set(gca,'FontSize',16,'FontWeight','Bold')
                axis([1 60 0.4 1])
                %
                % Make a legend telling how line color matches with outside weights
                for pp = 1:numel(uWK)
                    LC = round( size(colorsJ,1) * pp ./ numel(uWK) );
                    text(2, (1-0.02*pp), ['Wout = ',num2str(uWK(pp),2)],'Color',colorsJ(LC,:),'FontSize',16,'FontWeight','Bold')
                end
                
                % Label Curves by their within cluster weights
                for pp = 1:numel(uST)
                    y = find(SigW == uSW(i) & rmax == uRM(j) & C == uC(k) & Strng == uST(pp) );
                    
                    WinTot = 0;
                    for qq = 1:numel(y)
                        WinTot = WinTot + max(AUCm2{y(qq)});
                    end
                    WinTot = WinTot./numel(y);
                    text(60, WinTot, ['Win = ',num2str(uST(pp),2)],'Color','k','FontSize',16,'FontWeight','Bold')
                    
                end
                %
                % Save image
                saveGoodImg(h,[imgsDir,'AUC_v_T_Curve_SW',num2str(i),'_RM',num2str(j),'_C',num2str(k)],[0 0 1 1])
                close(h)  


            end  % loop over clusters
        end      % loop over rmax
    end          % loop over sigma w
    
    
    keyboard
    
end






%% Plot Win vs Wout with quality/speed/stability indicated by color.
if(0) % Plot Rmax & Cluster # combinations on different plots.
    
    
    for i = 1:numel(uSW)            % sigW - spread on natural frequency distribution
        for j = 1:numel(uRM)          % rmax - distance between oscillators overwhich influence can be felt
            for k = 1:numel(uC)        % number of clusters in the network
                for L = 1:numel(uPF)    % probability of far connection (beyond Rmax)
                
                    ind = find(SigW == uSW(i) & Rmax == uRM(j) & C == uC(k) & Pfar == uPF(L));
                    
                    if ~isempty(ind)

                        
                        
                        % For Kuramoto, Plot Win vs. Wout Segmentation Performance 
                        % (max AUC & % time spent at max) for parameter combination (SigW, Rmax, Pfar, C)
                        h=figure; subplot(8,1,1:7), hold on

                        for I = 1:numel(ind)

                            % position of scatter points on grid.
                            x = find( uWK == Wout(ind(I)) );
                            y = find( uST == Win(ind(I)) );

                            % use greyscale square to code for time to convergence
                            cindx = max( round(DTot.singleRuns.pctTimeOn_mn(ind(I)).*size(colorsB,1)), 1 );
                            scatter(x,y, 6500, 'Marker','s', 'MarkerEdgeColor','k', 'MarkerFaceColor',colorsB(cindx,:), 'LineWidth',1)
                            % use jet colormap on circle to code for quality of segmentation at convergence (AUC)
                            cindx = max( round(DTot.singleRuns.maxAUC_mn(ind(I)).*size(colorsB,1)), 1 );
                            scatter(x,y, 3500, 'Marker','o', 'MarkerEdgeColor','k', 'MarkerFaceColor',colorsJ(cindx,:), 'LineWidth',1)

                        end


                        % title and axes
                        set(gca, 'XTick',1:numel(uWK),'XTickLabel',uWKc,'YTick',1:numel(uST),'YTickLabel',uSTc,'FontSize',16,'FontWeight','Bold')
                        title(['Parameter Grid Search Results on ',num2str(runParams.Ndims(1)),'x',num2str(runParams.Ndims(2)),...
                            ' Oscillator Network -- ( \sigma_{\omega} = ',num2str(uSW(i)),' Hz) (#C = ',num2str(uC(k)),')',...
                            ' (Rmax = ',num2str(uRM(j)),') (pfar = ',num2str(uPF(L)*100),'%)'],'FontSize',16,'FontWeight','Bold')

                        xlabel(['Wout'],'FontSize',18,'FontWeight','Bold')
                        ylabel(['Win'],'FontSize',18,'FontWeight','Bold')
                        hcb=colorbar; set(hcb,'YTick',[0 .2 .4 .6 .8 1],'YTickLabel',{'0.5' '0.6','0.7','0.8','0.9','1'},'FontSize',16,'FontWeight','Bold')
                        ylabel(hcb,'max AUC (avg''d across 100 runs)','FontSize',18,'FontWeight','Bold')
                        axis([0.5 numel(uWK)+0.5 0.5 numel(uST)+0.5])
                        pbaspect([numel(uWK),numel(uST),1])
                        % Show a pseudo - colorbar for time to 98% of best AUC value
                        subplot(8,1,8), hold on,
                        for I = 1:64
                           cindx = max( I, 1 );
                           scatter(I./64,0.2, 350, 'Marker','s', 'MarkerEdgeColor','k', 'MarkerFaceColor',colorsB(cindx,:), 'LineWidth',1)
                        end
                        set(gca,'Ytick',[],'FontSize',16,'FontWeight','Bold')
                        xlabel('% time spent at max AUC solution','FontSize',16,'FontWeight','Bold')
                        pbaspect([(runParams.T-1),1,1])
                        xlim([0 runParams.Tsec])

                        % Save image
                        saveGoodImg(h,[imgsDir,'AUC_Kur_maxNpct_SW',num2str(i),'_RM',num2str(j),'_pfar',num2str(L),'_C',num2str(k)],[0 0 0.7 0.53])
                        close(h) 
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        % For Eigenvector, Plot Win vs. Wout Segmentation Performance (max AUC) for parameter combination (Rmax, Pfar, C)
                        if(i==1)
                            h=figure; subplot(8,1,1:7), hold on

                            for I = 1:numel(ind)

                                % position of scatter points on grid.
                                x = find( uWK == Wout(ind(I)) );
                                y = find( uST == Win(ind(I)) );
                                %white square, just so it looks same as Kuramoto other
                                scatter(x,y, 6500, 'Marker','s', 'MarkerEdgeColor','k', 'MarkerFaceColor',colorsB(1,:), 'LineWidth',1)
                                % use jet colormap on circle to code for quality of segmentation at convergence (AUC)
                                cindx = max( round(DTot.Eigen.AUCe(ind(I)).*size(colorsB,1)), 1 );
                                scatter(x,y, 3500, 'Marker','o', 'MarkerEdgeColor','k', 'MarkerFaceColor',colorsJ(cindx,:), 'LineWidth',1)

                            end


                            % title and axes
                            set(gca, 'XTick',1:numel(uWK),'XTickLabel',uWKc,'YTick',1:numel(uST),'YTickLabel',uSTc,'FontSize',16,'FontWeight','Bold')
                            title(['Parameter Grid Search Results on ',num2str(runParams.Ndims(1)),'x',num2str(runParams.Ndims(2)),...
                                ' (#C = ',num2str(uC(k)),')', '(Rmax = ',num2str(uRM(j)),') (pfar = ',num2str(uPF(L)*100),'%)'],'FontSize',16,'FontWeight','Bold')

                            xlabel(['Wout'],'FontSize',18,'FontWeight','Bold')
                            ylabel(['Win'],'FontSize',18,'FontWeight','Bold')
                            hcb=colorbar; set(hcb,'YTick',[0 .2 .4 .6 .8 1],'YTickLabel',{'0.5' '0.6','0.7','0.8','0.9','1'},'FontSize',16,'FontWeight','Bold')
                            ylabel(hcb,'max AUC (avg''d across 100 runs)','FontSize',18,'FontWeight','Bold')
                            axis([0.5 numel(uWK)+0.5 0.5 numel(uST)+0.5])
                            pbaspect([numel(uWK),numel(uST),1])
                            %
                            % Save image
                            saveGoodImg(h,[imgsDir,'AUC_Eig_max_RM',num2str(j),'_pfar',num2str(L),'_C',num2str(k)],[0 0 0.7 0.53])
                            close(h) 
                        end
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        % DONT KNOW IF ANY OF THIS REVEALS ANYTHING OR IS WORTH DOING.  MAYBE COME BACK TO IT.
                        if(0)
                            % What is variability of performance across Win & Wout range or dependence of performance on Win/Wout?
                            for a = 1:size(qFlat_W,2)
                                ind =  find(~isnan(qFlat_W(:,a)));
                                qFlat_mnAccWout(a) = mean(qFlat_W(ind,a));
                                qFlat_stdAccWout(a) = std(qFlat_W(ind,a));
                                %
                                tFlat_mnAccWout(a) = mean(tFlat_W(ind,a));
                                tFlat_stdAccWout(a) = std(tFlat_W(ind,a));
                                %
                                vFlat_mnAccWout(a) = mean(vFlat_W(ind,a));
                                vFlat_stdAccWout(a) = std(vFlat_W(ind,a));
                            end
                            %
                            disp('How much does each Wout value''s performance deviate from the mean performance using all other Wouts?')
                            devQ = abs( ( qFlat_W - repmat( qFlat_mnAccWout, numel(uWK), 1 ) ) ./ repmat( qFlat_stdAccWout, numel(uWK), 1 ) );
                            for b = 1:size(devQ,1)
                                ind =  find(~isnan(devQ(b,:)));
                                devQ_mn(b) = mean(devQ(b,ind));
                            end
                            [uWK; devQ_mn]
                            % 
                            disp('How much does each Wout value''s convergence time deviate from the mean conv. time using all other Wouts?')
                            devT = abs( ( tFlat_W - repmat( tFlat_mnAccWout, numel(uWK), 1 ) ) ./ repmat( tFlat_stdAccWout, numel(uWK), 1 ) );
                            for b = 1:size(devT,1)
                                ind =  find(~isnan(devT(b,:)));
                                devT_mn(b) = mean(devT(b,ind));
                            end
                            [uWK; devT_mn]
                            % 
                            disp('How much does each Wout value''s variability at conv. time deviate from the mean variability at conv. time using all other Wouts?')
                            devV = abs( ( vFlat_W - repmat( vFlat_mnAccWout, numel(uWK), 1 ) ) ./ repmat( vFlat_stdAccWout, numel(uWK), 1 ) );
                            for b = 1:size(devV,1)
                                ind =  find(~isnan(devV(b,:)));
                                devV_mn(b) = mean(devV(b,ind));
                            end
                            [uWK; devV_mn]
                            % 







                            disp('What is average performance using this Win (across all Wouts)?')
                            [uST; qFlat_mnAccWout; qFlat_stdAccWout]
                            %
                            disp('What is average convergence time using this Win (across all Wouts)?')
                            [uST; tFlat_mnAccWout; tFlat_stdAccWout]
                            %
                            disp('What is average volatility at convergence using this Win (across all Wouts)?')
                            [uST; vFlat_mnAccWout; vFlat_stdAccWout]
                            %






                            if(1)
                                g1=figure; 
                                subplot(211)
                                shadedErrorBar(1:4, qFlat_mnAccWout, qFlat_stdAccWout)
                                title('Marginalize out Wout')
                                xlabel('Win')
                                ylabel('Mean ROC AUC Performance')
                                set(gca,'XTick',1:4,'XTickLabel',uSTc)
                                subplot(212),
                                stem(devQ_mn)
                                title('Wout values with most deviant performance')
                                xlabel('Wout')
                                ylabel('Average Deviation from Mean Performance (units=STD)')
                                set(gca,'XTickLabel',uWKc)


                                g2=figure; 
                                subplot(211)
                                shadedErrorBar(1:4, tFlat_mnAccWout, tFlat_stdAccWout)
                                title('Marginalize out Wout')
                                xlabel('Win')
                                ylabel('Mean Convergence Time')
                                set(gca,'XTick',1:4,'XTickLabel',uSTc)
                                subplot(212),
                                stem(devT_mn)
                                title('Wout values with most deviant performance')
                                xlabel('Wout')
                                ylabel('Average Deviation from Mean Time (units=STD)')
                                set(gca,'XTickLabel',uWKc)

                                g3=figure; 
                                subplot(211)
                                shadedErrorBar(1:4, vFlat_mnAccWout, vFlat_stdAccWout)
                                title('Marginalize out Wout')
                                xlabel('Win')
                                ylabel('Mean Volatility of AUC at Convergence')
                                set(gca,'XTick',1:4,'XTickLabel',uSTc)
                                subplot(212),
                                stem(devV_mn)
                                title('Wout values with most deviant performance')
                                xlabel('Wout')
                                ylabel('Average Deviation from Mean Volatility (units=STD)')
                                set(gca,'XTickLabel',uWKc)
                            end

                            keyboard
                        end
                        
 
                    
                    end
                    
                end % Pfar
            end % C
        end % Rmax
    end % sigW
    
    
    
end





%% With Plots below and all plots in future, I ask a question ... and answer it.

% Q: Is there any clear dependence of 2 performance metrics (AUCmax & %T@max) on SigW alone?
% A: No, it is not clear from these plots.
% 
if(0) % Scatter Plot maxAUC vs.  %T@Max for different parameter settings
    figure, hold on
    for i = 1:numel(uSW)            % sigW - spread on natural frequency distribution
                
        subplot(1,numel(uSW),i)
        ind = find(SigW == uSW(i));
        if ~isempty(ind)
            scatter( DTot.singleRuns.maxAUC_mn(ind), DTot.singleRuns.pctTimeOn_mn(ind),100, 'MarkerFaceColor',colorScat(i,:), mrkrScat(i) )
        end
        title(['sigW=',num2str(uSW(i))])
        xlabel('max AUC')
        ylabel('% time @ max')
        
    end
end


%% (1). Kuramoto AUC vs. %T@max Plot Series over SigW

% Q:  How is segmentation performance affected by Win, Wout, Rmax, SigW?
% A:  Yes, It's all in this series of plots. Lots of structure !!
%
% Plot Rmax in columns, Win in rows & Wout colorcoded in each plot.
% Make Separate Plot for each SigW value.
%
if(1) % Scatter Plot maxAUC vs.  %T@Max for different parameter settings
    
    for BB=1:numel(uSW)
    
        SIG = BB;

        h=figure; hold on
        for i = 1:numel(uST)   
            for j = 1:numel(uWK)
                for k = 1:numel(uRM)

                    subplot(numel(uST), numel(uRM), k + (numel(uRM))*(i-1)), hold on
                    ind = find(SigW == uSW(SIG) & Win == uST(i) & Wout == uWK(j) & Rmax==uRM(k) & Pfar==uPF(1) ); % look at only 1 sigW (sigW = 0.6) and one Rmax (Rmax=inf)
                    if ~isempty(ind)
                        % scatter( DTot.singleRuns.maxAUC_mn(ind), DTot.singleRuns.pctTimeOn_mn(ind),10, 'Filled', colorScat(j,:), mrkrScat(1) )                                      % Plot result from each cluster size small
                        herrorbar(mean(DTot.singleRuns.maxAUC_mn(ind)), mean(DTot.singleRuns.pctTimeOn_mn(ind)), std(DTot.singleRuns.maxAUC_mn(ind)), colorScat2{j} )              % Plot error bars in AUCmax
                        errorbar(mean(DTot.singleRuns.maxAUC_mn(ind)), mean(DTot.singleRuns.pctTimeOn_mn(ind)), std(DTot.singleRuns.pctTimeOn_mn(ind)), colorScat2{j} )            % Plot error bars in %Time@max
                        scatter( mean(DTot.singleRuns.maxAUC_mn(ind)), mean(DTot.singleRuns.pctTimeOn_mn(ind)),50, 'Filled', mrkrScat(3), 'MarkerFaceColor',colorScat(j,:), 'MarkerEdgeColor','k' )   % Plot average across clusters

                        % Note: Should also scatter plot diamond on x axis for eigenvector AUC value
                        % DTot.Eigen.AUCe(ind)


                    end

                    if(k==1)
                        ylabel(['Win=',num2str(uST(i))],'FontSize',16,'FontWeight','Bold')
                        %xlabel('max AUC')
                        %ylabel('% time @ max')
                    end

                    if(i==1)
                        title(['Rmax=',num2str(uRM(k))],'FontSize',16,'FontWeight','Bold')
                    end

                    axis([0.5 1 0 1])
                    set(gca,'XTick',[],'YTick',[])

                end
            end
        end

        % label axes on one plot
        subplot(numel(uST), numel(uRM), numel(uST).*numel(uRM))
        xlabel('AUC max','FontSize',12,'FontWeight','Bold')
        ylabel('% T @ max','FontSize',12,'FontWeight','Bold')
        %
        % label the different colors for Wout
        text(0.51,0.95,['Wout=',num2str(uWK(1))],'Color',colorScat(1,:))
        text(0.51,0.82,['Wout=',num2str(uWK(2))],'Color',colorScat(2,:))
        text(0.51,0.70,['Wout=',num2str(uWK(3))],'Color',colorScat(3,:))
        text(0.51,0.58,['Wout=',num2str(uWK(4))],'Color',colorScat(4,:))
        text(0.51,0.45,['Wout=',num2str(uWK(5))],'Color',colorScat(5,:))
        text(0.51,0.32,['Wout=',num2str(uWK(6))],'Color',colorScat(6,:))
        text(0.51,0.18,['Wout=',num2str(uWK(7))],'Color',colorScat(7,:))
        text(0.51,0.05,['Wout=',num2str(uWK(8))],'Color',colorScat(8,:))
        set(gca,'XTick',[0.6 0.8 1],'YTick',[0 .5 1],'FontSize',12,'FontWeight','Bold')

        % title containing other information that is fixed.
        annotation('textbox', [0 0.9 1 0.1],'String', ...
                            ['Kuramoto Segmentation Performance : Avg Across Cluster Size (C=2 --> C=24). (N = ',num2str(DTot.params.N(1)),') ( \sigma_{\omega} = ',num2str(uSW(SIG)),' Hz) (Pfar = ',num2str(uPF(1)),')'],...
                            'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold')

                        
        %
        % Save image
        saveGoodImg(h,[imgsDir,'AUCmax_v_pctT_N',num2str(DTot.params.N(1)),'_Win_Wout_Rmax_SigW',num2str(BB),'_pfar0'],[0 0 1 1])
        close(h)
        
    end
    
end



%% (2). Eigenvector AUC (C Rmax ( Win Wout ) ) Plot Series over Pfar

% Note: Averaging performance across Cluster Size C may be bad because AUC value from Eigenvector segmentation is different for different C's
%
% Q: How does AUC obtained from Eigenvector Segmentation depend on Wout, Win, Rmax and C?
% A: Yes, there is a lot of structure in these 2 plots !!
%
% Plot Rmax in columns, C in rows, Win & Wout in each subplot. Colorcoding AUC in Greyscale.
%
if(1) % Scatter Plot maxAUC vs.  %T@Max for different parameter settings
    
    for BB=1:numel(uPF)
    
        h=figure; hold on
    
        for i = 1:numel(uC)
            for j = 1:numel(uRM)   
                
                subplot( numel(uC), numel(uRM), j + numel(uRM)*(i-1) ), hold on 
                
                ind = find( C == uC(i) & Rmax == uRM(j) & SigW == uSW(1) & Pfar == uPF(BB) ); % & Win == uST(k) & Wout == uWK(L)
                
                
                for I = 1:numel(ind)
                    % position of scatter points on grid.
                    x = find( uWK == Wout(ind(I)) );
                    y = find( uST == Win(ind(I)) );

                    % use greyscale square to code for time to convergence
                    cindx = max( round(    2*abs(DTot.Eigen.AUCe(ind(I))-0.5)    *size(colorsJ,1)), 1);
                    scatter(y,x, 60, 'Marker','s', 'MarkerEdgeColor','k', 'MarkerFaceColor',colorsJ(cindx,:), 'LineWidth',1)
                end
                
                axis([0.5 numel(uST)+0.5 0.5 numel(uWK)+0.5])
                pbaspect([numel(uST),numel(uWK),1])
                set(gca,'XTick',[],'YTick',[])
                
                
                
                 
                if(j==1)
                    ylabel(['C=',num2str(uC(i))],'FontSize',16,'FontWeight','Bold')
                    %xlabel('max AUC')
                    %ylabel('% time @ max')
                end
                
                if(i==1)
                    title(['Rmax=',num2str(uRM(j))],'FontSize',16,'FontWeight','Bold')
                end
               
                
                
            end
        end
        
        
        % label axes on one plot
        subplot(numel(uC), numel(uRM), numel(uC).*numel(uRM))
        ylabel('Wout','FontSize',12,'FontWeight','Bold')
        xlabel('Win','FontSize',12,'FontWeight','Bold')
        set(gca,'YTick',[1:numel(uWK)] ,'XTick',[1:numel(uST)], 'YTickLabel',uWKc, 'XTickLabel',uSTc)
        rotateXLabels(gca(),90)
        %
        
        % title containing other information that is fixed.
        annotation('textbox', [0 0.9 1 0.1],'String', ...
                            ['Eigenvector Segmentation Performance : (N = ',num2str(DTot.params.N(1)),') (Pfar = ',num2str(uPF(BB)),')'],...
                            'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold')

                        
        %
        % Save image
        saveGoodImg(h,[imgsDir,'AUCe_N',num2str(DTot.params.N(1)),'_Win_Wout_C_Rmax_pfar',num2str(BB)],[0 0 1 1])
        close(h)
        
%         keyboard
        
        
    end
    
end


%% (3). Kuramoto max AUC (C Rmax ( Win Wout ) ) Plot Series over SigW



% Q: How does AUC obtained from Kuramoto depend on Wout, Win, Rmax, C & SigW ?
% A: Yes, there is structure in these plots also !!
%
% Plot Rmax in columns, C in rows, Win & Wout in each subplot. Colorcoding AUC in Greyscale.
%
if(1) % Scatter Plot maxAUC vs.  %T@Max for different parameter settings
    
    for BB=1:1 % numel(uPF) - Does not seem to be effected by Pfar.
        for CC=1:numel(uSW)
    
            h=figure; hold on

            for i = 1:numel(uC)
                for j = 1:numel(uRM)   

                    subplot( numel(uC), numel(uRM), j + numel(uRM)*(i-1) ), hold on 

                    ind = find( C == uC(i) & Rmax == uRM(j) & SigW == uSW(CC) & Pfar == uPF(BB) ); % & Win == uST(k) & Wout == uWK(L)


                    for I = 1:numel(ind)
                        % position of scatter points on grid.
                        x = find( uWK == Wout(ind(I)) );
                        y = find( uST == Win(ind(I)) );

                        % use greyscale square to code for time to convergence
                        cindx = max( round(  2*abs(DTot.singleRuns.maxAUC_mn(ind(I))-0.5)    *size(colorsJ,1)), 1);
                        scatter(y,x, 60, 'Marker','s', 'MarkerEdgeColor','k', 'MarkerFaceColor',colorsJ(cindx,:), 'LineWidth',1)
                    end

                    axis([0.5 numel(uST)+0.5 0.5 numel(uWK)+0.5])
                    pbaspect([numel(uST),numel(uWK),1])
                    set(gca,'XTick',[],'YTick',[])




                    if(j==1)
                        ylabel(['C=',num2str(uC(i))],'FontSize',16,'FontWeight','Bold')
                        %xlabel('max AUC')
                        %ylabel('% time @ max')
                    end

                    if(i==1)
                        title(['Rmax=',num2str(uRM(j))],'FontSize',16,'FontWeight','Bold')
                    end



                end
            end


            % label axes on one plot
            subplot(numel(uC), numel(uRM), numel(uC).*numel(uRM))
            ylabel('Wout','FontSize',12,'FontWeight','Bold')
            xlabel('Win','FontSize',12,'FontWeight','Bold')
            set(gca,'YTick',[1:numel(uWK)] ,'XTick',[1:numel(uST)], 'YTickLabel',uWKc, 'XTickLabel',uSTc)
            rotateXLabels(gca(),90)
            %

            % title containing other information that is fixed.
            annotation('textbox', [0 0.9 1 0.1],'String', ...
                                ['Kuramoto Segmentation Performance : (N = ',num2str(DTot.params.N(1)),') (Pfar = ',num2str(uPF(BB)),') (SigW = ',num2str(uSW(CC)),')'],...
                                'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold')


            %
            % Save image
            saveGoodImg(h,[imgsDir,'AUCk_max_N',num2str(DTot.params.N(1)),'_Win_Wout_C_Rmax_pfar',num2str(BB),'_sigW',num2str(CC)],[0 0 1 1])
            close(h)

    %         keyboard
        
        
        end
    end

    
end





%% (4).  Kuramoto %T@max (C Rmax ( Win Wout ) ) Plot Series over SigW



% Q: How does time spent at good AUC by simulation depend on Wout, Win, Rmax, C, SigW?
% A: 
%
% Plot Rmax in columns, C in rows, Win & Wout in each subplot. Colorcoding AUC in Greyscale.
%
if(0) % Scatter Plot maxAUC vs.  %T@Max for different parameter settings
    
    for BB=1:1 % numel(uPF) - Does not seem to be effected by Pfar.
        for CC=1:numel(uSW)
    
            h=figure; hold on

            for i = 1:numel(uC)
                for j = 1:numel(uRM)   

                    subplot( numel(uC), numel(uRM), j + numel(uRM)*(i-1) ), hold on 

                    ind = find( C == uC(i) & Rmax == uRM(j) & SigW == uSW(CC) & Pfar == uPF(BB) ); % & Win == uST(k) & Wout == uWK(L)


                    for I = 1:numel(ind)
                        % position of scatter points on grid.
                        x = find( uWK == Wout(ind(I)) );
                        y = find( uST == Win(ind(I)) );

                        % use greyscale square to code for time to convergence
                        cindx = max( round(  2*abs(DTot.singleRuns.pctTimeOn_mn(ind(I))-0.5)    *size(colorsJ,1)), 1);
                        scatter(y,x, 60, 'Marker','s', 'MarkerEdgeColor','k', 'MarkerFaceColor',colorsJ(cindx,:), 'LineWidth',1)
                    end

                    axis([0.5 numel(uST)+0.5 0.5 numel(uWK)+0.5])
                    pbaspect([numel(uST),numel(uWK),1])
                    set(gca,'XTick',[],'YTick',[])




                    if(j==1)
                        ylabel(['C=',num2str(uC(i))],'FontSize',16,'FontWeight','Bold')
                        %xlabel('max AUC')
                        %ylabel('% time @ max')
                    end

                    if(i==1)
                        title(['Rmax=',num2str(uRM(j))],'FontSize',16,'FontWeight','Bold')
                    end



                end
            end


            % label axes on one plot
            subplot(numel(uC), numel(uRM), numel(uC).*numel(uRM))
            ylabel('Wout','FontSize',12,'FontWeight','Bold')
            xlabel('Win','FontSize',12,'FontWeight','Bold')
            set(gca,'YTick',[1:numel(uWK)] ,'XTick',[1:numel(uST)], 'YTickLabel',uWKc, 'XTickLabel',uSTc)
            rotateXLabels(gca(),90)
            %

            % title containing other information that is fixed.
            annotation('textbox', [0 0.9 1 0.1],'String', ...
                                ['Kuramoto % Time @ max AUC : (N = ',num2str(DTot.params.N(1)),') (Pfar = ',num2str(uPF(BB)),') (SigW = ',num2str(uSW(CC)),')'],...
                                'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold')


            %
            % Save image
            saveGoodImg(h,[imgsDir,'AUCk_pctT_N',num2str(DTot.params.N(1)),'_Win_Wout_C_Rmax_pfar',num2str(BB),'_sigW',num2str(CC)],[0 0 1 1])
            close(h)

    %         keyboard
        
        
        end
    end

    
end






%% (5).  Kuramoto DurationOn (C Rmax ( Win Wout ) ) Plot Series over SigW



% Q: Does time that each oscillation spends at good segmentation depend on Wout, Win, Rmax, C, sigW?
% A: 
%
% Plot Rmax in columns, C in rows, Win & Wout in each subplot. Colorcoding AUC in Greyscale.
%
if(0) % Scatter Plot maxAUC vs.  %T@Max for different parameter settings
    
    for BB=1:1 % numel(uPF) - Does not seem to be effected by Pfar.
        for CC=1:numel(uSW)
    
            h=figure; hold on

            for i = 1:numel(uC)
                for j = 1:numel(uRM)   

                    subplot( numel(uC), numel(uRM), j + numel(uRM)*(i-1) ), hold on 

                    ind = find( C == uC(i) & Rmax == uRM(j) & SigW == uSW(CC) & Pfar == uPF(BB) ); % & Win == uST(k) & Wout == uWK(L)


                    for I = 1:numel(ind)
                        % position of scatter points on grid.
                        x = find( uWK == Wout(ind(I)) );
                        y = find( uST == Win(ind(I)) );

                        % use greyscale square to code for time to convergence
                        cindx = max( round(  DTot.singleRuns.meanOnDuration(ind(I))./(runParams.T./runParams.spp)    *size(colorsJ,1)), 1);
                        scatter(y,x, 60, 'Marker','s', 'MarkerEdgeColor','k', 'MarkerFaceColor',colorsJ(cindx,:), 'LineWidth',1)
                    end

                    axis([0.5 numel(uST)+0.5 0.5 numel(uWK)+0.5])
                    pbaspect([numel(uST),numel(uWK),1])
                    set(gca,'XTick',[],'YTick',[])




                    if(j==1)
                        ylabel(['C=',num2str(uC(i))],'FontSize',16,'FontWeight','Bold')
                        %xlabel('max AUC')
                        %ylabel('% time @ max')
                    end

                    if(i==1)
                        title(['Rmax=',num2str(uRM(j))],'FontSize',16,'FontWeight','Bold')
                    end



                end
            end


            % label axes on one plot
            subplot(numel(uC), numel(uRM), numel(uC).*numel(uRM))
            ylabel('Wout','FontSize',12,'FontWeight','Bold')
            xlabel('Win','FontSize',12,'FontWeight','Bold')
            set(gca,'YTick',[1:numel(uWK)] ,'XTick',[1:numel(uST)], 'YTickLabel',uWKc, 'XTickLabel',uSTc)
            rotateXLabels(gca(),90)
            %

            % title containing other information that is fixed.
            annotation('textbox', [0 0.9 1 0.1],'String', ...
                                ['Kuramoto % Time @ max AUC : (N = ',num2str(DTot.params.N(1)),') (Pfar = ',num2str(uPF(BB)),') (SigW = ',num2str(uSW(CC)),')'],...
                                'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold')


            %
            % Save image
            saveGoodImg(h,[imgsDir,'AUCk_OnDur_N',num2str(DTot.params.N(1)),'_Win_Wout_C_Rmax_pfar',num2str(BB),'_sigW',num2str(CC)],[0 0 1 1])
            close(h)

    %         keyboard
        
        
        end
    end

    
end



keyboard




















% What is the ratio of Wout/Win?  This plot shows it.
figure, surf( repmat(uWK',1,numel(uST)).*repmat(uST,numel(uWK),1) )