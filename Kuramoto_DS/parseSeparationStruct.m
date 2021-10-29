
dirPre = onCluster;

% Directory to look for data output from Kuramoto main to analyze here
dataDir = [dirPre,'output/Kuramoto/HandCookedNetwork/data/MetaCluster/'];

matfile = 'SeparationStruct_phase_wNegs';
load([dataDir,matfile,'.mat'])

SeparationStruct.tp = SeparationStruct.tp'; % flip this so combining code just below works right.


SStr = SeparationStruct; % Creating this structure to combine negative and non neg ones.


matfile = 'SeparationStruct_phase_noNeg';
load([dataDir,matfile,'.mat'])

SeparationStruct.tp = SeparationStruct.tp'; % flip this so combining code just below works right.




% %Combine SStr and SeparationStruct
% vars = fieldnames(SeparationStruct);
% for i=1:numel(vars)
%     eval(['SeparationStruct.',vars{i},'=[SeparationStruct.',vars{i},',SStr.',vars{i},'];'])
% end





% Directory to save output plots into...
imgsDir = [dirPre,'output/Kuramoto/HandCookedNetwork/imgs/MetaCluster/',matfile,'/'];
if ~exist(imgsDir,'dir') 
    mkdir(imgsDir)
end









colors = 'kbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmy';

colorsJ = colormap('jet');
colorsB = colormap('bone');

cmapRWB = rd_plotColorbar('redwhiteblue',256);
cmapWR = rd_plotColorbar('whitered',256);
% cmapWR = [[0 0 0]; cmapWR];      % Set zero = black because

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



%% Make Plots of Slices through parameter Space that show interesting Trends

% (0). extract parameters from Separation structure
SigW = SeparationStruct.sigW;
Win  = SeparationStruct.Win;
Wout = SeparationStruct.Wout;
C    = SeparationStruct.C;
Rmax = SeparationStruct.rmax;
meanCSep = SeparationStruct.meanOverRuns_CSep;
stdCSep = SeparationStruct.stdOverRuns_CSep;
medianCSep = SeparationStruct.medianOverRuns_CSep;
maxCSep = SeparationStruct.maxOverRuns_CSep;
%
tp = SeparationStruct.tp;




% Do some analysis on time that mean Separation measure goes positive.
% Find probabilty that Separation will go positive and mean time it crosses when it does.
T = 600;  % this is number of time points (I prob wanna save this from mCV)
runs = size(tp,2);
for i = 1:size(tp,1)
    pos = find(tp(i,:) ~= T);
    tp_pos_prob(i) = numel(pos) ./ runs;
    tp_mean(i) = mean(tp(i,pos));
end



uSW = unique(SigW); % find numbers of unique parameter settings to loop over below.
uC = unique(C);
uST = unique(Win);
uWK = unique(Wout);
uRM = unique(Rmax);


% to Label axis of plot
for i = 1:numel(uWK)
    uWKc{i} = num2str(uWK(i));
end
%
for i = 1:numel(uST)
    uSTc{i} = num2str(uST(i));
end
%
for i = 1:numel(uSW)
    uSWc{i} = num2str(uSW(i));
end
%
for i = 1:numel(uRM)
    uRMc{i} = num2str(uRM(i));
end



%% Plot a 2D Histogram of meanSeparation vs. stdSeparation (for different ind choices, subsets of conditions)
if(0)

    h=figure;
    
    for k = 1:numel(uRM)
    
%         k=1;
%         ind = find( SigW == uSW(k) ); % C == uC(i) & Rmax == uRM(j) & 
        ind = find( Rmax == uRM(k) );
%         ind = 1:numel(SeparationStruct.meanOverRuns_CSep);
        
        subplot(1,numel(uRM),k)
        
        [x,y] = hist3([SeparationStruct.meanOverRuns_CSep(ind);SeparationStruct.stdOverRuns_CSep(ind)]');
         hold on
        imagesc(y{2},y{1},x)
        caxis([0 300])
        plot([y{2}(1) y{2}(end)],[y{2}(1) y{2}(end)],'w--','LineWidth',2)
%         hcb=colorbar;
        %set(hcb,'Ylabel','7500')
        if(k==1)
        ylabel(['Mean of Separation'],'FontSize',18,'FontWeight','Bold')
        end
        if(k==2)
        xlabel(['STD of Separation (',num2str(size(SeparationStruct.tp,1)),' runs)'],'FontSize',18,'FontWeight','Bold')
        end
        set(gca,'FontSize',16,'FontWeight','Bold')
        title(['Rmax = ',num2str(uRM(k))],'FontSize',20,'FontWeight','Bold')
%         title(['Separation Statistics (',num2str(numel(SeparationStruct.C)),' Parameter Settings)'],'FontSize',20,'FontWeight','Bold')
        axis tight
        
    end
    
%     % title containing other information that is fixed.
%         annotation('textbox', [0 0.9 1 0.1],'String', ['Separation Stats'],...
%             'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')
    
    %
    % Save image
    saveGoodImg(h,[imgsDir,'SeparationStatisticsHist2D'],[0 0 1 0.4]) % [0 0 0.7 1] if 1.  [0 0 0.7 1] if 3.
    close(h)
    
end




%% Plot a 1D histogram of meanSep - stdSep.
if(0)

    h=figure;
    
    for k = 1 %1:numel(uRM)
    
%         ind = find( SigW == uSW(k) ); % C == uC(i) & Rmax == uRM(j) & 
%         ind = find( Rmax == uRM(k) );
        ind = 1:numel(SeparationStruct.meanOverRuns_CSep);
        
        [x,y] = hist([SeparationStruct.meanOverRuns_CSep(ind) - SeparationStruct.stdOverRuns_CSep(ind)]');
        hold on
        plot(y,x./sum(x),colors(k),'LineWidth',2)
        
        
        
        
%         ind = 1:numel(SStr.meanOverRuns_CSep);
%         [x,y] = hist([SStr.meanOverRuns_CSep(ind) - SStr.stdOverRuns_CSep(ind)]');
%         hold on
%         plot(y,x./sum(x),colors(2),'LineWidth',2)
        
        
        
        if(k==1)
            ylabel(['% Parameter Settings'],'FontSize',18,'FontWeight','Bold')
            xlabel(['Mean - STD of Separation (',num2str(size(SeparationStruct.tp,1)),' runs)'],'FontSize',18,'FontWeight','Bold')
        end
        set(gca,'FontSize',16,'FontWeight','Bold')

        title(['Separation Mean - STD'],'FontSize',20,'FontWeight','Bold')
        
    end
    
    plot([0 0], [0 max(x./sum(x))],'r--','LineWidth',2)

%     legend('R=1','R=4','R=\infty')  %('\sigma_\omega=0','\sigma_\omega=0.3','\sigma_\omega=0.6','\sigma_\omega=0.9','\sigma_\omega=1.2') %  % ('Weak Wout (+)', 'Repulsive Wout (-)')
    
    %
    % Save image
    saveGoodImg(h,[imgsDir,'SeparationStatisticsHist1D'],[0 0 1 0.4]) % [0 0 0.7 1] if 1.  [0 0 0.7 1] if 3.
    close(h)
    
end






%% #1B.  For Rmax=4, plot Win vs. Sep for different C in different colors. Do this because
%   we see vertical strips in the ClusterSep images.
if(0) 
    
    data2plot = SeparationStruct.meanOverRuns_CSep - SeparationStruct.stdOverRuns_CSep;
    %           this can be mean, median, max, etc.
    datamax = max(data2plot);
    datamin = min(data2plot);
    datalarger = max( datamax, abs(datamin) );
    
    for k = 1:numel(uSW)

        for i = 1:numel(uC)
            for j = 1:numel(uRM)  
                for L = 1:numel(uST)

                    ind = find( C == uC(i) & Rmax == uRM(j) & SigW == uSW(k) & Win == uST(L));

                    datamean(i,j,L) = mean(data2plot(ind));
                    datastd(i,j,L) = std(data2plot(ind));
                
                end

            end
        end

        
        
        % Note: Dim 2 is Rmax.  2 means Rmax = 4.
        h=figure; hold on
        plot(squeeze(datamean(:,2,:))','LineWidth',2)
        legend('C=2','3','4','6','8'),
        errorbar( squeeze(datamean(:,2,:))', squeeze(datastd(:,2,:))' ,'LineWidth',2) 
        plot([0 numel(uST)],[0 0],'k--','LineWidth',2)
        axis([1 numel(uST) -0.5 1])
        title(['Rmax=',num2str(uRM(2)),' \sigma_\omega=',num2str(uSW(3))],'FontSize',20,'FontWeight','Bold')
        xlabel('Win Strength','FontSize',18,'FontWeight','Bold')
        ylabel('Separation Metric (Mean - STD)','FontSize',18,'FontWeight','Bold')
        set(gca,'FontSize',16,'FontWeight','Bold')
        

        %
        % Save image
        saveGoodImg(h,[imgsDir,'iuno',num2str(k)],[0 0 0.5 0.8])
        close(h)
        

    end
    
end



%% #1.  Plot Median (over runs) of mean (over clusters) Separation Statistic on a 
% grid of Win & Wout inside a grid of C & Rmax.  Each figure will be SigW value.
if(0) 
    
    data2plot = SeparationStruct.meanOverRuns_CSep - SeparationStruct.stdOverRuns_CSep;
    %           this can be mean, median, max, etc.
    datamax = max(data2plot);
    datamin = min(data2plot);
    datalarger = max( datamax, abs(datamin) );
    
    for k = 1:numel(uSW)
    
        h=figure; hold on

        for i = 1:numel(uC)
            for j = 1:numel(uRM)   

                subplot( numel(uC), numel(uRM), j + numel(uRM)*(i-1) ), hold on 

                ind = find( C == uC(i) & Rmax == uRM(j) & SigW == uSW(k) );

                keyboard

                for I = 1:numel(ind)
                    % position of scatter points on grid.
                    x = find( uWK == Wout(ind(I)) );
                    y = find( uST == Win(ind(I)) );

                    % use colorscale to indicate Separation
                    cindx = max( round( ( data2plot(ind(I)) ./ datamax ) .* size(cmapRWB,1) ), 1);
                    scatter(y,x, 60, 'Marker','s', 'MarkerEdgeColor','k', 'MarkerFaceColor',cmapRWB(cindx,:), 'LineWidth',1)
                end

                axis([0.5 numel(uST)+0.5 0.5 numel(uWK)+0.5])
                pbaspect([numel(uST),numel(uWK),1])
                set(gca,'XTick',[],'YTick',[])

                if(j==1)
                    ylabel(['C=',num2str(uC(i))],'FontSize',16,'FontWeight','Bold')
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
                            ['Cluster Separation Measure \sigma_{NF} = ',num2str(uSW(k)), '  (extremes = ['...
                            '\color{blue}{',num2str(datamin,3),'},\color{red}{',num2str(datamax,3),'}\color{black}{])}'],...
                            'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold')

        %
        % Save image
        saveGoodImg(h,[imgsDir,'ClusterSepMnStd_SigW',num2str(k)],[0 0 0.5 0.8])
        close(h)
        

    end
    
end




%% #2.  Plot Probability that Separation measure goes positive for parameter set
%       and plot average time when it goes positive for the first time.

if(0) 
    
    data2plot = tp_mean;
    %           this can be tp_prob_pos or tp_mean.
    datamax = max(data2plot);
    
    for k = 1:numel(uSW)
    
        h=figure; hold on

        for i = 1:numel(uC)
            for j = 1:numel(uRM)   

                subplot( numel(uC), numel(uRM), j + numel(uRM)*(i-1) ), hold on 

                ind = find( C == uC(i) & Rmax == uRM(j) & SigW == uSW(k) );


                for I = 1:numel(ind)
                    % position of scatter points on grid.
                    x = find( uWK == Wout(ind(I)) );
                    y = find( uST == Win(ind(I)) );

                    % use colorscale to indicate Separation
                    cindx = max( round( ( data2plot(ind(I)) ./ datamax ) .* size(cmapWR,1) ), 1);
                    scatter(y,x, 60, 'Marker','s', 'MarkerEdgeColor','k', 'MarkerFaceColor',cmapWR(cindx,:), 'LineWidth',1)
                end

                axis([0.5 numel(uST)+0.5 0.5 numel(uWK)+0.5])
                pbaspect([numel(uST),numel(uWK),1])
                set(gca,'XTick',[],'YTick',[])

                if(j==1)
                    ylabel(['C=',num2str(uC(i))],'FontSize',16,'FontWeight','Bold')
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
                            ['Cluster Separation Measure \sigma_{NF} = ',num2str(uSW(k))],...
                            'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold')

        %
        % Save image
        saveGoodImg(h,[imgsDir,'ClusterSep_TPmean_SigW',num2str(k)],[0 0 0.5 0.8])
        close(h)
        

    end
    
end































%% HERE IS CODE LOOKING INTO STATISTICS WITH POSITIVE WOUT.  DOES NONZERO SIGW IMPROVE?
% load('/Users/world7one/Desktop/Grad_School/Berkeley/Work/Fritz_Work/Projects/output/Kuramoto/HandCookedNetwork/data/MetaCluster/SeparationStruct_phase_noNeg.mat')
%
% NOTE: I have to comment out the section about combining SeparationStruct and SStr above.



if(1)
    
    data2plot = meanCSep - stdCSep;

    % Look at all with postive Kout
    figure, subplot(3,1,1)
    [x,y] = hist(data2plot);
    plot(y,x,'LineWidth',2);
    title('Mean - Std Separation for For all +Kout')
    xlim([-0.5 0.5])

    % Split them up by SigW - Does adding desynchronization by larger sigW help in cases of +Kout?  It should)
    subplot(3,1,2), hold on
    for k = 1:numel(uSW)
        ind = find( SigW == uSW(k) );
        [x,y] = hist(data2plot(ind));
        plot(y,x,'Color',colorScat2{k},'LineWidth',2);
    end
    xlim([-0.5 0.5])
    legend(uSWc)
    title(['varying sigW'])

    % Split them up by Rmax - Smaller neighborhood does better because of over synchronization.
    subplot(3,1,3), hold on
    for k = 1:numel(uRM)
        ind = find( Rmax == uRM(k) );
        [x,y] = hist(data2plot(ind));
        plot(y,x,'Color',colorScat2{k},'LineWidth',2);
    end
    xlim([-0.5 0.5])
    legend(uRMc)
    title(['varying Rmax'])
    
    % Split them up by SigW with Rmax = 1)
    figure
    subplot(3,1,1), hold on
    for k = 1:numel(uSW)
        ind = find( SigW == uSW(k) & Rmax == uRM(1));
        [x,y] = hist(data2plot(ind));
        plot(y,x,'Color',colorScat2{k},'LineWidth',2);
    end
    xlim([-0.5 0.5])
    legend(uSWc)
    title(['Mean - Std Separation varying sigW with Rmax = ',num2str(uRM(1))])
    
    % Split them up by SigW with Rmax = 4)
    subplot(3,1,2), hold on
    for k = 1:numel(uSW)
        ind = find( SigW == uSW(k) & Rmax == uRM(2));
        [x,y] = hist(data2plot(ind));
        plot(y,x,'Color',colorScat2{k},'LineWidth',2);
    end
    xlim([-0.5 0.5])
    legend(uSWc)
    title(['varying sigW with Rmax = ',num2str(uRM(2))])
    
    % Split them up by SigW with Rmax = inf)
    subplot(3,1,3), hold on
    for k = 1:numel(uSW)
        ind = find( SigW == uSW(k) & Rmax == uRM(3));
        [x,y] = hist(data2plot(ind));
        plot(y,x,'Color',colorScat2{k},'LineWidth',2);
    end
    xlim([-0.5 0.5])
    legend(uSWc)
    title(['varying sigW with Rmax = ',num2str(uRM(3))])

end


















keyboard





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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