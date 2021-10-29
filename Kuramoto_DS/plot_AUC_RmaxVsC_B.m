
pltRmaxVsC_flg = 0; % plot Rmax vs AUC for different C values below
pltRmax_avgC_flg = 0; % plot Rmax vs <AUC> (where <.> means average across all C values)

% Directory to look for data output from Kuramoto main to analyze here
dataDir = ['/Users/world7one/Desktop/Grad_School/Berkeley/Work/Fritz_Work/',...
    'Projects/output/Kuramoto/HandCookedNetwork/data/'];

load([dataDir,'AUC_ROC_stats_N1x48_RM1-10_C2-24_pfar0-10.mat'])



% not sure why this happens, maybe look back and figure out.  Ignore for now.
qFlat(isnan(qFlat)) = 0.5;



% Directory to save output plots into...
imgsDir = ['/Users/world7one/Desktop/Grad_School/Berkeley/Work/Fritz_Work/',...
    'Projects/output/Kuramoto/HandCookedNetwork/imgs/RmaxVsC_explore/'];
if ~exist(imgsDir,'dir') 
    mkdir(imgsDir)
end

colorsJ = colormap('jet');
colorsB = flipud(colormap('bone'));

tau = runParams.tau;
tstep = runParams.tau*runParams.spp; 

%% (0). extract parameters from ROCparams structure
SigW = ROCparams.sigW;
Strng = ROCparams.Strng;
Weak = ROCparams.Weak;
C = ROCparams.C;
rmax = ROCparams.rmax;
pfar = ROCparams.pfar;


uSW = unique(SigW); % find numbers of unique parameter settings to loop over below.
uST = unique(Strng);
uWK = unique(Weak);
uC = unique(C);
uRM = unique(rmax);
uPF = unique(pfar);


% to Label axis of plot
for i = 1:numel(uWK)
    uWKc{i} = num2str(uWK(i),2);
end
%
for i = 1:numel(uST)
    uSTc{i} = num2str(uST(i),2);
end
%
for i = 1:numel(uC)
    uCc{i} = num2str(uC(i),2);
    uCLc{i} = num2str(runParams.N./uC(i),2);
    leguC{i} = ['#C=',uCc{i},' : L=',num2str(48/uC(i))];
end
%
for i = 1:numel(uRM)
    uRMc{i} = num2str(uRM(i),2);
end


% uSW
% uC
% uST
% uWK
% uRM
% uPF



%% (1).  VERY USEFUL: Plot AUC vs. T curves for specific parameter settings.  
% Can use this to look further into runs that look interesting or dont make 
% sense in the grid scatter plot that shows AUC qFlat and tFlat for certain 
% param settings.
%
if(0)
    % which2Examine = 1:numel(tFlat); % look at all curves for all parameters settings
    which2Examine = find( SigW == uSW(4) & Strng == uST(4) & C == uC(3) & rmax == uRM(1) & pfar == uPF(2) ); % & Weak == uWK(1)

    disp(['Looking at ',num2str(numel(which2Examine)),' AUC vs. T curves.'])

    for i = which2Examine

        runParamsTag = ['NF ',num2str(ROCparams.sigW(i)),' | Win ',num2str(ROCparams.Strng(i)),' | Wout ',num2str(ROCparams.Weak(i)), ...
                        ' | C ',num2str(ROCparams.C(i)),' | Rmax ',num2str(ROCparams.rmax(i)),' | Pfar ',num2str(ROCparams.pfar(i))];


         % Plot mean AUC vs T curve for each run with Best & Flattening values on there
         figure, hold on,
         shadedErrorBar(1:numel(AUCm2{i}), AUCm2{i}, AUCs2{i})
         plot(AUCm2{i},'LineWidth',2),  
         %
         scatter(tFlat(i),AUCm2{i}(tFlat(i)),'rx','LineWidth',2)
         text(tFlat(i),0.99*AUCm2{i}(tFlat(i)),'Flat','VerticalAlignment','Top','HorizontalAlignment','Left','Color','r','FontWeight','Bold','FontSize',16)
         %
         scatter(bestT(i),AUCm2{i}(bestT(i)),'go','LineWidth',2)
         text(bestT(i),1.01*AUCm2{i}(bestT(i)),'Max','VerticalAlignment','Bottom','HorizontalAlignment','Left','Color','g','FontWeight','Bold','FontSize',16)
         %
         scatter(fast(i),AUCm2{i}(fast(i)),'cs','LineWidth',2)
         text(fast(i),1.01*AUCm2{i}(fast(i)),'t_{98%}','VerticalAlignment','Bottom','HorizontalAlignment','Left','Color','c','FontWeight','Bold','FontSize',16)
         %
         if( AUCm2{i}(round(end/2)) > 0.75)
             text(round(numel(AUCs2{i})/2), AUCm2{i}(round(end/2)) - AUCs2{i}(round(end/2)), 'Variability over 100 trials','VerticalAlignment','Top')
         else
             text(round(numel(AUCs2{i})/2), AUCm2{i}(round(end/2)) + AUCs2{i}(round(end/2)), 'Variability over 100 trials','VerticalAlignment','Bottom')
         end
         %
         xlabel('Time - periods of 60Hz signal','FontWeight','Bold','FontSize',16)
         ylabel('Area Underneath ROC Curve','FontWeight','Bold','FontSize',16)
         title([num2str(i),' :: ',runParamsTag],'FontWeight','Bold','FontSize',16)
         %
         axis([0 numel(AUCm2{i}) 0.5 1])
         pause
         close

         if (tFlat(i) > 60)
             keyboard
         end

    end

end








% Just to show that these two measures of Clustering Quality (q/tFlat & Best/t98%)
% are not ever very different, but the times they occur at are very different.
% Flattening out measure is better because tFlat has more meaning than bestT.
if(0)
    h=figure; 
    subplot(311), hist(best-qFlat,100)  % How much does qFlat underestimate quality?
    title(['Best larger than Flat by'],'FontWeight','Bold','FontSize',20)
    set(gca,'FontSize',16,'FontWeight','Bold')
    xlabel('AUC Quality Measure ranges from 0.5 (useless) to 1 (perfect)','FontSize',16,'FontWeight','Bold')
    subplot(312), hist(bestT-tFlat,100) % How much longer do I have to wait for best quality?
    title(['Flat quicker than Best by'],'FontWeight','Bold','FontSize',20)
    set(gca,'FontSize',16,'FontWeight','Bold')
    xlabel('Periods of 60Hz osc.','FontSize',16,'FontWeight','Bold')
    subplot(313), hist(fast-tFlat,100) % Relationship between t_{98%} and tFlat?
    title(['Flat quicker than t_{98%} by'],'FontWeight','Bold','FontSize',20)
    set(gca,'FontSize',16,'FontWeight','Bold')
    xlabel('Periods of 60Hz osc.','FontSize',16,'FontWeight','Bold')
    %
    % Save image
    saveGoodImg(h,[imgsDir,'qFlat_vs_qBest_Comparison'],[0 0 1 1])
    close(h) 
end









% keyboard







%% (2a). Plot AUC (averaged over different C values) vs Rmax for fixed values of Win, Wout.
%        Data generated with sigNF = 1.2Hz and Pfar = 0;
%

AUC_meanAcrossC_vs_rmax = zeros( numel(uST), numel(uWK), numel(uRM) );
AUC_stdAcrossC_vs_rmax = zeros( numel(uST), numel(uWK), numel(uRM) );

for j = 1:numel(uST)
    for k = 1:numel(uWK)
        for L = 1:numel(uRM)
            
            ind = find( SigW == uSW(1) & pfar == uPF(1) & Strng == uST(j) & Weak == uWK(k) & rmax == uRM(L) );
            [X,I] = sort(C(ind),'ascend'); 
        
        
%             [SigW(ind(I)); Strng(ind(I)); Weak(ind(I)); pfar(ind(I))]
%             [rmax(ind(I)); C(ind(I)); qFlat(ind(I))]
%             [mean(qFlat(ind)), std(qFlat(ind))]
            
            
            AUC_meanAcrossC_vs_rmax(j,k,L) = mean(best(ind));
            AUC_stdAcrossC_vs_rmax(j,k,L) = std(best(ind));
            %
            % Can also grab tFlat to look at speed of convergence dependce on Rmax.

        end
    end
end






% Plot separately for each Rmax/C pair the AUC.  You can see a trend with
% line sloping up or down or being flat for different Rmax values. Plot also
% the mean and std AUC for a single Rmax value across different size clusters.
if(pltRmax_avgC_flg)
    for j = 1:numel(uST)
        for k = 1:numel(uWK)

            runParamsTag = ['pfar ',num2str(uPF(1)),' | NF ',num2str(uSW(1)),' | Win ',num2str(uST(j)),' | Wout ',num2str(uWK(k))];

            % Plot Rmax vs. AUC (averaged over C values) - REINSTATE THIS. MOST META.
            h1=figure; hold on
            sp1=subplot(1,10,1:7); hold on
            sp2=subplot(1,10,8:10); hold on

            for L = 1:numel(uRM)

                ind = find( SigW == uSW(1) & pfar == uPF(1) & Strng == uST(j) & Weak == uWK(k) & rmax == uRM(L) );
                [X,I] = sort(C(ind),'ascend'); 

                subplot(sp1);
                plot(1:numel(I), best(ind(I)),'Color',colorsJ( round(size(colorsJ,1)*L/numel(uRM)), : ),'Marker','x', 'LineWidth',2)
                %
                subplot(sp2);
                errorbar(L, mean(best(ind)), std(best(ind)),'Color',colorsJ( round(size(colorsJ,1)*L/numel(uRM)), : ),'Marker','x', 'LineWidth',2)

            end


            figure(h1);
            subplot(sp1),
            axis([1 numel(uC) 0.5 1]);
            ylabel(['Segmentation Quality - (AUC ROC)'],'FontSize',18,'FontWeight','Bold')
            xlabel('Size of Clusters','FontSize',18,'FontWeight','Bold')
            title(['Cluster Extent vs. Rmax :: [',runParamsTag,']'],'FontSize',20,'FontWeight','Bold')
            set(gca,'XTick',1:numel(uC),'XTickLabel',uCLc,'FontSize',16,'FontWeight','Bold')
            %
            subplot(sp2),
            axis([0.5 numel(uRM)+0.5 0.5 1]);
            xlabel(['Rmax'],'FontSize',18,'FontWeight','Bold')
            set(gca,'XTick',1:numel(uRM),'XTickLabel',uRMc,'YTickLabel',[],'FontSize',16,'FontWeight','Bold')
            title(['Mean & Std across Cluser Size'],'FontSize',20,'FontWeight','Bold')
            %
            % Save image
            saveGoodImg(h1,[imgsDir,'AUC_avgAcrossC_NF',num2str(1),'_Win',num2str(j),'_Wout',num2str(k),'_Pfar',num2str(1)],[0 0 1 1])
            close(h1) 

        end % Loop over Wout (k)
    end     % Loop over Win  (j)
end







%% (2bB). Plot the same Grid Win vs Wout for Rmax = 1,2,3 with Pfar = 0.01, 0.05, 0.10. And for different sig NF.
%
%
for A = 1:numel(uSW)
    for B = 1:numel(uPF)

        AUC_meanAcrossC_vs_rmax_pfar{A,B} = zeros( numel(uST), numel(uWK), 3 );
        AUC_stdAcrossC_vs_rmax_pfar{A,B} = zeros( numel(uST), numel(uWK), 3 );

        for j = 1:numel(uST)
            for k = 1:numel(uWK)
                for L = 1:3      % Only got data for Rmax = 1,2,3.

                    ind = find( SigW == uSW(A) & pfar == uPF(B) & Strng == uST(j) & Weak == uWK(k) & rmax == uRM(L) & ...
                        (C == 2 | C == 3 | C == 4 | C == 6 | C == 8) ); % Need to specify cluster number for equal comparison after adding pfar below.
                    AUC_meanAcrossC_vs_rmax_pfar{A,B}(j,k,L) = mean(best(ind));
                    AUC_stdAcrossC_vs_rmax_pfar{A,B}(j,k,L) = std(best(ind));

                end % Loop over Rmax (L)
            end    % Loop over Wout (k)
        end     % Loop over Win  (j)

    end    % Loop over Pfar (0.01, 0.05, 0.10) %
end     % Loop over SigW (1.2, 1.8, 2.4, 3.0) Hz




% Find maximum value of std of AUC so colorbars are all consistent.
maxStd = max(AUC_stdAcrossC_vs_rmax(:));
for A = 1:numel(uSW)
    for B = 2:numel(uPF)   
        maxStd = max( maxStd, max(AUC_stdAcrossC_vs_rmax_pfar{A,B}(:)));
    end
end



%% (2b). Plot a grid Win vs. Wout (different one for each rmax) - with mean & std
%            
%
if(0)
    for i = 1:numel(uRM)

        h=figure; subplot(8,1,1:7),hold on

        for j = 1:numel(uST)
            for k = 1:numel(uWK)

                % use greyscale square to code for time to convergence
                cindx = max( round( AUC_stdAcrossC_vs_rmax(j,k,i)./max(AUC_stdAcrossC_vs_rmax(:)) * size(colorsB,1)), 1 );
                scatter(k,j, 6500, 'Marker','s', 'MarkerEdgeColor','k', 'MarkerFaceColor',colorsB(cindx,:), 'LineWidth',1)
                % use jet colormap on circle to code for quality of segmentation at convergence (AUC)
                cindx = max( round(2*( AUC_meanAcrossC_vs_rmax(j,k,i)-0.5)*size(colorsJ,1)), 1 );
                scatter(k,j, 3500, 'Marker','o', 'MarkerEdgeColor','k', 'MarkerFaceColor',colorsJ(cindx,:), 'LineWidth',1)

            end
        end
        %
        % title and axes
        set(gca, 'XTick',1:numel(uWK),'XTickLabel',uWKc,'YTick',1:numel(uST),'YTickLabel',uSTc,'FontSize',16,'FontWeight','Bold')
        title(['Averaging Across Cluster Size: on ',num2str(runParams.Ndims(1)),'x',num2str(runParams.Ndims(2)),...
            ' Oscillator Network -- ( \sigma_{\omega} = ',num2str(uSW(1)),' Hz)',...
            ' (Rmax = ',num2str(uRM(i)),') (pfar = ',num2str(uPF(1),2),')'],'FontSize',16,'FontWeight','Bold')

        xlabel(['Wout'],'FontSize',18,'FontWeight','Bold')
        ylabel(['Win'],'FontSize',18,'FontWeight','Bold')
        hcb=colorbar; set(hcb,'YTick',[0 .2 .4 .6 .8 1],'YTickLabel',{'0.5' '0.6','0.7','0.8','0.9','1'},'FontSize',16,'FontWeight','Bold')
        ylabel(hcb,'AUC ROC (mean across C)','FontSize',18,'FontWeight','Bold')
        axis([0.5 numel(uWK)+0.5 0.5 numel(uST)+0.5])
        pbaspect([numel(uWK),numel(uST),1])
        %
        % Show a pseudo - colorbar for std across C of AUC segmentation quality measure
        subplot(8,1,8), hold on,
        for I = 1:64
           cindx = max( round(I/64*size(colorsB,1)), 1 );
           scatter(I/64*maxStd, 0.2, 350, 'Marker','s', 'MarkerEdgeColor','k', 'MarkerFaceColor',colorsB(cindx,:), 'LineWidth',1)
        end
        set(gca,'Ytick',[],'FontSize',16,'FontWeight','Bold')
        xlabel('AUC ROC (std across C)','FontSize',16,'FontWeight','Bold')
        pbaspect([64,1,1])
        xlim([0 maxStd])
        %
        % Save image
        saveGoodImg(h,[imgsDir,'AUC_MeanAndStd_acrossC_N',num2str(runParams.Ndims(1)),'x',num2str(runParams.Ndims(2)),'_NF1_Pfar1_RM',num2str(i)],[0 0 0.7 0.53])
        close(h) 


    end
end
        
        




%% Plot and Save Grid for different combinations of SigW and Pfar.
if(0)    
    for A = 1:numel(uSW)
        for B = 2:numel(uPF)            
            for i = 1:3

                h=figure; subplot(8,1,1:7),hold on

                for j = 1:numel(uST)
                    for k = 1:numel(uWK)

                        % use greyscale square to code for time to convergence
                        cindx = max( round( AUC_stdAcrossC_vs_rmax_pfar{A,B}(j,k,i)./max(AUC_stdAcrossC_vs_rmax_pfar{A,B}(:)) * size(colorsB,1)), 1 );
                        scatter(k,j, 6500, 'Marker','s', 'MarkerEdgeColor','k', 'MarkerFaceColor',colorsB(cindx,:), 'LineWidth',1)
                        % use jet colormap on circle to code for quality of segmentation at convergence (AUC)
                        cindx = max( round(2*( AUC_meanAcrossC_vs_rmax_pfar{A,B}(j,k,i)-0.5)*size(colorsJ,1)), 1 );
                        scatter(k,j, 3500, 'Marker','o', 'MarkerEdgeColor','k', 'MarkerFaceColor',colorsJ(cindx,:), 'LineWidth',1)

                    end
                end
                %
                % title and axes
                set(gca, 'XTick',1:numel(uWK),'XTickLabel',uWKc,'YTick',1:numel(uST),'YTickLabel',uSTc,'FontSize',16,'FontWeight','Bold')
                title(['Averaging Across Cluster Size: on ',num2str(runParams.Ndims(1)),'x',num2str(runParams.Ndims(2)),...
                    ' Oscillator Network -- ( \sigma_{\omega} = ',num2str(uSW(A)),' Hz)',...
                    ' (Rmax = ',num2str(uRM(i)),') (pfar = ',num2str(uPF(B),2),')'],'FontSize',16,'FontWeight','Bold')

                xlabel(['Wout'],'FontSize',18,'FontWeight','Bold')
                ylabel(['Win'],'FontSize',18,'FontWeight','Bold')
                hcb=colorbar; set(hcb,'YTick',[0 .2 .4 .6 .8 1],'YTickLabel',{'0.5' '0.6','0.7','0.8','0.9','1'},'FontSize',16,'FontWeight','Bold')
                ylabel(hcb,'AUC ROC (mean across C)','FontSize',18,'FontWeight','Bold')
                axis([0.5 numel(uWK)+0.5 0.5 numel(uST)+0.5])
                pbaspect([numel(uWK),numel(uST),1])
                %
                % Show a pseudo - colorbar for std across C of AUC segmentation quality measure
                subplot(8,1,8), hold on,
                for I = 1:64
                   cindx = max( round(I/64*size(colorsB,1)), 1 );
                   scatter(I/64*maxStd, 0.2, 350, 'Marker','s', 'MarkerEdgeColor','k', 'MarkerFaceColor',colorsB(cindx,:), 'LineWidth',1)
                end
                set(gca,'Ytick',[],'FontSize',16,'FontWeight','Bold')
                xlabel('AUC ROC (std across C)','FontSize',16,'FontWeight','Bold')
                pbaspect([64,1,1])
                xlim([0 maxStd])
                %
                % Save image
                saveGoodImg(h,[imgsDir,'AUC_MeanAndStd_acrossC_N',num2str(runParams.Ndims(1)),'x',num2str(runParams.Ndims(2)),'_NF',num2str(A),'_Pfar',num2str(B),'_RM',num2str(i)],[0 0 0.7 0.53])
                close(h) 

            end % Loop over Rmax (1,2,3) pix
        end   % Loop over Pfar = (0.01, 0.05, 0.10)
    end   % Loop over SigW = (1.2, 1.8, 2.4, 3.0)
    
end % if statement














%% XXX. Compare performance for different Pfar values vs performance when Pfar = 0.
%
% May not be entirely fair comparison because one when pfar=0 calculated
% using larger variety of cluster sizes (nonzero pfar one does not use
% cases with many small clusters and these are typically where small Rmax
% performs very well).  Fix it.  

for i = 1:3 % Pfar = 0.01, 0.05, 0.10
    for j = 1:3 % Rmax = 1,2,3
        rat_PF = AUC_meanAcrossC_vs_rmax_pfar{1,1+i}(:,:,j) ./ AUC_meanAcrossC_vs_rmax_pfar{1,1}(:,:,j); 
        pfar_ratio_meanAcrossW(i,j) = mean(rat_PF(:));
        pfar_ratio_stdAcrossW(i,j) = std(rat_PF(:));
    end
end

h=figure; hold on
colorr = 'rgb';
for i = 3:-1:1
    errorbar(1:3, pfar_ratio_meanAcrossW(i,:), pfar_ratio_stdAcrossW(i,:),'Color',colorr(i),'LineWidth',2)
end
%
plot([0.5 3.5], [1 1],'k--')
text( 2.0, 1.1, {['\sigma_W = ',num2str(uSW(1))],['N1x48'],['#C = 2,3,4,6,8']},'FontSize',16,'FontWeight','Bold')
legend('Pfar = 10%','Pfar = 5%','Pfar = 1%')
title(['Improvement by Adding Few Random Far Connections?'],'FontSize',20,'FontWeight','Bold')
xlabel(['Rmax'],'FontSize',18,'FontWeight','Bold')
ylabel(['Performance Improvement Ratio'],'FontSize',18,'FontWeight','Bold')
set(gca,'FontSize',16,'FontWeight','Bold','XTick',1:3)
%
% Save image
saveGoodImg(h,[imgsDir,'Pfar_Improvement',num2str(runParams.Ndims(1)),'x',num2str(runParams.Ndims(2)),'_NF',num2str(1)],[0 0 0.7 0.53])
close(h) 



% Should think of how to look at other NF's without generating more data for Pfar = 0.  For larger values of NF, it
% is not really clear what is going on - whether including / increasing pfar helps performance.




















keyboard
















%% (2d). Plot AUC / speed of convergence vs Rmax.  Fix Win & Wout, SigW, C & Pfar. 
%
% THIS PLOT OF RMAX VS AUC IS NOT REALLY CLEAR RIGHT NOW.  NEED TO CLEAN UP. PROBABLY SOMETHING HERE THO.
%
for h = 1:numel(uPF)          % Loop over different Probability of Far Connection
    for i = 1:numel(uSW)      % Loop over different spread on Natural Frequency Distribution

        bestR_CvR_rmax = zeros( numel(uST), numel(uWK), numel(uC) );
        bestR_CvR_Auc = zeros( numel(uST), numel(uWK), numel(uC) );

        for j = 1:numel(uST)
            for k = 1:numel(uWK)

                runParamsTag = ['NF ',num2str(uSW(i)),' | Win ',num2str(uST(j)),' | Wout ',num2str(uWK(k))];

                % Plot Rmax vs. AUC (for different C values) - REINSTATE THIS. MOST META.
                if(pltRmaxVsC_flg)
                    h1=figure; hold on
                end
        
                
                
                for L = 1:numel(uC)
                    ind = find(SigW == uSW(i) & Strng == uST(j) & Weak == uWK(k) & C == uC(L));
                    [X,I] = sort(rmax(ind),'ascend'); 
                    %
                    if(pltRmaxVsC_flg)
                        figure(h1);
                        %subplot(211), hold on, 
                        plot(1:numel(I), qFlat(ind(I)),'Color',colorsJ( round(size(colorsJ,1)*L/numel(uC)), : ),'Marker','x', 'LineWidth',2)
                        % This assumes that files processed in order.  If any are present, none before them are missing. Maybe valid
                        %subplot(212), hold on, plot(1:numel(I), fast(ind(I)),'Color',colorsJ( round(size(colorsJ,1)*L/numel(uC)), : ),'Marker','x', 'LineWidth',2)
                        % This assumes that files processed in order.  If any are present, none before them are missing. Maybe valid
                    end


                    % Find Rmax value that yields the best AUC for each cluster size (given sigW, Win & Wout)
                    perf_CvR = qFlat(ind(I));
                    bestR_CvR_Auc(j,k,L) = max(perf_CvR);        % Highest AUC value given Win,Wout,C & sigW.  Vary Rmax.
                    rbest = find( perf_CvR == max(perf_CvR) );
                    bestR_CvR_rmax(j,k,L) = rbest(1);            % Rmax value when AUC was highest given set Win,Wout,C & sigW.

                end
                
                
                
                
                if(pltRmaxVsC_flg)
                    figure(h1);
                    %subplot(211),
                    ylabel(['Seg. Quality (AUC ROC)'],'FontSize',18,'FontWeight','Bold')
                    %xlabel([runParamsTag],'FontSize',18,'FontWeight','Bold')
                    title(['Cluster Extent / # Clusters vs. Rmax :: ',runParamsTag],'FontSize',20,'FontWeight','Bold')
                    set(gca,'XTick',1:numel(X),'XTickLabel',uRMc,'FontSize',16,'FontWeight','Bold')
                    axis([1 numel(I) 0.5 1]);
                    legend(leguC,'Location','Best')
                    %
                    %subplot(212),ylabel(['Seg. Speed'],'FontSize',18,'FontWeight','Bold')
                    xlabel(['Rmax (Coupling Interaction Extent)'],'FontSize',18,'FontWeight','Bold')
                    %set(gca,'XTick',1:numel(X),'XTickLabel',uRMc,'FontSize',16,'FontWeight','Bold')
                end










%                 keyboard




                % Plot mean AUC vs Time for different Rmax values.  Fixing C, sigW, Win, Wout.
                if(0)
                    for L = 1:numel(uC)

                        ind = find(SigW == uSW(i) & Strng == uST(j) & Weak == uWK(k) & C == uC(L));
                        [X,I] = sort(rmax(ind),'ascend'); 
                        figure, hold on,
                        for pp = 1:numel(ind)
                            plot(AUCm2{ind(I(pp))},'Color',colorsJ( round(size(colorsJ,1)*pp/numel(ind)), : ),'Marker','x', 'LineWidth',2)
                            % text()
                        end
                        legend(uRMc)

                        title([runParamsTag,' C = ',num2str(uC(L))],'FontSize',20,'FontWeight','Bold')
                        xlabel(['Time'],'FontSize',18,'FontWeight','Bold')
                        ylabel(['AUC'],'FontSize',18,'FontWeight','Bold')

                    end
                end



            end % Loop over Wout (k)
        end     % Loop over Win  (j)
        
        
        
        % If I want to do something here with the matrix of Win vs. Wout
        % vs. C I made I need to do it here.  Because it gets written over.
        
        mean(bestR_CvR_Auc,3)
        std(bestR_CvR_Auc,3)
        keyboard
        
        
        
        
    end         % Loop over SigW (i)
end             % Loop over Pfar (h)





















%% (3). Count up number of times that each value of Rmax yielded best performance
% for each different value of cluster size (C) and imagesc plot it C extent vs Rmax
if(1)
    for i = 1:numel(uC)
        x = bestR_CvR_rmax(:,:,i);
        y(:,i) = hist(x(:), numel(uRM));
    end
    %
    h=figure; imagesc(y'./40),
    hcb=colorbar;
    set(gca,'YTick',1:numel(uC),'YTickLabel',leguC,'XTick',1:numel(uRM),'XTickLabel',uRMc,'FontSize',16,'FontWeight','Bold')
    xlabel(['Rmax Value (Phase Interaction Spatial Extent)'],'FontSize',16,'FontWeight','Bold')
    ylabel(['Embedded Clusters'' Geometry'],'FontSize',16,'FontWeight','Bold')
    title(['Clustering Performance: Rmax vs. Cluster Extent (Avg''d across 4-Win 10-Wout Combinations)'],'FontSize',20,'FontWeight','Bold')
    ylabel(hcb,'Fraction of times Rmax Yielded Best Clustering (/40 each row)')
    %
    % Save image
    saveGoodImg(h,[imgsDir,'Rmax_vs_Cextent_ignoringWinWout'],[0 0 1 1])
    close(h) 
end






%% Plot Rmax vs Cluster Extent or different combinations of Win & Wout ??
if(1)
    
    for i = 1:numel(uSW)
        for j = 1:numel(uST)
            for k = 1:numel(uWK)
        
            
            runParamsTag = ['NF ',num2str(uSW(i)),' | Win ',num2str(uST(j)),' | Wout ',num2str(uWK(k))];
            
                for L = 1:numel(uC)
                    ind = find(SigW == uSW(i) & Strng == uST(j) & Weak == uWK(k) & C == uC(L));
                    
                    % Plot Rmax vs. C with quality as a colored circle and time as a grayscale box behind (like plot_ROC_AUC_stats)
                    
                    
                end
            end
        end
    end
    
    
    
    
end