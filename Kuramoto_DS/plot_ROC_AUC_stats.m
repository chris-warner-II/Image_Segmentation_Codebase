% Directory to look for data output from Kuramoto main to analyze here
dataDir = ['/Users/world7one/Desktop/Grad_School/Berkeley/Work/Fritz_Work/',...
    'Projects/output/Kuramoto/HandCookedNetwork/data/ROC_KurEig_Distil/'];

matfile = 'Stitch_KurEig_ROC_Distil_1_112';
load([dataDir,matfile,'.mat'])




% Directory to save output plots into...
imgsDir = ['/Users/world7one/Desktop/Grad_School/Berkeley/Work/Fritz_Work/',...
    'Projects/output/Kuramoto/HandCookedNetwork/imgs/SegStatsParamSearch/',matfile,'/'];
if ~exist(imgsDir,'dir') 
    mkdir(imgsDir)
end

colorsJ = colormap('jet');
colorsB = flipud(colormap('bone'));


tstep = DTot.params.runParams{1}.tau*DTot.params.runParams{1}.spp; 
% just indexing into first run since tau & spp didnt change.







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
    uWKc{i} = num2str(uWK(i),2);
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
if(1) % Plot Rmax & Cluster # combinations on different plots.
    
    
    for i = 1:numel(uSW)            % sigW - spread on natural frequency distribution
        for j = 1:numel(uRM)          % rmax - distance between oscillators overwhich influence can be felt
            for k = 1:numel(uC)        % number of clusters in the network
                for L = 1:numel(uPF)    % probability of far connection (beyond Rmax)
                
                    ind = find(SigW == uSW(i) & Rmax == uRM(j) & C == uC(k) & Pfar == uPF(L));
                    
                    if ~isempty(ind)

                        h=figure; subplot(8,1,1:7),hold on
                        
%                         tFlat_W = nan( numel(uWK), numel(uST) );
%                         qFlat_W = nan( numel(uWK), numel(uST) );
%                         vFlat_W = nan( numel(uWK), numel(uST) );

                        for I = 1:numel(ind)

                            % position of scatter points on grid.
                            x = find( uWK == Wout(ind(I)) );
                            y = find( uST == Win(ind(I)) );

                            % use greyscale square to code for time to convergence
                            cindx = max( round(tFlat(ind(I))./(runParams.T/runParams.spp-1)*size(colorsB,1)), 1 );
                            scatter(x,y, 6500, 'Marker','s', 'MarkerEdgeColor','k', 'MarkerFaceColor',colorsB(cindx,:), 'LineWidth',1)
                            % use jet colormap on circle to code for quality of segmentation at convergence (AUC)
                            cindx = max( round(2*(qFlat(ind(I))-0.5)*size(colorsJ,1)), 1 );
                            scatter(x,y, 3500, 'Marker','o', 'MarkerEdgeColor','k', 'MarkerFaceColor',colorsJ(cindx,:), 'LineWidth',1)
                            
%                             % Make matrices of Win vs. Wout for tFlat & qFlat so we can look at mean/std across Win or Wout
%                             tFlat_W(x,y) = tFlat(ind(I));
%                             qFlat_W(x,y) = qFlat(ind(I));
%                             vFlat_W(x,y) = vFlat(ind(I));




                        end
                        
                        
                        

                        

                        % title and axes
                        set(gca, 'XTick',1:numel(uWK),'XTickLabel',uWKc,'YTick',1:numel(uST),'YTickLabel',uSTc,'FontSize',16,'FontWeight','Bold')
                        title(['Parameter Grid Search Results on ',num2str(runParams.Ndims(1)),'x',num2str(runParams.Ndims(2)),...
                            ' Oscillator Network -- ( \sigma_{\omega} = ',num2str(uSW(i)),' Hz) (#C = ',num2str(uC(k)),')',...
                            ' (Rmax = ',num2str(uRM(j)),') (pfar = ',num2str(uPF(L),2),')'],'FontSize',16,'FontWeight','Bold')

                        % Put a Title over all subplots on the figure
        %                 annotation('textbox', [0 0.85 1 0.1],'String', ...
        %                     ['Parameter Grid Search Results on 12x1 Oscillator Network -- ( \sigma_{\omega} = ',num2str(uSW(i)),' Hz) (Rmax = ',...
        %                     num2str(uRM(j)),') (#C = ',num2str(uC(k)),')'],'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')

                        xlabel(['Wout'],'FontSize',18,'FontWeight','Bold')
                        ylabel(['Win'],'FontSize',18,'FontWeight','Bold')
                        hcb=colorbar; set(hcb,'YTick',[0 .2 .4 .6 .8 1],'YTickLabel',{'0.5' '0.6','0.7','0.8','0.9','1'},'FontSize',16,'FontWeight','Bold')
                        ylabel(hcb,'AUC ROC','FontSize',18,'FontWeight','Bold')
                        axis([0.5 numel(uWK)+0.5 0.5 numel(uST)+0.5])
                        pbaspect([numel(uWK),numel(uST),1])
                        % Show a pseudo - colorbar for time to 98% of best AUC value
                        subplot(8,1,8), hold on,
                        for I = 1:(runParams.T/runParams.spp-1)
                           cindx = max( round(I/(runParams.T/runParams.spp-1)*size(colorsB,1)), 1 );
                           scatter(I*tstep,0.2, 350, 'Marker','s', 'MarkerEdgeColor','k', 'MarkerFaceColor',colorsB(cindx,:), 'LineWidth',1)
                        end
                        set(gca,'Ytick',[],'FontSize',16,'FontWeight','Bold')
                        xlabel('Time to convergence (in seconds)','FontSize',16,'FontWeight','Bold')
                        pbaspect([(runParams.T-1),1,1])
                        xlim([0 runParams.Tsec])

                        % Save image
                        saveGoodImg(h,[imgsDir,'AUC_grid_SW',num2str(i),'_RM',num2str(j),'_pfar',num2str(L),'_C',num2str(k)],[0 0 0.7 0.53])
                        close(h) 
                        
                        
                        
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
                    
                end
            end
        end
    end
    
    
    
end