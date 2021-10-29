function [H, D, DivMarg] =  scatter_plot_DistPairs(MethodMan,StrawMan,varsRange,dirMeth,dirStraw,fileGeneral,fileSize,segMethod,netMethod,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,GtNumSegs,pltFlg)


        % Syntax: [H, DistSep, DivMarg] =  scatter_plot_DistPairs(MethodMan,StrawMan,fileGeneral,fileSize,segMethod,netMethod,rM,sD,sP,sW,Ts,Ks)
        %
        % This function takes in 2 x #{patches & groundTruth pairs} array
        % of measures of average pairwise distance within-cluster (row #1)
        % and average pairwise distance across-cluster (row #2) for a given
        % method (Kur, Eig1, ...) in METHODMAN variable and for the raw
        % image pixels in STRAWMAN variable.  FILEGENERAL is something like
        % 'BSDS_patch'. FILESIZE is something like '51x51_ds1' - containing
        % size of the image patch and the downsample amount.  METHODNAME
        % will be the name of the network construction method: either
        % 'Mod_SKHAdj' or 'GLnrm' or whatever. The others are values of the
        % different parameters for the network construction {rM, sD, sP} and
        % Kuramoto simulation {sW, Ts, Ks}.  They are strings already and
        % should look something like {rM,sD,sP,sW,Ts,Ks} = {'7','Inf','0p2','0','1','300'}
        %
        % It outputs a handle to the plot (H) and subtractive differences
        % between STRAWMAN & METHODMAN for divisive margin (in DIVMARG) and
        % keeping distances in-cluster and across-cluster separate (in DISTSEP)
        %
        % I am adding the ability to plot image patch and results of
        % segmentation (either Kuramoto phases or Eigenvector or
        % Eigenvector with Nonlinear Visualization)
        
        disp(['Scatter Plot of ',fileGeneral,' ',fileSize,' ',segMethod])

        SegRes =  (MethodMan(1,:)./MethodMan(2,:));
        StrawRs = (StrawMan(1,:)./StrawMan(2,:));

        DivMarg = ( StrawRs - SegRes );        
        
%         % check sizes of veriables because I am getting an errror.
%         size(StrawRs)
%         size(StrawMan(2,:))
%         size(GtNumSegs)
%         size(SegRes)
%         size(MethodMan(2,:))
%         size(varsRange)
        
        
        Pmtric_S = (1 - StrawRs).*StrawMan(2,:).*(GtNumSegs-1); % varsRange for strawman = 1;
        Pmtric_M = (1 - SegRes).*MethodMan(2,:).*(GtNumSegs-1)./varsRange;
        

        D = Pmtric_M - Pmtric_S; % ( DistSep(1,:) - DistSep(2,:) );
        
        
        
        
        % For plotting scatter points with color showing metric value.
        cb = colormap('jet');
        D_bounds = linspace(min(D),max(D),size(cb,1)+1);
        DM_bounds = linspace(min(DivMarg),max(DivMarg),size(cb,1)+1);
        
        Seg_bounds = linspace(min(SegRes),max(SegRes),size(cb,1)+1);
        Straw_bounds = linspace(min(StrawRs),max(StrawRs),size(cb,1)+1);
        
        maxX = 1;
        maxY = 1;


%         % Maybe a useful diagnostic visualization.
%         k=1000; [(StrawMan(1,1:k) - MethodMan(1,1:k)); (StrawMan(2,1:k) - MethodMan(2,1:k)); nan(1,k); (StrawMan(1,1:k)./StrawMan(2,1:k)); (MethodMan(1,1:k)./MethodMan(2,1:k))]
        


        % Set which metric we are looking at primarily, and which is secondary
        metricA = DivMarg;
        metA_bnds = DM_bounds;
        %
        metricB = D;
        metB_bnds = D_bounds;


        
        
        % If MethodMan or StrawMan are all Nan's, BreakOut
        if ( isempty(find(~isnan(MethodMan))) | isempty(find(~isnan(StrawMan))) | ~pltFlg )
            H=0;
            return
        else
            H=figure; 
        end
        
        
        
        % I am no longer plotting this because it is not so informative.
        if(0)
            % Keep Within-Cluster and Across-Cluster Distance Points Separated
            subplot(221)
            hold on
            scatter( StrawMan(1,:) , MethodMan(1,:) , 3 ,'b' ), 
            scatter( StrawMan(2,:) , MethodMan(2,:) , 3 ,'r' ),
            %
            scatter( mean(StrawMan(1,:)) , mean(MethodMan(1,:)) , 100 ,'g', 'LineWidth', 2 ),
            scatter( mean(StrawMan(2,:)) , mean(MethodMan(2,:)) , 100 ,'k', 'LineWidth', 2 ),
            %
            errorbar( mean(StrawMan(1,:)) , mean(MethodMan(1,:)) , std(MethodMan(1,:)) , 'g', 'LineWidth', 1.5 )
            herrorbar( mean(StrawMan(1,:)) , mean(MethodMan(1,:)) , std(StrawMan(1,:)), 'g' )
            %
            errorbar( mean(StrawMan(2,:)) , mean(MethodMan(2,:)) , std(MethodMan(2,:)) , 'k', 'LineWidth', 1.5 )
            herrorbar( mean(StrawMan(2,:)) , mean(MethodMan(2,:)) , std(StrawMan(2,:)), 'k' )
            %
            %legend('In-Cluster Pairs','Out-Cluster Pairs','Location','EastOutside')
            axis square
            set(gca,'FontSize',16,'FontWeight','Bold'); %,'XTick',[0 pi/2 pi],'XTickLabel',{'0','\pi/2','\pi'},'YTick',[0 pi/2 pi],'YTickLabel',{'0','\pi/2','\pi'},'Interpreter','Tex')
            ylabel(['DistsPW ',MethodName],'FontSize',18,'FontWeight','Bold')
            xlabel('DistsPW Img','FontSize',18,'FontWeight','Bold')
            %
            text(1.05*maxX,0.8*maxY,'\color{blue}{witin-cluster}','FontSize',14,'FontWeight','Bold','HorizontalAlignment','left')
            text(1.05*maxX,0.9*maxY,'\color{red}{across-cluster}','FontSize',14,'FontWeight','Bold','HorizontalAlignment','left')
            %


            plot([0 max(maxX,maxY)], [0 max(maxX,maxY)],'k--','LineWidth',1.5)
            axis([0 maxX 0 maxY])
            %
            text(maxX,0.1*maxY,'\color{black}{clustered}','FontSize',18,'FontWeight','Bold','HorizontalAlignment','right')
            text(0,0.95*maxY,'\color{black}{separated}','FontSize',18,'FontWeight','Bold','HorizontalAlignment','left')
            %text(0.95*maxX,0.95*maxY,'\color{black}{unchanged}','FontSize',18,'FontWeight','Bold','HorizontalAlignment','center')
            %text(maxX,0.65*maxY, ['\color{black}{chng In=',num2str(mean(DistSep(1)),3),'}'],'FontSize',18,'FontWeight','Bold','HorizontalAlignment','right')
            %text(maxX,0.50*maxY,['\color{black}{chng Out=',num2str(mean(DistSep(2),3)),'}'],'FontSize',18,'FontWeight','Bold','HorizontalAlignment','right')
            % text(maxX,0.35*maxY,['\color{black}{chng Dif=',num2str(mean(DistSep(1)-DistSep(2)),3),'}'],'FontSize',18,'FontWeight','Bold','HorizontalAlignment','right')
        end
        
        
        
        
        
        
        % BELOW, I WANNA PLOT 2D scatter space an histogram of dist from diagonal in that space in these plots too.
        [xx,yy] = hist(D,100);

        
        Lo = [-1 -1]; % min(DistSep.Kur,[],2);
        Hi = [1 1]; % max(DistSep.Kur,[],2);

        subplot(221), hold on
        
        
        % scatter(DistSep(1,:),DistSep(2,:),3)
        
        % OR
        for i = 1:size(cb,1)            
            ind = find( metricA>=metA_bnds(i) & metricA<=metA_bnds(i+1) );
            scatter(DistSep(1,ind),DistSep(2,ind),3,cb(i,:))
        end

        %
        errorbar( mean(DistSep(1,:)), mean(DistSep(2,:)), std(DistSep(1,:)), 'k', 'LineWidth', 1.5 )
        herrorbar( mean(DistSep(1,:)), mean(DistSep(2,:)), std(DistSep(2,:)), 'k')
        %
        
        plot([-1 1],[-1 1],'k--','LineWidth',2)
        plot([-1 1], [0 0],'k--')
        plot([0 0], [-1 1],'k--')
        
        
        scatter( mean(DistSep(1,:)), mean(DistSep(2,:)), 100 ,'k', 'LineWidth', 2 ),
        text( mean(DistSep(1,:)), mean(DistSep(2,:)), ['Mean = ',num2str(mean(D(~isnan(D))),3)],'FontWeight','Bold','FontSize',16,'VerticalAlignment','Top','HorizontalAlignment','Left');

        text(1, -0.95,'\color{green}{better}','FontSize',18,'FontWeight','Bold','HorizontalAlignment','right')
        text(-1, 0.95, '\color{red}{worse}','FontSize',18,'FontWeight','Bold','HorizontalAlignment','left')

        text(1, 0.95,'\color{black}{separate all}','FontSize',18,'FontWeight','Bold','HorizontalAlignment','right')
        text(-1,-0.95, '\color{black}{cluster all}','FontSize',18,'FontWeight','Bold','HorizontalAlignment','left')
        
        %axis([Lo(2) Hi(2) Lo(1) Hi(1)])
        axis([-1 1 -1 1])
        axis square
        set(gca,'XTick',[-1 -1/2 0 1/2 1],'YTick',[-1 -1/2 0 1/2 1],'FontWeight','Bold','FontSize',18)

        ylabel('D_{out} = <D_{ij}>_s^{out} - <D_{ij}>_k^{out}','FontWeight','Bold','FontSize',18)
        xlabel('D_{in} = <D_{ij}>_s^{in} - <D_{ij}>_k^{in}','FontWeight','Bold','FontSize',18)
        title('Keeping separate','FontWeight','Bold','FontSize',20)

        %text(1.2*Lo(2),0.85*Hi(1),{'\color{red}{Out-Cluster}','\color{red}{Pts Clustered}'},'HorizontalAlignment','Center','VerticalAlignment','Middle','FontWeight','Bold','FontSize',18)
        %text(1.2*Lo(2),0.85*Lo(1),{'\color{green}{Out-Cluster}','\color{green}{Pts Spread}'},'HorizontalAlignment','center','VerticalAlignment','Middle','FontWeight','Bold','FontSize',18)

        %text(0.85*Hi(2), 1.2*Lo(1),{'\color{green}{In-Cluster}','\color{green}{Pts Clustered}'},'HorizontalAlignment','Center','VerticalAlignment','Middle','FontWeight','Bold','FontSize',18)
        %text(0.85*Lo(2), 1.2*Lo(1),{'\color{red}{In-Cluster}','\color{red}{Pts Spread}'},'HorizontalAlignment','Center','VerticalAlignment','Middle','FontWeight','Bold','FontSize',18)

        
        
        annotation('textbox', [0 0.9 1 0.1],'String', ...
                [fileGeneral,' ',fileSize,' - \{rM',rM,',sD',sD,',sP',sP,',sW',sW,',Ts',Ts,',Ks',Ks,'\} - #patches=',num2str(size(MethodMan,2))], ...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')
        
        
        %
        %
        % Show Divisive Margin Metric (Within-Cluster Average Pairwise Distance / Across-Cluster Average Pairwise Distance
        subplot(222), hold on
        
        %scatter( StrawRs , SegRes , 3 ,'b' ), 
        
        % OR
        for i = 1:size(cb,1)            
            ind = find( metricA>=metA_bnds(i) & metricA<=metA_bnds(i+1) );
            scatter( StrawRs(ind), SegRes(ind) , 3 ,cb(i,:))
        end
        %
        errorbar( mean(StrawRs) , mean(SegRes) , std(SegRes) , 'k', 'LineWidth', 1.5 )
        herrorbar( mean(StrawRs) , mean(SegRes) , std(StrawRs), 'k' )
        %
        scatter( mean(StrawRs) , mean(SegRes) , 100 ,'k', 'LineWidth', 2 ),
        text(mean(StrawRs) , mean(SegRes), ['Mean = ',num2str(mean(DivMarg(~isnan(DivMarg))),3)],'FontWeight','Bold','FontSize',16,'VerticalAlignment','Top','HorizontalAlignment','Left');
        
        %
        plot([0 1], [0 1],'k--','LineWidth',1.5)
        axis([0 1 0 1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold');
        [hx,hy] = format_ticks(gca,{'0','1'},{'0','1'},[0 1],[0 1]);
        ylabel(['DivMarg ',segMethod],'FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        %
        text(1, 0.05,'\color{green}{better}','FontSize',18,'FontWeight','Bold','HorizontalAlignment','right')
        text(0, 0.95, '\color{red}{worse}','FontSize',18,'FontWeight','Bold','HorizontalAlignment','left')
        %text(0.9,0.9,'\color{black}{neutral}','FontSize',18,'FontWeight','Bold','HorizontalAlignment','center')
        %text(0.1,0.85,['\color{black}{% chng DM=',num2str(mean(DivMarg),3),'}'],'FontSize',18,'FontWeight','Bold','HorizontalAlignment','left')
        % 
        title(['Divisive Margin'],'FontSize',20,'FontWeight','Bold')
        
        

        
        % Plot Histogram of Subtractive Difference Within vs. Without Metric
        if(0)
            subplot(625)
            hist(D,100);
            hold on
            plot([0 0], [0 0.1*numel(D)],'w--')
            xlabel('Diagonal Distance from Diagonal (D_{in}=D_{out}) Line','FontWeight','Bold','FontSize',18)
            ylabel('# Counts','FontWeight','Bold','FontSize',18)
            axis([-1 1 0 0.1*numel(D)]) % -.7071 .7071
            %
            a=arrowline([-1/2 -.7071],[0.05*numel(D) 0.05*numel(D)]);
            set(a,'Color','red','LineWidth',2)
            text(-.7071, 0.05*numel(D), '\color{red}{worse}','VerticalAlignment','Top','HorizontalAlignment','Left','FontWeight','Bold','FontSize',18)
            %
            a=arrowline([1/2 .7071],[0.05*numel(D) 0.05*numel(D)]);
            set(a,'Color','green','LineWidth',2)
            text(.7071, 0.05*numel(D), '\color{green}{better}','VerticalAlignment','Top','HorizontalAlignment','Right','FontWeight','Bold','FontSize',18)


            text(1/8, 0.08*numel(D), ['Mean = ',num2str(mean(D(~isnan(D))),3)],'FontWeight','Bold','FontSize',18);
            set(gca,'XTick',[-.7071 0 .7071],'FontWeight','Bold','FontSize',18)
        end

        
        
        
        
        
        % Plot Histogram of Subtractive Difference of Divisive Margin Metric
        if(0)
            subplot(626)
            hist(DivMarg,100), 
            hold on
            xlabel('Distance from DivMarg Diagonal Line','FontWeight','Bold','FontSize',18)
            ylabel('# Counts','FontWeight','Bold','FontSize',18)
            plot([0 0], [0 0.1*numel(DivMarg)],'w--')
            axis([-1 1 0 0.1*numel(DivMarg)]) % -.7071 .7071

            text(1/8, 0.08*numel(DivMarg), ['Mean = ',num2str(mean(DivMarg(~isnan(DivMarg))),3)],'FontWeight','Bold','FontSize',18);
            set(gca,'XTick',[-.7071 0 .7071],'FontWeight','Bold','FontSize',18)
            %
            a=arrowline([-1/2 -.7071],[0.05*numel(DivMarg) 0.05*numel(DivMarg)]);
            set(a,'Color','red','LineWidth',2)
            text(-.7071, 0.05*numel(DivMarg), '\color{red}{worse}','VerticalAlignment','Top','HorizontalAlignment','Left','FontWeight','Bold','FontSize',18)
            %
            a=arrowline([1/2 .7071],[0.05*numel(DivMarg) 0.05*numel(DivMarg)]);
            set(a,'Color','green','LineWidth',2)
            text(.7071, 0.05*numel(DivMarg), '\color{green}{better}','VerticalAlignment','Top','HorizontalAlignment','Right','FontWeight','Bold','FontSize',18)
        end
        
        
        
         % Plot Image Patches, Segmentations & Ground Truths on Scatter plot figure for sampling of performance values.

        % Distance from Diagonal Line in Left plots is our metric for now.
        % Could be Divisive Margin too.
       
        [Y,I] = sort(metricA);

        figure(H),

        numIms = 0; %12;
        ind = round(linspace(1,numel(find(~isnan(Y))),numIms));
        indCB = round(linspace(1,size(cb,1),numIms));
        
        
        
        
        % add x's onto scatter plot to indicate image patches shown
        subplot(222), scatter(StrawRs(I(ind)), SegRes(I(ind)), 30, 'kx')
        
        for i = 1:numIms            
            
            %plot Image Patch - Strawman
            ptch = load([dirStraw, ImgPtchID{I(ind(i))}]);
            subplot(7,numIms,numIms*3+i), imagesc(ptch.im), colormap('bone')
            axis square
            set(gca,'Xtick',[],'Ytick',[])
            freezeColors
            if(i==1)
                ylabel('Img','FontSize',18,'FontWeight','Bold')
                xlabel(['D_i=',num2str(  StrawRs(I(ind(i))) , '%.2f' )],'FontSize',12,'FontWeight','Bold')
            else
                xlabel([num2str(  StrawRs(I(ind(i))) , '%.2f' )],'FontSize',12,'FontWeight','Bold')
                
            end
            
            
            %plot Segmentation Results
            subplot(7,numIms,numIms*4+i)
            
            
            
            switch segMethod
                
                case{'Kur'}
                    
                    Kur = load([dirMeth,'KurMC_',ImgPtchID{I(ind(i))},'_rM',rM,'_sD',sD,'_sP',sP,'_NF_60_',sW,'_kscale',Ks,'_tscale',Ts,'_runs1.mat']);
                    ph = visKurPhase_inHSV( ptch.im, reshape(Kur.metaCluster.phaseAtClk(:,end), Kur.netParams.Ndims(1), Kur.netParams.Ndims(2) ) );
                    imagesc(ph), colormap('hsv')
                    
                case{'Eig1'}
                    
                    Eig = load([dirMeth,'Evecs_',ImgPtchID{I(ind(i))},'_rM',rM,'_sD',sD,'_sP',sP,'.mat']);
                    ph = reshape(Eig.EVecsML(:,1), Eig.netParams.Ndims(1), Eig.netParams.Ndims(2) );
                    imagesc(ph), colormap('jet')
                
                case{'Eig1v'}
                    
                    Eig = load([dirMeth,'Evecs_',ImgPtchID{I(ind(i))},'_rM',rM,'_sD',sD,'_sP',sP,'.mat']);
                    ph = EvecVizF( reshape(Eig.EVecsML(:,1), Eig.netParams.Ndims(1), Eig.netParams.Ndims(2) ), Eig.MC.vizNonLin );
                    imagesc(ph), colormap('jet')
                    
                case{'Eig2'}
                    
                    Eig = load([dirMeth,'Evecs_',ImgPtchID{I(ind(i))},'_rM',rM,'_sD',sD,'_sP',sP,'.mat']);
                    ph = reshape(Eig.EVecsML(:,2), Eig.netParams.Ndims(1), Eig.netParams.Ndims(2) );
                    imagesc(ph), colormap('jet')
                    
                case{'Eig2v'}
                    
                    Eig = load([dirMeth,'Evecs_',ImgPtchID{I(ind(i))},'_rM',rM,'_sD',sD,'_sP',sP,'.mat']);
                    ph = EvecVizF( reshape(Eig.EVecsML(:,2), Eig.netParams.Ndims(1), Eig.netParams.Ndims(2) ), Eig.MC.vizNonLin );
                    imagesc(ph), colormap('jet')
                
                case{'Eig3'}
                    
                    Eig = load([dirMeth,'Evecs_',ImgPtchID{I(ind(i))},'_rM',rM,'_sD',sD,'_sP',sP,'.mat']);
                    ph = reshape(Eig.EVecsML(:,3), Eig.netParams.Ndims(1), Eig.netParams.Ndims(2) );
                    imagesc(ph), colormap('jet')
                
                case{'Eig3v'}
                    
                    Eig = load([dirMeth,'Evecs_',ImgPtchID{I(ind(i))},'_rM',rM,'_sD',sD,'_sP',sP,'.mat']);
                    ph = EvecVizF( reshape(Eig.EVecsML(:,3), Eig.netParams.Ndims(1), Eig.netParams.Ndims(2) ), Eig.MC.vizNonLin );
                    imagesc(ph), colormap('jet')
                    
                case{'Eig1-2'}
                    
                    Eig = load([dirMeth,'Evecs_',ImgPtchID{I(ind(i))},'_rM',rM,'_sD',sD,'_sP',sP,'.mat']);
                    ph = reshape(sum(Eig.EVecsML(:,1:2),2), Eig.netParams.Ndims(1), Eig.netParams.Ndims(2) );
                    imagesc(ph), colormap('jet')

                case{'Eig1-2v'}
                    
                    Eig = load([dirMeth,'Evecs_',ImgPtchID{I(ind(i))},'_rM',rM,'_sD',sD,'_sP',sP,'.mat']);
                    ph = EvecVizF( reshape(sum(Eig.EVecsML(:,1:2),2), Eig.netParams.Ndims(1), Eig.netParams.Ndims(2)) , Eig.MC.vizNonLin );
                    imagesc(ph), colormap('jet')

                case{'Eig1-3'}
                    
                    Eig = load([dirMeth,'Evecs_',ImgPtchID{I(ind(i))},'_rM',rM,'_sD',sD,'_sP',sP,'.mat']);
                    ph = reshape(sum(Eig.EVecsML(:,1:3),2), Eig.netParams.Ndims(1), Eig.netParams.Ndims(2) );
                    imagesc(ph), colormap('jet')

                case{'Eig1-3v'}
                    
                    Eig = load([dirMeth,'Evecs_',ImgPtchID{I(ind(i))},'_rM',rM,'_sD',sD,'_sP',sP,'.mat']);
                    ph = EvecVizF( reshape(sum(Eig.EVecsML(:,1:3),2), Eig.netParams.Ndims(1), Eig.netParams.Ndims(2)) , Eig.MC.vizNonLin );
                    imagesc(ph), colormap('jet')
                 
            end
            
            
            
            
            axis square
            set(gca,'Xtick',[],'Ytick',[])
            freezeColors
            caxis([0 2*pi])
            if(i==1)
                ylabel('Seg','FontSize',18,'FontWeight','Bold')
                xlabel(['D_s=',num2str(  SegRes(I(ind(i))) , '%.2f' )],'FontSize',12,'FontWeight','Bold')
            else
                xlabel([num2str(  SegRes(I(ind(i))) , '%.2f' )],'FontSize',12,'FontWeight','Bold')
            end

            
            %plot ground truth
            subplot(7,numIms,numIms*5+i)
            try
            imagesc(ptch.gT{ImgGtID(ind(i))}), colormap('jet')
            catch
                
            end
            axis square
            set(gca,'Xtick',[],'Ytick',[])
            freezeColors
            if(i==1)
                ylabel('GT','FontSize',18,'FontWeight','Bold')
                xlabel({['D_A=',num2str(metricA(I(ind(i))),'%.2f')],['D_B=',num2str(metricB(I(ind(i))),'%.2f')]},'Color',cb(indCB(i),:),'FontSize',12,'FontWeight','Bold')
            else
                xlabel({[num2str(metricA(I(ind(i))),'%.2f')],[num2str(metricB(I(ind(i))),'%.2f')]},'Color',cb(indCB(i),:),'FontSize',12,'FontWeight','Bold')
            end
        end
        
        %keyboard
        
        
        % Plot D value vs DivMarg value
        if(1)
            
            [Ds,Is] = sort(metricA);
            
            subplot(7,1,7), plot(metricB(Is)./max(abs(metricB)),'color',[0.5,0.5,0.5]); % plot unsorted one in gray.
            
            hold on, 
            
            ind = round(linspace(1,size(metricA,2),size(cb,1)+1));
            
            for i = 1:size(cb,1)            
                plot([ind(i):ind(i+1)], Ds(ind(i):ind(i+1))./max(abs(Ds)),'color',cb(i,:)); % plot sorted one in with colors as in scatter plots.
            end
            
            plot([1 size(metricA,2)],[0 0],'k--')
            
            axis off
            
            
            % keyboard
            
        end
    
end
    
    
    
        
        
        
        
        
