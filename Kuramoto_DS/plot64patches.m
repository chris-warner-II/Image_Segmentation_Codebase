function [H] =  plot64patches(MethodMan,StrawMan,StrawMan2sortBy,dirMeth,dirStraw,fileGeneral,fileSize,segMethod,netMethod,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,vert,horz)


        % Syntax: [H, DistSep, DivMarg] =  scatter_plot_DistPairs(MethodMan,StrawMan,fileGeneral,fileSize,segMethod,netMethod,rM,sD,sP,sW,Ts,Ks)
        %
        % This function takes in 2 x #{patches & groundTruth pairs} array
        % of measures of average pairwise distance within-cluster (row #1)
        % and average pairwise distance across-cluster (row #2) for a given
        % method (Kur, Eig1, ...) in METHODMAN variable and for the raw
        % image pixels in STRAWMAN variable.  FILEGENERAL is something like
        % 'BSDS_patch'. FILESIZE is something like '51x51_ds1' - containing
        % size of the image patch and the downsample amount.  SEGMETHOD
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
        
        [dirPre,sizeGoodIm] = onCluster;
        
        
        disp(['Plot 64 Patches and Results ',fileGeneral,' ',fileSize,' ',netMethod,' ',segMethod])
        
        
        
        SegRes =  (MethodMan(1,:)./MethodMan(2,:));
        StrawRs = (StrawMan(1,:)./StrawMan(2,:));

        DivMarg = ( StrawRs - SegRes );

        DistSep(1,:) = ( StrawMan(1,:) - MethodMan(1,:) ); % In-Cluster:  Best if +  
        DistSep(2,:) = ( StrawMan(2,:) - MethodMan(2,:) ); % Out-Cluster: Best if -

        D = ( DistSep(1,:) - DistSep(2,:) );
        
%         % For plotting scatter points with color showing metric value.
%         cb = colormap('jet');
%         D_bounds = linspace(min(D),max(D),size(cb,1)+1);
%         DM_bounds = linspace(min(DivMarg),max(DivMarg),size(cb,1)+1);
%         
%         Seg_bounds = linspace(min(SegRes),max(SegRes),size(cb,1)+1);
%         Straw_bounds = linspace(min(StrawRs),max(StrawRs),size(cb,1)+1);
%         
%         maxX = 1;
%         maxY = 1;


%         % Maybe a useful diagnostic visualization.
%         k=1000; [(StrawMan(1,1:k) - MethodMan(1,1:k)); (StrawMan(2,1:k) - MethodMan(2,1:k)); nan(1,k); (StrawMan(1,1:k)./StrawMan(2,1:k)); (MethodMan(1,1:k)./MethodMan(2,1:k))]
        


        % Set which metric we are looking at primarily, and which is secondary
        metricA = DivMarg;
        %metA_bnds = DM_bounds;
        %
        metricB = D;
        %metB_bnds = D_bounds;

        
        
        
        
%         vert=5;  % number of patches to put on figure in vertical axis
%         horz=12; % number of patches to put on figure in horizontal axis
        
        
        
        
        
        H = figure;
        Ha = tight_subplot((vert+1)*3, horz*3, 0, 0.05, 0.01); % [0.001, 0]
        
        for j = 1:numel(Ha)
            axes(Ha(j)),axis off
        end
        
        % Find ordering of image patches in terms of Pixel Divisive Margin.  How complicated are they expected to be to segment?
        DMSM = (StrawMan2sortBy(1,:)./StrawMan2sortBy(2,:));
        [Dsm,Ism] = sort(DMSM);

        ind = round(linspace(1,numel(ImgPtchID),vert*horz-1));
        %
        
        axes(Ha(1)), text(1,1,'Seg','HorizontalAlignment','Center','VerticalAlignment','Middle','FontSize',16,'FontWeight','Bold'), axis off, axis([0 2 0 2])
        axes(Ha(2)), text(1,1,'Img','HorizontalAlignment','Center','VerticalAlignment','Middle','FontSize',16,'FontWeight','Bold'), axis off, axis([0 2 0 2])
        axes(Ha(1+horz*3)),text(1,1,'gT','HorizontalAlignment','Center','VerticalAlignment','Middle','FontSize',16,'FontWeight','Bold'), axis off, axis([0 2 0 2])
        axes(Ha(2+horz*3)), axis off, axis([0 2 0 3])
        text(1,2,['\color{blue}{DM}'],'HorizontalAlignment','Center','VerticalAlignment','Middle','FontSize',14,'FontWeight','Bold') %,'color',cb())
        text(1,1,['\color{green}{DS}'],'HorizontalAlignment','Center','VerticalAlignment','Middle','FontSize',14,'FontWeight','Bold') %,'color',cb())
        
        J=0; % offset to make patches skip a row in between each row set.
        
        %
        for i = 1:numel(ind)

            ptch = load([dirStraw, ImgPtchID{Ism(ind(i))}]);
            
            % make sure to jump a row in between rows of patches
            if (mod(i,horz)==0)
                J=J+2*3*horz;
            end
    
            
            
            
            switch segMethod
                
                case{'Kur'}
                    
                    Kur = load([dirMeth,'KurMC_',ImgPtchID{Ism(ind(i))},'_rM',rM,'_sD',sD,'_sP',sP,'_NF_60_',sW,'_kscale',Ks,'_tscale',Ts,'_runs1.mat']);
                    ph = visKurPhase_inHSV( ptch.im, reshape(Kur.metaCluster.phaseAtClk(:,end), Kur.netParams.Ndims(1), Kur.netParams.Ndims(2) ) );
                    axes(Ha(1+3*i+J)), imagesc(ph), colormap('hsv'), axis off, freezeColors
                    
                case{'Eig1'}
                    
                    Eig = load([dirMeth,'Evecs_',ImgPtchID{Ism(ind(i))},'_rM',rM,'_sD',sD,'_sP',sP,'.mat']);
                    ph = reshape(Eig.EVecsML(:,1), Eig.netParams.Ndims(1), Eig.netParams.Ndims(2) );
                    axes(Ha(1+3*i+J)), imagesc(ph), colormap('jet'), axis off, freezeColors
                
                case{'Eig1v'}
                    
                    Eig = load([dirMeth,'Evecs_',ImgPtchID{Ism(ind(i))},'_rM',rM,'_sD',sD,'_sP',sP,'.mat']);
                    ph = EvecVizF( reshape(Eig.EVecsML(:,1), Eig.netParams.Ndims(1), Eig.netParams.Ndims(2) ), Eig.MC.vizNonLin );
                    axes(Ha(1+3*i+J)), imagesc(ph), colormap('jet'), axis off, freezeColors
                    
                case{'Eig2'}
                    
                    Eig = load([dirMeth,'Evecs_',ImgPtchID{Ism(ind(i))},'_rM',rM,'_sD',sD,'_sP',sP,'.mat']);
                    ph = reshape(Eig.EVecsML(:,2), Eig.netParams.Ndims(1), Eig.netParams.Ndims(2) );
                    axes(Ha(1+3*i+J)), imagesc(ph), colormap('jet'), axis off, freezeColors
                    
                case{'Eig2v'}
                    
                    Eig = load([dirMeth,'Evecs_',ImgPtchID{Ism(ind(i))},'_rM',rM,'_sD',sD,'_sP',sP,'.mat']);
                    ph = EvecVizF( reshape(Eig.EVecsML(:,2), Eig.netParams.Ndims(1), Eig.netParams.Ndims(2) ), Eig.MC.vizNonLin );
                    axes(Ha(1+3*i+J)), imagesc(ph), colormap('jet'), axis off, freezeColors
                
                case{'Eig3'}
                    
                    Eig = load([dirMeth,'Evecs_',ImgPtchID{Ism(ind(i))},'_rM',rM,'_sD',sD,'_sP',sP,'.mat']);
                    ph = reshape(Eig.EVecsML(:,3), Eig.netParams.Ndims(1), Eig.netParams.Ndims(2) );
                    axes(Ha(1+3*i+J)), imagesc(ph), colormap('jet'), axis off, freezeColors
                
                case{'Eig3v'}
                    
                    Eig = load([dirMeth,'Evecs_',ImgPtchID{Ism(ind(i))},'_rM',rM,'_sD',sD,'_sP',sP,'.mat']);
                    ph = EvecVizF( reshape(Eig.EVecsML(:,3), Eig.netParams.Ndims(1), Eig.netParams.Ndims(2) ), Eig.MC.vizNonLin );
                    axes(Ha(1+3*i+J)), imagesc(ph), colormap('jet'), axis off, freezeColors
                    
                case{'Eig1-2'}
                    
                    Eig = load([dirMeth,'Evecs_',ImgPtchID{Ism(ind(i))},'_rM',rM,'_sD',sD,'_sP',sP,'.mat']);
                    ph = reshape(sum(Eig.EVecsML(:,1:2),2), Eig.netParams.Ndims(1), Eig.netParams.Ndims(2) );
                    axes(Ha(1+3*i+J)), imagesc(ph), colormap('jet'), axis off, freezeColors

                case{'Eig1-2v'}
                    
                    Eig = load([dirMeth,'Evecs_',ImgPtchID{Ism(ind(i))},'_rM',rM,'_sD',sD,'_sP',sP,'.mat']);
                    ph = EvecVizF( reshape(sum(Eig.EVecsML(:,1:2),2), Eig.netParams.Ndims(1), Eig.netParams.Ndims(2)) , Eig.MC.vizNonLin );
                    axes(Ha(1+3*i+J)), imagesc(ph), colormap('jet'), axis off, freezeColors

                case{'Eig1-3'}
                    
                    Eig = load([dirMeth,'Evecs_',ImgPtchID{Ism(ind(i))},'_rM',rM,'_sD',sD,'_sP',sP,'.mat']);
                    ph = reshape(sum(Eig.EVecsML(:,1:3),2), Eig.netParams.Ndims(1), Eig.netParams.Ndims(2) );
                    axes(Ha(1+3*i+J)), imagesc(ph), colormap('jet'), axis off, freezeColors

                case{'Eig1-3v'}
                    
                    Eig = load([dirMeth,'Evecs_',ImgPtchID{Ism(ind(i))},'_rM',rM,'_sD',sD,'_sP',sP,'.mat']);
                    ph = EvecVizF( reshape(sum(Eig.EVecsML(:,1:3),2), Eig.netParams.Ndims(1), Eig.netParams.Ndims(2)) , Eig.MC.vizNonLin );
                    axes(Ha(1+3*i+J)), imagesc(ph), colormap('jet'), axis off, freezeColors
                 
            end
            
            axes(Ha(2+3*i+J)), imagesc(ptch.im), colormap('bone'), axis off, freezeColors
            axes(Ha(3+3*i+J)), text(1,1,num2str(i),'HorizontalAlignment','Left','VerticalAlignment','Middle','FontSize',14,'FontWeight','Bold'), axis off, axis([0.8 2 0 2])
            axes(Ha(1+3*i+horz*3+J)),
            try
                imagesc(ptch.gT{ImgGtID(ind(i))}), colormap('jet'), axis off, freezeColors
            catch
                text(1,1,'N/A','HorizontalAlignment','Center','VerticalAlignment','Middle'), axis off, axis([0 2 0 2])
            end
            
            axes(Ha(2+3*i+horz*3+J)), axis off, axis([0 2 0 3])
            text(1,2,[num2str(metricA(Ism(ind(i))),'%.2f')],'HorizontalAlignment','Center','VerticalAlignment','Middle','FontSize',14,'FontWeight','Bold') %,'color',cb())
            text(1,1,[num2str(metricB(Ism(ind(i))),'%.2f')],'HorizontalAlignment','Center','VerticalAlignment','Middle','FontSize',14,'FontWeight','Bold') %,'color',cb())

        end
        
        
        % Plot Da & Db metric values sorted by D_i.  Sorted by how easy the image is to segment.
        subplot(vert+1,1,vert+1), hold on
        scatter(1:numel(Ism),metricA(Ism),3,'b'); 
        scatter(1:numel(Ism),metricB(Ism),3,'g');
        plot(ind, Dsm(ind),'r--','LineWidth',2)
        % scatter(ind, Dsm(ind),50,'kx')
        text(1,-2,'\color{red}{easier}','HorizontalAlignment','Left','VerticalAlignment','Middle','FontSize',12,'FontWeight','Bold')
        text(numel(Ism),-2,'\color{red}{harder}','HorizontalAlignment','Right','VerticalAlignment','Middle','FontSize',12,'FontWeight','Bold')
       
        
        for i = 1:numel(ind)
            text(ind(i),-0.8,num2str(i))
        end
        
        plot([1 numel(Ism)],[0 0],'k--','LineWidth',2)
        plot([0 0],[-1 1],'k-','LineWidth',2)
        text(-1,1,['+'],'HorizontalAlignment','Right','VerticalAlignment','Top','FontSize',12,'FontWeight','Bold')
        text(-1,-1,['-'],'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',12,'FontWeight','Bold')
        axis([-1 numel(Ism) -2 2])
        axis off
        
        text(numel(Ism),1,['\color{red}{DM_{sm}}'],'HorizontalAlignment','Left','VerticalAlignment','Middle','FontSize',12,'FontWeight','Bold')
        text(numel(Ism),0.5,['\color{blue}{DivMarg} (mn=',num2str(mean(DivMarg),'%.2f'),')'],'HorizontalAlignment','Left','VerticalAlignment','Middle','FontSize',12,'FontWeight','Bold')
        text(numel(Ism),0,['\color{green}{DistSep} (mn=',num2str(mean(D),'%.2f'),')'],'HorizontalAlignment','Left','VerticalAlignment','Middle','FontSize',12,'FontWeight','Bold')
        
        %legend({'\Delta DivMarg','DistSep','DM_{sm}'},'Location','EastOutside')
       
        
        annotation('textbox', [0 0.9 1 0.1],'String', ...
                [fileGeneral,' ',fileSize,' \color{blue}{',netMethod,'} \color{red}{',segMethod,'} \color{black}{ - \{rM}\color{green}{',rM,'}\color{black}{,sD',sD,',sP',sP,',sW',sW,',Ts',Ts,',Ks',Ks,...
                '\} - #patches=',num2str(size(MethodMan,2)),'} '],'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',...
                20,'FontWeight','Bold')
            
            
        % add HSV circular colorbar in the bottom left of the figure.
        if strncmp(segMethod,'Kur',3)
            subplot('position',[0 0.05 0.12 0.12]),imshow([dirPre,'images/HSV_colorwheel.jpeg'])
        end
        
        
        %% add DistSep Scatter Plot in top Left.
%         subplot('position',[0 0.05 0.12 0.12])   (floor(vert/2),floor(horz/2),1), hold on
%         
%         
%         % scatter(DistSep(1,:),DistSep(2,:),3)
%         
%         % OR
%         for i = 1:size(cb,1)            
%             ind = find( metricA>=metA_bnds(i) & metricA<=metA_bnds(i+1) );
%             scatter(DistSep(1,ind),DistSep(2,ind),3,cb(i,:))
%         end
% 
%         %
%         errorbar( mean(DistSep(1,:)), mean(DistSep(2,:)), std(DistSep(1,:)), 'k', 'LineWidth', 1.5 )
%         herrorbar( mean(DistSep(1,:)), mean(DistSep(2,:)), std(DistSep(2,:)), 'k')
%         %
%         
%         plot([-1 1],[-1 1],'k--','LineWidth',2)
%         plot([-1 1], [0 0],'k--')
%         plot([0 0], [-1 1],'k--')
%         
%         
%         scatter( mean(DistSep(1,:)), mean(DistSep(2,:)), 100 ,'k', 'LineWidth', 2 ),
%         text( mean(DistSep(1,:)), mean(DistSep(2,:)), ['Mean = ',num2str(mean(D(~isnan(D))),3)],'FontWeight','Bold','FontSize',16,'VerticalAlignment','Top','HorizontalAlignment','Left');
% 
%         text(1, -0.95,'\color{green}{better}','FontSize',18,'FontWeight','Bold','HorizontalAlignment','right')
%         text(-1, 0.95, '\color{red}{worse}','FontSize',18,'FontWeight','Bold','HorizontalAlignment','left')
% 
%         text(1, 0.95,'\color{black}{separate all}','FontSize',18,'FontWeight','Bold','HorizontalAlignment','right')
%         text(-1,-0.95, '\color{black}{cluster all}','FontSize',18,'FontWeight','Bold','HorizontalAlignment','left')
%         
%         %axis([Lo(2) Hi(2) Lo(1) Hi(1)])
%         axis([-1 1 -1 1])
%         axis square
%         set(gca,'XTick',[-1 -1/2 0 1/2 1],'YTick',[-1 -1/2 0 1/2 1],'FontWeight','Bold','FontSize',18)
% 
%         ylabel('D_{out} = <D_{ij}>_s^{out} - <D_{ij}>_k^{out}','FontWeight','Bold','FontSize',18)
%         xlabel('D_{in} = <D_{ij}>_s^{in} - <D_{ij}>_k^{in}','FontWeight','Bold','FontSize',18)
%         title('Keeping separate','FontWeight','Bold','FontSize',20)
% 
%         %text(1.2*Lo(2),0.85*Hi(1),{'\color{red}{Out-Cluster}','\color{red}{Pts Clustered}'},'HorizontalAlignment','Center','VerticalAlignment','Middle','FontWeight','Bold','FontSize',18)
%         %text(1.2*Lo(2),0.85*Lo(1),{'\color{green}{Out-Cluster}','\color{green}{Pts Spread}'},'HorizontalAlignment','center','VerticalAlignment','Middle','FontWeight','Bold','FontSize',18)
% 
%         %text(0.85*Hi(2), 1.2*Lo(1),{'\color{green}{In-Cluster}','\color{green}{Pts Clustered}'},'HorizontalAlignment','Center','VerticalAlignment','Middle','FontWeight','Bold','FontSize',18)
%         %text(0.85*Lo(2), 1.2*Lo(1),{'\color{red}{In-Cluster}','\color{red}{Pts Spread}'},'HorizontalAlignment','Center','VerticalAlignment','Middle','FontWeight','Bold','FontSize',18)
% 
%         
%         
%         annotation('textbox', [0 0.9 1 0.1],'String', ...
%                 [fileGeneral,' ',fileSize,' - \{rM',rM,',sD',sD,',sP',sP,',sW',sW,',Ts',Ts,',Ks',Ks,'\} - #patches=',num2str(size(MethodMan,2))], ...
%                 'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')
        
        
        
        %% add DivMarg Scatter Plot in top Right
        
        
        
%         keyboard
            