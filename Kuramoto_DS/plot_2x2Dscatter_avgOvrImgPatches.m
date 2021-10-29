% This script will loop through mat files that were produced using
% exp_SepVPar_acgOvrImgPatches.m function.  It will produce 2x2Dscatter
% plots of Divisive Margin for Kuramoto_vs_Strawman and
% Eigenvector_vs_Strawman.  
%

[dirPre,sizeGoodIm] = onCluster;
dirScat = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/imgs/Kur_PIF_Fourier1/Mod_SKHAdj/AvgAcrossImgs/'];

files = dir([dirScat,'*tscale1*sP0p2*rM4*sDInf*sW0_*.mat']);

for F = 1:numel(files)
    
    % Load in Data from DivMarg mat files
    load([dirScat,files(F).name])
    
    % for parameters, turn numbers into strings for image filename.
    sPs = num2str(sigP);
    sPs(sPs=='.')='p';
    rMs = num2str(Rmax);
    rMs(rMs=='.')='p';
    sDs = num2str(sigD);
    sDs(sDs=='.')='p';
    sWs = num2str(sigW);
    sWs(sWs=='.')='p';
    Kss = num2str(Kscale);
    Kss(Kss=='.')='p';
    
    % if 2x2D Scatter Plot 'XXX' already exists, do not reproduce it.
    if 0 & exist([dirScat,'DivMarg_Scatter_2D_tscale',Tscale,'_sP',sPs,'_rM',rMs,'_sD',sDs,'_sW',sWs,'_Kscale',Kss,'_XXX.jpg'],'file')
        disp(['Plot DivMarg_Scatter_2D_tscale',Tscale,'_sP',sPs,'_rM',rMs,'_sD',sDs,'_sW',sWs,'_Kscale',Kss,'_XXX.jpg already exists. Not Reproducing it.'])
        continue
    end
    
    
    % If this misaligned data is not in the mat file, it is old and I can not use it.
    if ~exist('misalignedID','var')
        disp(['Old File: sP= ',num2str(sigP),', Rm: ',num2str(Rmax),', sD= ',num2str(sigD),', sW= ',num2str(sigW), ', Ks: ',num2str(Kscale),', Ts: ',Tscale])
        continue
    end
    

    % go through misaligned patches and find their indecies in data vectors.
    misaligned_indx = []; 
    for i = 1:numel(misalignedID)
        misaligned_indx = [misaligned_indx, find(strcmp(misalignedID{i},ImgPtchID))];
    end
    %
    good_indx = setxor(1:numel(ImgPtchID),misaligned_indx);


    % Calculate avg. distance of point from diagonal or some other metric on this 2D scatter point cloud.
    metric_Kur = mean( SMDivMarg(good_indx) - KurDivMarg(good_indx) );
    metric_Eig = mean( SMDivMarg(good_indx) - EigDivMarg(good_indx) );
    

    % Plot
    if(1)
        hSc = figure;   
        %
        hKS = subplot(121); 
        scatter( SMDivMarg(good_indx) , KurDivMarg(good_indx) ), hold on  
        %scatter( SMDivMarg(misaligned_indx) , KurDivMarg(misaligned_indx) ,'rx') 
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg of Coupled Oscillator Model','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg of Pixel Contrast Model','FontSize',18,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Kur,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        hES = subplot(122); 
        scatter( SMDivMarg(good_indx) , EigDivMarg(good_indx) ), hold on
        %scatter( SMDivMarg(misaligned_indx) , EigDivMarg(misaligned_indx) ,'rx')
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg of Eigenvector Computation','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg of Pixel Contrast Model','FontSize',18,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Eig,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        annotation('textbox', [0 0.9 1 0.1],'String', ...
            ['# Pts= ',num2str(numel(ImgPtchID)-numel(misaligned_indx)),' / ',num2str(numel(ImgPtchID)),', \sigmaP= ',num2str(sigP),', Rmax: ',num2str(Rmax),...
            ', \sigmaD= ',num2str(sigD),', \sigmaW= ',num2str(sigW), ', Kscale: ',num2str(Kscale),', Tscale: ',Tscale,' '], ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')
        %
        saveGoodImg(hSc,[dirScat,'DivMarg_Scatter_2D_tscale',Tscale,'_sP',sPs,'_rM',rMs,'_sD',sDs,'_sW',sWs,'_Kscale',Kss,'_XXX'],sizeGoodIm)
        close(hSc)
    end
    
    
    % store 2x2D scatter metric values into a vector
    scatter_Metric.sigP(F) = sigP;
    scatter_Metric.Rmax(F) = Rmax;
    scatter_Metric.sigD(F) = sigD;
    scatter_Metric.sigW(F) = sigW;
    scatter_Metric.Kscale(F) = Kscale;
    scatter_Metric.Tscale{F} = Tscale;
    %
    scatter_Metric.Kur_cloud_metric(F) = metric_Kur;
    scatter_Metric.Eig_cloud_metric(F) = metric_Eig;
    %
    scatter_Metric.num_good(F) = numel(good_indx);
    scatter_Metric.num_tot(F) = numel(ImgPtchID);
end






% Here save mat file with scatter_Metric data structure in it.
if(0)
    
    
end







% Here plot a scatter file with Rmax on x-axis and scatter_Metric on y-axis
% Scatter point colors code for sigP.  Scatter point shapes (x,o,s,...)
% code for Tscale. What codes for sigW & sigD?  Different plots, maybe.
if(0)
    
    scatter_Metric.Rmax(scatter_Metric.Rmax==inf) = 5; % Hard code this for now ...

    % find different values for each parameter
    sPu = unique(scatter_Metric.sigP);
    rMu = unique(scatter_Metric.Rmax);
    sDu = unique(scatter_Metric.sigD);
    sWu = unique(scatter_Metric.sigW);
    Ksu = unique(scatter_Metric.Kscale);
    Tsu = unique(scatter_Metric.Tscale);
    
    scat_colors = 'rgbk'; % colors encode Tscale -------- sigmaPix (0.1, 0.2, 0.3, 0.4)
    scat_shapes = '^v'; % shapes to encode sigDist  --------- 'oxsd';  % shapes encode Tscale (0.1, 0.5, 1)
    
    H=figure; hold on
    
    for i = 1:numel(sPu)
        for j = 1:numel(Tsu)

            ind = find( scatter_Metric.sigP==sPu(i) & strcmp(Tsu(j),scatter_Metric.Tscale) );
            
            indSDinf = find( scatter_Metric.sigD(ind) == inf ); 
            indSDfin = find( scatter_Metric.sigD(ind) ~= inf );
            
            

            scatter( scatter_Metric.Rmax(ind(indSDinf)) + i/(numel(sPu)+1) - 0.5, scatter_Metric.Kur_cloud_metric(ind(indSDinf)), 300, [scat_colors(j),scat_shapes(1)], 'LineWidth',3 )
            scatter( scatter_Metric.Rmax(ind(indSDfin)) + i/(numel(sPu)+1) - 0.5, scatter_Metric.Kur_cloud_metric(ind(indSDfin)), 300, [scat_colors(j),scat_shapes(2)], 'LineWidth',3 )

            % [scatter_Metric.sigP(ind); scatter_Metric.Rmax(ind); scatter_Metric.sigD(ind); scatter_Metric.sigW(ind); scatter_Metric.Kur_cloud_metric(ind)]
            
        end
    end
    %
    % draw lines to connect results with same Rmax, Tscale, sigD & sigW
    for i = 1:numel(rMu)
        for j = 1:numel(Tsu)
            for k = 1:numel(sDu)
                for L = 1:numel(sWu)
                   
                    ind = find( scatter_Metric.Rmax==rMu(i) & strcmp(Tsu(j),scatter_Metric.Tscale) & scatter_Metric.sigD==sDu(k)& scatter_Metric.sigW==sWu(L) );
                    
                    indSW0 = find( scatter_Metric.sigW(ind) == 0 );
                    indSWpos = find( scatter_Metric.sigW(ind) ~= 0 );
                    
                    plot( scatter_Metric.Rmax(ind(indSW0)) + 10*scatter_Metric.sigP(ind(indSW0))/(numel(sPu)+1) - 0.5, scatter_Metric.Kur_cloud_metric(ind(indSW0)), [scat_colors(j),'--'], 'LineWidth',1 )
                    plot( scatter_Metric.Rmax(ind(indSWpos)) + 10*scatter_Metric.sigP(ind(indSWpos))/(numel(sPu)+1) - 0.5, scatter_Metric.Kur_cloud_metric(ind(indSWpos)), [scat_colors(j),'-'], 'LineWidth',1 )
                    
                end
            end
        end
    end
    
    title('Evaluating Segmentation Performance over Parameter Space','FontSize',20,'FontWeight','Bold')
    xlabel('Rmax','FontSize',18,'FontWeight','Bold')
    ylabel('2D Scatter Point Cloud Metric','FontSize',18,'FontWeight','Bold')
    %
    plot([0.5,5.5],[0,0],'k--')
    plot([0.5,0.5],[min(scatter_Metric.Kur_cloud_metric), max(scatter_Metric.Kur_cloud_metric)],'k--')
    plot([1.5,1.5],[min(scatter_Metric.Kur_cloud_metric), max(scatter_Metric.Kur_cloud_metric)],'k--')
    plot([2.5,2.5],[min(scatter_Metric.Kur_cloud_metric), max(scatter_Metric.Kur_cloud_metric)],'k--')
    plot([3.5,3.5],[min(scatter_Metric.Kur_cloud_metric), max(scatter_Metric.Kur_cloud_metric)],'k--')
    plot([4.5,4.5],[min(scatter_Metric.Kur_cloud_metric), max(scatter_Metric.Kur_cloud_metric)],'k--','LineWidth',3)
    plot([5.5,5.5],[min(scatter_Metric.Kur_cloud_metric), max(scatter_Metric.Kur_cloud_metric)],'k--')
    %
    set(gca,'XTick',[1,2,3,4,5],'XTickLabel',{'1','2','3','4','\infty'},'FontSize',16,'FontWeight','Bold')
    %
    % put a 2nd x axis for sigma_pix using the text command.
    text( (1 + 2.5/(numel(sPu)+1) - 0.5), 1.05*min(scatter_Metric.Kur_cloud_metric), '\sigma_P','FontSize',16,'FontWeight','Bold','HorizontalAlignment','Center')
    text( (1 + 1/(numel(sPu)+1) - 0.5), min(scatter_Metric.Kur_cloud_metric), num2str(sPu(1),1),'FontSize',16,'FontWeight','Bold','HorizontalAlignment','Center')
    text( (1 + 2/(numel(sPu)+1) - 0.5), min(scatter_Metric.Kur_cloud_metric), num2str(sPu(2),1),'FontSize',16,'FontWeight','Bold','HorizontalAlignment','Center')
    text( (1 + 3/(numel(sPu)+1) - 0.5), min(scatter_Metric.Kur_cloud_metric), num2str(sPu(3),1),'FontSize',16,'FontWeight','Bold','HorizontalAlignment','Center')
    text( (1 + 4/(numel(sPu)+1) - 0.5), min(scatter_Metric.Kur_cloud_metric), num2str(sPu(4),1),'FontSize',16,'FontWeight','Bold','HorizontalAlignment','Center')
    
    
    % Save the plot.
    saveGoodImg(H,'Scatter_Point_Cloud_Results',sizeGoodIm)
    close(H)
    
end

