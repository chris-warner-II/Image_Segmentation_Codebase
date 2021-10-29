function visualize_pts_from_2x2D_scatter(fileType, fileSize, fileSpecifier, phaseMovie)

% 'DivMarg_Scatter_2D_tscale0p1_sP0p1_rM1_sD0p25_sW0_Kscale300.mat'

% This script/function will take in a mat file saved from
% exp_SepVPar_avgOverImgPatches.  That file contains the Divisive Margin
% results for Kuramoto, Eigenvector & Strawman segmentations. We devise a
% metric on the scatter points and sort points in order of this metric.
% Then we visualize the image patch and different segmentationsvis
% corresponding to certain points.
%
% fileType = 'BSDS_patch';
% fileSize = '51x51_ds1';
% fileSpecifier = '_sP0p2_rM4_sDInf_sW0_Ks300_Ts1' or something like.
%
% phaseMovie = 1; flag to make avi movie of oscillator phase evolution during Kuramoto simulation




% 'KurDivMarg','EigDivMarg','SMDivMarg','ImgPtchID','ImgGtID','fileType','fileSize','Tscale','sigP','Rmax','sigD','sigW','Kscale'



%% Load in mat file and set up directories to different data.
[dirPre,sizeGoodIm] = onCluster;



% for now, delete this. Saving these variables in the explore_compare_DivMarg code.
sP = '0p2';
rM = '4';
sD = 'Inf';
Ks = '300';
Ts = '1';
sW = '0';



paramsStr = fileSpecifier;
paramsStr(paramsStr=='_')=' ';

fileTypeStr = fileType;
fileTypeStr(fileTypeStr=='_')=' ';

fileSizeStr = fileSize;
fileSizeStr(fileSizeStr=='_')=' ';

% mat file saved from exp_SepVPar_avgOverImgPatches.
matFilesDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/Kur_PIF_Fourier1/Mod_SKHAdj/'];

matFiles = dir([matFilesDir,'DivMarg_Results_allPatches*',fileSpecifier,'*.mat']); % <-- Can add more specifiers here to choose only certain parameter value combinations.


% directory to save output image files into (includes a 2D Histogram of the DivMarg Scatter Plots for Kur_Vs_Straw & Eig1vis_Vs_Straw)
imgKurDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/imgs/Kur_PIF_Fourier1/Mod_SKHAdj/'];
%
imgEigDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/imgs/spectral/Mod_SKHAdj/'];


for j = 1:numel(matFiles)

    disp(['File # ',num2str(j),' / ',num2str(numel(matFiles))])
    load([matFilesDir,matFiles(j).name]);
    
    
    % go through misaligned patches and find their indecies in data vectors.
    misaligned_indx = []; 
    for i = 1:numel(misalignedID)
        misaligned_indx = [misaligned_indx, find(strcmp(misalignedID{i},ImgPtchID))];
    end
    %
    good_indx = setxor(1:numel(ImgPtchID),misaligned_indx);

    
    
    
    % Keep only data points that are not misaligned (for some, Kuramoto looks at different patch out of image than it should).
    KurDivMarg = KurDivMarg(good_indx)';
    EigDivMarg = Eig1vDivMarg_L(good_indx);
    SMDivMarg_L = SMDivMarg_L(good_indx);
    SMDivMarg_C = SMDivMarg_C(good_indx);
    ImgPtchID = ImgPtchID(good_indx); % this IS the correct way to index into this cell array.
    ImgGtID = ImgGtID(good_indx);
    % still seems to have misaligned patch in top performers even tho I get rid of them here...why?



    % directory to image patches extracted from BSDS.
    patchesDir = [dirPre,'images/',fileType,'/',fileSize,'/'];

    % directory to KurMC files and to Kur_metaSummary files (think I need the KurMC ones)
    dirKurMats = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/Kur_PIF_Fourier1/Mod_SKHAdj/'];

    % directory to the Evecs files
    dirEigMats = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/spectral/Mod_SKHAdj/'];


    
    
%     % turn Parameter value numbers into strings for image filename.
%     sPs = num2str(sigP);
%     sPs(sPs=='.')='p';
%     rMs = num2str(Rmax);
%     rMs(rMs=='.')='p';
%     sDs = num2str(sigD);
%     sDs(sDs=='.')='p';
%     sWs = num2str(sigW);
%     sWs(sWs=='.')='p';
%     Kss = num2str(Kscale);
%     Kss(Kss=='.')='p';
    
    
    
    
    % Calculate a 'goodness of Kuramoto/Eigen-comp over Strawman' metric in 2D scatter space
    metric_Kur = SMDivMarg_C - KurDivMarg;
    metric_Eig = SMDivMarg_L - EigDivMarg; % distance above (+) or below (-) diagonal line.  Good if below the line.
    [Y,Itop] = sort(metric_Kur,'descend');                          % best performing patches
    [Y,Ibot] = sort(metric_Kur,'ascend');                           % worst performing patches
    [Y,Imn] = sort(abs(metric_Kur - mean(metric_Kur) ),'ascend');   % patches closest to mean performance
    TBM = {'Itop'}; % 'Ibot', ,'Imn'

    disp('Mean Metric Value')
    mean(metric_Kur)
    




    %% Plot 2x2D scatter.  Strawman vs. Kuramoto.  And Strawman vs. Eigenvector.
    if(1)
        
        disp('Plotting and Saving 2D Histogram of DivMarg Scatter Plots')
        tic

        Xaxis = linspace(0,1,20);
        Yaxis = linspace(0,1,20);
        
        % Plot
        hSc = figure;   
        %
        hKS = subplot(121); 
        %scatter( SMDivMarg_C , KurDivMarg ), hold on 
        Hist2 = hist2d([SMDivMarg_C ; KurDivMarg],20,20,[0 1],[0 1]); close
        imagesc(Xaxis, Yaxis, Hist2./sum(Hist2(:))), hold on, axis xy
        plot([0 1], [0 1],'w--','LineWidth',1.5)
        axis([0 1 0 1])
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg of Coupled Oscillator Model','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg of Pixel Contrast Model','FontSize',18,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{+}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{-}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.1,'\color{black}{0}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,0.95,['\color{red}{<metric>=',num2str(mean(metric_Kur),2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        cb=colorbar('NorthOutside');
        set(cb,'FontSize',16,'FontWeight','Bold');
        %
        hES = subplot(122); 
        %scatter( SMDivMarg_L , EigDivMarg  ), hold on
        Hist2 = hist2d([SMDivMarg_L ; EigDivMarg],20,20,[0 1],[0 1]); close
        imagesc(Xaxis, Yaxis, Hist2./sum(Hist2(:))), hold on, axis xy
        plot([0 1], [0 1],'w--','LineWidth',1.5)
        axis([0 1 0 1])
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg of Evec1 Viz Computation','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg of Pixel Contrast Model','FontSize',18,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{+}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{-}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.1,'\color{black}{0}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,0.95,['\color{red}{<metric>=',num2str(mean(metric_Eig),2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        cb=colorbar('NorthOutside');
        set(cb,'FontSize',16,'FontWeight','Bold');
        %
        annotation('textbox', [0 0.9 1 0.1],'String', ...
            [fileTypeStr,' ',fileSizeStr,' : ',paramsStr,' : ','# Pts= ',num2str(numel(ImgPtchID))], ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')
        %
        disp(['Saving ',imgKurDir,'DivMarg_Hist2d_',fileSpecifier])
        saveGoodImg(hSc,[imgKurDir,'DivMarg_Hist2d_',fileSpecifier],sizeGoodIm)
        close(hSc)
        
        toc
        
    end


    %% Sort points according to metric on 2D Kuramoto Strawman Divisive Margin Space.
    
    % if mean(metric_Kur) > some Threshold Value then do this stuff below.
    if( mean(metric_Kur) > 0.08 ) 
        
        num2viz = 3;  % number of points to visualize.
    
        for KK = 1:numel(TBM) % types of points to visualize (top, bottom, mean)
        
            % Keep track of which image patches I visualize in the num2viz loop below. Then I can avoid visualizing the same 
            % one several times when it occurs for a slightly different ground truth.
            alreadyVizzed = {};  % cell array to hold names of image patches we already visualized.
            vizNum = 1;          % counter to index into alreadyVizzed array


            % index into top or bottom or mean scatter points / image patches
            eval(['I=',TBM{KK},';'])


            for i = 1:num2viz

                yy = find(strcmp(alreadyVizzed,ImgPtchID{I(i)}));

                if ~isempty(yy) % if we havent already visualized this image patch

                    disp(['Already visualized ',ImgPtchID{I(i)}])

                else
                    
                    % File names for image and movie looking at particular examples on single image patches
                    figKurName = [imgKurDir,'DivMarg_Scatter_2D_Ts',Ts,'_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_',TBM{KK},'_',num2str(num2viz)];
                    figEigName = [imgEigDir,'DivMarg_Scatter_2D_Ts',Ts,'_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_',TBM{KK},'_',num2str(num2viz)];
                    
                    
                    movName = [imgKurDir,'MovPhaseEvoL_tscale',Ts,'_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Kscale',Ks,'_',TBM{KK},'_',num2str(num2viz)];

                    if ( exist([figKurName,'.jpg'],'file') &  exist([figEigName,'.jpg'],'file') & ( ( phaseMovie & exist([movName,'.avi'],'file') ) | ~phaseMovie ) )
                        
                        disp(['File Already Exists: ',num2str(i),' / ',num2str(num2viz),' : ',figName])
                        % put ImgPtchID into the alreadyVizzed array after we made and saved plot.
                        alreadyVizzed{vizNum} = ImgPtchID{I(i)};
                        vizNum = vizNum+1;
                        
                    else

                        disp('Loading Mats')

                        tic

                        % (1). Load in Image Patch to plot it and ground truth.
                        ImPtch = load([patchesDir,ImgPtchID{I(i)},'.mat']);


                        % (2). Load in KurMC file and plot a.) ending solution and 
                        %      b.) Kur DivMarg time evolution along with Strawman DivMarg.
                        KurMC = load([dirKurMats,'KurMC_',ImgPtchID{I(i)},'_rM',rM,'_sD',sD,'_sP',sP,'_NF_60_',sW,'_kscale',Ks,'_tscale',Ts,'_runs1.mat']);
                        %
                        % DivMarg for Kuramoto at all times of simulation.
                        DivMarg_K = KurMC.MC{1}.DistAvgPW(:,1,ImgGtID(I(i))) ./ KurMC.MC{1}.DistAvgPW(:,2,ImgGtID(I(i)));


                        % (3). Load in Evecs file to plot Eigenvector solution and Eig DivMarg.
                        % EigDivMarg(I(i)) % Note: Not dealing with Eig right now...
                        EigMC = load([dirEigMats,'Evecs_',ImgPtchID{I(i)},'_rM',rM,'_sD',sD,'_sP',sP,'.mat']);
                        %
                        % DivMarg for Kuramoto at all times of simulation.
                        DivMarg_E1v = EigMC.MC.ev1v_l{ImgGtID(I(i))}.DistAvgPW(1) ./ EigMC.MC.ev1v_l{ImgGtID(I(i))}.DistAvgPW(2);
                        DivMarg_E1 = EigMC.MC.ev1_l{ImgGtID(I(i))}.DistAvgPW(1) ./ EigMC.MC.ev1_l{ImgGtID(I(i))}.DistAvgPW(2);
                        
                        
                        
                        toc
                        
                        
                        keyboard
                        
                        
                        

                        % (4). Plot all that shit.
                        if ~exist([figKurName,'.jpg'],'file')
                            disp('Plotting')

                            tic
                            
                            % KurPhaseInit = visKurPhase_inBone(KurMC.netParams.im, reshape(KurMC.metaCluster.phaseAtClk(:,1),KurMC.netParams.Ndims(1),KurMC.netParams.Ndims(2)));
                            
                            

                            H = figure; 
                            subplot(121), 
                            %scatter( SMDivMarg_C , KurDivMarg ), hold on 
                            Hist2 = hist2d([SMDivMarg_C ; KurDivMarg],20,20,[0 1],[0 1]); close
                            imagesc(Xaxis, Yaxis, Hist2./sum(Hist2(:))), hold on, axis xy
                            freezeColors
                            plot([0 1], [0 1],'w--','LineWidth',1.5)
                            axis([0 1 0 1])
                            set(gca,'FontSize',16,'FontWeight','Bold')
                            ylabel('DivMarg of Coupled Oscillator Model','FontSize',18,'FontWeight','Bold')
                            xlabel('DivMarg of Pixel Contrast Model','FontSize',18,'FontWeight','Bold')
                            text(0.9,0.1,'\color{green}{+}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
                            text(0.1,0.9,'\color{red}{-}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
                            text(0.1,0.1,'\color{black}{0}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
                            text(0.5,0.95,['\color{red}{<metric>=',num2str(mean(metric_Kur),2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
                            cb=colorbar('NorthOutside');
                            set(cb,'FontSize',16,'FontWeight','Bold');
                            cbfreeze
                            
                            % mean(Hist2(:))./sum(Hist2(:))
                            
                            
                            subplot(num2viz,4,3), % plot img
                            imagesc(KurMC.netParams.im)  % Kuramoto Seg
                            axis square
                            colormap('bone'), freezeColors
                            title('Img','FontSize',20,'FontWeight','Bold')
                            set(gca,'Xtick',[],'Ytick',[],'FontSize',16,'FontWeight','Bold')
%                             ylabel(['DM = ',num2str(DivMarg_K(end),3)],'FontSize',18,'FontWeight','Bold')
%                             colormap('hsv'), caxis([0 2*pi]), freezeColors
%                             cb=colorbar('SouthOutside'); 
%                             set(cb,'XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','$\pi$/2','$\pi$','3$\pi$/2','2$\pi$'})
%                             cbfreeze
                            
                            
                            subplot(num2viz,4,4), % plot kur
                            KurPhaseFinl = visKurPhase_inBone(KurMC.netParams.im, reshape(KurMC.metaCluster.phaseAtClk(:,end),KurMC.netParams.Ndims(1),KurMC.netParams.Ndims(2)));
                            imagesc(KurPhaseFinl)  % Kuramoto Seg
                            axis square
                            colormap('bone'), freezeColors
                            title('Kur','FontSize',20,'FontWeight','Bold')
                            set(gca,'Xtick',[],'Ytick',[],'FontSize',16,'FontWeight','Bold')
%                             ylabel(['DM = ',num2str(DivMarg_K(end),3)],'FontSize',18,'FontWeight','Bold')
%                             colormap('hsv'), caxis([0 2*pi]), freezeColors
%                             cb=colorbar('SouthOutside'); 
%                             set(cb,'XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','$\pi$/2','$\pi$','3$\pi$/2','2$\pi$'})
%                             cbfreeze
                            


















                            subplot(242), imagesc(ImPtch.im)                                           % Image Patch
                            axis square
                            title('Image Patch / Strawman','FontSize',20,'FontWeight','Bold')
                            set(gca,'Xtick',[],'Ytick',[],'FontSize',16,'FontWeight','Bold')
                            %ylabel(['DM = ',num2str(DivMarg_S,3)],'FontSize',18,'FontWeight','Bold')
                            colormap('jet'), freezeColors
                            colorbar('SouthOutside'), cbfreeze


                            subplot(243), imagesc(ImPtch.gT{ImgGtID(I(i))})                            % Ground Truth
                            axis square
                            title(['Ground Truth #',num2str(ImgGtID(I(i)))],'FontSize',20,'FontWeight','Bold')
                            set(gca,'Xtick',[],'Ytick',[],'FontSize',16,'FontWeight','Bold')
                            colormap('bone'), freezeColors
                            cb=colorbar('SouthOutside'); set(cb,'Visible', 'off'), cbfreeze


                            subplot(244), imagesc(ImPtch.imFull), hold on                              % Full Image
                            %axis square
                            plot([ImPtch.pach.xpbeg, ImPtch.pach.xpbeg],[ImPtch.pach.ypbeg, ImPtch.pach.ypfin],'r')
                            plot([ImPtch.pach.xpfin, ImPtch.pach.xpfin],[ImPtch.pach.ypbeg, ImPtch.pach.ypfin],'r')
                            plot([ImPtch.pach.xpbeg, ImPtch.pach.xpfin],[ImPtch.pach.ypbeg, ImPtch.pach.ypbeg],'r')
                            plot([ImPtch.pach.xpbeg, ImPtch.pach.xpfin],[ImPtch.pach.ypfin, ImPtch.pach.ypfin],'r') 
                            x = ImgPtchID{I(i)};
                            x(x=='_')=' ';
                            title({'Full Image:',x},'FontSize',20,'FontWeight','Bold')
                            set(gca,'Xtick',[],'Ytick',[],'FontSize',16,'FontWeight','Bold')
                            colormap('bone'), freezeColors


                            subplot(234), hold on
                            axis square
                            scatter( SMDivMarg_C,KurDivMarg, 'k.' )                    % all data points
                            % find other points in the 2x2D scatter plot belonging to same image patch (different ground truth or ...)
                            xx = find(strcmp(ImgPtchID,ImgPtchID{I(i)}));
                            % ImgPtchID{xx}
                            % ImgGtID(xx)
                            scatter( SMDivMarg_C(xx),KurDivMarg(xx), 200, 'g.' )       % other data points from this same img patch
                            scatter( SMDivMarg_C(I(i)),KurDivMarg(I(i)), 300, 'r.' )   % this data point (pertaining to kuramoto phase image)
                            plot([0 1.2], [0 1.2],'r--','LineWidth',1.5)
                            axis ([0 1.2 0 1.2])
                            plot([0 1], [1 1],'r--','LineWidth',1.5)
                            plot([1 1], [0 1],'r--','LineWidth',1.5)
                            set(gca,'FontSize',16,'FontWeight','Bold')
                            ylabel('Kuramoto DM','FontSize',18,'FontWeight','Bold')
                            xlabel('Strawman DM','FontSize',18,'FontWeight','Bold')
                            %text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','HorizontalAlignment','Right')
                            %text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold')
                            %text(0.2,0.25,'\color{cyan}{Neutral}','FontSize',16,'FontWeight','Bold','HorizontalAlignment','Center')
                            title(['metric: (this pt = ',num2str(metric_Kur(I(i)),2),') & (avg = ',num2str(mean(metric_Kur),2),')'],'FontSize',20,'FontWeight','Bold')


                            subplot(2,3,[5:6]), hold on                                                % DivMarg Kur & Straw
                            plot(KurMC.kurParams.tau*(1:numel(DivMarg_K)),DivMarg_K,'b','LineWidth',2)
                            %plot(KurMC.kurParams.tau*(1:numel(DivMarg_K)),repmat(DivMarg_S,size(DivMarg_K)),'r','LineWidth',2)
                            xlabel('Time (sec)','FontSize',18,'FontWeight','Bold')
                            ylabel('Divisive Margin','FontSize',18,'FontWeight','Bold')
                            legend('Kuramoto Osc','Pixel Strawman','FontSize',16,'FontWeight','Bold')
                            title(['Parameters: Tscale=',Ts,' | \sigma_P=',sP,' | Rmax=',rM,' | \sigma_D=',sD,' | \sigma_W=',sW,' | Kscale=',Ks],'FontSize',20,'FontWeight','Bold')
                            set(gca,'FontSize',16,'FontWeight','Bold')
                            
                            toc

                            disp('Saving')
                            tic

                            disp(['Saving ',num2str(i),' / ',num2str(num2viz),' : ',figName])
                            saveGoodImg(H,figName,sizeGoodIm)
                            close(H)
                            
                            toc
                            
                        end
                        
                        
                        % Make a video of phase evolution of oscillators
                        if(phaseMovie & ~exist([movName,'.avi'],'file') )
                            
                            disp('Making Movie File:')
                            tic
                        
                            makeMovie_phaseEvolution( movName, KurMC, ImPtch, ImgPtchID{I(i)}, ImgGtID(I(i)) )
                            
                            toc
                        
                        end
                        

                        % put ImgPtchID into the alreadyVizzed array after we made and saved plot.
                        alreadyVizzed{vizNum} = ImgPtchID{I(i)};
                        vizNum = vizNum+1;

                        
                    
                    end % if the jpg file has been (or hasnt) generated already.

                end % if we havent already visualized this patch

            end % loop through number of patches to visualize
        
        end % types of points to visualize (top, bottom, mean)

    end % if the 2D scatter metric exceeds a threshold

end % loop over all mat files (all parameter value combinations)