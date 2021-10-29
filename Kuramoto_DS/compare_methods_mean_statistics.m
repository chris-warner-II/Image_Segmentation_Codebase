function compare_methods_mean_statistics(fileType, fileSize, methodType, rMx,sDx,sPx,sWx,Tsx,Ksx)

    % syntax:  compare_methods_mean_statistics(fileType, fileSize, methodType, rM,sD,sP,sW,Ts,Ks)
    %
    %          This function will take in the DistsPW_Results_allPatches*.mat
    %          files from the DistPairs directory. It will output a plot of
    %          mean (across all image patches) Divisive Margin & DistSep for
    %          scatter plots that can be made using scatter_plot_DistPairs.m.
    %
    %       What I want you to do with it right now:
    %
    %           (1). Input single methodType 
    %                Input array of rM values
    %                Input array of sP value
    %                Input single sD, sW, Ts, Ks.
    %                This will plot rM vs. sP vs. mean(D) where D is DivMarg
    %                and DistSep in different subplots.
    %
    %
    % fileType = 'BSDS_patch';
    % fileSize = '51x51_ds1';
    % methodType = {'Mod_SKHAdj'};
    %
    % rMx = (1:10,inf)
    % sDx = inf
    % sPx = (0.1,0.2,0.3,0.4)
    % sWx = 0
    % Tsx = 1
    % Ksx = 3
    %




    %% Load in mat file and set up directories to different data.
    [dirPre,sizeGoodIm] = onCluster;




    for K = 1:numel(methodType)
    

        % output directory to save the RDdata mat files into.
        outDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/RateDistPlots/',methodType{K},'/'];
        if ~exist(outDir)
            mkdir(outDir)
        end

        % mat file saved from exp_SepVPar_avgOverImgPatches.
        matFilesDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/DistPairs/',methodType{K},'/'];
        if ~exist(matFilesDir)
            mkdir(matFilesDir)
        end


        % 
        imgFilesDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/imgs/DistPairs/',methodType{K},'/'];
        if ~exist(imgFilesDir)
            mkdir(imgFilesDir)
        end


        % directory to image patches extracted from BSDS.
        patchesDir = [dirPre,'images/',fileType,'/',fileSize,'/'];

        % directory to KurMC files and to Kur_metaSummary files (think I need the KurMC ones)
        dirKurMats = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/Kur_PIF_Fourier1/',methodType{K},'/'];

        % directory to the Evecs files
        dirEigMats = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/spectral/',methodType{K},'/'];


        % For Plot title displays and whatnot later...
        fileTypeStr = fileType;
        fileTypeStr(fileTypeStr=='_')=' ';

        fileSizeStr = fileSize;
        fileSizeStr(fileSizeStr=='_')=' ';

        methodTypeStr = methodType{K};
        methodTypeStr(methodTypeStr=='_')=' ';



        % Making Labels and Legends for plots.
        for i = 1:numel(rMx)
            rM_lgnd{i} =  num2str(rMx(i));
        end
        %
        for j = 1:numel(sPx)
            sP_lgnd{j} =  num2str(sPx(j));
        end



         % Preallocate memory for array of mean values.
         segMethods = {'kur';'im';'ev1';'ev2o';'ev3o'}; % ;'ev2';'ev3';'ev2w';'ev3w'};

    %     
    %     %
    %     for k = 1:numel(segMethods)
    %         DivMargMean.(segMethods{k}) = nan(numel(rMx),numel(sPx));
    %         PmtricMean.(segMethods{k}) = nan(numel(rMx),numel(sPx));
    %         %
    %         DivMargRat.(segMethods{k}) = nan(numel(rMx),numel(sPx));
    %         PmtricRat.(segMethods{k}) = nan(numel(rMx),numel(sPx));
    %     end



        %% Loop through parameter values (preset by user) to extract mean 
        % (across all image patches processed) DivMarg or DistSub.
        for j = 1:numel(sPx)

            for i = 1:numel(rMx)


                % Turn parameter values into strings  
                rM = num2str(rMx(i)); % '1';
                sD = num2str(sDx); % 'Inf';
                sP = num2str(sPx(j)); % '0p2';
                sP(sP=='.')='p';
                Ts = num2str(Tsx); % '1';
                Ks = num2str(Ksx); % '300';
                sW = num2str(sWx); % '0';

                fileSpecifier = ['_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts];

                paramsStr = fileSpecifier;
                paramsStr(paramsStr=='_')=' ';




                files = dir([matFilesDir,'Clustering_Results_allPatches',fileSpecifier,'_files*.mat']);
                try
                    load([matFilesDir,files(end).name]);
                catch
                    disp(['File not there : Clustering_Results_allPatches',fileSpecifier])
                    continue
                end



                if ~exist('RateDistSig','var')
                    RateDistSig = [0,0.01,0.1,0.3,0.5,0.7,0.9];
                end





                for S = 1:numel(segMethods) % number of methods (kur, im, ev1, ...)

                    for R = 1:7 % number of rate distortion values.

                        x = d_prime.(segMethods{S})(:,R);
                        %
                        infs = find(isinf(x));
                        nans = find(isnan(x));
                        fins = find(~isinf(x)&~isnan(x));
                        %
                        numNans(S,R) = numel(nans);
                        numInf(S,R) = numel(infs);
                        meanFins(S,R) = mean(x(fins));
                        
                        y = d_prime_wtd.(segMethods{S})(:,R);
                        %
                        infs = find(isinf(y));
                        nans = find(isnan(y));
                        fins = find(~isinf(y)&~isnan(y));
                        %
                        numNansWt(S,R) = numel(nans);
                        numInfWt(S,R) = numel(infs);
                        meanFinsWt(S,R) = mean(x(fins));
                        

                    end

                end




                % [numInf;meanFins]


    %             figure, plot( RateDistSig(2:end),  meanFins(:,2:end)' ,'LineWidth',2)
    %             legend(segMethods)
    %             title([fileTypeStr,' ',fileSizeStr,' ',methodTypeStr,' ',paramsStr])
    %             xlabel('Rate Distortion \sigma')
    %             ylabel('d''')



                save([outDir,'RDdata_',methodType{K},fileSpecifier,'.mat'],'RateDistSig','meanFins','numInf','numNans','meanFinsWt','numInfWt','numNansWt','segMethods')





    %             % Can show scatter plots and histograms of performance across image patches for diagnosis.
    %             if (0)
    %                 [H, DS, DM] = scatter_plot_DistPairs(KurDistsPW,SMDistsPW_C,dirKurMats,patchesDir,'BSDS_patch','51x51_ds1','Kur',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID,ImgGtID,1);
    %                 % hmm, maybe have to change this SMDistsPW_C to be SMDistsPW_L.
    %             end







    %             % Fill arrays with statistics (across images patches) for different parameter combinations.
    %             %segMethods = fieldnames(DivMarg);
    %             for k = 1:numel(segMethods)
    %                 DivMargMean.(segMethods{k})(i,j) = mean( DivMarg.(segMethods{k}) );
    %                 PmtricMean.(segMethods{k})(i,j)  = mean(  Pmtric.(segMethods{k}) );
    %                 %
    %                 better = numel(find(DivMarg.(segMethods{k})>0)); 
    %                 worse = numel(find(DivMarg.(segMethods{k})<0));
    %                 DivMargRat.(segMethods{k})(i,j) = better ./ (better+worse);
    %                 %
    %                 better = numel(find( Pmtric.(segMethods{k}) > 0 ) );
    %                 worse = numel(find( Pmtric.(segMethods{k}) < 0 ) );
    %                 PmtricRat.(segMethods{k})(i,j) = better ./ (better+worse);
    %             end




            end  % loop over rM

        end% loop over sP


    
    
    end
    
    
    
%     % Set explicitly anything weighted by eigenvalues to Nan because I am
%     % not doing that right now.  Have to go back and do it at some point.
%     DivMargMean.Eig2w = nan(numel(rMx),numel(sPx));
%     %DivMargMean.Eig2wv = nan(numel(rMx),numel(sPx));
%     DivMargMean.Eig3w = nan(numel(rMx),numel(sPx));
%     %DivMargMean.Eig3wv = nan(numel(rMx),numel(sPx));
%     %
%     PmtricMean.Eig2w = nan(numel(rMx),numel(sPx));
%     %PmtricMean.Eig2wv = nan(numel(rMx),numel(sPx));
%     PmtricMean.Eig3w = nan(numel(rMx),numel(sPx));
%     %PmtricMean.Eig3wv = nan(numel(rMx),numel(sPx));
%     %
%     DivMargRat.Eig2w = nan(numel(rMx),numel(sPx));
%     %DivMargRat.Eig2wv = nan(numel(rMx),numel(sPx));
%     DivMargRat.Eig3w = nan(numel(rMx),numel(sPx));
%     %DivMargRat.Eig3wv = nan(numel(rMx),numel(sPx));
%     %
%     PmtricRat.Eig2w = nan(numel(rMx),numel(sPx));
%     %PmtricRat.Eig2wv = nan(numel(rMx),numel(sPx));
%     PmtricRat.Eig3w = nan(numel(rMx),numel(sPx));
%     %PmtricRat.Eig3wv = nan(numel(rMx),numel(sPx));
    
    
    
    
    
    
    
    
    
    
%     %% Loop thru different segMethods to find max D and parameter values that lead to it.
%     
%     %segMethods = fieldnames(DivMargMean);
%     
%     for k = 1:numel(segMethods)
%         [M,x,y] = find_max( DivMargMean.(segMethods{k}) );
%         try
%             maxDMM.(segMethods{k})(1:3) = [M,rMx(x),sPx(y)];
%         catch
%             maxDMM.(segMethods{k})(1:3) = [M,x,y];
%         end
%         %
%         [M,x,y] = find_max( PmtricMean.(segMethods{k}) );
%         try
%             maxDSM.(segMethods{k})(1:3) = [M,rMx(x),sPx(y)];
%         catch
%             maxDSM.(segMethods{k})(1:3) = [M,x,y];
%         end
%         %
%         [M,x,y] = find_max( DivMargRat.(segMethods{k}) );
%         try
%             maxDMR.(segMethods{k})(1:3) = [M,rMx(x),sPx(y)];
%         catch
%             maxDMR.(segMethods{k})(1:3) = [M,x,y];
%         end
%         %
%         [M,x,y] = find_max( PmtricRat.(segMethods{k}) );
%         try
%             maxDSR.(segMethods{k})(1:3) = [M,rMx(x),sPx(y)];
%         catch
%             maxDSR.(segMethods{k})(1:3) = [M,x,y];
%         end
%     end
    
    
    
    

    
    
    
    
%     %% DO SOME PLOTS.
%     if(1)
%         
%         %segMethods = fieldnames(DivMargMean);
%     
%         for k = 1:numel(segMethods)
%         
%             % Plot mean D vs Rmax & sig_P
%             if(1) % ( max(DivMargMean.(segMethods{k})(:)) > 0 ) | ( max(PmtricMean.(segMethods{k})(:)) > 0 ) )
%                 Hf = plot_D_vs_rM_vs_sP(DivMargMean.(segMethods{k}),PmtricMean.(segMethods{k}),segMethods{k},methodType{K},sP_lgnd,rM_lgnd);
%                 saveGoodImg(Hf,[imgFilesDir,'Dmean_vs_rM_vs_sP_',methodType{K},'_',segMethods{k}],sizeGoodIm)
%                 close(Hf)
%             end
%             
%             % Plot ratio of D vs Rmax & sig_P
%             if(1)
%                 Hf = plot_D_vs_rM_vs_sP(DivMargRat.(segMethods{k}),PmtricRat.(segMethods{k}),segMethods{k},methodType{K},sP_lgnd,rM_lgnd);
%                 saveGoodImg(Hf,[imgFilesDir,'Drat_vs_rM_vs_sP_',methodType{K},'_',segMethods{k}],sizeGoodIm)
%                 close(Hf)
%             end
%         end
% 
%     end
    

    
    % How to include max DMM, DSM, DMR, DSR ?   Save them into mat files.
    % Then, plot netMethod vs segMethod vs max performance noting rM & sP values.
    %save([matFilesDir,'best_perf_values_rM_vs_sP'],'maxDMM','maxDSM','maxDMR','maxDSR','DivMargMean',...
    %    'PmtricMean','DivMargRat','PmtricRat','segMethods', 'rMx','sDx','sPx','sWx','Tsx','Ksx')
    
    
    
    
    
    
    
    
    
    
    
    
    
    disp('Function Successfully Completed.')
    clock

end % end main function




% % Make the actual plot of rM vs. sP vs. D
% function [Hf] = plot_D_vs_rM_vs_sP(DM,DS,segMeth,netMeth,sP_lgnd,rM_lgnd)
%     Hf = figure;
%     subplot(211), hold on
%     plot(DM,'LineWidth',2)
%     plot([1,numel(rM_lgnd)],[0, 0],'k--')
%     axis([1 numel(rM_lgnd) -0.5 1])
%     set(gca,'XTick',[],'FontSize',18,'FontWeight','Bold')
%     ylabel('\Delta DivMarg','FontSize',18,'FontWeight','Bold')
%     hl=legend(sP_lgnd,'Location','Best');
%     set(get(hl,'Title'),'String','\sigma_p','FontSize',18,'FontWeight','Bold')
%     title([netMeth,' :: ',segMeth],'FontSize',20,'FontWeight','Bold')
%     %
%     [M,x,y] = find_max(DM);
%     scatter(x(1),M,100,'ro','LineWidth',2)
% 
%     subplot(212), hold on
%     plot(DS,'LineWidth',2)
%     plot([1,numel(rM_lgnd)],[0, 0],'k--')
%     axis([1 numel(rM_lgnd) -0.5 1])
%     set(gca,'XTick',[1:numel(rM_lgnd)],'XTickLabel',rM_lgnd,'FontSize',18,'FontWeight','Bold')
%     xlabel('Rmax','FontSize',18,'FontWeight','Bold')
%     ylabel('\Delta P-metric','FontSize',18,'FontWeight','Bold')
%     %
%     [M,x,y] = find_max(DS);
%     scatter(x(1),M,100,'ro','LineWidth',2)
% end
% 
% % Find maximum value and location of it in a 2D array.
% function [M,x,y] = find_max(input)
% 
%     if( isempty(find(~isnan(input))) )
%         M=nan; x = nan; y = nan;
%     else
%         M = max(input(:));
%         [x,y] = find( input == M);
%         x = x(1); y = y(1); 
%     end
% end

