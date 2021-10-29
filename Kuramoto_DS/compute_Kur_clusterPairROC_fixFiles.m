function compute_Kur_clusterPairROC_fixFiles(InputImages,MethodType,PTCH,RM,SP)


% This function will loop through the Evecs files and will calculate mean
% and standard deviation of units within each cluster.  This will be used
% later for the d' sensitivity metric.


    [dirPre, sizeGoodIm] = onCluster;
    filesInDir = [dirPre,'output/Kuramoto/NetsFromImgs/',InputImages,'_101x101_ds1/data/Kur_PIF_Fourier1/',MethodType,'/'];
    files = dir([filesInDir,'KurMC*ptch',PTCH,'_rM',RM,'*_sP',SP,'*.mat']);
    
    
    filesOutDir = [dirPre,'output/Kuramoto/NetsFromImgs/',InputImages,'_101x101_ds1/data/Kur_PIF_Fourier1_fixed/',MethodType,'/'];
    if ~exist(filesOutDir,'dir')
        mkdir(filesOutDir)
    end

    disp([num2str(numel(files)),' Files To Process.'])

    for F = 1:numel(files)

        disp([num2str(F),' / ',num2str(numel(files)),' - ',files(F).name])
        try
            load([filesInDir,files(F).name])
        catch
            disp('File Corrupted or Missing.  I dunno.  Moving on.')
            continue
        end
        
        if exist('AUC_ROC_1D','var')
            disp('This file already updated.  Moving it directly to the new directory...')
            [filesInDir,files(F).name]
            movefile([filesInDir,files(F).name],filesOutDir)
            clear idCluster sizeCluster AUC_ROC_1D AUC_ROC_GS
            continue
        end

        im = netParams.im;
        kur = reshape(metaCluster.phaseAtClk(:,end), netParams.Ndims);
        gT = netParams.gT;
        
        
        [AUC_ROC_1D.kur, AUC_ROC_GS.kur] = calc_ClusterPairROC_B(gT,kur,1);
        [AUC_ROC_1D.im]                  = calc_ClusterPairROC_B(gT,im,0);
        
        
        
        
        for i = 1:numel(gT) % loop thru ground truths
            C = unique(gT{i});
            for j = 1:numel(C)
                idCluster{i}(j) = C(j);
                sizeCluster{i}(j) = numel(find(gT{i}==C(j)));
            end
        end
        
        
        
        

        % save Evecs file right on top of the old one...
        save([filesOutDir,files(F).name],'metaCluster', 'kurParams', 'kurflags','netParams','netflags', 'idCluster', 'sizeCluster', 'AUC_ROC_1D', 'AUC_ROC_GS')
        
        
        clear idCluster sizeCluster AUC_ROC_1D AUC_ROC_GS
        delete([filesInDir,files(F).name])
        

    end % loop over files
    
    disp('This function completed successfully.')
    clock

end % main function end
