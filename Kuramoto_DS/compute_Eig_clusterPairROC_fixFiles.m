function compute_Eig_clusterPairROC_fixFiles(InputImages,MethodType,PTCH,RM,SP)


% This function will loop through the Evecs files and will calculate mean
% and standard deviation of units within each cluster.  This will be used
% later for the d' sensitivity metric.


    [dirPre, sizeGoodIm] = onCluster;
    filesInDir = [dirPre,'output/Kuramoto/NetsFromImgs/',InputImages,'_51x51_ds1/data/spectral/',MethodType,'/'];
    files = dir([filesInDir,'Evecs*ptch',PTCH,'_rM',RM,'*_sP',SP,'.mat']);
    
    
    filesOutDir = [dirPre,'output/Kuramoto/NetsFromImgs/',InputImages,'_51x51_ds1/data/spectral_fixed/',MethodType,'/'];
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
            clear idCluster sizeCluster AUC_ROC_1D
            continue
        end

        im = netParams.im;
        ev1 = Normlze(reshape(EVecsML(:,1), netParams.Ndims));
        ev2o = Normlze(reshape(EVecsML(:,2), netParams.Ndims));
        ev3o = Normlze(reshape(EVecsML(:,3), netParams.Ndims));

        ev2 = Normlze(reshape(sum(EVecsML(:,1:2),2), netParams.Ndims));
        ev3 = Normlze(reshape(sum(EVecsML(:,1:3),2), netParams.Ndims));
        ev2w = Normlze(reshape(EVecsML(:,1).*EValsML(1) + EVecsML(:,2).*EValsML(2), netParams.Ndims));
        ev3w = Normlze(reshape(EVecsML(:,1).*EValsML(1) + EVecsML(:,2).*EValsML(2) + EVecsML(:,3).*EValsML(3), netParams.Ndims));

        gT = netParams.gT;

        % This function will compute mean cluster centers & standard deviations as well as 95% confidence intervals of each
        disp('Computing ROC Curve and Area Underneath for {im,ev1,ev2o,ev3o,ev2,ev3,ev2w,ev3w}')                            
        [AUC_ROC_1D.im] = calc_ClusterPairROC_B(gT,im,0);
        [AUC_ROC_1D.ev1] = calc_ClusterPairROC_B(gT,ev1,0);
        [AUC_ROC_1D.ev2o] = calc_ClusterPairROC_B(gT,ev2o,0);
        [AUC_ROC_1D.ev3o] = calc_ClusterPairROC_B(gT,ev3o,0);

        [AUC_ROC_1D.ev2] = calc_ClusterPairROC_B(gT,ev2,0);
        [AUC_ROC_1D.ev3] = calc_ClusterPairROC_B(gT,ev3,0);
        [AUC_ROC_1D.ev2w] = calc_ClusterPairROC_B(gT,ev2w,0);
        [AUC_ROC_1D.ev3w] = calc_ClusterPairROC_B(gT,ev3w,0);

        for i = 1:numel(gT) % loop thru ground truths
            C = unique(gT{i});
            for j = 1:numel(C)
                idCluster{i}(j) = C(j);
                sizeCluster{i}(j) = numel(find(gT{i}==C(j)));
            end
        end


        % save new fixed Evecs file ...
        save([filesOutDir,files(F).name],'EValsML', 'EVecPM', 'EVecsML','netParams','netflags','idCluster', 'sizeCluster', 'AUC_ROC_1D')
            
        
        clear idCluster sizeCluster AUC_ROC_1D
        delete([filesInDir,files(F).name]) % delete the old one
        

    end % loop over files
    
    disp('THis function completed successfully.')
    clock


end % main function end



%% 

% Renormalize Eigenvector to be between 0 & 1 (same as I did for image
% patch) so that I know its range for further processing.
function [x] = Normlze(x)

    x = x - min(x(:));
    x = x ./max(x(:));

end