function compute_Kur_bDgradMetric_fixFiles(InputImages,PatchSize,MethodType,RM,KS)


% This function will loop through the Evecs files and will calculate mean
% and standard deviation of units within each cluster.  This will be used
% later for the d' sensitivity metric.


    [dirPre, sizeGoodIm] = onCluster;
    filesInDir = [dirPre,'output/Kuramoto/NetsFromImgs/',InputImages,'_',PatchSize,'/data/Kur_PIF_Fourier1/',MethodType,'/'];
    files = dir([filesInDir,'KurMC*','_rM',RM,'_*_ks',KS,'*.mat']);
    
    
    filesOutDir = [dirPre,'output/Kuramoto/NetsFromImgs/',InputImages,'_',PatchSize,'/data/Kur_PIF_Fourier1_fixed/',MethodType,'/'];
    if ~exist(filesOutDir,'dir')
        mkdir(filesOutDir)
    end

    disp([num2str(numel(files)),' Files To Process.'])

    for F = 1:numel(files)

        disp([num2str(F),' / ',num2str(numel(files)),' - ',files(F).name])
        
        
        % Try to load in the file
        try
            load([filesInDir,files(F).name])
        catch
            disp('File Corrupted or Missing.  I dunno.  Moving on.')
            continue
        end
        
        
        
        % If variables I am planning to create already exist, do no more.
        if(0)% NOT CHECKING THIS NOW.. exist('MC','var')
            disp('This file already updated.  Moving it directly to the new directory...')
            [filesInDir,files(F).name]
            movefile([filesInDir,files(F).name],filesOutDir)
            clear MC
            continue
        end
        
        
        
        try
            bD = netParams.bD;
        catch
            
            % If bD doesnt exist in netParams, open ground truth file in images directory and add bD to netParams.
            ind = strfind(files(F).name,'_rM')-1;
            GTFile = load([dirPre,'images/',InputImages,'/',PatchSize,'/',files(F).name(7:ind)]);
            netParams.bD = GTFile.bD;
            bD = netParams.bD;
        end


        % Blur ground truth boundary lines for analysis
        MC.bDc_blurs_info = {'(bD>0)','(bD>0)blur_binom2','(bD>0)blur_binom3','(bD>1)','(bD>1)blur_binom2','(bD>1)blur_binom3','bD'};
        bDc_blurs(:,:,1) = single(bD>0);
        bDc_blurs(:,:,2) = single(logical( blur( single(bD>0),1,'binom2' ) ) );
        bDc_blurs(:,:,3) = single(logical( blur( single(bD>0),1,'binom3' ) ) );
        bDc_blurs(:,:,4) = single(bD>1);
        bDc_blurs(:,:,5) = single(logical( blur( single(bD>1),1,'binom2' ) ) );
        bDc_blurs(:,:,6) = single(logical( blur( single(bD>1),1,'binom3' ) ) ); % note: 'power' increases with increased blurring
        bDc_blurs(:,:,7) = bD;
        
        % Compute boundary based Gradient Metric
        X = visKurPhase_inHSV(netParams.im, reshape(metaCluster.phaseAtClk(:,end),netParams.Ndims));
        [MC.F, MC.M, MC.S, MC.D] = compute_BoundaryGradientMetric(X, bDc_blurs, numel(netParams.gT), 1, MethodType);
        


        % save Evecs file right on top of the old one...
        save([filesOutDir,files(F).name],'metaCluster', 'kurParams', 'kurflags','netParams','netflags','MC')
        
        
        clear MC
        delete([filesInDir,files(F).name])
        

    end % loop over files
    
    disp('This function completed successfully.')
    clock

end % main function end
