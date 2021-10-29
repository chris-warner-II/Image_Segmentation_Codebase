function compute_Eig_bDgradMetric_fixFiles(InputImages,PatchSize,MethodType,RM,SP)


% This function will loop through the Evecs files and will calculate mean
% and standard deviation of units within each cluster.  This will be used
% later for the d' sensitivity metric.


    [dirPre, sizeGoodIm] = onCluster;
    filesInDir = [dirPre,'output/Kuramoto/NetsFromImgs/',InputImages,'_',PatchSize,'/data/spectral/',MethodType,'/'];
    
    
    % Note: IsoDiff files named somewhat differently because they dont have SP parameter.
    if strmatch(MethodType,'IsoDiff')
        files = dir([filesInDir,'Evecs*rM',RM,'.mat']);
    else
        files = dir([filesInDir,'Evecs*rM',RM,'*_sP',SP,'.mat']);
    end
    
    
    filesOutDir = [dirPre,'output/Kuramoto/NetsFromImgs/',InputImages,'_',PatchSize,'/data/spectral_fixed/',MethodType,'/'];
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

        if(0) % NOT CHECKING FOR EXISTING VARIABLE NOW... exist('MC','var')
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

        
        im = netParams.im;
        ev1 = reshape(EVecsML(:,1), netParams.Ndims);
        ev2o = reshape(EVecsML(:,2), netParams.Ndims);
        ev3o = reshape(EVecsML(:,3), netParams.Ndims);

        ev2 = reshape(sum(EVecsML(:,1:2),2), netParams.Ndims);
        ev3 = reshape(sum(EVecsML(:,1:3),2), netParams.Ndims);
        ev2w = reshape(EVecsML(:,1).*EValsML(1) + EVecsML(:,2).*EValsML(2), netParams.Ndims);
        ev3w = reshape(EVecsML(:,1).*EValsML(1) + EVecsML(:,2).*EValsML(2) + EVecsML(:,3).*EValsML(3), netParams.Ndims);

        
        

        % Compute boundary based Gradient Metric
        [MC.im.F, MC.im.M, MC.im.S, MC.im.D] = compute_BoundaryGradientMetric(im, bDc_blurs, numel(netParams.gT), 0, 'Im Pix');
        

        [MC.ev1.F, MC.ev1.M, MC.ev1.S, MC.ev1.D] = compute_BoundaryGradientMetric(ev1, bDc_blurs, numel(netParams.gT), 0, 'Evec 1');


        [MC.ev2o.F, MC.ev2o.M, MC.ev2o.S, MC.ev2o.D] = compute_BoundaryGradientMetric(ev2o, bDc_blurs, numel(netParams.gT), 0, 'Evec 2'); 

        
        [MC.ev3o.F, MC.ev3o.M, MC.ev3o.S, MC.ev3o.D] = compute_BoundaryGradientMetric(ev3o, bDc_blurs, numel(netParams.gT), 0, 'Evec 3');
 

        [MC.ev2.F, MC.ev2.M, MC.ev2.S, MC.ev2.D] = compute_BoundaryGradientMetric(ev2, bDc_blurs, numel(netParams.gT), 0, 'Evecs 1-2');


        [MC.ev3.F, MC.ev3.M, MC.ev3.S, MC.ev3.D] = compute_BoundaryGradientMetric(ev3, bDc_blurs, numel(netParams.gT), 0, 'Evecs 1-3');


        [MC.ev2w.F, MC.ev2w.M, MC.ev2w.S, MC.ev2w.D] = compute_BoundaryGradientMetric(ev2w, bDc_blurs, numel(netParams.gT), 0, 'Evecs 1-2w');


        [MC.ev3w.F, MC.ev3w.M, MC.ev3w.S, MC.ev3w.D] = compute_BoundaryGradientMetric(ev3w, bDc_blurs, numel(netParams.gT), 0, 'Evecs 1-3w');


        % % %
        
        
        [MC.imP.F, MC.imP.M, MC.imP.S, MC.imP.D] = compute_BoundaryGradientMetric(im.*pi, bDc_blurs, numel(netParams.gT), 0, 'Im Pix');


        % save new fixed Evecs file ...
        save([filesOutDir,files(F).name],'EValsML', 'EVecPM', 'EVecsML','netParams','netflags','MC')
        
        clear EVecsML EValsML EVecPM MC
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