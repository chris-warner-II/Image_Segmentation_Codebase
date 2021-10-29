function plot_netMethods_vs_segMethods_best_perf(fileType, fileSize, methodType)

    % syntax:  compare_methods_mean_statistics(fileType, fileSize, methodType, rM,sD,sP,sW,Ts,Ks)
    %
    %          This function will take in the best_perf_values_rM_vs_sP.mat
    %          files that sit in BSDS_patch_size/data/DistPairs/netMethod
    %          directory. That file contains maximum value of 4 performance
    %          metrics {(1). DivMargMean, (2). DistSubMean, (3). DivMargRat,
    %          (4). DistSubRat} and the {rM,sP} parameter values at which
    %          they occurred.  
    %
    %          This function will generate 4 separate bar plots of those 
    %          maximum performance values for different combinations of
    %          {netMethod, segMethod} parameters.


    %% Load in mat file and set up directories to different data.
    [dirPre,sizeGoodIm] = onCluster;

    % For Plot title displays and whatnot later...
    fileTypeStr = fileType;
    fileTypeStr(fileTypeStr=='_')=' ';

    fileSizeStr = fileSize;
    fileSizeStr(fileSizeStr=='_')=' ';


    % Directory to save output images into.
    imgFilesDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/imgs/compare_best_perf/'];
    if ~exist(imgFilesDir,'dir')
        mkdir(imgFilesDir)
    end
    
    
    for k = 1:numel(methodType)
    
    
        % mat file saved from exp_SepVPar_avgOverImgPatches.
        matFilesDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/DistPairs/',methodType{k},'/'];
    
    
        load([matFilesDir,'best_perf_values_rM_vs_sP.mat']);
        
        for j = 1:numel(segMethods)
        
            % compile arrays {netMethods x segMethods} of max performance values
            comp_max_DMM(k,j) = maxDMM.(segMethods{j})(1); 
            comp_max_DMM_rM(k,j) = maxDMM.(segMethods{j})(2); 
            comp_max_DMM_sP(k,j) = maxDMM.(segMethods{j})(3); 
            %
            comp_max_DSM(k,j) = maxDSM.(segMethods{j})(1);
            comp_max_DSM_rM(k,j) = maxDSM.(segMethods{j})(2); 
            comp_max_DSM_sP(k,j) = maxDSM.(segMethods{j})(3); 
            %
            comp_max_DMR(k,j) = maxDMR.(segMethods{j})(1);
            %
            comp_max_DSR(k,j) = maxDSR.(segMethods{j})(1); 

        end

    end

    
    H = makeBarPlots(comp_max_DMM, comp_max_DMM_rM, comp_max_DMM_sP, 'Max DivMarg', methodType, segMethods);
    saveGoodImg(H,[imgFilesDir,'Bar_max_DMM_netMethod_vs_segMethod'],sizeGoodIm)
    close(H)
    %
    H = makeBarPlots(comp_max_DSM, comp_max_DSM_rM, comp_max_DSM_sP, 'Max P-metric', methodType, segMethods);
    saveGoodImg(H,[imgFilesDir,'Bar_max_PMM_netMethod_vs_segMethod'],sizeGoodIm)
    close(H)
    %
%     H = makeBarPlots(comp_max_DMR, 'Ratio DivMarg', methodType, segMethods);
%     saveGoodImg(H,[imgFilesDir,'Bar_max_DMR_netMethod_vs_segMethod'],sizeGoodIm)
%     close(H)
%     %
%     H = makeBarPlots(comp_max_DSR, 'Ratio DistSub', methodType, segMethods);
%     saveGoodImg(H,[imgFilesDir,'Bar_max_DSR_netMethod_vs_segMethod'],sizeGoodIm)
%     close(H)
    

    disp('Function Finished Properly')
    clock



    
end
    
    
% Make Bar Plots.
function [H] = makeBarPlots(array,rM,sP,metric_label,netMeth,segMeth)
    
    cb = colormap('jet');
    [x,y] = find(array>0);
    
    H=figure; 
    subplot(211), hold on
    bar(array)
    for k=1:numel(x)
        text(x(k)-0.45+0.8*y(k)./size(array,2), array(x(k),y(k)), ['(',num2str(rM(x(k),y(k))),' , ',num2str(sP(x(k),y(k))),')'],...
            'VerticalAlignment','Middle','HorizontalAlignment','Left','FontSize',12,'Color',cb(round(y(k)./size(array,2).*size(cb,1)),:),'Rotation',90)
    end
    text(1-0.5, array(1,1), ['(r_M , \sigma_P)'],'VerticalAlignment','Middle','HorizontalAlignment','Left','FontSize',14,'FontWeight','Bold','Rotation',90)
    plot([0 numel(netMeth)+1],[array(1,1), array(1,1)],'k--')
    legend(segMeth,'Location','EastOutside')
    set(gca,'XTick',[1:numel(netMeth)],'XTickLabel',netMeth,'FontSize',16,'FontWeight','Bold')
    ylabel(metric_label,'FontSize',18,'FontWeight','Bold')
    axis([0 size(array,1)+1 -0.1 max(array(:))+0.15])
    
    subplot(212), hold on
    bar(array')
    for k=1:numel(x)
        text(y(k)-0.45+0.8*x(k)./size(array,1), array(x(k),y(k)), ['(',num2str(rM(x(k),y(k))),' , ',num2str(sP(x(k),y(k))),')'],...
            'VerticalAlignment','Middle','HorizontalAlignment','Left','FontSize',12,'Color',cb(round(x(k)./size(array,1).*size(cb,1)),:),'Rotation',90)
    end
    text(1-0.5, array(1,1), ['(r_M , \sigma_P)'],'VerticalAlignment','Middle','HorizontalAlignment','Left','FontSize',14,'FontWeight','Bold','Rotation',90)
    plot([0 numel(segMeth)+1],[array(1,1), array(1,1)],'k--')
    legend(netMeth,'Location','EastOutside')
    set(gca,'XTick',[1:numel(segMeth)],'XTickLabel',segMeth,'FontSize',16,'FontWeight','Bold')
    rotateXLabels( gca(), 45 )
    ylabel(metric_label,'FontSize',18,'FontWeight','Bold')
    axis([0 size(array,2)+1 -0.1 max(array(:))+0.15])
    
end
    