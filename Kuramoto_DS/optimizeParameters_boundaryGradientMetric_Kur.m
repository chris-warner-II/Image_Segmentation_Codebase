function [] = optimizeParameters_boundaryGradientMetric_Kur(method,which_bD)


% which_bD = 7;


[dirPre,sizeGoodIm] = onCluster;

dirImgSave = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/imgs/Kur_PIF_Fourier1/allParams/'];
if ~exist(dirImgSave,'dir')
    mkdir(dirImgSave);
end


kur_data_dir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/',method,'/'];

fname = dir([kur_data_dir,'../ImPix/BoundaryGradientMetricKur_rM_KS_files*.mat']);
ImPix = load([kur_data_dir,'../ImPix/',fname(1).name]);

%
disp('rM1')
fname = dir([kur_data_dir,'BoundaryGradientMetricKur_rM1_KSsml_files*.mat']);
rM1_KSsml = load([kur_data_dir,fname(1).name]);
rM1_KSsml = rmfield(rM1_KSsml,'netParams'); % doing this so I dont overload memory with Mod_N&G ( can get rid of this)
[x,ind_rM1_KSsml,indPix_rM1_KSsml] = intersect(rM1_KSsml.ImgPtchID,ImPix.ImgPtchID);
fname = dir([kur_data_dir,'BoundaryGradientMetricKur_rM1_KSmid_files*.mat']);
rM1_KSmid = load([kur_data_dir,fname(1).name]);
rM1_KSmid = rmfield(rM1_KSmid,'netParams'); % doing this so I dont overload memory with Mod_N&G ( can get rid of this)
[x,ind_rM1_KSmid,indPix_rM1_KSmid] = intersect(rM1_KSmid.ImgPtchID,ImPix.ImgPtchID);
fname = dir([kur_data_dir,'BoundaryGradientMetricKur_rM1_KSlrg_files*.mat']);
rM1_KSlrg = load([kur_data_dir,fname(1).name]);
rM1_KSlrg = rmfield(rM1_KSlrg,'netParams'); % doing this so I dont overload memory with Mod_N&G ( can get rid of this)
[x,ind_rM1_KSlrg,indPix_rM1_KSlrg] = intersect(rM1_KSlrg.ImgPtchID,ImPix.ImgPtchID);

%
disp('rM3')
fname = dir([kur_data_dir,'BoundaryGradientMetricKur_rM3_KSsml_files*.mat']);
rM3_KSsml = load([kur_data_dir,fname(1).name]);
rM3_KSsml = rmfield(rM3_KSsml,'netParams'); % doing this so I dont overload memory with Mod_N&G ( can get rid of this)
[x,ind_rM3_KSsml,indPix_rM3_KSsml] = intersect(rM3_KSsml.ImgPtchID,ImPix.ImgPtchID);
fname = dir([kur_data_dir,'BoundaryGradientMetricKur_rM3_KSmid_files*.mat']);
rM3_KSmid = load([kur_data_dir,fname(1).name]);
rM3_KSmid = rmfield(rM3_KSmid,'netParams'); % doing this so I dont overload memory with Mod_N&G ( can get rid of this)
[x,ind_rM3_KSmid,indPix_rM3_KSmid] = intersect(rM3_KSsml.ImgPtchID,ImPix.ImgPtchID);
fname = dir([kur_data_dir,'BoundaryGradientMetricKur_rM3_KSlrg_files*.mat']);
rM3_KSlrg = load([kur_data_dir,fname(1).name]);
rM3_KSlrg = rmfield(rM3_KSlrg,'netParams'); % doing this so I dont overload memory with Mod_N&G ( can get rid of this)
[x,ind_rM3_KSlrg,indPix_rM3_KSlrg] = intersect(rM3_KSlrg.ImgPtchID,ImPix.ImgPtchID);

%
disp('rM5')
fname = dir([kur_data_dir,'BoundaryGradientMetricKur_rM5_KSsml_files*.mat']);
rM5_KSsml = load([kur_data_dir,fname(1).name]);
rM5_KSsml = rmfield(rM5_KSsml,'netParams'); % doing this so I dont overload memory with Mod_N&G ( can get rid of this)
[x,ind_rM5_KSsml,indPix_rM5_KSsml] = intersect(rM5_KSsml.ImgPtchID,ImPix.ImgPtchID);
fname = dir([kur_data_dir,'BoundaryGradientMetricKur_rM5_KSmid_files*.mat']);
rM5_KSmid = load([kur_data_dir,fname(1).name]);
rM5_KSmid = rmfield(rM5_KSmid,'netParams'); % doing this so I dont overload memory with Mod_N&G ( can get rid of this)
[x,ind_rM5_KSmid,indPix_rM5_KSmid] = intersect(rM5_KSmid.ImgPtchID,ImPix.ImgPtchID);
fname = dir([kur_data_dir,'BoundaryGradientMetricKur_rM5_KSlrg_files*.mat']);
rM5_KSlrg = load([kur_data_dir,fname(1).name]);
rM5_KSlrg = rmfield(rM5_KSlrg,'netParams'); % doing this so I dont overload memory with Mod_N&G ( can get rid of this)
[x,ind_rM5_KSlrg,indPix_rM5_KSlrg] = intersect(rM5_KSlrg.ImgPtchID,ImPix.ImgPtchID);

%
disp('rM10')
fname = dir([kur_data_dir,'BoundaryGradientMetricKur_rM10_KSsml_files*.mat']);
rM10_KSsml = load([kur_data_dir,fname(1).name]);
rM10_KSsml = rmfield(rM10_KSsml,'netParams'); % doing this so I dont overload memory with Mod_N&G ( can get rid of this)
[x,ind_rM10_KSsml,indPix_rM10_KSsml] = intersect(rM10_KSsml.ImgPtchID,ImPix.ImgPtchID);
fname = dir([kur_data_dir,'BoundaryGradientMetricKur_rM10_KSmid_files*.mat']);
rM10_KSmid = load([kur_data_dir,fname(1).name]);
rM10_KSmid = rmfield(rM10_KSmid,'netParams'); % doing this so I dont overload memory with Mod_N&G ( can get rid of this)
[x,ind_rM10_KSmid,indPix_rM10_KSmid] = intersect(rM10_KSmid.ImgPtchID,ImPix.ImgPtchID);
fname = dir([kur_data_dir,'BoundaryGradientMetricKur_rM10_KSlrg_files*.mat']);
rM10_KSlrg = load([kur_data_dir,fname(1).name]);
rM10_KSlrg = rmfield(rM10_KSlrg,'netParams'); % doing this so I dont overload memory with Mod_N&G ( can get rid of this)
[x,ind_rM10_KSlrg,indPix_rM10_KSlrg] = intersect(rM10_KSlrg.ImgPtchID,ImPix.ImgPtchID);




clear fname x








%%  Butterfly Plots
Hc = figure; 
%
subplot(3,4,1),
[mn_del_dp(1,1),std_del_dp(1,1)] = butterfly(ImPix.boundaryDiscriminability(which_bD,indPix_rM1_KSsml),rM1_KSsml.boundaryDiscriminability(which_bD,ind_rM1_KSsml));
title(['rM=1'],'FontSize',18,'FontWeight','Bold')
ylabel(['Ks=sml'],'FontSize',18,'FontWeight','Bold')
subplot(3,4,2),
[mn_del_dp(1,1),std_del_dp(1,1)] = butterfly(ImPix.boundaryDiscriminability(which_bD,indPix_rM3_KSsml),rM3_KSsml.boundaryDiscriminability(which_bD,ind_rM3_KSsml));
title(['rM=3'],'FontSize',18,'FontWeight','Bold')
subplot(3,4,3), 
[mn_del_dp(1,1),std_del_dp(1,1)] = butterfly(ImPix.boundaryDiscriminability(which_bD,indPix_rM5_KSsml),rM5_KSsml.boundaryDiscriminability(which_bD,ind_rM5_KSsml));
title(['rM=5'],'FontSize',18,'FontWeight','Bold')
subplot(3,4,4),
[mn_del_dp(1,1),std_del_dp(1,1)] = butterfly(ImPix.boundaryDiscriminability(which_bD,indPix_rM10_KSsml),rM10_KSsml.boundaryDiscriminability(which_bD,ind_rM10_KSsml));
title(['rM=10'],'FontSize',18,'FontWeight','Bold')
%
subplot(3,4,5),
[mn_del_dp(1,1),std_del_dp(1,1)] = butterfly(ImPix.boundaryDiscriminability(which_bD,indPix_rM1_KSmid),rM1_KSmid.boundaryDiscriminability(which_bD,ind_rM1_KSmid));
ylabel(['Ks=mid'],'FontSize',18,'FontWeight','Bold')
subplot(3,4,6),
[mn_del_dp(1,1),std_del_dp(1,1)] = butterfly(ImPix.boundaryDiscriminability(which_bD,indPix_rM3_KSmid),rM3_KSmid.boundaryDiscriminability(which_bD,ind_rM3_KSmid));
subplot(3,4,7),
[mn_del_dp(1,1),std_del_dp(1,1)] = butterfly(ImPix.boundaryDiscriminability(which_bD,indPix_rM5_KSmid),rM5_KSmid.boundaryDiscriminability(which_bD,ind_rM5_KSmid));
subplot(3,4,8),
[mn_del_dp(1,1),std_del_dp(1,1)] = butterfly(ImPix.boundaryDiscriminability(which_bD,indPix_rM10_KSmid),rM10_KSmid.boundaryDiscriminability(which_bD,ind_rM10_KSmid));
%
subplot(3,4,9),
[mn_del_dp(1,1),std_del_dp(1,1)] = butterfly(ImPix.boundaryDiscriminability(which_bD,indPix_rM1_KSlrg),rM1_KSlrg.boundaryDiscriminability(which_bD,ind_rM1_KSlrg));
ylabel(['Ks=lrg'],'FontSize',18,'FontWeight','Bold')
subplot(3,4,10),
[mn_del_dp(1,1),std_del_dp(1,1)] = butterfly(ImPix.boundaryDiscriminability(which_bD,indPix_rM3_KSlrg),rM3_KSlrg.boundaryDiscriminability(which_bD,ind_rM3_KSlrg));
subplot(3,4,11),
[mn_del_dp(1,1),std_del_dp(1,1)] = butterfly(ImPix.boundaryDiscriminability(which_bD,indPix_rM5_KSlrg),rM5_KSlrg.boundaryDiscriminability(which_bD,ind_rM5_KSlrg));
subplot(3,4,12),
[mn_del_dp(1,1),std_del_dp(1,1)] = butterfly(ImPix.boundaryDiscriminability(which_bD,indPix_rM10_KSlrg),rM10_KSlrg.boundaryDiscriminability(which_bD,ind_rM10_KSlrg));
%
annotation('textbox', [0 0.9 1 0.1],'String',[method],'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold')
saveGoodImg(Hc,[dirImgSave,'/Butterfly_',method,'_bD',num2str(which_bD)],sizeGoodIm)
close(Hc)















%% Plot 1. Consensus Boundaries vs Any Boundaries (Is there an improvement by throwing away single-segmenter boundaries? Not really.)
if(0)
    
    % note: can loop thru different rM and KS values and compare location of mean or something
    anyBound = rM10_KSmid.boundaryDiscriminability(1,:);
    consensus = rM10_KSmid.boundaryDiscriminability(4,:);


    figure, hold on,
    scatter(anyBound,consensus,40,'b.')
    scatter(mean(anyBound),mean(consensus),40,'r.')
    errorbar( mean(anyBound), mean(consensus), std(consensus), 'r')
    herrorbar( mean(anyBound), mean(consensus), std(anyBound), 'r')
    axis([0,2,0,2])
    plot([0 2],[0 2],'k--')

    title('Improvement using consensus boundaries (bD>1) rather than any boundaries (bD>0) ?')
    xlabel('d'' any Boundary')
    ylabel('d'' consensus Boundaries')

end







%% Plot 2. Compare d' across different parameter settings.


% NOTE: These only work if number of image patch files processed in each combine file are the same.
if(0)
    % Matrix of d' values for different image patches and different parameter
    % values where GT boundaries were any boundary drawn (bD>0)
    Xany = [...
            rM1_KSsml.boundaryDiscriminability(1,:) - ImPix.boundaryDiscriminability(1,:); ...
            rM1_KSmid.boundaryDiscriminability(1,:) - ImPix.boundaryDiscriminability(1,:); ...
            rM1_KSlrg.boundaryDiscriminability(1,:) - ImPix.boundaryDiscriminability(1,:); ...
            rM3_KSsml.boundaryDiscriminability(1,:) - ImPix.boundaryDiscriminability(1,:); ...
            rM3_KSmid.boundaryDiscriminability(1,:) - ImPix.boundaryDiscriminability(1,:); ...
            rM3_KSlrg.boundaryDiscriminability(1,:) - ImPix.boundaryDiscriminability(1,:); ...
            rM5_KSsml.boundaryDiscriminability(1,:) - ImPix.boundaryDiscriminability(1,:); ...
            rM5_KSmid.boundaryDiscriminability(1,:) - ImPix.boundaryDiscriminability(1,:); ...
            rM5_KSlrg.boundaryDiscriminability(1,:) - ImPix.boundaryDiscriminability(1,:); ...
            rM10_KSsml.boundaryDiscriminability(1,:) - ImPix.boundaryDiscriminability(1,:); ...
            rM10_KSmid.boundaryDiscriminability(1,:) - ImPix.boundaryDiscriminability(1,:); ...
            rM10_KSlrg.boundaryDiscriminability(1,:) - ImPix.boundaryDiscriminability(1,:) ...
            ];

    % Matrix of d' values for different image patches and different parameter
    % values where GT boundaries were consensus boundaries (bD>1)
    Xcon = [...
            rM1_KSsml.boundaryDiscriminability(4,:) - ImPix.boundaryDiscriminability(4,:); ...
            rM1_KSmid.boundaryDiscriminability(4,:) - ImPix.boundaryDiscriminability(4,:); ...
            rM1_KSlrg.boundaryDiscriminability(4,:) - ImPix.boundaryDiscriminability(4,:); ...
            rM3_KSsml.boundaryDiscriminability(4,:) - ImPix.boundaryDiscriminability(4,:); ...
            rM3_KSmid.boundaryDiscriminability(4,:) - ImPix.boundaryDiscriminability(4,:); ...
            rM3_KSlrg.boundaryDiscriminability(4,:) - ImPix.boundaryDiscriminability(4,:); ...
            rM5_KSsml.boundaryDiscriminability(4,:) - ImPix.boundaryDiscriminability(4,:); ...
            rM5_KSmid.boundaryDiscriminability(4,:) - ImPix.boundaryDiscriminability(4,:); ...
            rM5_KSlrg.boundaryDiscriminability(4,:) - ImPix.boundaryDiscriminability(4,:); ...
            rM10_KSsml.boundaryDiscriminability(4,:) - ImPix.boundaryDiscriminability(4,:); ...
            rM10_KSmid.boundaryDiscriminability(4,:) - ImPix.boundaryDiscriminability(4,:); ...
            rM10_KSlrg.boundaryDiscriminability(4,:) - ImPix.boundaryDiscriminability(4,:) ...
            ];   
        
        
    ParamsLabel = {'{1,sml}','{1,mid}','{1,lrg}','{3,sml}','{3,mid}','{3,lrg}','{5,sml}','{5,mid}','{5,lrg}','{10,sml}','{10,mid}','{10,lrg}'};

    
    figure, hold on,
    errorbar([1:12], mean(Xany'), std(Xany'), 'bx')
    errorbar([1:12]+0.1, mean(Xcon'), std(Xcon'), 'rx')
    set(gca,'XTick',[1:12],'XTickLabel',ParamsLabel,'FontSize',16,'FontWeight','Bold')
    %rotateXLabels( gca(), 45)
    legend({'any','consensus'})
    title([method,' : Finding Best Parameters'],'FontSize',20,'FontWeight','Bold')
    xlabel('Parameters \{rM,KS\}','FontSize',18,'FontWeight','Bold')
    ylabel(['\Deltad'' (',method,' - ImPix)'],'FontSize',18,'FontWeight','Bold')
    grid on
        
        
        
        
end




end




%% Function to do single pane in subplot figure of Butterfly Scatter Plot of Method vs ImPix d'.
function [mn_del_dp,std_del_dp] = butterfly(sigX,sigY)

    minn = -1; %min( [sigX, sigY] );
    maxx = 3;  %max( [sigX, sigY] );
    scatter(sigX, sigY, 'b.'), hold on, 
    plot([minn maxx], [minn maxx], 'k--')
    plot([0 0], [minn maxx], 'k--')
    plot([minn maxx], [0 0], 'k--')
    axis([minn maxx minn maxx])
    axis square, 
    set(gca,'FontSize',16,'FontWeight','Bold','XTick',linspace(floor(minn),ceil(maxx),5),'YTick',linspace(floor(minn),ceil(maxx),5) )
    %
    mn_del_dp = mean(sigY-sigX);
    std_del_dp = std(sigY-sigX);
    text(maxx,minn,['\Deltad''=',num2str(mn_del_dp,2)],'VerticalAlignment','Bottom','HorizontalAlignment','Right','FontSize',20,'FontWeight','Bold')

end

