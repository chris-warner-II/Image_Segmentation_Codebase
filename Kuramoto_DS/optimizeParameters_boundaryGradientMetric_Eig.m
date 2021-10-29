function [] = optimizeParameters_boundaryGradientMetric_Eig(method,which_bD)


% which_bD = 7;


[dirPre,sizeGoodIm] = onCluster;

dirImgSave = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/imgs/spectral/allParams/'];
if ~exist(dirImgSave,'dir')
    mkdir(dirImgSave);
end


spec_data_dir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/spectral/allParams/',method,'/'];

fname = dir([spec_data_dir,'BoundaryGradientMetricEig_rM1_files*.mat']);
rM1 = load([spec_data_dir,fname(1).name]);
%
fname = dir([spec_data_dir,'BoundaryGradientMetricEig_rM3_files*.mat']);
rM3 = load([spec_data_dir,fname(1).name]);
%
fname = dir([spec_data_dir,'BoundaryGradientMetricEig_rM5_files*.mat']);
rM5 = load([spec_data_dir,fname(1).name]);
%
fname = dir([spec_data_dir,'BoundaryGradientMetricEig_rM10_files*.mat']);
rM10 = load([spec_data_dir,fname(1).name]);

clear fname

% NOTE:  I do not have to "register" ImPix to rM values because ImPix is
% calculated individually in each file (so numbers match).




%% Butterfly Plots for different parameters and different eigenvector combinations.
H=figure;



subplot(3,4,1), 
[mn_del_dp(1,1),std_del_dp(1,1)] = butterfly( rM1.boundaryDiscriminability.im(which_bD,:), rM1.boundaryDiscriminability.ev1(which_bD,:) ); 
ylabel('ev1','FontSize',20,'FontWeight','Bold')
title('rM=1','FontSize',20,'FontWeight','Bold')
subplot(3,4,5), 
[mn_del_dp(1,2),std_del_dp(1,2)] = butterfly( rM1.boundaryDiscriminability.im(which_bD,:), rM1.boundaryDiscriminability.ev2o(which_bD,:) ); 
ylabel('ev2','FontSize',20,'FontWeight','Bold')
subplot(3,4,9), 
[mn_del_dp(1,3),std_del_dp(1,3)] = butterfly( rM1.boundaryDiscriminability.im(which_bD,:), rM1.boundaryDiscriminability.ev3o(which_bD,:) ); 
ylabel('ev3','FontSize',20,'FontWeight','Bold')
%
subplot(3,4,2),
[mn_del_dp(1,1),std_del_dp(1,1)] = butterfly( rM3.boundaryDiscriminability.im(which_bD,:), rM3.boundaryDiscriminability.ev1(which_bD,:) ); 
title('rM=3','FontSize',20,'FontWeight','Bold')
subplot(3,4,6), 
[mn_del_dp(1,2),std_del_dp(1,2)] = butterfly( rM3.boundaryDiscriminability.im(which_bD,:), rM3.boundaryDiscriminability.ev2o(which_bD,:) ); 
%title('ev2o','FontSize',20,'FontWeight','Bold')
subplot(3,4,10), 
[mn_del_dp(1,3),std_del_dp(1,3)] = butterfly( rM3.boundaryDiscriminability.im(which_bD,:), rM3.boundaryDiscriminability.ev3o(which_bD,:) ); 
%title('ev3o','FontSize',20,'FontWeight','Bold')
%
subplot(3,4,3),
[mn_del_dp(1,1),std_del_dp(1,1)] = butterfly( rM5.boundaryDiscriminability.im(which_bD,:), rM5.boundaryDiscriminability.ev1(which_bD,:) ); 
title('rM=5','FontSize',20,'FontWeight','Bold')
subplot(3,4,7), 
[mn_del_dp(1,2),std_del_dp(1,2)] = butterfly( rM5.boundaryDiscriminability.im(which_bD,:), rM5.boundaryDiscriminability.ev2o(which_bD,:) ); 
%title('ev2o','FontSize',20,'FontWeight','Bold')
subplot(3,4,11), 
[mn_del_dp(1,3),std_del_dp(1,3)] = butterfly( rM5.boundaryDiscriminability.im(which_bD,:), rM5.boundaryDiscriminability.ev3o(which_bD,:) ); 
%title('ev3o','FontSize',20,'FontWeight','Bold')
%
subplot(3,4,4),
[mn_del_dp(1,1),std_del_dp(1,1)] = butterfly( rM10.boundaryDiscriminability.im(which_bD,:), rM10.boundaryDiscriminability.ev1(which_bD,:) ); 
title('rM=10','FontSize',20,'FontWeight','Bold')
subplot(3,4,8), 
[mn_del_dp(1,2),std_del_dp(1,2)] = butterfly( rM10.boundaryDiscriminability.im(which_bD,:), rM10.boundaryDiscriminability.ev2o(which_bD,:) ); 
%title('ev2o','FontSize',20,'FontWeight','Bold')
subplot(3,4,12), 
[mn_del_dp(1,3),std_del_dp(1,3)] = butterfly( rM10.boundaryDiscriminability.im(which_bD,:), rM10.boundaryDiscriminability.ev3o(which_bD,:) ); 
%title('ev3o','FontSize',20,'FontWeight','Bold')



annotation('textbox', [0 0.9 1 0.1],'String',[method],'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold')







saveGoodImg(H,[dirImgSave,'/Butterfly_',method,'_bD',num2str(which_bD)],sizeGoodIm)
close(H)







% H2=figure;
% subplot(2,4,1), 
% [mn_del_dp(2,1),std_del_dp(2,1)] = butterfly( rM3.boundaryDiscriminability.im(which_bD,:), rM3.boundaryDiscriminability.ev1(which_bD,:) ); 
% title('ev1','FontSize',20,'FontWeight','Bold')
% subplot(2,4,2), 
% [mn_del_dp(2,2),std_del_dp(2,2)] = butterfly( rM3.boundaryDiscriminability.im(which_bD,:), rM3.boundaryDiscriminability.ev2o(which_bD,:) ); 
% title('ev2o','FontSize',20,'FontWeight','Bold')
% subplot(2,4,3), 
% [mn_del_dp(2,3),std_del_dp(2,3)] = butterfly( rM3.boundaryDiscriminability.im(which_bD,:), rM3.boundaryDiscriminability.ev3o(which_bD,:) ); 
% title('ev3o','FontSize',20,'FontWeight','Bold')
% subplot(2,4,4), text(0.5,0.5,{method,'rM=3'},'FontSize',20,'FontWeight','Bold'), axis off
% subplot(2,4,5), 
% [mn_del_dp(2,4),std_del_dp(2,4)] = butterfly( rM3.boundaryDiscriminability.im(which_bD,:), rM3.boundaryDiscriminability.ev2(which_bD,:) ); 
% title('ev1-2','FontSize',20,'FontWeight','Bold')
% subplot(2,4,6), 
% [mn_del_dp(2,5),std_del_dp(2,5)] = butterfly( rM3.boundaryDiscriminability.im(which_bD,:), rM3.boundaryDiscriminability.ev3(which_bD,:) ); 
% title('ev1-3','FontSize',20,'FontWeight','Bold')
% subplot(2,4,7), 
% [mn_del_dp(2,6),std_del_dp(2,6)] = butterfly( rM3.boundaryDiscriminability.im(which_bD,:), rM3.boundaryDiscriminability.ev2w(which_bD,:) ); 
% title('ev1-2w','FontSize',20,'FontWeight','Bold')
% subplot(2,4,8), 
% [mn_del_dp(2,7),std_del_dp(2,7)] = butterfly( rM3.boundaryDiscriminability.im(which_bD,:), rM3.boundaryDiscriminability.ev3w(which_bD,:) ); 
% title('ev1-3w','FontSize',20,'FontWeight','Bold')
% saveGoodImg(H2,[dirImgSave,'/Butterfly_',method,'_rM3_bD',num2str(which_bD)],sizeGoodIm)
% close(H2)
% 
% 
% H3=figure;
% subplot(2,4,1), 
% [mn_del_dp(3,1),std_del_dp(3,1)] = butterfly( rM5.boundaryDiscriminability.im(which_bD,:), rM5.boundaryDiscriminability.ev1(which_bD,:) ); 
% title('ev1','FontSize',20,'FontWeight','Bold')
% subplot(2,4,2), 
% [mn_del_dp(3,2),std_del_dp(3,2)] = butterfly( rM5.boundaryDiscriminability.im(which_bD,:), rM5.boundaryDiscriminability.ev2o(which_bD,:) ); 
% title('ev2o','FontSize',20,'FontWeight','Bold')
% subplot(2,4,3), 
% [mn_del_dp(3,3),std_del_dp(3,3)] = butterfly( rM5.boundaryDiscriminability.im(which_bD,:), rM5.boundaryDiscriminability.ev3o(which_bD,:) ); 
% title('ev3o','FontSize',20,'FontWeight','Bold') 
% subplot(2,4,4), text(0.5,0.5,{method,'rM=5'},'FontSize',20,'FontWeight','Bold'), axis off
% subplot(2,4,5), 
% [mn_del_dp(3,4),std_del_dp(3,4)] = butterfly( rM5.boundaryDiscriminability.im(which_bD,:), rM5.boundaryDiscriminability.ev2(which_bD,:) ); 
% title('ev1-2','FontSize',20,'FontWeight','Bold')
% subplot(2,4,6), 
% [mn_del_dp(3,5),std_del_dp(3,5)] = butterfly( rM5.boundaryDiscriminability.im(which_bD,:), rM5.boundaryDiscriminability.ev3(which_bD,:) ); 
% title('ev1-3','FontSize',20,'FontWeight','Bold')
% subplot(2,4,7), 
% [mn_del_dp(3,6),std_del_dp(3,6)] = butterfly( rM5.boundaryDiscriminability.im(which_bD,:), rM5.boundaryDiscriminability.ev2w(which_bD,:) ); 
% title('ev1-2w','FontSize',20,'FontWeight','Bold')
% subplot(2,4,8), 
% [mn_del_dp(3,7),std_del_dp(3,7)] = butterfly( rM5.boundaryDiscriminability.im(which_bD,:), rM5.boundaryDiscriminability.ev3w(which_bD,:) ); 
% title('ev1-3w','FontSize',20,'FontWeight','Bold')
% saveGoodImg(H3,[dirImgSave,'/Butterfly_',method,'_rM5_bD',num2str(which_bD)],sizeGoodIm)
% close(H3)
% 
% 
% H4=figure;
% subplot(2,4,1), 
% [mn_del_dp(4,1),std_del_dp(4,1)] = butterfly( rM10.boundaryDiscriminability.im(which_bD,:), rM10.boundaryDiscriminability.ev1(which_bD,:) ); 
% title('ev1','FontSize',20,'FontWeight','Bold')
% subplot(2,4,2), 
% [mn_del_dp(4,2),std_del_dp(4,2)] = butterfly( rM10.boundaryDiscriminability.im(which_bD,:), rM10.boundaryDiscriminability.ev2o(which_bD,:) ); 
% title('ev2o','FontSize',20,'FontWeight','Bold')
% subplot(2,4,3), 
% [mn_del_dp(4,3),std_del_dp(4,3)] = butterfly( rM10.boundaryDiscriminability.im(which_bD,:), rM10.boundaryDiscriminability.ev3o(which_bD,:) ); 
% title('ev3o','FontSize',20,'FontWeight','Bold')
% subplot(2,4,4), text(0.5,0.5,{method,'rM=10'},'FontSize',20,'FontWeight','Bold'), axis off
% subplot(2,4,5), 
% [mn_del_dp(4,4),std_del_dp(4,4)] = butterfly( rM10.boundaryDiscriminability.im(which_bD,:), rM10.boundaryDiscriminability.ev2(which_bD,:) ); 
% title('ev1-2','FontSize',20,'FontWeight','Bold')
% subplot(2,4,6), 
% [mn_del_dp(4,5),std_del_dp(4,5)] = butterfly( rM10.boundaryDiscriminability.im(which_bD,:), rM10.boundaryDiscriminability.ev3(which_bD,:) );
% title('ev1-3','FontSize',20,'FontWeight','Bold')
% subplot(2,4,7),
% [mn_del_dp(4,6),std_del_dp(4,6)] = butterfly( rM10.boundaryDiscriminability.im(which_bD,:), rM10.boundaryDiscriminability.ev2w(which_bD,:) ); 
% title('ev1-2w','FontSize',20,'FontWeight','Bold')
% subplot(2,4,8), 
% [mn_del_dp(4,7),std_del_dp(4,7)] = butterfly( rM10.boundaryDiscriminability.im(which_bD,:), rM10.boundaryDiscriminability.ev3w(which_bD,:) ); 
% title('ev1-3w','FontSize',20,'FontWeight','Bold')
% saveGoodImg(H4,[dirImgSave,'/Butterfly_',method,'_rM10_bD',num2str(which_bD)],sizeGoodIm)
% close(H4)





end



%% Function to do single pane in subplot figure of Butterfly Scatter Plot of Method vs ImPix d'.
function [mn_del_dp,std_del_dp] = butterfly(sigX,sigY)

    minn = min( [sigX, sigY] );
    maxx = max( [sigX, sigY] );
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
