           
GL1 = load('/Users/world7one/Desktop/Grad_School/Berkeley/Work/Fritz_Work/Projects/output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/RateDistPlots/RDdata_GLnrm_sP0p2_rM1_sDInf_sW0_Ks300_Ts1.mat');
GL7 = load('/Users/world7one/Desktop/Grad_School/Berkeley/Work/Fritz_Work/Projects/output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/RateDistPlots/RDdata_GLnrm_sP0p2_rM7_sDInf_sW0_Ks300_Ts1.mat');
GLi = load('/Users/world7one/Desktop/Grad_School/Berkeley/Work/Fritz_Work/Projects/output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/RateDistPlots/RDdata_GLnrm_sP0p2_rMInf_sDInf_sW0_Ks300_Ts1.mat');

NG1 = load('/Users/world7one/Desktop/Grad_School/Berkeley/Work/Fritz_Work/Projects/output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/RateDistPlots/RDdata_Mod_N&G_sP0p2_rM1_sDInf_sW0_Ks300_Ts1.mat');
NG7 = load('/Users/world7one/Desktop/Grad_School/Berkeley/Work/Fritz_Work/Projects/output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/RateDistPlots/RDdata_Mod_N&G_sP0p2_rM7_sDInf_sW0_Ks300_Ts1.mat');
NGi = load('/Users/world7one/Desktop/Grad_School/Berkeley/Work/Fritz_Work/Projects/output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/RateDistPlots/RDdata_Mod_N&G_sP0p2_rMInf_sDInf_sW0_Ks300_Ts1.mat');

SKH1 = load('/Users/world7one/Desktop/Grad_School/Berkeley/Work/Fritz_Work/Projects/output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/RateDistPlots/RDdata_Mod_SKHAdj_sP0p2_rM1_sDInf_sW0_Ks300_Ts1.mat');
SKH7 = load('/Users/world7one/Desktop/Grad_School/Berkeley/Work/Fritz_Work/Projects/output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/RateDistPlots/RDdata_Mod_SKHAdj_sP0p2_rM7_sDInf_sW0_Ks300_Ts1.mat');
SKHi = load('/Users/world7one/Desktop/Grad_School/Berkeley/Work/Fritz_Work/Projects/output/Kuramoto/NetsFromImgs/BSDS_patch_51x51_ds1/data/RateDistPlots/RDdata_Mod_SKHAdj_sP0p2_rMInf_sDInf_sW0_Ks300_Ts1.mat');


RateDistSig = GL1.RateDistSig;
segMethods = GL1.segMethods;

x1 = GL1.meanFins(:,2:end);
x2 = GL7.meanFins(:,2:end);
x3 = NG1.meanFins(:,2:end);
x4 = NG7.meanFins(:,2:end);
x5 = SKH1.meanFins(:,2:end);
x6 = SKH7.meanFins(:,2:end);
%
x7 = GLi.meanFins(:,2:end);
x8 = NGi.meanFins(:,2:end);
x9 = SKHi.meanFins(:,2:end);

ymax = max( [max(x1(:)), max(x2(:)), max(x3(:)), max(x4(:)), max(x5(:)), max(x6(:)) ] ); % max(x7(:)), max(x8(:)), max(x8(:))


figure, 


subplot(2,3,1),plot( RateDistSig(2:end),  x1' ,'LineWidth',3), hold on
area(RateDistSig(2:end),x4(2,:)','FaceColor',[0.8 0.8 0.8]), alpha(0.5)
title(['Graph Laplacian (R=1)'],'FontSize',20,'FontWeight','Bold')
xlabel('\sigma_{RD} ','FontSize',18,'FontWeight','Bold')
ylabel('d''','FontSize',18,'FontWeight','Bold')
axis([RateDistSig(2) RateDistSig(end) 0 ymax])
set(gca,'FontSize',16,'FontWeight','Bold')

subplot(2,3,4), plot( RateDistSig(2:end),  x2' ,'LineWidth',3), hold on
area(RateDistSig(2:end),x4(2,:)','FaceColor',[0.8 0.8 0.8]), alpha(0.5)
title(['(R=7)'],'FontSize',20,'FontWeight','Bold')
xlabel('\sigma_{RD} ','FontSize',18,'FontWeight','Bold')
ylabel('d''','FontSize',18,'FontWeight','Bold')
axis([RateDistSig(2) RateDistSig(end) 0 ymax])
set(gca,'FontSize',16,'FontWeight','Bold')


% subplot(3,3,7), plot( RateDistSig(2:end),  x7' ,'LineWidth',3), hold on
% area(RateDistSig(2:end),x4(2,:)','FaceColor',[0.8 0.8 0.8]), alpha(0.5)
% title(['(R=\infty)'],'FontSize',20,'FontWeight','Bold')
% xlabel('\sigma_{RD} ','FontSize',18,'FontWeight','Bold')
% ylabel('d''','FontSize',18,'FontWeight','Bold')
% axis([RateDistSig(2) RateDistSig(end) 0 ymax])
% set(gca,'FontSize',16,'FontWeight','Bold')



% % %

subplot(2,3,2),plot( RateDistSig(2:end),  x5' ,'LineWidth',3), hold on
area(RateDistSig(2:end),x4(2,:)','FaceColor',[0.8 0.8 0.8]), alpha(0.5)
title(['Topographic Modularity (R=1)'],'FontSize',20,'FontWeight','Bold')
xlabel('\sigma_{RD} ','FontSize',18,'FontWeight','Bold')
ylabel('d''','FontSize',18,'FontWeight','Bold')
axis([RateDistSig(2) RateDistSig(end) 0 ymax])
set(gca,'FontSize',16,'FontWeight','Bold')

subplot(2,3,5), plot( RateDistSig(2:end),  x6' ,'LineWidth',3), hold on
area(RateDistSig(2:end),x4(2,:)','FaceColor',[0.8 0.8 0.8]), alpha(0.5)
title(['(R=7)'],'FontSize',20,'FontWeight','Bold')
xlabel('\sigma_{RD} ','FontSize',18,'FontWeight','Bold')
ylabel('d''','FontSize',18,'FontWeight','Bold')
axis([RateDistSig(2) RateDistSig(end) 0 ymax])
set(gca,'FontSize',16,'FontWeight','Bold')

% subplot(3,3,8), plot( RateDistSig(2:end),  x9' ,'LineWidth',3), hold on
% area(RateDistSig(2:end),x4(2,:)','FaceColor',[0.8 0.8 0.8]), alpha(0.5)
% title(['(R=\infty)'],'FontSize',20,'FontWeight','Bold')
% xlabel('\sigma_{RD} ','FontSize',18,'FontWeight','Bold')
% ylabel('d''','FontSize',18,'FontWeight','Bold')
% axis([RateDistSig(2) RateDistSig(end) 0 ymax])
% set(gca,'FontSize',16,'FontWeight','Bold')



% % %

subplot(2,3,3), plot( RateDistSig(2:end),  x3' ,'LineWidth',3), hold on
area(RateDistSig(2:end),x4(2,:)','FaceColor',[0.8 0.8 0.8]), alpha(0.5)
title(['Nontopographic Modularity (R=1)'],'FontSize',20,'FontWeight','Bold')
xlabel('\sigma_{RD} ','FontSize',18,'FontWeight','Bold')
ylabel('d''','FontSize',18,'FontWeight','Bold')
axis([RateDistSig(2) RateDistSig(end) 0 ymax])
set(gca,'FontSize',16,'FontWeight','Bold')

subplot(2,3,6), plot( RateDistSig(2:end),  x4' ,'LineWidth',3), hold on
area(RateDistSig(2:end),x4(2,:)','FaceColor',[0.8 0.8 0.8]), alpha(0.5)
legend(segMethods)
title(['(R=7)'],'FontSize',20,'FontWeight','Bold')
xlabel('\sigma_{RD} ','FontSize',18,'FontWeight','Bold')
ylabel('d''','FontSize',18,'FontWeight','Bold')
set(gca,'FontSize',16,'FontWeight','Bold')
axis([RateDistSig(2) RateDistSig(end) 0 ymax])

% subplot(3,3,9), plot( RateDistSig(2:end),  x8' ,'LineWidth',3), hold on
% area(RateDistSig(2:end),x4(2,:)','FaceColor',[0.8 0.8 0.8]), alpha(0.5)
% legend(segMethods)
% title(['(R=\infty)'],'FontSize',20,'FontWeight','Bold')
% xlabel('\sigma_{RD} ','FontSize',18,'FontWeight','Bold')
% ylabel('d''','FontSize',18,'FontWeight','Bold')
% set(gca,'FontSize',16,'FontWeight','Bold')
% axis([RateDistSig(2) RateDistSig(end) 0 ymax])


