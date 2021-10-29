 
hIm=figure; subplot(121), imagesc(imFull),colormap('bone'), hold on, axis square
title('Full Image','FontSize',20,'FontWeight','Bold')
set(gca,'FontSize',16,'FontWeight','Bold')

% plot a box around the patch
plot([pach.xpbeg, pach.xpbeg],[pach.ypbeg, pach.ypfin],'r')
plot([pach.xpfin, pach.xpfin],[pach.ypbeg, pach.ypfin],'r')
plot([pach.xpbeg, pach.xpfin],[pach.ypbeg, pach.ypbeg],'r')
plot([pach.xpbeg, pach.xpfin],[pach.ypfin, pach.ypfin],'r')


subplot(2,6,[4,5]), imagesc(im), colormap('bone'),  axis square off % freezeColors,
title(['Patch ',num2str(pach.xpat),'x',num2str(pach.ypat)],'FontSize',20,'FontWeight','Bold')



% plot Ground Truth Segmentations.
for b = 1:numel(gT)
    
    subplot(numel(gT), 6, b*6), imagesc(gT{b}), axis square
    set(gca,'XTick',[],'YTick',[])
    colormap('jet'), freezeColors
    
    if(b==1)
        title({'Gnd Truth','Segmentations'},'FontSize',20,'FontWeight','Bold')
    end
    
    % xlabel(num2str(Avg_AUC_ROC(b),2))
    ylabel(['C=',num2str(numel(unique(gT{b})))],'FontSize',18,'FontWeight','Bold')
    
end

% plot Ground Truth Boundaries.
subplot(2,6,[10,11]),
imagesc(max(bD(:)) - bD), colormap('bone'), %freezeColors
title('Gnd Truth Boundaries','FontSize',20,'FontWeight','Bold')
set(gca,'XTick',[],'YTick',[])
axis square