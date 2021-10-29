% This script is taken from and now called from explore_Separation_vs_Parameters .
% It will plot the time evolution of the DivMarg metric across the 
% 100 simulation runs.  It also plots the distribution of DivMarg at the
% end of simulations.  



pctBetterThan = cumsum(cntNd)./sum(cntNd); % normalized cumulative sum of histogram

hDM=figure;
subplot(1,6,1:3), plot(linspace(0,kurParams.Tsec,kurParams.T./kurParams.spp+1), DivMarg'), 
xlim([0 kurParams.Tsec]), ylim([0 1.1])
hold on, plot([0 kurParams.Tsec], [0.385 0.385], 'k--','LineWidth',2)
title([netflags.imageIn,' ',fileSubset,' ',netflags.method,' ',kurflags.KurTitleTag,' Kscale ',num2str(kurParams.Kscale)],'FontSize',14,'FontWeight','Bold')
xlabel('Time')
ylabel({'Divisive Relative Margin','(Within Cluster / Across Cluster)'})
%
subplot(1,6,4), barh(HistBins,cntNd./sum(cntNd)), ylim([0 1.1])
title('at End T')
xlabel('%')
set(gca,'YTick',[])



for g = 1:numBins
    if(cntNd(numBins-g+1)~=0)

        subplot(numBins,6,5 + (g-1)*6)
        imagesc(reshape(metaCluster(runInd_histNd(numBins-g+1)).phaseAtClk(:,end-1), netParams.Ndims(1), netParams.Ndims(2))), colormap('hsv'), 
        set(gca,'XTick',[],'YTick',[])
        axis square
        if(g==1)
            title({'Exemplar','Phase Img'})
        end

        subplot(numBins,6,6 + (g-1)*6),
        plot_phaseAtClk(metaCluster,MC,netParams,kurParams,runInd_histNd(numBins-g+1),1)
        set(gca,'XTick',[],'YTick',[])
        if(g==numBins)
            xlabel('t')
            ylabel('\phi')
        end
        if(g==1)
            title('Phase Evolution')
        end

     end

end

for g = 1:numBins
    if(cntNd(numBins-g+1)~=0)
        annotation('textbox', [0.90 0.08*(numBins-g+1) 0.1 0.1],'String', [num2str(pctBetterThan(numBins-g+1))],'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',12,'FontWeight','Bold')
    end
end

annotation('textbox', [0.8 0.9 0.1 0.1],'String', ['% Runs As Good Or Better'],'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',12,'FontWeight','Bold')

axes('Position',[.71 .05 .025 .045]), box on, imagesc(colorwheel), axis off

saveGoodImg(hDM,[imgKur,'DivMargDist_',fImg],[0 0 1 1])
close(hDM)