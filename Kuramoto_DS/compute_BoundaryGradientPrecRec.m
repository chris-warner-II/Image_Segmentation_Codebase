function [precision,recall,f_measure] = compute_BoundaryGradientPrecRec(segResult, boundaryGT, numTH, methodName)


% This function takes in a segmentation results (either phase or pixels or
% eigenvector) and the segmentation ground truth in the form of boundaries
% hilighted from human segmenters and computes quality metric


[Fx,Fy] = gradient(segResult);
F = sqrt(Fx.^2 + Fy.^2);        % gradient magnitude of segmentation result

% if strmatch(methodName,'Modularity SK')
%     keyboard
% end

% set thresholds between min & max gradient value to compute precision and recall
TH = linspace(min(F(:)),max(F(:)),numTH);

for i = 1:size(boundaryGT,3)

    bD = boundaryGT(:,:,i);
    
    for j = 1:numTH
        
       G = (F>TH(j));
       
       out = Evaluate(G(:),bD);
       
       precision(i,j) = out(4);
       recall(i,j) = out(5);
       f_measure(i,j) = out(6);
        
    end

end



% Plot Precision vs Recall and F-measure vs. Threshold on Gradient
if(1)

    figure, 
    subplot(3,1,[1:2])
    plot(recall', precision','LineWidth',2)
    xlabel('Recall','FontSize',20,'FontWeight','Bold')
    ylabel('Precision','FontSize',20,'FontWeight','Bold')
    axis square
    legend({'blur1','blur2','blur3','blur4','blur5'})
    set(gca,'FontSize',16,'FontWeight','Bold')
    title(methodName)
    %
    subplot(3,10,[21:28])
    plot(TH,f_measure,'LineWidth',2)
    xlabel('threshold on gradient','FontSize',20,'FontWeight','Bold')
    ylabel('F-measure','FontSize',20,'FontWeight','Bold')
    set(gca,'FontSize',16,'FontWeight','Bold')
    axis([min(TH) max(TH) 0 1])
    text(TH(end-2), 0.9, ['\color{magenta}{',num2str(max(f_measure(5,:)), 2),'}'])
    text(TH(end-2), 0.8, ['\color{cyan}{',num2str(max(f_measure(4,:)), 2),'}'])
    text(TH(end-2), 0.7, ['\color{red}{',num2str(max(f_measure(3,:)), 2),'}'])
    text(TH(end-2), 0.6, ['\color{green}{',num2str(max(f_measure(2,:)), 2),'}'])
    text(TH(end-2), 0.5, ['\color{blue}{',num2str(max(f_measure(1,:)), 2),'}'])
    %
    subplot(3,10,29)
    imagesc(F), axis square
    set(gca,'XTick',[],'YTick',[]),
    title('Gradient','FontSize',20,'FontWeight','Bold'),
    colormap('Jet')
    freezeColors
    %
    subplot(3,10,30)
    imagesc(bD), axis square
    set(gca,'XTick',[],'YTick',[]),
    title('Gnd Truth','FontSize',20,'FontWeight','Bold'),
    colormap('Bone')
    freezeColors
    
end


