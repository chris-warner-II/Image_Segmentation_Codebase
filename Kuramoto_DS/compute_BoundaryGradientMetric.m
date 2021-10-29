function [F,M,S,D] = compute_BoundaryGradientMetric(segResult, boundaryGT, numGTs, circ, methodName)


% This function takes in a segmentation results (either phase or pixels or
% eigenvector) and the segmentation ground truth in the form of boundaries
% hilighted from human segmenters and computes quality metric


[Fx,Fy] = gradientB(segResult,circ); % this is a function I wrote to consider circular variables (CW)
                                     % Two differences between it and gradient()::
                                     % (1). It returns abs of gradients because I dont care whether they are negative.
                                     % (2). For circular variables, gradient is computed with shortest distance around circle.

%F = sqrt(Fx.^2 + Fy.^2);        % gradient magnitude of segmentation result

%F = segResult; % get rid of this!! (was useful when inputting boundary GT as segmentation)

% F = mean(F2,[],3); ??

% I could do these other things to compute F from Fx & Fy, but just max is best I think.
F2(:,:,1) = Fx;
F2(:,:,2) = Fy;
F = max(F2,[],3);


% Troubleshoot Diagnostic Plot.
if(0)
   figure,
   subplot(3,3,2), imagesc(segResult), caxis([0 2*pi]), axis square off, colormap('hsv'), freezeColors, colorbar, cbfreeze, title('Phase')
   subplot(3,3,4), imagesc(abs(Fx)), caxis([0 pi]), axis square off, colormap('jet'), freezeColors, colorbar, cbfreeze, title('Horiz Gradients')
   subplot(3,3,5), imagesc(1-boundaryGT(:,:,1)), axis square off, colormap('bone'), freezeColors, title('GT Boundaries')
   subplot(3,3,6), imagesc(abs(Fy)), caxis([0 pi]), axis square off, colormap('jet'), freezeColors, colorbar, cbfreeze, title('Vertical Gradients')
   subplot(3,3,8), imagesc(F), caxis([0 pi]), axis square off, colormap('jet'), freezeColors, colorbar, cbfreeze, title(' Gradients') 
end





for i = 1:size(boundaryGT,3)

    bD = boundaryGT(:,:,i);

    
    
    
    % NOTE: REWRITE THE LINES in the else statement TO HANDLE WEIGHTED BD MATRIX OR GENERAL CASE AT BEST.
    if(i==7)
        
        on = bD;
        PQ = on.*F;

        mean_on = sum(PQ(:)) ./ sum(on(:));
        mid_on = (F-mean_on).^2.*on;
        std_on = sqrt( sum(mid_on(:)) ./ sum(on(:))  ) ;

        M(i,:) = [mean_on, std_on];

        off = numGTs - bD;
        PR = off.*F;

        mean_off = sum(PR(:)) ./ sum(off(:));
        mid_off = (F-mean_off).^2.*off;
        std_off = sqrt( sum(mid_off(:)) ./ sum(off(:))) ;

        S(i,:) = [mean_off, std_off];


        if(0)
            figure, 
            subplot(231), imagesc(F), axis square, title('Gradients in Phase Field'), colorbar
            subplot(232), imagesc(on), axis square, title('Ground Truth Boundaries (on)'), colorbar
            subplot(233), imagesc(PQ), axis square, title('Alignment of Gradients and Boundaries'), colorbar
            %
            subplot(234),  axis square, off
            subplot(235), imagesc(off), axis square, title('Ground Truth Non-Boundaries (off)'), colorbar
            subplot(236), imagesc(PR), axis square, title('Alignment of Gradients and Non-Boundaries'), colorbar
        end
        
    else
        
        
        
        
        % NOTE: THESE LINES OF CODE WORK WHEN BD IS LOGICAL (1 or 0)
        PQ = bD.*F;

        % M quantifies how well gradients in phase field line up with
        % boundaries in ground truth (and how steep they are).
        M(i,:) = [mean(PQ(bD>0)), std(PQ(bD>0))];


        % quantify how many spurrious gradients (ones not on boundaries)
        % exist in phase fields and what statistics on them are.
        spur = (~bD).*F;
        S(i,:) = [mean(spur(~bD)),std(spur(~bD))];
    
    
        
        

    end
    
    
    
    
    
        
    if(0) % strmatch(methodName,'Modularity SK')
        
        cmap = colormap('jet');
        cmap(1,:) = [1,1,1];
        colormap(cmap);

        figure, 
        subplot(231), imagesc(1-bD), colormap('bone'), freezeColors, set(gca,'XTick',[],'YTick',[]), axis square
        title(['GT Boundaries : Blur = ',num2str(i)],'FontSize',18,'FontWeight','Bold')
        
        subplot(232), imagesc(spur - bD), colormap(cmap), colorbar, freezeColors, set(gca,'XTick',[],'YTick',[]), cbfreeze, axis square
        title({'Gradient Energy','\color{red}{Off} \color{black}{Boundaries}'},'FontSize',18,'FontWeight','Bold')
        
        bins=linspace(0,pi,100);
        [noff] = hist(spur(~bD),bins);
        [non] = hist(grad(bD),bins);
        
        subplot(233), imagesc(bD.*F - ~bD), colormap(cmap), colorbar, freezeColors, cbfreeze, set(gca,'XTick',[],'YTick',[]), axis square
        title({'Gradient Energy','\color{blue}{On} \color{black}{Boundaries}'},'FontSize',18,'FontWeight','Bold')
        
        subplot(212), hold on,
        plot(bins,noff./sum(noff),'r','LineWidth',2), 
        plot(bins,non./sum(non),'b','LineWidth',2), 
        %
        herrorbar(S(i,1),0.3,S(i,2),'rx')
        herrorbar(M(i,1),0.25,M(i,2),'bx')
        
        
        title('Phase Gradient Energy Distribution','FontSize',20,'FontWeight','Bold')
        axis([0 pi 0 0.4])
        xlabel('Phase Gradient','FontSize',18,'FontWeight','Bold')
        ylabel('Counts','FontSize',18,'FontWeight','Bold')
        legend({'Off Boundary','On Boundary'})
        set(gca,'FontSize',16,'FontWeight','Bold')
    end
    
    
    
    


end


% Compute Boundary Discriminability (d') metric from Mean & STD of Gradients On vs Off Boundaries.
D = ( M(:,1) - S(:,1) ) ./ sqrt( 0.5.*(M(:,2).^2 + S(:,2).^2 ) ); 
% BoundaryDiscriminability was not right because previously, I was squaring the standard deviations.
% Actually, no it WAS right !!! The standard deviations SHOULD be squared.
% Variance = STD^2.





if(0) % strmatch(methodName,'Modularity SK')
    keyboard
end
