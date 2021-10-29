function [F] = compute_Spatial_Gradient(segResult, circ)


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

