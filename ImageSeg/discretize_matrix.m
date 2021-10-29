function M = discretize_matrix(M,bins)

% syntax: Mout = discretize_matrix(Min,bins);
%
% This function takes in a continuous, positive real-valued matrix of numbers
% and outputs a discretized version of it that has been thresholded into a 
% user prescribed number of bins.  For Hillar's idea of max Entropy
% segmentation.

M = M ./ max(max(M)); % normalize matrix
M = floor(M.*bins);   % discretize
