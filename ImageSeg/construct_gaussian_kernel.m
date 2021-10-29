function [kern] = construct_gaussian_kernel(sig)

% This function inputs the standard deviation of the gaussian kernel.  It 
% will then choose the smallest square gaussian kernel of odd dimension 
% (3x3 or 5x5 or 17x17) that has the the values drop to zero at the edges
% ie.{ sum(kern(1,:)) < 1e-6 }.

tol = 1e-6; % I can change this to be even smaller.
sz = 3;

kern = fspecial('gaussian',[sz sz], sig);

while sum(kern(1,:)) >= tol
    
    sz = sz+2;
    kern = fspecial('gaussian',[sz sz], sig);
    
end

