function [kern] = construct_DoG_kernel(sigC,sigS,Krat)

center = construct_gaussian_kernel(sigC);

surround = Krat.*construct_gaussian_kernel(sigS);

C = padarray(center,[(size(surround,1) - size(center,1))/2,(size(surround,1) - size(center,1))/2]);

kern = C - surround;

% %



if(0)
    x=1:size(surround); y=1:size(surround);
    %
    disp('Integrated volume under center Gaussian is:')
    trapz(y,trapz(x,C,2),1)
    %
    disp('Integrated volume under surround Gaussian is:')
    trapz(y,trapz(x,surround,2),1)
end