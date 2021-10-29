function [THopt,Iopt] = optimize_threshold(img,truth,N)

% Syntax: THopt = optimize_threshold(img,truth,N)
%
% This function takes in an image (or equivalently, its eigenvector) and
% determines the threshold that will maximize the information contained in
% a segmentation of that image at that threshold given the ground truth
% entered also.  It finds the threshold by doing a linesearch and checking
% a number of thresholds between max and min pixel values. One can specify 
% the number of thresholds to check with N.
%
% Note: a better version of this is the optimize_threshold2 function.


if ~exist('N','var')
    N=100;
end
info = zeros(1,N);
Thr = linspace(min(img(:)),max(img(:)),N);

for i = 1:N
    seg = img>Thr(i);
    info(i) = calc_infoB(seg,truth);
end

ind = find(info==max(info));
%
THopt = Thr(ind);
THopt = THopt(1);
%
Iopt = info(ind);
Iopt = Iopt(1);

if(0)
    seg = (img>THopt);
    figure, subplot(131), plot(Thr,info),xlabel('Threshold'), ylabel('Information')
	subplot(132), imagesc(img), colorbar, title('Image / Eigenvector'),set(gca,'xtick',[],'ytick',[])
    subplot(133), imagesc(seg),title(['Segmentation at Optimal Threshold ',num2str(THopt)]), set(gca,'xtick',[],'ytick',[])
    keyboard
end

