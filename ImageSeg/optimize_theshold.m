function THopt = optimize_threshold(img,truth,N)

% Syntax: THopt = optimize_threshold(img,truth,N)
%
% This function takes in an image (or equivalently, its eigenvector) and
% determines the threshold that will maximize the information contained in
% a segmentation of that image at that threshold given the ground truth
% entered also.  One can specify the number of thresholds to check with N.
%



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
THopt = Thr(ind);

if(0)
    seg = img>THopt;
    figure, subplot(131), plot(Thr,info),xlabel('Threshold'), ylabel('Information')
	subplot(132), imagesc(img), colorbar, title('Image / Eigenvector'),set(gca,'xtick',[],'ytick',[])
    subplot(133), imagesc(seg),title('Segmentation at Optimal Threshold'), set(gca,'xtick',[],'ytick',[])
end