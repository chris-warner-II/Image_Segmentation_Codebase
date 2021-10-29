% This was just a script to explore the interplay between size of gaussian
% blurring kernel and the resulting blurred image.  I find that the
% resulting blurred images are the same (or close enough to the same) as 
% long as the gaussian kernel falls close enough to zero at its edges.
%
% I also learned that I want the gaussian kernel size to be odd numbers
% because they will shift the pixels in the blurred image around relative
% to the original image.





sig = 2;
sz = [3,5,7,11,15,21];
% NOTE: Having sz be an even number shifts all pixels in the pb image
% relative to original image (I think)


kern1 = fspecial('gaussian',[sz(1) sz(1)], sig);
pb(:,:,1) = make_gaussian_blur_pb_png_img_patches(sz(1),sig);

kern2 = fspecial('gaussian',[sz(2) sz(2)], sig);
pb(:,:,2) = make_gaussian_blur_pb_png_img_patches(sz(2),sig);

kern3 = fspecial('gaussian',[sz(3) sz(3)], sig);
pb(:,:,3) = make_gaussian_blur_pb_png_img_patches(sz(3),sig);

kern4 = fspecial('gaussian',[sz(4) sz(4)], sig);
pb(:,:,4) = make_gaussian_blur_pb_png_img_patches(sz(4),sig);

kern5 = fspecial('gaussian',[sz(5) sz(5)], sig);
pb(:,:,5) = make_gaussian_blur_pb_png_img_patches(sz(5),sig);

kern6 = fspecial('gaussian',[sz(6) sz(6)], sig);
pb(:,:,6) = make_gaussian_blur_pb_png_img_patches(sz(6),sig);



% This is a measure of how the kernel has dropped to zero at the edges
sum(kern5(1,:))


% Hardcoding what image to load in
load('/Users/world7one/Desktop/Grad_School/Berkeley/Work/Fritz_Work/Projects/images/BSDS_patch/101x101_ds1/100007_ptch1.mat')
pb(:,:,7) = MC.F./max(MC.F(:));


% Make some plots just to visualize.
if(0)
    figure
    subplot(121), imagesc(kern1), axis square, colormap(bone), colorbar('Location','SouthOutside')
    title(['Sz = ',num2str(sz(1))])
    subplot(122), imagesc(pb(:,:,1)), axis square, colormap(bone), colorbar('Location','SouthOutside')


    figure
    subplot(121), imagesc(kern2), axis square, colormap(bone), colorbar('Location','SouthOutside')
    title(['Sz = ',num2str(sz(2))])
    subplot(122), imagesc(pb(:,:,2)), axis square, colormap(bone), colorbar('Location','SouthOutside')
end




% What is the largest maximal difference between the pb image produced from
% different size gaussian kernels with same sigma value?
for i = 1:size(pb,3)
    for j = i+1:size(pb,3)
        res = pb(:,:,i) - pb(:,:,j);
        AA(i,j) = max(abs(res(:)));
    end
end
