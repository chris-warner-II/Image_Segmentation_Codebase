function pb = make_gaussian_blur_pb_png_img_patches(sig)

% This function will loop through all the BSDS image patch mat files in the
% matInDir.  It will blur them by the gaussian kernel filter built using
% the prescribed filter.  The filter kernel is build to be the correct size
% given the constraint on sigma.  This function will place all the blurred
% png image patches in a directory called blur_sz#_sig#/pb_png within the
% matInDir

% sig = 0.5; % sigma or standard deviation of gaussian blur filter kernel - default = 0.5

[kern] = construct_gaussian_kernel(sig);
sz = size(kern,1);


[dirPre,sizeGoodIm] = onCluster;



%% (1). This is for Image file (just spatial gradients in image).

disp(['Creating PB_png image files for spatial gradients in actual image.'])

sig_tag = num2str(sig);
sig_tag(sig_tag=='.')='p';

% location of mat files with extracted image patches
matInDir = [dirPre,'images/BSDS_patch/101x101_ds1']; 
% matInDir = [dirPre,'images/BSDS_full/ds1/']; % for full images, not patches.


% location to write out the probabalistic boundaries png image file (input to benchmark code).
pbOutDir = [matInDir,'blur_sz',num2str(sz),'_sig',sig_tag,'/pb_png/'];
if ~exist(pbOutDir,'dir')
    mkdir(pbOutDir)
end


% 
files = dir([matInDir,'*.mat']);

disp(['Creating Probabalistic Boundary PB image files for Image Spatial Gradients'])

% Loop through all the files to convert F matrix (of spatial gradients) to png pb image file.
for i = 1:numel(files)
    
    fname = files(i).name;
    
    % do not rerun this if the pb_png files already exist.
    if exist([pbOutDir,fname(1:end-4),'.png'],'file')
        disp([num2str(i),' Probabalistic Boundary PNG File: ',fname(1:end-4),'.png',' :already exists. Moving on.'])
        pb = double(imread([pbOutDir,fname(1:end-4),'.png']))/255;
        continue
    end
    
    disp([num2str(i),' / ',num2str(numel(files)),' : ', fname(1:end-4),'.png'])
    
    
    % load image patch mat file
    load([matInDir,fname])
    
    
    imB = imfilter(im,kern,'symmetric');


    [Fx,Fy] = gradientB(imB,0); % this is a function I wrote to consider circular variables (CW)
                                         % Two differences between it and gradient()::
                                         % (1). It returns abs of gradients because I dont care whether they are negative.
                                         % (2). For circular variables, gradient is computed with shortest distance around circle.

    % I could do these other things to compute F from Fx & Fy, but just max is best I think.
    F2(:,:,1) = Fx;
    F2(:,:,2) = Fy;
    FimB = max(F2,[],3);
    
    
    
    if(0)
        [Fx,Fy] = gradientB(im,0); % this is a function I wrote to consider circular variables (CW)
        % I could do these other things to compute F from Fx & Fy, but just max is best I think.
        F2(:,:,1) = Fx;
        F2(:,:,2) = Fy;
        Fim = max(F2,[],3);
    
        figure, 
        subplot(221), imagesc(im), colormap(bone), colorbar, title('original'), axis square
        subplot(222), imagesc(imB), colormap(bone), colorbar, title('gaussian blur'), axis square
        subplot(223), imagesc(Fim./max(Fim(:))), colormap(bone), colorbar, title('original grads'), axis square
        subplot(224), imagesc(FimB./max(FimB(:))), colormap(bone), colorbar, title('blured grads'), axis square
    end
    
    
    
    
    pb = FimB ./ max(FimB(:)); % turn spatial gradients into "probabalistic boundaries"
    
    imwrite(pb,[pbOutDir,fname(1:end-4),'.png'],'PNG','BitDepth',8)
    

    clear F2 FimB
    
end









