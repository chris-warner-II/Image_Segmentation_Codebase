% This script is going to do some initial diagnostic plots for eigenvector
% analysis.

[dirPre,sizeGoodIm] = onCluster;


% Flag to use a gaussian kernel (sig=1) to preblur image before running network Kuramoto computation.
blur_flg=1;
if(blur_flg)
    blur_tag_M = '_blur_sig1';
    blur_tit = 'w/ Pre-Blurring (\sigma=1)';
    blur_tag_I = 'blur_sz13_sig1/';
else
    blur_tag_M = ''; % if we are not blurring.
    blur_tit = 'w/ No Pre-Blurring.';
    blur_tag_I = '/';
end







method = 'Mod_SKHAdj';

rM = '5';

Evs_data_dir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/spectral/',method,'/'];
%Kur_data_dir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/Kur_PIF_Fourier1/',method,'/'];

Evs_pb_dir = [Evs_data_dir,'pb_png/rM',rM,'/sDInf/sP0p2/'];
%Kur_pb_dir = 0;



Evs_benchmark_results = [Evs_data_dir,'benchmark_results/rM',rM,'/sDInf/sP0p2/']; % Have to run these ...
%Kur_benchmark_results = 0;



evFiles = dir([Evs_data_dir,'*_rM',rM,'_*.mat']);





% must add ground truth dir to calculate F-measure.
gTdir = [dirPre,'images/BSDS_patch/101x101_ds1/groundTruth/'];





for i = 1:numel(evFiles)
    
    % load EV mat file
    load([Evs_data_dir,evFiles(i).name])
    
    % gT file for computing F-measure below
    gTfile = [gTdir,netflags.fname,'.mat'];

    % pull out first 3 eigenvectors individually.
    ev1 = reshape(EVecsML(:,1),netParams.Ndims);
    ev2 = reshape(EVecsML(:,2),netParams.Ndims);
    ev3 = reshape(EVecsML(:,3),netParams.Ndims);

    
    % calculate pb images from eigenvectors directly in mat file
    ev1_pb = compute_Spatial_Gradient(ev1, 0);
    ev1_pb = ev1_pb./max(ev1_pb(:));
    ev2_pb = compute_Spatial_Gradient(ev2, 0);
    ev2_pb = ev2_pb./max(ev2_pb(:));
    ev3_pb = compute_Spatial_Gradient(ev3, 0);
    ev3_pb = ev3_pb./max(ev3_pb(:));
    %
    im_pb = compute_Spatial_Gradient(netParams.im, 0);
    im_pb = im_pb./max(im_pb(:));
    

    % load in pb image files from directory for comparison.
    pbFile1 = [Evs_pb_dir,'/ev1/',netflags.fname,'.png'];
    ev1b_pb = double(imread(pbFile1))/255;
    pbFile2 = [Evs_pb_dir,'/ev2o/',netflags.fname,'.png'];
    ev2b_pb = double(imread(pbFile2))/255;
    pbFile3 = [Evs_pb_dir,'/ev3o/',netflags.fname,'.png'];
    ev3b_pb = double(imread(pbFile3))/255;
    
    
    
    
    % Compare pb images calculated now vs ones gotten from images
    % previously. Only draw attention if different.
    ev1_pb_comp = abs(ev1_pb - ev1b_pb);
    ev2_pb_comp = abs(ev2_pb - ev2b_pb); 
    ev3_pb_comp = abs(ev3_pb - ev3b_pb); 
    %
    tol = 0.01;
    if(max(ev1_pb_comp(:)) > tol)
        disp('Ev1 pb different')
        figure,
        subplot(121), imagesc(ev1_pb), colormap(bone), title('computed from ev1')
        subplot(122), imagesc(ev1b_pb), colormap(bone), title('from pb img')
    end
    if(max(ev2_pb_comp(:)) > tol)
        disp('Ev2 pb different')
        figure,
        subplot(121), imagesc(ev2_pb), colormap(bone), title('computed from ev2')
        subplot(122), imagesc(ev2b_pb), colormap(bone), title('from pb img')
    end
    if(max(ev3_pb_comp(:)) > tol)
        disp('Ev3 pb different')
        figure,
        subplot(121), imagesc(ev3_pb), colormap(bone), title('computed from ev3')
        subplot(122), imagesc(ev3b_pb), colormap(bone), title('from pb img')
    end
    

    
    
    
    
    
    % write these to pb_png images to be used by evaluation_bdry_imageB function to compute F-measure.
    imwrite(ev1_pb,'./temp_ev1_pb.png','PNG','BitDepth',8);
    imwrite(ev2_pb,'./temp_ev2_pb.png','PNG','BitDepth',8);
    imwrite(ev3_pb,'./temp_ev3_pb.png','PNG','BitDepth',8);
    
    
    % compute R,P,F and save them in temporary text files
    nthr = 10;
    dt = 2; % in pixels
    xxx = evaluation_bdry_imageB('./temp_ev1_pb.png',gTfile,'temporary_ev1_prd_file',nthr,dt);
    xxx = evaluation_bdry_imageB('./temp_ev2_pb.png',gTfile,'temporary_ev2_prd_file',nthr,dt);
    xxx = evaluation_bdry_imageB('./temp_ev3_pb.png',gTfile,'temporary_ev3_prd_file',nthr,dt);



    % read in from those temporary text files saved just above (the contents is:) 
    %          1        2     3      4        5      6      7        8      9      10      11      12       13 
    vars = ['thresh', 'dt', 'mnR', 'stdR', 'Rmax', 'mnP', 'stdP', 'Pmax', 'mnF', 'stdF', 'Fmax', '#GT', 'BestGT'];
    for j = 1:4
        eval(['Aev1(:,:,j) = dlmread(''./temporary_ev1_pd',num2str(j),'_rd_file'');'])
        eval(['Aev2(:,:,j) = dlmread(''./temporary_ev2_pd',num2str(j),'_rd_file'');'])
        eval(['Aev3(:,:,j) = dlmread(''./temporary_ev3_pd',num2str(j),'_rd_file'');'])
    end



    % Find best Fmax & Fmean for each {ks,dt} combination
    for j = 1:4
        best_Fmax_ev1(j) = max(Aev1(:,11,j));
        best_Fmean_ev1(j) = max(Aev1(:,9,j));
        %
        best_Fmax_ev2(j) = max(Aev2(:,11,j));
        best_Fmean_ev2(j) = max(Aev2(:,9,j));
        %
        best_Fmax_ev3(j) = max(Aev3(:,11,j));
        best_Fmean_ev3(j) = max(Aev3(:,9,j));
    end
        
    
    
    
    
    
    
    
    
    
    
    
    % Plot image & eigenvectors and pb_imgs
    figure
    subplot(3,3,1), imagesc(netParams.im), colormap(bone), freezeColors
    axis square, set(gca,'XTick',[],'YTick',[]), title('Image')
    subplot(3,3,4), imagesc(im_pb), colormap(bone), freezeColors
    axis square, set(gca,'XTick',[],'YTick',[]), title('pb im') 
    subplot(3,3,7), imagesc(netParams.bD), colormap(bone), freezeColors
    axis square, set(gca,'XTick',[],'YTick',[]), title('Ground Truth')
    %
    subplot(3,3,2), imagesc(ev1), colormap(bone), freezeColors
    axis square, set(gca,'XTick',[],'YTick',[]), title('Ev1')
    subplot(3,3,5), imagesc(ev2), colormap(bone), freezeColors
    axis square, set(gca,'XTick',[],'YTick',[]), title('Ev2')
    subplot(3,3,8), imagesc(ev3), colormap(bone), freezeColors
    axis square, set(gca,'XTick',[],'YTick',[]), title('Ev3')
    %
    subplot(3,3,3), imagesc(ev1_pb), colormap(bone), freezeColors, title('pb')
    axis square, set(gca,'XTick',[],'YTick',[]),
    subplot(3,3,6), imagesc(ev2_pb), colormap(bone), freezeColors
    axis square, set(gca,'XTick',[],'YTick',[]),
    subplot(3,3,9), imagesc(ev3_pb), colormap(bone), freezeColors
    axis square, set(gca,'XTick',[],'YTick',[]),
    %

    
    
    
    
    
    
    keyboard
    

end