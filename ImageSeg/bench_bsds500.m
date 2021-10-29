function bench_bsds500(kurOrSpec, method, rM, sD, sP, NF, ks, blur_flg)

% This function will take in a pb_png image file sitting in the pbDir
% (note: the directory structure for the pbDir or outDir is different 
% depending on whether you are doing (1). Kur or (2). Spec or (3). ImPix or
% ImBlur/ImDoG. There are several if statements in the function to handle those
% different cases.


% kurOrSpec = 'Kur_PIF_Fourier1';
% method = 'Mod_SKHEuc' %{'IsoDiff','AAnrm','GLnrm','Mod_SKHEuc','Mod_SKHAdj','Mod_N&G'}; 'ImPix', 'ImBlur', 'ImDoG' also.
% 
% rM = 'rM5' %{'rM1','rM3','rM5','rM10'};
% sD = 'sDInf';
% sP = 'sP0p2';
% NF = 'NF_60_0';
% ks = 'kssml' %{'kssml','ksmid','kslrg'};
% blur_flg = 1
%
%
%
% kurOrSpec = 'spectral';
% method = {'AAnrm','GLnrm','Mod_SKHAdj','Mod_N&G'};
% 
% rM = {'rM1','rM3','rM5','rM10'};
% sD = 'sDInf';
% sP = 'sP0p2';
% NF = 'na';
% ks = 'na';



if(blur_flg==2) % this means DoG filter
    blur_tag_M = '_blur_sigC1_S8_Kr0p01';        
    blur_tag_I = 'blur_sigC1_S8_Kr0p01/'; 
elseif(blur_flg==1) % this means Gaussian filter.    
    blur_tag_M = '_blur_sig1';      % '_blur_sig1';         
    blur_tag_I = 'blur_sz13_sig1/'; % 'blur_sz13_sig1/';
else
    blur_tag_M = ''; % if we are not blurring.
    blur_tag_I = '';
end


thinpb = false; % always false. Fvck boundary thinning operation!
maxDist = 2; % (units = pixels) 2 hardcoded into evaluation_bdry_imageB for now.

[dirPre,sizeGoodIm] = onCluster;
addpath([dirPre,'images/BSDS_images/BSR/bench/benchmarks/'])




% paths changed to run from Projects directory
imgDir = [dirPre,'images/BSDS_patch/101x101_ds1/'];
gtDir = [dirPre,'images/BSDS_patch/101x101_ds1/groundTruth/'];
nthresh = 10;






% If we are looking at Spatial Gradients in Phase Maps resulting from Kuramoto Coupled Oscillator Simulations, using this...
if strcmp(kurOrSpec,'Kur_PIF_Fourier1')

    if strcmp(method,'IsoDiff')
        pbDir{1} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/pb_png/',rM,'/',NF,'/',ks,'/'];
        outDir{1} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/benchmark_results/',rM,'/',NF,'/',ks,'/'];
    else
        pbDir{1} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/pb_png/',rM,'/',sD,'/',sP,'/',NF,'/',ks,'/'];
        outDir{1} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/benchmark_results/',rM,'/',sD,'/',sP,'/',NF,'/',ks,'/'];
    end
    
end



% If we are looking at Spatial Gradients in Eigenvectors, using this... (have to treat each eigenvector combination individually)
if strcmp(kurOrSpec,'spectral')    
    
    if strcmp(method,'IsoDiff')
        pbDir{1} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/pb_png/',rM,'/ev1/'];
        outDir{1} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/benchmark_results/',rM,'/ev1/'];
        %
        pbDir{2} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/pb_png/',rM,'/ev2o/'];
        outDir{2} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/benchmark_results/',rM,'/ev2o/'];
        %
        pbDir{3} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/pb_png/',rM,'/ev3o/'];
        outDir{3} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/benchmark_results/',rM,'/ev3o/'];
        %
        pbDir{4} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/pb_png/',rM,'/ev2/'];
        outDir{4} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/benchmark_results/',rM,'/ev2/'];
        %
        pbDir{5} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/pb_png/',rM,'/ev3/'];
        outDir{5} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/benchmark_results/',rM,'/ev3/'];
        %
        pbDir{6} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/pb_png/',rM,'/ev2w/'];
        outDir{6} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/benchmark_results/',rM,'/ev2w/'];
        %
        pbDir{7} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/pb_png/',rM,'/ev3w/'];
        outDir{7} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/benchmark_results/',rM,'/ev3w/'];
    else
        pbDir{1} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/pb_png/',rM,'/',sD,'/',sP,'/ev1/'];
        outDir{1} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/benchmark_results/',rM,'/',sD,'/',sP,'/ev1/'];
        %
        pbDir{2} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/pb_png/',rM,'/',sD,'/',sP,'/ev2o/'];
        outDir{2} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/benchmark_results/',rM,'/',sD,'/',sP,'/ev2o/'];
        %
        pbDir{3} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/pb_png/',rM,'/',sD,'/',sP,'/ev3o/'];
        outDir{3} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/benchmark_results/',rM,'/',sD,'/',sP,'/ev3o/'];
        %
        pbDir{4} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/pb_png/',rM,'/',sD,'/',sP,'/ev2/'];
        outDir{4} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/benchmark_results/',rM,'/',sD,'/',sP,'/ev2/'];
        %
        pbDir{5} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/pb_png/',rM,'/',sD,'/',sP,'/ev3/'];
        outDir{5} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/benchmark_results/',rM,'/',sD,'/',sP,'/ev3/'];
        %
        pbDir{6} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/pb_png/',rM,'/',sD,'/',sP,'/ev2w/'];
        outDir{6} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/benchmark_results/',rM,'/',sD,'/',sP,'/ev2w/'];
        %
        pbDir{7} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/pb_png/',rM,'/',sD,'/',sP,'/ev3w/'];
        outDir{7} = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/',kurOrSpec,'/',method,'/benchmark_results/',rM,'/',sD,'/',sP,'/ev3w/'];
        %
    end
    
end
   

% If we are looking at Blurred Image Pixel Spatial Gradients, using this...
if ( strcmp(method,'ImBlur') || strcmp(method,'ImDoG') )
    pbDir{1} = [dirPre,'images/BSDS_patch/101x101_ds1/',blur_tag_I,'pb_png/'];
    outDir{1} = [dirPre,'images/BSDS_patch/101x101_ds1/',blur_tag_I,'benchmark_results/'];
end


% If we are looking at Raw Pixel Spatial Gradients, using this... (Should be last - after other strcmp if checker statements)
if strcmp(method,'ImPix')
    pbDir{1} = [dirPre,'images/BSDS_patch/101x101_ds1/pb_png/'];
    outDir{1} = [dirPre,'images/BSDS_patch/101x101_ds1/benchmark_results/'];
end




    
for i = 1:numel(outDir) % this will loop thru different eigenvectors if spectral. Will be only 1 if Kur.

    if ~exist(outDir{i},'dir')
        mkdir(outDir{i})
    end


    tic;
    boundaryBench(imgDir, gtDir, pbDir{i}, outDir{i}, nthresh, maxDist, thinpb);
    toc;

    % plot_eval(outDir{i});


end
