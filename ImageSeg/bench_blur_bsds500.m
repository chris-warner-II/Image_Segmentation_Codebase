function bench_blur_bsds500(sz,sig)

% This function will loop through all the *.png files in the
% blur_sz#_sig#/pb_png directory and will run them through the benchmark
% code to output resultant *.txt files  in the
% blur_sz#_sig#/benchmark_results directory.


sig_tag = num2str(sig);
sig_tag(sig_tag=='.')='p';


thinpb = false;
maxDist = 1e-3; % 0.0075; %

[dirPre,sizeGoodIm] = onCluster;
addpath([dirPre,'images/BSDS_images/BSR/bench/benchmarks/'])




% paths changed to run from Projects directory
imgDir = [dirPre,'images/BSDS_patch/101x101_ds1/'];
gtDir = [dirPre,'images/BSDS_patch/101x101_ds1/groundTruth/'];
nthresh = 10;



% Directory to input pb.png images and to where benchmark results will be stored.
pbDir{1} = [dirPre,'images/BSDS_patch/101x101_ds1/blur_sz',num2str(sz),'_sig',sig_tag,'/pb_png/'];
outDir{1} = [dirPre,'images/BSDS_patch/101x101_ds1/blur_sz',num2str(sz),'_sig',sig_tag,'/benchmark_results/'];




    
for i = 1:numel(outDir) % this will loop thru different eigenvectors if spectral. Will be only 1 if Kur.

    if ~exist(outDir{i},'dir')
        mkdir(outDir{i})
    end


    tic;
    boundaryBench(imgDir, gtDir, pbDir{i}, outDir{i}, nthresh, maxDist, thinpb);
    toc;

    % plot_eval(outDir{i});


end
