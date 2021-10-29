% This script/function written by CW 11/15 to plot multiple precision
% recall curves on the same isoF plot.  We can use this to compare
% different methods or different parameter settings for the same method.


[dirPre,sizeGoodIm] = onCluster;
addpath([dirPre,'images/BSDS_images/BSR/bench/benchmarks/'])


method = {'GLnrm','AAnrm','Mod_SKHAdj','Mod_N&G','IsoDiff'};

method_short = {'GL','AA','TM','M','ISO'};
method_color = {'g','y','r','b','k'};


rM = {'1','3','5','10'};
shapes = {'d-','o-','s-','^-'};
         % {'r','m','b','c'};

ev = {'ev1','ev2o','ev3o','ev2','ev3','ev2w','ev3w'}; 
ev_long = {'ev 1 only','ev 2 only','ev 3 only','ev 1&2','ev 1-3','ev 1&2 weighted','ev 1-3 weighted'};

colors = {'r','m','b','c'}; % {'r','g','b','c','m','y','k'};

plot_max = 1; % set to 1 to plot average of Precision & Recall leading to largest F, regardless of threshold
              % set to 0 to plot average of Precision & Recall at each threshold value for each image

num_cPdist = 4; 
cPd_colors = {'r','g','b','k'};

which_errbars = 'sem'; % either 'std' for standard deviation or 'sem' for standard error.    

switch which_errbars
    case 'sem'
        plot_name = 'SE';
    case'std'
        plot_name = '\sigma';
end

          
              
% Flag to use a gaussian kernel (sig=1) to preblur image before running network Kuramoto computation.
blur_flg=1;
if(blur_flg)
    blur_tag_M = '_blur_sig1';
    blur_tit = 'w/ Pre-Blurring (\sigma=1)';
else
    blur_tag_M = ''; % if we are not blurring.
    blur_tit = 'w/ No Pre-Blurring.';
end
blur_tag_I = 'blur_sz13_sig1/';





imPixBenchmarkDir = [dirPre,'images/BSDS_patch/101x101_ds1/benchmark_results/'];
imBlurBenchmarkDir = [dirPre,'images/BSDS_patch/101x101_ds1/',blur_tag_I,'benchmark_results/'];




outImgDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/imgs/spectral/'];
%
if ~exist(outImgDir,'dir')
    mkdir(outImgDir)
end




% For just the method - the absolute value. Not comparing to any strawman or null model.
%
% (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
justMethod.maxF_newMean_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
justMethod.maxF_newMean_semAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
justMethod.maxF_newStd_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
%
justMethod.maxR_newMean_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
justMethod.maxR_newMean_semAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
justMethod.maxR_newStd_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
%
justMethod.maxP_newMean_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
justMethod.maxP_newMean_semAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
justMethod.maxP_newStd_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
%
% (4). New way to compute P-R-F - (Max value using best single ground truther.). 
justMethod.maxF_newMax_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
justMethod.maxF_newMax_semAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
%
justMethod.maxR_newMax_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
justMethod.maxR_newMax_semAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
%
justMethod.maxP_newMax_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
justMethod.maxP_newMax_semAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));







% justImPix - mean & sem performance across 500 image patches using blurred image pixel spatial gradients.
%
% (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
justImPix.maxF_newMean_meanAccImgs = zeros(1,num_cPdist);
justImPix.maxF_newMean_semAccImgs = zeros(1,num_cPdist);
justImPix.maxF_newStd_meanAccImgs = zeros(1,num_cPdist);
%
justImPix.maxR_newMean_meanAccImgs = zeros(1,num_cPdist);
justImPix.maxR_newMean_semAccImgs = zeros(1,num_cPdist);
justImPix.maxR_newStd_meanAccImgs = zeros(1,num_cPdist);
%
justImPix.maxP_newMean_meanAccImgs = zeros(1,num_cPdist);
justImPix.maxP_newMean_semAccImgs = zeros(1,num_cPdist);
justImPix.maxP_newStd_meanAccImgs = zeros(1,num_cPdist);
%
% (4). New way to compute P-R-F - (Max value using best single ground truther.). 
justImPix.maxF_newMax_meanAccImgs = zeros(1,num_cPdist);
justImPix.maxF_newMax_semAccImgs = zeros(1,num_cPdist);
%
justImPix.maxR_newMax_meanAccImgs = zeros(1,num_cPdist);
justImPix.maxR_newMax_semAccImgs = zeros(1,num_cPdist);
%
justImPix.maxP_newMax_meanAccImgs = zeros(1,num_cPdist);
justImPix.maxP_newMax_semAccImgs = zeros(1,num_cPdist);







% relImPix - values of P/R/F of method relative to just using Raw Image Pixels.
%
% (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
relImPix.maxF_newMean_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
relImPix.maxF_newMean_semAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
relImPix.maxF_newStd_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
%
relImPix.maxR_newMean_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
relImPix.maxR_newMean_semAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
relImPix.maxR_newStd_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
%
relImPix.maxP_newMean_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
relImPix.maxP_newMean_semAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
relImPix.maxP_newStd_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));


% (4). New way to compute P-R-F - (Max value using best single ground truther.). 
relImPix.maxF_newMax_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
relImPix.maxF_newMax_semAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
%
relImPix.maxR_newMax_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
relImPix.maxR_newMax_semAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
%
relImPix.maxP_newMax_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
relImPix.maxP_newMax_semAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));




if(blur_flg)
    
    
    
    % justImBlur - mean & sem performance across 500 image patches using blurred image pixel spatial gradients.
    %
    % (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
    justImBlur.maxF_newMean_meanAccImgs = zeros(1,num_cPdist);
    justImBlur.maxF_newMean_semAccImgs = zeros(1,num_cPdist);
    justImBlur.maxF_newStd_meanAccImgs = zeros(1,num_cPdist);
    %
    justImBlur.maxR_newMean_meanAccImgs = zeros(1,num_cPdist);
    justImBlur.maxR_newMean_semAccImgs = zeros(1,num_cPdist);
    justImBlur.maxR_newStd_meanAccImgs = zeros(1,num_cPdist);
    %
    justImBlur.maxP_newMean_meanAccImgs = zeros(1,num_cPdist);
    justImBlur.maxP_newMean_semAccImgs = zeros(1,num_cPdist);
    justImBlur.maxP_newStd_meanAccImgs = zeros(1,num_cPdist);
    %
    % (4). New way to compute P-R-F - (Max value using best single ground truther.). 
    justImBlur.maxF_newMax_meanAccImgs = zeros(1,num_cPdist);
    justImBlur.maxF_newMax_semAccImgs = zeros(1,num_cPdist);
    %
    justImBlur.maxR_newMax_meanAccImgs = zeros(1,num_cPdist);
    justImBlur.maxR_newMax_semAccImgs = zeros(1,num_cPdist);
    %
    justImBlur.maxP_newMax_meanAccImgs = zeros(1,num_cPdist);
    justImBlur.maxP_newMax_semAccImgs = zeros(1,num_cPdist);

    
    
    
    % relImBlur - values of P/R/F of method relative to optimal Gaussian Blurring (sig=1).
    %
    % (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
    relImBlur.maxF_newMean_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
    relImBlur.maxF_newMean_semAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
    relImBlur.maxF_newStd_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
    %
    relImBlur.maxR_newMean_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
    relImBlur.maxR_newMean_semAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
    relImBlur.maxR_newStd_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
    %
    relImBlur.maxP_newMean_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
    relImBlur.maxP_newMean_semAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
    relImBlur.maxP_newStd_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
    %
    % (4). New way to compute P-R-F - (Max value using best single ground truther.). 
    relImBlur.maxF_newMax_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
    relImBlur.maxF_newMax_semAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
    %
    relImBlur.maxR_newMax_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
    relImBlur.maxR_newMax_semAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
    %
    relImBlur.maxP_newMax_meanAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
    relImBlur.maxP_newMax_semAccImgs = zeros(numel(rM),numel(ev),num_cPdist,numel(method));
end

numFiles = zeros(numel(rM),numel(ev),numel(method));




% matrix of parameter values for nice labeling later.
param_matrix = cell(numel(rM),numel(ev));
for i = 1:numel(rM)
    for j = 1:numel(ev)
        param_matrix{i,j} = ['rM',rM{i},':',ev{j}];
    end   
end             
              
              
        

plot_IsoF_flag = 0; % no longer plotting this for now.

              
              
              
for A = 1:numel(method)              
              
    method{A}
    
    for j = 1:numel(ev)            % loop over eigenvector configuration
        for i = 1:numel(rM)        % loop over max distance
            for d = 1:num_cPdist   % loop over correspondPixels distance allowance.

                disp( ['ev: ',num2str(j),' / ',num2str(numel(ev)),' , rM: ',num2str(i),' / ',num2str(numel(rM)),' , cP: ',num2str(d),' / ',num2str(num_cPdist)] )
                
                % Set up evalDir which points to all the benchmark results (ev2.txt files)
                if strcmp(method{A},'IsoDiff')
                    evalDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/spectral/',method{A},'/benchmark_results/rM',rM{i},'/',ev{j},'/'];
                else
                    evalDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/spectral/',method{A},'/benchmark_results/rM',rM{i},'/sDInf/sP0p2/',ev{j},'/'];
                end


                % Loop through each image patch and grab maxF. So later I can compute their mean and std.
                files = dir([evalDir,'*_d',num2str(d),'_ev2.txt']);
                %
                maxF_new_max = zeros(1,numel(files));    % For these next 3, we compute P/R/F for each groundtruth associated with an
                maxR_new_max = zeros(1,numel(files));    % image patch.  Here we just grab the max value.  What is the algorithm's
                maxP_new_max = zeros(1,numel(files));    % best match to any single human?
                maxF_bestGT  = zeros(1,numel(files));
                thr_new_max  = zeros(1,numel(files)); 
                %
                maxF_new_mean = zeros(1,numel(files));   % For a single image patch, the average P/R/F across the different groundtruths
                maxR_new_mean = zeros(1,numel(files));
                maxP_new_mean = zeros(1,numel(files));
                thr_new_mean  = zeros(1,numel(files)); 
                %
                maxF_new_std = zeros(1,numel(files));    % Standard Deviation of P/R/F across different groundtruths for single image ptch
                maxR_new_std = zeros(1,numel(files));
                maxP_new_std = zeros(1,numel(files));
                %
                numFiles(i,j,A) = numel(files);              % number of image patches processed for a given method or blur
                %
                for k = 1:numel(files) 
                    filename = fullfile(evalDir,files(k).name);
                    AA  = dlmread(filename);
                    thr = AA(:, 1); 
                    cP_d = AA(:, 2);
                    %
                    Rmean = AA(:, 3);
                    Rstd = AA(:, 4);
                    Rmax = AA(:, 5);
                    %
                    Pmean = AA(:, 6);
                    Pstd = AA(:, 7);
                    Pmax = AA(:, 8);
                    %
                    Fmean = AA(:, 9);
                    Fstd = AA(:, 10);
                    Fmax = AA(:, 11);
                    %
                    num_gTs = AA(:, 12);
                    Fmax_whichGTs = AA(:, 13);
                    %
                    maxF_new_max(k) = max(Fmax);
                    ind = find(Fmax==max(Fmax));
                    ind = ind(1);
                    maxR_new_max(k) = Rmax(ind);
                    maxP_new_max(k) = Pmax(ind);
                    maxF_bestGT(k) = Fmax_whichGTs(ind);
                    thr_new_max(k) = thr(ind);
                    %
                    maxF_new_mean(k) = max(Fmean);
                    ind = find(Fmean==max(Fmean));
                    ind = ind(1);
                    maxR_new_mean(k) = Rmean(ind);
                    maxP_new_mean(k) = Pmean(ind);
                    maxF_new_std(k) = Fstd(ind);   % std across groundtruths as threshold that gave max avgF value.
                    thr_new_mean(k) = thr(ind);
                    %
                    imgPtchName{k} = files(k).name(1:end-8);
                    %
                    %k
                end % loop over image patches
                
                
                
                
                
                
                
                
                % Compute statistics of this method & parameters relative to imPix.
                %
                % (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
                justMethod.maxF_newMean_meanAccImgs(i,j,d,A) = mean(maxF_new_mean);
                justMethod.maxF_newMean_semAccImgs(i,j,d,A) = sem(maxF_new_mean,which_errbars);
                justMethod.maxF_newStd_meanAccImgs(i,j,d,A) = mean(maxF_new_std);
                %
                justMethod.maxR_newMean_meanAccImgs(i,j,d,A) = mean(maxR_new_mean);
                justMethod.maxR_newMean_semAccImgs(i,j,d,A) = sem(maxR_new_mean,which_errbars);
                justMethod.maxR_newStd_meanAccImgs(i,j,d,A) = mean(maxR_new_std);
                %
                justMethod.maxP_newMean_meanAccImgs(i,j,d,A) = mean(maxP_new_mean);
                justMethod.maxP_newMean_semAccImgs(i,j,d,A) = sem(maxP_new_mean,which_errbars);
                justMethod.maxP_newStd_meanAccImgs(i,j,d,A) = mean(maxP_new_std);
                %
                % (4). New way to compute P-R-F - (Max value using best single ground truther.). 
                justMethod.maxF_newMax_meanAccImgs(i,j,d,A) = mean(maxF_new_max);
                justMethod.maxF_newMax_semAccImgs(i,j,d,A) = sem(maxF_new_max,which_errbars);
                %
                justMethod.maxR_newMax_meanAccImgs(i,j,d,A) = mean(maxR_new_max);
                justMethod.maxR_newMax_semAccImgs(i,j,d,A) = sem(maxR_new_max,which_errbars);
                %
                justMethod.maxP_newMax_meanAccImgs(i,j,d,A) = mean(maxP_new_max);
                justMethod.maxP_newMax_semAccImgs(i,j,d,A) = sem(maxP_new_max,which_errbars);
                
                
                
                
                
                
                


                % Using just raw image pixels.
                maxF_new_max_imPix = zeros(1,numel(files));    % For these next 3, we compute P/R/F for each groundtruth associated with an
                maxR_new_max_imPix = zeros(1,numel(files));    % image patch.  Here we just grab the max value.  What is the algorithm's
                maxP_new_max_imPix = zeros(1,numel(files));    % best match to any single human?
                maxF_bestGT_imPix  =  zeros(1,numel(files));
                thr_new_max_imPix  = zeros(1,numel(files));
                %
                maxF_new_mean_imPix = zeros(1,numel(files));   % For a single image patch, the average P/R/F across the different groundtruths
                maxR_new_mean_imPix = zeros(1,numel(files));
                maxP_new_mean_imPix = zeros(1,numel(files));
                thr_new_mean_imPix  = zeros(1,numel(files));
                %
                maxF_new_std_imPix = zeros(1,numel(files));    % Standard Deviation of P/R/F across different groundtruths for single image ptch
                maxR_new_std_imPix = zeros(1,numel(files));
                maxP_new_std_imPix = zeros(1,numel(files));
                %
                for k = 1:numel(files)       
                    filename = fullfile(imPixBenchmarkDir,files(k).name);
                    AA  = dlmread(filename);
                    thr = AA(:, 1); 
                    cP_d = AA(:, 2);
                    %
                    Rmean = AA(:, 3);
                    Rstd = AA(:, 4);
                    Rmax = AA(:, 5);
                    %
                    Pmean = AA(:, 6);
                    Pstd = AA(:, 7);
                    Pmax = AA(:, 8);
                    %
                    Fmean = AA(:, 9);
                    Fstd = AA(:, 10);
                    Fmax = AA(:, 11);
                    %
                    num_gTs = AA(:, 12);
                    Fmax_whichGTs = AA(:, 13);
                    %
                    maxF_new_max_imPix(k) = max(Fmax);
                    ind = find(Fmax==max(Fmax));
                    ind = ind(1);
                    maxR_new_max_imPix(k) = Rmax(ind);
                    maxP_new_max_imPix(k) = Pmax(ind);
                    maxF_bestGT_imPix(k) = Fmax_whichGTs(ind);
                    thr_new_max_imPix(k) = thr(ind);
                    %
                    maxF_new_mean_imPix(k) = max(Fmean);
                    ind = find(Fmean==max(Fmean));
                    ind = ind(1);
                    maxR_new_mean_imPix(k) = Rmean(ind);
                    maxP_new_mean_imPix(k) = Pmean(ind);
                    maxF_new_std_imPix(k) = Fstd(ind);   % std across groundtruths as threshold that gave max avgF value.
                    thr_new_mean_imPix(k) = thr(ind);
                    %
                    %k
                end
                %
                
                
                
                
                
                
                % justImPix - mean & sem performance across 500 image patches using just image pixel spatial gradients.
                %
                % (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
                if(A==1 && i==i && j==1)
                    justImPix.maxF_newMean_meanAccImgs(d) = mean(maxF_new_mean_imPix);
                    justImPix.maxF_newMean_semAccImgs(d) = sem(maxF_new_mean_imPix,which_errbars);
                    justImPix.maxF_newStd_meanAccImgs(d) = mean(maxF_new_std_imPix);
                    %
                    justImPix.maxR_newMean_meanAccImgs(d) = mean(maxR_new_mean_imPix);
                    justImPix.maxR_newMean_semAccImgs(d) = sem(maxR_new_mean_imPix,which_errbars);
                    justImPix.maxR_newStd_meanAccImgs(d) = mean(maxR_new_std_imPix);
                    %
                    justImPix.maxP_newMean_meanAccImgs(d) = mean(maxP_new_mean_imPix);
                    justImPix.maxP_newMean_semAccImgs(d) = sem(maxP_new_mean_imPix,which_errbars);
                    justImPix.maxP_newStd_meanAccImgs(d) = mean(maxP_new_std_imPix);
                    %
                    % (4). New way to compute P-R-F - (Max value using best single ground truther.). 
                    justImPix.maxF_newMax_meanAccImgs(d) = mean(maxF_new_max_imPix);
                    justImPix.maxF_newMax_semAccImgs(d) = sem(maxF_new_max_imPix,which_errbars);
                    %
                    justImPix.maxR_newMax_meanAccImgs(d) = mean(maxR_new_max_imPix);
                    justImPix.maxR_newMax_semAccImgs(d) = sem(maxR_new_max_imPix,which_errbars);
                    %
                    justImPix.maxP_newMax_meanAccImgs(d) = mean(maxP_new_max_imPix);
                    justImPix.maxP_newMax_semAccImgs(d) = sem(maxP_new_max_imPix,which_errbars);
                end
                
                
                
                
                
                % Compute statistics of this method & parameters relative to imPix.
                %
                % (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
                relImPix.maxF_newMean_meanAccImgs(i,j,d,A) = mean(maxF_new_mean-maxF_new_mean_imPix);
                relImPix.maxF_newMean_semAccImgs(i,j,d,A) = sem(maxF_new_mean-maxF_new_mean_imPix,which_errbars);
                relImPix.maxF_newStd_meanAccImgs(i,j,d,A) = mean(maxF_new_std-maxF_new_std_imPix);
                %
                relImPix.maxR_newMean_meanAccImgs(i,j,d,A) = mean(maxR_new_mean-maxF_new_mean_imPix);
                relImPix.maxR_newMean_semAccImgs(i,j,d,A) = sem(maxR_new_mean-maxR_new_mean_imPix,which_errbars);
                relImPix.maxR_newStd_meanAccImgs(i,j,d,A) = mean(maxR_new_std-maxR_new_std_imPix);
                %
                relImPix.maxP_newMean_meanAccImgs(i,j,d,A) = mean(maxP_new_mean-maxP_new_mean_imPix);
                relImPix.maxP_newMean_semAccImgs(i,j,d,A) = sem(maxP_new_mean-maxP_new_mean_imPix,which_errbars);
                relImPix.maxP_newStd_meanAccImgs(i,j,d,A) = mean(maxP_new_std-maxP_new_std_imPix);
                %
                % (4). New way to compute P-R-F - (Max value using best single ground truther.). 
                relImPix.maxF_newMax_meanAccImgs(i,j,d,A) = mean(maxF_new_max-maxF_new_max_imPix);
                relImPix.maxF_newMax_semAccImgs(i,j,d,A) = sem(maxF_new_max-maxF_new_max_imPix,which_errbars);
                %
                relImPix.maxR_newMax_meanAccImgs(i,j,d,A) = mean(maxR_new_max-maxR_new_max_imPix);
                relImPix.maxR_newMax_semAccImgs(i,j,d,A) = sem(maxR_new_max-maxR_new_max_imPix,which_errbars);
                %
                relImPix.maxP_newMax_meanAccImgs(i,j,d,A) = mean(maxP_new_max-maxP_new_max_imPix);
                relImPix.maxP_newMax_semAccImgs(i,j,d,A) = sem(maxP_new_max-maxP_new_max_imPix,which_errbars);
                
                
                
                
                
                % Using optimal gaussian blurring.
                maxF_new_max_imBlur = zeros(1,numel(files));    % For these next 3, we compute P/R/F for each groundtruth associated with an
                maxR_new_max_imBlur = zeros(1,numel(files));    % image patch.  Here we just grab the max value.  What is the algorithm's
                maxP_new_max_imBlur = zeros(1,numel(files));    % best match to any single human?
                maxF_bestGT_imBlur  =  zeros(1,numel(files));
                thr_new_max_imBlur  = zeros(1,numel(files)); 
                %
                maxF_new_mean_imBlur = zeros(1,numel(files));   % For a single image patch, the average P/R/F across the different groundtruths
                maxR_new_mean_imBlur = zeros(1,numel(files));
                maxP_new_mean_imBlur = zeros(1,numel(files));
                thr_new_mean_imBlur  = zeros(1,numel(files)); 
                %
                maxF_new_std_imBlur = zeros(1,numel(files));    % Standard Deviation of P/R/F across different groundtruths for single image ptch
                maxR_new_std_imBlur = zeros(1,numel(files));
                maxP_new_std_imBlur = zeros(1,numel(files));
                %
                for k = 1:numel(files)       
                    filename = fullfile(imBlurBenchmarkDir,files(k).name);
                    AA  = dlmread(filename);
                    thr = AA(:, 1); 
                    cP_d = AA(:, 2);
                    %
                    Rmean = AA(:, 3);
                    Rstd = AA(:, 4);
                    Rmax = AA(:, 5);
                    %
                    Pmean = AA(:, 6);
                    Pstd = AA(:, 7);
                    Pmax = AA(:, 8);
                    %
                    Fmean = AA(:, 9);
                    Fstd = AA(:, 10);
                    Fmax = AA(:, 11);
                    %
                    num_gTs = AA(:, 12);
                    Fmax_whichGTs = AA(:, 13);
                    %
                    maxF_new_max_imBlur(k) = max(Fmax);
                    ind = find(Fmax==max(Fmax));
                    ind = ind(1);
                    maxR_new_max_imBlur(k) = Rmax(ind);
                    maxP_new_max_imBlur(k) = Pmax(ind);
                    maxF_bestGT_imBlur(k) = Fmax_whichGTs(ind);
                    thr_new_max_imBlur(k) = thr(ind);
                    %
                    maxF_new_mean_imBlur(k) = max(Fmean);
                    ind = find(Fmean==max(Fmean));
                    ind = ind(1);
                    maxR_new_mean_imBlur(k) = Rmean(ind);
                    maxP_new_mean_imBlur(k) = Pmean(ind);
                    maxF_new_std_imBlur(k) = Fstd(ind);   % std across groundtruths as threshold that gave max avgF value.
                    thr_new_mean_imBlur(k) = thr(ind);
                    %
                    %k
                end
                
                
                
                
                
                % justImBlur - mean & sem performance across 500 image patches using blurred image pixel spatial gradients.
                %
                % (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
                if(A==1 && i==i && j==1)
                    justImBlur.maxF_newMean_meanAccImgs(d) = mean(maxF_new_mean_imBlur);
                    justImBlur.maxF_newMean_semAccImgs(d) = sem(maxF_new_mean_imBlur,which_errbars);
                    justImBlur.maxF_newStd_meanAccImgs(d) = mean(maxF_new_std_imBlur);
                    %
                    justImBlur.maxR_newMean_meanAccImgs(d) = mean(maxR_new_mean_imBlur);
                    justImBlur.maxR_newMean_semAccImgs(d) = sem(maxR_new_mean_imBlur,which_errbars);
                    justImBlur.maxR_newStd_meanAccImgs(d) = mean(maxR_new_std_imBlur);
                    %
                    justImBlur.maxP_newMean_meanAccImgs(d) = mean(maxP_new_mean_imBlur);
                    justImBlur.maxP_newMean_semAccImgs(d) = sem(maxP_new_mean_imBlur,which_errbars);
                    justImBlur.maxP_newStd_meanAccImgs(d) = mean(maxP_new_std_imBlur);
                    %
                    % (4). New way to compute P-R-F - (Max value using best single ground truther.). 
                    justImBlur.maxF_newMax_meanAccImgs(d) = mean(maxF_new_max_imBlur);
                    justImBlur.maxF_newMax_semAccImgs(d) = sem(maxF_new_max_imBlur,which_errbars);
                    %
                    justImBlur.maxR_newMax_meanAccImgs(d) = mean(maxR_new_max_imBlur);
                    justImBlur.maxR_newMax_semAccImgs(d) = sem(maxR_new_max_imBlur,which_errbars);
                    %
                    justImBlur.maxP_newMax_meanAccImgs(d) = mean(maxP_new_max_imBlur);
                    justImBlur.maxP_newMax_semAccImgs(d) = sem(maxP_new_max_imBlur,which_errbars);
                end

                
                
                
                
                
                % Compute statistics of this method & parameters relative to imBlur.
                %
                % (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
                relImBlur.maxF_newMean_meanAccImgs(i,j,d,A) = mean(maxF_new_mean-maxF_new_mean_imBlur);
                relImBlur.maxF_newMean_semAccImgs(i,j,d,A) = sem(maxF_new_mean-maxF_new_mean_imBlur,which_errbars);
                relImBlur.maxF_newStd_meanAccImgs(i,j,d,A) = mean(maxF_new_std-maxF_new_std_imBlur);
                %
                relImBlur.maxR_newMean_meanAccImgs(i,j,d,A) = mean(maxR_new_mean-maxR_new_mean_imBlur);
                relImBlur.maxR_newMean_semAccImgs(i,j,d,A) = sem(maxR_new_mean-maxR_new_mean_imBlur,which_errbars);
                relImBlur.maxR_newStd_meanAccImgs(i,j,d,A) = mean(maxR_new_std-maxR_new_std_imBlur);
                %
                relImBlur.maxP_newMean_meanAccImgs(i,j,d,A) = mean(maxP_new_mean-maxP_new_mean_imBlur);
                relImBlur.maxP_newMean_semAccImgs(i,j,d,A) = sem(maxP_new_mean-maxP_new_mean_imBlur,which_errbars);
                relImBlur.maxP_newStd_meanAccImgs(i,j,d,A) = mean(maxP_new_std-maxP_new_std_imBlur);
                %
                % (4). New way to compute P-R-F - (Max value using best single ground truther.). 
                relImBlur.maxF_newMax_meanAccImgs(i,j,d,A) = mean(maxF_new_max-maxF_new_max_imBlur);
                relImBlur.maxF_newMax_semAccImgs(i,j,d,A) = sem(maxF_new_max-maxF_new_max_imBlur,which_errbars);
                %
                relImBlur.maxR_newMax_meanAccImgs(i,j,d,A) = mean(maxR_new_max-maxR_new_max_imBlur);
                relImBlur.maxR_newMax_semAccImgs(i,j,d,A) = sem(maxR_new_max-maxR_new_max_imBlur,which_errbars);
                %
                relImBlur.maxP_newMax_meanAccImgs(i,j,d,A) = mean(maxP_new_max-maxP_new_max_imBlur);
                relImBlur.maxP_newMax_semAccImgs(i,j,d,A) = sem(maxP_new_max-maxP_new_max_imBlur,which_errbars);

                %[d,i,j]
    
        	end % loop over d of cpD
        end % loop over rM parameter
    end % loop over ev parameter
    
    
end % loop over method

which_eig_combo = zeros(1,numel(method));
for m = 1:numel(method)
    which_eig_combo(m) =  find( mean(squeeze(justMethod.maxF_newMax_meanAccImgs(4,:,:,m)),2) == max(mean(squeeze(justMethod.maxF_newMax_meanAccImgs(4,:,:,m)),2)) );
    best_legend{m} = [method_short{m},' (',ev_long{which_eig_combo(m)},')'];
end



% Save into a mat file all the things I plot in 4 & 5 above so I can plot them with Kuramoto stuff.
save(['Eig_plot_results',blur_tag_M,'_',which_errbars], 'justImPix','justImBlur','justMethod','relImPix','relImBlur','method_short','ev',...
    'which_errbars','blur_tag_M','blur_tit','which_eig_combo','best_legend')





% Free up some memory space
clear maxF_bestGT  maxF_bestGT_imBlur maxF_bestGT_imPix ...
      maxF_new_max maxF_new_max_imBlur maxF_new_max_imPix ...
      maxF_new_mean maxF_new_mean_imBlur maxF_new_mean_imPix ...
      maxF_new_std maxF_new_std_imBlur maxF_new_std_imPix ...
      maxP_new_max maxP_new_max_imBlur maxP_new_max_imPix ...
      maxP_new_mean maxP_new_mean_imBlur maxP_new_mean_imPix ...
      maxP_new_std maxP_new_std_imBlur maxP_new_std_imPix ...
      maxR_new_max maxR_new_max_imBlur maxR_new_max_imPix ... 
      maxR_new_mean maxR_new_mean_imBlur maxR_new_mean_imPix ...
      maxR_new_std maxR_new_std_imBlur maxR_new_std_imPix ...
      files thr_new_max thr_new_max_imBlur thr_new_max_imPix ...
      thr_new_mean thr_new_mean_imBlur thr_new_mean_imPix 



  
  
  
  
  
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  

%% Now, make some plots. Look in plot_eval_compare_bsdsBench_KurB.


% (1). Plot F-measure (mean & sem errorbars) for justMethod for all methods, all Evecs, all methods, all dt, all rM.
if(0)
    H1=figure; 
    ha = tight_subplot(numel(method), num_cPdist, [0.03 0.01], [0.07 0.05], [0.05 0.01]);

    for m = 1:numel(method) % {GLnrm AAnrm Mod_SKHAdj Mod_N&G IsoDiff}

        maxY = max(max(max( justMethod.maxF_newMax_meanAccImgs(:,:,:,m)+justMethod.maxF_newMax_semAccImgs(:,:,:,m)   )));
        minY = min(min(min( justMethod.maxF_newMax_meanAccImgs(:,:,:,m)-justMethod.maxF_newMax_semAccImgs(:,:,:,m)   )));

        for d = 1:num_cPdist

            subplot(ha(d + (m-1)*num_cPdist)), hold on
            %errorbar(justMethod.maxF_newMean_meanAccImgs(:,:,d,m), justMethod.maxF_newMean_semAccImgs(:,:,d,m), 'LineWidth',2,'LineStyle','--')
            errorbar(justMethod.maxF_newMax_meanAccImgs(:,:,d,m), justMethod.maxF_newMax_semAccImgs(:,:,d,m), 'LineWidth',2,'LineStyle','-')
            %
            if(m==1)
                title(['(d_t = ',num2str(d),')'],'FontSize',16,'FontWeight','Bold')
            end
            %
            if(d==1)
                if(m==5)
                    legend(ev)
                    ylabel('F','FontSize',16,'FontWeight','Bold')
                    xlabel('rM','FontSize',16,'FontWeight','Bold')
                end
                text(1,0.9*maxY,[method_short{m}],'FontSize',20,'FontWeight','Bold','VerticalAlignment','Top')
                % set(gca,'YTick',1:5,'YTickLabel',linspace(minY,maxY,5))
            end
            %
            grid on
            ylim([minY maxY])
            %
            if(m==5)
                set(gca,'XTick',1:numel(rM),'XTickLabel',rM)
            else
                set(gca,'XTick',1:numel(rM),'XTickLabel',[])
            end
            %
            set(gca,'YTick',linspace(minY,maxY,4))
            n=get(gca,'YTick');
            if(d==1)
                set(gca,'Yticklabel',sprintf('%.2f |',n'));
            else
                set(gca,'Yticklabel',[]);
            end
            %
            set(gca,'FontSize',14,'FontWeight','Bold')


        end

    end
    %
    annotation('textbox', [0 0.9 1 0.1],'String',[ '[Just Method - Spectral ',blur_tit,']' ], ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',18,'FontWeight','Bold')
    %
    disp([outImgDir,'compareSpectralF_everything_justMethod_errbars_',which_errbars,blur_tag_M,'.svg'])
    plot2svg([outImgDir,'compareSpectralF_everything_justMethod_errbars_',which_errbars,blur_tag_M,'.svg'],H1)
    %
    saveGoodImg(H1,[outImgDir,'compareSpectralF_everything_justMethod_errbars_',which_errbars,blur_tag_M,'.jpg'],sizeGoodIm)
    close(H1)
end






% (2). Plot F-measure (mean & sem errorbars) for relImPix for all methods, all Evecs, all methods, all dt, all rM.
if(0)
    H1=figure; 
    ha = tight_subplot(numel(method), num_cPdist, [0.03 0.01], [0.07 0.05], [0.05 0.01]);

    for m = 1:numel(method) % {GLnrm AAnrm Mod_SKHAdj Mod_N&G IsoDiff}

        maxY = max(max(max( relImPix.maxF_newMax_meanAccImgs(:,:,:,m)+relImPix.maxF_newMax_semAccImgs(:,:,:,m)   )));
        minY = min(min(min( relImPix.maxF_newMax_meanAccImgs(:,:,:,m)-relImPix.maxF_newMax_semAccImgs(:,:,:,m)   )));

        for d = 1:num_cPdist

            subplot(ha(d + (m-1)*num_cPdist)), hold on
            %errorbar(relImPix.maxF_newMean_meanAccImgs(:,:,d,m), relImPix.maxF_newMean_semAccImgs(:,:,d,m), 'LineWidth',2,'LineStyle','--')
            errorbar(relImPix.maxF_newMax_meanAccImgs(:,:,d,m), relImPix.maxF_newMax_semAccImgs(:,:,d,m), 'LineWidth',2,'LineStyle','-')
            %
            if(m==1)
                title(['(d_t = ',num2str(d),')'],'FontSize',16,'FontWeight','Bold')
            end
            %
            if(d==1)
                if(m==5)
                    legend(ev)
                    ylabel('\Delta F_i','FontSize',16,'FontWeight','Bold')
                    xlabel('rM','FontSize',16,'FontWeight','Bold')
                end
                text(1,0.9*maxY,[method_short{m}],'FontSize',20,'FontWeight','Bold','VerticalAlignment','Top')
                % set(gca,'YTick',1:5,'YTickLabel',linspace(minY,maxY,5))
            end
            %
            grid on
            ylim([minY maxY])
            %
            if(m==5)
                set(gca,'XTick',1:numel(rM),'XTickLabel',rM)
            else
                set(gca,'XTick',1:numel(rM),'XTickLabel',[])
            end
            %
            set(gca,'YTick',linspace(minY,maxY,4))
            n=get(gca,'YTick');
            if(d==1)
                set(gca,'Yticklabel',sprintf('%.2f |',n'));
            else
                set(gca,'Yticklabel',[]);
            end
            %
            set(gca,'FontSize',14,'FontWeight','Bold')


        end

    end
    %
    annotation('textbox', [0 0.9 1 0.1],'String',[ '[rel Im Pix - Spectral ',blur_tit,']' ], ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',18,'FontWeight','Bold')
    %
    disp([outImgDir,'compareSpectralF_everything_relImPix_errbars_',which_errbars,blur_tag_M,'.svg'])
    plot2svg([outImgDir,'compareSpectralF_everything_relImPix_errbars_',which_errbars,blur_tag_M,'.svg'],H1)
    %
    saveGoodImg(H1,[outImgDir,'compareSpectralF_everything_relImPix_errbars_',which_errbars,blur_tag_M,'.jpg'],sizeGoodIm)
    close(H1)
end



% (3). Plot F-measure (mean & sem errorbars) for relImBlur for all methods, all Evecs, all methods, all dt, all rM.
if(0)    
    H1=figure; 
    ha = tight_subplot(numel(method), num_cPdist, [0.03 0.01], [0.07 0.05], [0.05 0.01]);

    for m = 1:numel(method) % {GLnrm AAnrm Mod_SKHAdj Mod_N&G IsoDiff}

        maxY = max(max(max( relImBlur.maxF_newMax_meanAccImgs(:,:,:,m)+relImBlur.maxF_newMax_semAccImgs(:,:,:,m)   )));
        minY = min(min(min( relImBlur.maxF_newMax_meanAccImgs(:,:,:,m)-relImBlur.maxF_newMax_semAccImgs(:,:,:,m)   )));

        for d = 1:num_cPdist

            subplot(ha(d + (m-1)*num_cPdist)), hold on
            %errorbar(relImBlur.maxF_newMean_meanAccImgs(:,:,d,m), relImBlur.maxF_newMean_semAccImgs(:,:,d,m), 'LineWidth',2,'LineStyle','--')
            errorbar(relImBlur.maxF_newMax_meanAccImgs(:,:,d,m), relImBlur.maxF_newMax_semAccImgs(:,:,d,m), 'LineWidth',2,'LineStyle','-')
            %
            if(m==1)
                title(['(d_t = ',num2str(d),')'],'FontSize',16,'FontWeight','Bold')
            end
            %
            if(d==1)
                if(m==5)
                    legend(ev)
                    ylabel('\Delta F_b','FontSize',16,'FontWeight','Bold')
                    xlabel('rM','FontSize',16,'FontWeight','Bold')
                end
                text(1,0.9*maxY,[method_short{m}],'FontSize',20,'FontWeight','Bold','VerticalAlignment','Top')
                % set(gca,'YTick',1:5,'YTickLabel',linspace(minY,maxY,5))
            end
            %
            grid on
            ylim([minY maxY])
            %
            if(m==5)
                set(gca,'XTick',1:numel(rM),'XTickLabel',rM)
            else
                set(gca,'XTick',1:numel(rM),'XTickLabel',[])
            end
            %
            set(gca,'YTick',linspace(minY,maxY,4))
            n=get(gca,'YTick');
            if(d==1)
                set(gca,'Yticklabel',sprintf('%.2f |',n'));
            else
                set(gca,'Yticklabel',[]);
            end
            %
            set(gca,'FontSize',14,'FontWeight','Bold')


        end

    end
    %
    annotation('textbox', [0 0.9 1 0.1],'String',[ '[rel Im Blur - Spectral ',blur_tit,']' ], ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',18,'FontWeight','Bold')
    %
    disp([outImgDir,'compareSpectralF_everything_relImBlur_errbars_',which_errbars,blur_tag_M,'.svg'])
    plot2svg([outImgDir,'compareSpectralF_everything_relImBlur_errbars_',which_errbars,blur_tag_M,'.svg'],H1)
    %
    saveGoodImg(H1,[outImgDir,'compareSpectralF_everything_relImBlur_errbars_',which_errbars,blur_tag_M,'.jpg'],sizeGoodIm)
    close(H1)
end


% % % % % %



%% COMBINE 4 & 5.


% (4).  Plot F-measure (mean & sem errbars) for [justMethod (x5), justImPix, justImBlur]
%       for each Evec combination separately. Take best rM. Across all dt or cpDist values.
if(1)
    H1=figure; 
    ha = tight_subplot(2,2, [0.03 0.01], [0.07 0.05], [0.05 0.01]); % make one subplot for each eigenvector combination.

    switch which_errbars
        case 'sem'
            maxY = 0.4;
            minY = 0;
        case 'std'
            maxY = 0.55;
            minY = 0;
    end

    for v = 1:3

        subplot(ha(v)), hold on

        % performance using just imPix
        errorbar([1:4]+0.1*(-1), justImPix.maxF_newMax_meanAccImgs, justImPix.maxF_newMax_semAccImgs, 'LineWidth',3,'LineStyle','-','Color','magenta')
        %
        % performance using just imBlur
        errorbar([1:4]+0.1*(-2), justImBlur.maxF_newMax_meanAccImgs, justImBlur.maxF_newMax_semAccImgs, 'LineWidth',3,'LineStyle','-','Color','cyan')
        %
        % Show maxGT for methods, imPix and imBlur
        for m = 1:numel(method) % {GLnrm AAnrm Mod_SKHAdj Mod_N&G IsoDiff}
            errorbar([1:4]+0.1*(m-1), justMethod.maxF_newMax_meanAccImgs(4,v,:,m), justMethod.maxF_newMax_semAccImgs(4,v,:,m), 'LineWidth',3,'LineStyle','-','Color',method_color{m})
        end
        %
        % if show meanGT in addition to maxGT results.
        if(1)
            % performance using just imPix
            errorbar([1:4]+0.1*(-1), justImPix.maxF_newMean_meanAccImgs, justImPix.maxF_newMean_semAccImgs, 'LineWidth',3,'LineStyle','--','Color','magenta')
            %
            % performance using just imBlur
            errorbar([1:4]+0.1*(-2), justImBlur.maxF_newMean_meanAccImgs, justImBlur.maxF_newMean_semAccImgs, 'LineWidth',3,'LineStyle','--','Color','cyan')
            %
            for m = 1:numel(method) % {GLnrm AAnrm Mod_SKHAdj Mod_N&G IsoDiff}
                errorbar([1:4]+0.1*(m-1), justMethod.maxF_newMean_meanAccImgs(4,v,:,m), justMethod.maxF_newMean_semAccImgs(4,v,:,m), 'LineWidth',3,'LineStyle','--','Color',method_color{m})
            end
        end

        if(v==1)
           xlabel('d_t','FontSize',18,'FontWeight','Bold')
           ylabel('F','FontSize',18,'FontWeight','Bold')
        end

        title(ev_long{v},'FontSize',18,'FontWeight','Bold')
        ylim([minY maxY])
        grid on
        set(gca,'XTick',1:4,'XTickLabel',{'0','1','1.4','2'},'YTick',0:0.1:0.4,'FontSize',16,'FontWeight','Bold')

        axis square

    end

    %legend('Pix','Blur',method_short{:},'Location','SouthEast')
    %
    %annotation('textbox', [0 0.9 1 0.1],'String',['Spectral ',blur_tit,' and rM = ',rM{4}], ...
    %       'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')
    %
    % disp([outImgDir,'compareSpectralF_Evecs123_errbars_',which_errbars,blur_tag_M,'.svg'])
    % plot2svg([outImgDir,'compareSpectralF_Evecs123_errbars_',which_errbars,blur_tag_M,'.svg'],H1)
    % %
    % saveGoodImg(H1,[outImgDir,'compareSpectralF_Evecs123_errbars_',which_errbars,blur_tag_M,'.jpg'],sizeGoodIm)
    % close(H1)


    % (5).  Plot F-measure (mean & sem errbars) for [justMethod (x5), justImPix, justImBlur]
    %       for best Evec combination for each method. Take best rM. Across all dt or cpDist values.

    switch which_errbars
        case 'sem'
            maxY = 0.4;
            minY = 0;
        case 'std'
            maxY = 0.55;
            minY = 0;
    end



    subplot(ha(4)), hold on

    % performance using just imPix
    errorbar([1:4]+0.1*(-1), justImPix.maxF_newMax_meanAccImgs, justImPix.maxF_newMax_semAccImgs, 'LineWidth',3,'LineStyle','-','Color','magenta')
    %
    % performance using just imBlur
    errorbar([1:4]+0.1*(-2), justImBlur.maxF_newMax_meanAccImgs, justImBlur.maxF_newMax_semAccImgs, 'LineWidth',3,'LineStyle','-','Color','cyan')
    %
    % Show maxGT for methods, imPix and imBlur
    for m = 1:numel(method) % {GLnrm AAnrm Mod_SKHAdj Mod_N&G IsoDiff}
        errorbar([1:4]+0.1*(m-1), justMethod.maxF_newMax_meanAccImgs(4,which_eig_combo(m),:,m), justMethod.maxF_newMax_semAccImgs(4,which_eig_combo(m),:,m), 'LineWidth',3,'LineStyle','-','Color',method_color{m})
    end


    % if show meanGT in addition to maxGT results.
    if(1)
        % performance using just imPix
        errorbar([1:4]+0.1*(-1), justImPix.maxF_newMean_meanAccImgs, justImPix.maxF_newMean_semAccImgs, 'LineWidth',3,'LineStyle','--','Color','magenta')
        %
        % performance using just imBlur
        errorbar([1:4]+0.1*(-2), justImBlur.maxF_newMean_meanAccImgs, justImBlur.maxF_newMean_semAccImgs, 'LineWidth',3,'LineStyle','--','Color','cyan')
        %
        for m = 1:numel(method) % {GLnrm AAnrm Mod_SKHAdj Mod_N&G IsoDiff}
            errorbar([1:4]+0.1*(m-1), justMethod.maxF_newMean_meanAccImgs(4,v,:,m), justMethod.maxF_newMean_semAccImgs(4,v,:,m), 'LineWidth',3,'LineStyle','--','Color',method_color{m})
        end
        %
    end

    legend('Pix','Blur',best_legend{:},'Location','SouthEast')

    %xlabel('d_t','FontSize',18,'FontWeight','Bold')
    %ylabel('F','FontSize',18,'FontWeight','Bold')

    title(['Best ev. combo'],'FontSize',18,'FontWeight','Bold')
    ylim([minY maxY])
    grid on
    set(gca,'XTick',1:4,'XTickLabel',{'0','1','1.4','2'},'YTick',0:0.1:0.4,'FontSize',16,'FontWeight','Bold')
    axis square
    %
    disp([outImgDir,'compareSpectralF_bestEvecs_errbars_',which_errbars,blur_tag_M,'.svg'])
    plot2svg([outImgDir,'compareSpectralF_bestEvecs_errbars_',which_errbars,blur_tag_M,'.svg'],H1)
    %
    saveGoodImg(H1,[outImgDir,'compareSpectralF_bestEvecs_errbars_',which_errbars,blur_tag_M,'.jpg'],sizeGoodIm)
    close(H1)
end












keyboard