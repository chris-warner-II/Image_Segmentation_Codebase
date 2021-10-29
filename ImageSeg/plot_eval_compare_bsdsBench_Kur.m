% This script/function written by CW 11/15 to plot multiple precision
% recall curves on the same isoF plot.  We can use this to compare
% different methods or different parameter settings for the same method.


[dirPre,sizeGoodIm] = onCluster;
addpath([dirPre,'images/BSDS_images/BSR/bench/benchmarks/'])


method = {'IsoDiff','AAnrm','GLnrm','Mod_SKHAdj','Mod_N&G'}; % 'IsoDiff','AAnrm','GLnrm','Mod_SKHAdj','Mod_N&G';

% labels for plotting later.
which_F_computation = {'benchmarkF','benchmarkF','benchmarkF','meanGT','meanGT','meanGT','maxGT','maxGT','maxGT','unionGT','unionGT','unionGT'};
relative_to_what = {'justMethod','relImPix','relImBlur','justMethod','relImPix','relImBlur','justMethod','relImPix','relImBlur','justMethod','relImPix','relImBlur'};
    


rM = {'1','3','5','10'};
colors =  {'r','b','c','m'};

ks =  {'sml','mid','lrg'}; 
shapes = {'v-','d-','^-'};

plot_max = 1; % set to 1 to plot average of Precision & Recall leading to largest F, regardless of threshold
              % set to 0 to plot average of Precision & Recall at each threshold value for each image
              

% Flag to use a gaussian kernel (sig=2) to preblur image before running network Kuramoto computation.
blur_flg = 1;
if(blur_flg)
    blur_tag_M = '_blur_sig2';
    blur_tit = 'w/ Pre-Blurring (\sigma=2)';
else
    blur_tag_M = ''; % if we are not blurring.
    blur_tit = 'w/ No Pre-Blurring.';
end
blur_tag_I = 'blur_sz21_sig2/';


% directory to ground truth mat files.
gTdir = [dirPre,'images/BSDS_patch/101x101_ds1/groundTruth/'];

% directory to image patches mat files.
imDir = [dirPre,'images/BSDS_patch/101x101_ds1/'];

% directory to probabalistic boundary png images for imPix & Blurring
pbDirBlur = [dirPre,'images/BSDS_patch/101x101_ds1/',blur_tag_I,'pb_png/'];
pbDirIm = [dirPre,'images/BSDS_patch/101x101_ds1/pb_png/'];





%% Preallocate memory for arrays to hold statistics (across image patches) for each different blur.


% % Way it was being done previously.    
F_all = zeros(numel(rM),numel(ks),numel(method));
% F_max = zeros(numel(method)+2,1);
% 
% mean_maxF = zeros(numel(rM),numel(ks),numel(method));
% std_maxF = zeros(numel(rM),numel(ks),numel(method));
% 
% mean_maxF_relBlur = zeros(numel(rM),numel(ks),numel(method));
% std_maxF_relBlur = zeros(numel(rM),numel(ks),numel(method));
% %
% mean_maxF_relPix = zeros(numel(rM),numel(ks),numel(method));
% std_maxF_relPix = zeros(numel(rM),numel(ks),numel(method));
% 
% numFiles = zeros(numel(rM),numel(ks),numel(method));
%
% % data structure that will take F-measure for method, for method-ImBlur, and for method-ImPix
% maxF_struct(i,j,A).method
% maxF_struct(i,j,A).method_relBlur
% maxF_struct(i,j,A).method_relPix
%
% % Now I am doing this below...


% For just the method - the absolute value. Not comparing to any strawman or null model.

% (1). The way it was originally implemented in the benchmark code - Not Self Consistent.
justMethod.maxF_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
justMethod.maxF_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
justMethod.maxF_skewAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
justMethod.maxR_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
justMethod.maxR_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
justMethod.maxP_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
justMethod.maxP_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));


% (2). Our reiteration of the original implementation (should be same as 1)
justMethod.maxF_old_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
justMethod.maxF_old_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
justMethod.maxF_old_skewAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
justMethod.maxR_old_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
justMethod.maxR_old_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
justMethod.maxP_old_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
justMethod.maxP_old_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));


% (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
justMethod.maxF_newMean_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
justMethod.maxF_newMean_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
justMethod.maxF_newMean_skewAccImgs = zeros(numel(rM),numel(ks),numel(method));
justMethod.maxF_newStd_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
justMethod.maxR_newMean_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
justMethod.maxR_newMean_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
justMethod.maxR_newStd_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
justMethod.maxP_newMean_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
justMethod.maxP_newMean_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
justMethod.maxP_newStd_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));


% (4). New way to compute P-R-F - (Max value using best single ground truther.). 
justMethod.maxF_newMax_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
justMethod.maxF_newMax_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
justMethod.maxF_newMax_skewAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
justMethod.maxR_newMax_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
justMethod.maxR_newMax_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
justMethod.maxP_newMax_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
justMethod.maxP_newMax_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));


% (5). New way to compute R & F - Kinda matching what benchmark was
% doing to compute Precision originally. Unioning all Ground Truthers
% into one GroundTruth and computing P & R & F using that unioned GT.
justMethod.maxF_unionGT_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
justMethod.maxF_unionGT_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
justMethod.maxF_unionGT_skewAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
justMethod.maxR_unionGT_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
justMethod.maxR_unionGT_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
justMethod.maxP_unionGT_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
justMethod.maxP_unionGT_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));





% relImPix - values of P/R/F of method relative to just using Raw Image Pixels.

% (1). The way it was originally implemented in the benchmark code - Not Self Consistent.
relImPix.maxF_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImPix.maxF_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImPix.maxF_skewAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
relImPix.maxR_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImPix.maxR_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
relImPix.maxP_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImPix.maxP_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));


% (2). Our reiteration of the original implementation (should be same as 1)
relImPix.maxF_old_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImPix.maxF_old_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImPix.maxF_old_skewAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
relImPix.maxR_old_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImPix.maxR_old_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
relImPix.maxP_old_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImPix.maxP_old_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));


% (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
relImPix.maxF_newMean_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImPix.maxF_newMean_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImPix.maxF_newMean_skewAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImPix.maxF_newStd_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
relImPix.maxR_newMean_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImPix.maxR_newMean_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImPix.maxR_newStd_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
relImPix.maxP_newMean_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImPix.maxP_newMean_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImPix.maxP_newStd_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));


% (4). New way to compute P-R-F - (Max value using best single ground truther.). 
relImPix.maxF_newMax_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImPix.maxF_newMax_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImPix.maxF_newMax_skewAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
relImPix.maxR_newMax_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImPix.maxR_newMax_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
relImPix.maxP_newMax_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImPix.maxP_newMax_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));


% (5). New way to compute R & F - Kinda matching what benchmark was
% doing to compute Precision originally. Unioning all Ground Truthers
% into one GroundTruth and computing P & R & F using that unioned GT.
relImPix.maxF_unionGT_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImPix.maxF_unionGT_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImPix.maxF_unionGT_skewAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
relImPix.maxR_unionGT_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImPix.maxR_unionGT_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
relImPix.maxP_unionGT_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImPix.maxP_unionGT_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));




% relImBlur - values of P/R/F of method relative to optimal Gaussian Blurring (sig=2).

% (1). The way it was originally implemented in the benchmark code - Not Self Consistent.
relImBlur.maxF_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImBlur.maxF_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImBlur.maxF_skewAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
relImBlur.maxR_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImBlur.maxR_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
relImBlur.maxP_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImBlur.maxP_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));


% (2). Our reiteration of the original implementation (should be same as 1)
relImBlur.maxF_old_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImBlur.maxF_old_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImBlur.maxF_old_skewAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
relImBlur.maxR_old_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImBlur.maxR_old_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
relImBlur.maxP_old_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImBlur.maxP_old_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));


% (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
relImBlur.maxF_newMean_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImBlur.maxF_newMean_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImBlur.maxF_newMean_skewAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImBlur.maxF_newStd_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
relImBlur.maxR_newMean_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImBlur.maxR_newMean_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImBlur.maxR_newStd_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
relImBlur.maxP_newMean_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImBlur.maxP_newMean_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImBlur.maxP_newStd_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));


% (4). New way to compute P-R-F - (Max value using best single ground truther.). 
relImBlur.maxF_newMax_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImBlur.maxF_newMax_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImBlur.maxF_newMax_skewAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
relImBlur.maxR_newMax_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImBlur.maxR_newMax_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
relImBlur.maxP_newMax_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImBlur.maxP_newMax_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));


% (5). New way to compute R & F - Kinda matching what benchmark was
% doing to compute Precision originally. Unioning all Ground Truthers
% into one GroundTruth and computing P & R & F using that unioned GT.
relImBlur.maxF_unionGT_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImBlur.maxF_unionGT_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImBlur.maxF_unionGT_skewAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
relImBlur.maxR_unionGT_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImBlur.maxR_unionGT_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));
%
relImBlur.maxP_unionGT_meanAccImgs = zeros(numel(rM),numel(ks),numel(method));
relImBlur.maxP_unionGT_stdAccImgs = zeros(numel(rM),numel(ks),numel(method));




numFiles = zeros(numel(rM),numel(ks),numel(method));






% Preallocate memory to store mean & std & skew of each method with optimal
% parameter settings for different ways of computing F-measure.
maxF_meanAccImgs = zeros(numel(method),numel(which_F_computation)); 
maxR_meanAccImgs = zeros(numel(method),numel(which_F_computation)); 
maxP_meanAccImgs = zeros(numel(method),numel(which_F_computation)); 
%
maxF_stdAccImgs = zeros(numel(method),numel(which_F_computation)); 
maxR_stdAccImgs = zeros(numel(method),numel(which_F_computation)); 
maxP_stdAccImgs = zeros(numel(method),numel(which_F_computation)); 
%
maxF_skewAccImgs = zeros(numel(method),numel(which_F_computation)); 
    


gT_agreement = cell(numel(method),numel(which_F_computation));
gT_agreement2 = cell(numel(method),numel(which_F_computation));

%% Loop through Different Methods & Parameter settings to calculate statistics across all available image patches on
% Precision/Recall/F-measure calculated in a couple different ways and on those things relative to ImPix and Optimal ImBlur.


% matrix of parameter values for nice labeling later.
param_matrix = cell(numel(rM),numel(ks));
for i = 1:numel(rM)
    for j = 1:numel(ks)
        param_matrix{i,j} = ['rM',rM{i},':ks',ks{j}];
    end   
end



plot_IsoF_flag = 0; % no longer plotting this for now.


              
for A = 1:numel(method)

    cntr = 0;

    H=open('isoF.fig');


    for i = 1:numel(rM)
        for j = 1:numel(ks)
            
            cntr=cntr+1;


            if strcmp(method{A},'IsoDiff')
                evalDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/Kur_PIF_Fourier1/',method{A},'/benchmark_results/rM',rM{i},'/NF_60_0/ks',ks{j},'/'];
            else
                evalDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/Kur_PIF_Fourier1/',method{A},'/benchmark_results/rM',rM{i},'/sDInf/sP0p2/NF_60_0/ks',ks{j},'/'];
            end


            [H, evalRes] = plot_eval(evalDir,[colors{i},shapes{j}],H,plot_max); % 

            % Note: evalRes = [bestT,bestR,bestP,bestF,R_max,P_max,F_max,Area_PR]
            disp(['rM:',rM{i},' ks:',ks{j}])
            disp(['      Same TH -- R:',num2str(evalRes(2),2),' P:',num2str(evalRes(3),2),' F:',num2str(evalRes(4),2)])
            disp(['       Any TH -- R:',num2str(evalRes(5),2),' P:',num2str(evalRes(6),2),' F:',num2str(evalRes(7),2)])

            F_all(i,j,A) = evalRes(7);
            
            

            % Loop through each image patch and grab maxF. So later I can compute their mean and std.
            files = dir([evalDir,'*_ev1.txt']);

            maxF = zeros(1,numel(files));            % original benchmark way of computing P/R/F - not self consistent
            maxR = zeros(1,numel(files)); 
            maxP = zeros(1,numel(files)); 
            thr1 = zeros(1,numel(files)); 

            maxF_old = zeros(1,numel(files));        % Our reimplementation of benchmark P/R/F (should be same as above)
            maxR_old = zeros(1,numel(files));
            maxP_old = zeros(1,numel(files));
            thr_old  = zeros(1,numel(files)); 

            maxF_unionGT = zeros(1,numel(files));    % Take union of all K groundtruths as single groundtruth before computing P/R/F
            maxR_unionGT = zeros(1,numel(files));
            maxP_unionGT = zeros(1,numel(files));
            thr_unionGT  = zeros(1,numel(files)); 

            maxF_new_max = zeros(1,numel(files));    % For these next 3, we compute P/R/F for each groundtruth associated with an
            maxR_new_max = zeros(1,numel(files));    % image patch.  Here we just grab the max value.  What is the algorithm's
            maxP_new_max = zeros(1,numel(files));    % best match to any single human?
            maxF_bestGT  = zeros(1,numel(files));
            thr_new_max  = zeros(1,numel(files)); 

            maxF_new_mean = zeros(1,numel(files));   % For a single image patch, the average P/R/F across the different groundtruths
            maxR_new_mean = zeros(1,numel(files));
            maxP_new_mean = zeros(1,numel(files));
            thr_new_mean  = zeros(1,numel(files)); 
            %
            maxF_new_std = zeros(1,numel(files));    % Standard Deviation of P/R/F across different groundtruths for single image ptch
            maxR_new_std = zeros(1,numel(files));
            maxP_new_std = zeros(1,numel(files));

            numFiles(i,j,A) = numel(files);              % number of image patches processed for a given method or blur
            
            for k = 1:numel(files)

                filename = fullfile(evalDir,files(k).name);
                AA  = dlmread(filename);
                thr = AA(:, 1); 
                %cntR = AA(:, 2);
                sumR = AA(:, 3);
                %cntP = AA(:, 4);
                sumP = AA(:, 5);
                %
                cntR0 = AA(:, 6);                      % CW: added these new ways to compute pixel correspondence in evaluation_bdry_image.
                cntP0 = AA(:, 7);

                sumR0 = AA(:, 8);

                R_old = AA(:, 9);
                P_old = AA(:, 10);

                R_new_mean = AA(:, 11);
                R_new_std = AA(:, 12);
                R_new_max = AA(:, 13);

                P_new_mean = AA(:, 14);
                P_new_std = AA(:, 15);
                P_new_max = AA(:, 16);
                
                F_new_meanA = AA(:, 17);
                F_new_stdA = AA(:, 18);
                F_new_maxA = AA(:, 19);

                R_unionGT = AA(:, 20);

                num_gTs = AA(:, 21);
                Fmax_whichGTs = AA(:, 22);

                R = cntR0 ./ (sumR + (sumR==0)); % Note: this (sumR==0) in the denominator is just to ensure we do not
                P = cntP0 ./ (sumP + (sumP==0)); %       divide by zero. If sumR/P == 0, then we divide by 1. (CW).

                F = 2*P.*R./(P+R+((P+R)==0));    % F-measure.

                F_old = 2*P_old.*R_old./(P_old+R_old+((P_old+R_old)==0));

                F_new_max = 2*P_new_max.*R_new_max./(P_new_max+R_new_max+((P_new_max+R_new_max)==0)); 

                F_new_mean = 2*P_new_mean.*R_new_mean./(P_new_mean+R_new_mean+((P_new_mean+R_new_mean)==0)); 

                F_new_std = 2*P_new_std.*R_new_std./(P_new_std+R_new_std+((P_new_std+R_new_std)==0)); 

                F_unionGT = 2*P.*R_unionGT./(P+R_unionGT+((P+R_unionGT)==0));



                % Collect max F (computed different ways) vs threshold into long vector.
                maxF(k) = max(F);
                ind = find(F==max(F));
                ind = ind(1);
                maxR(k) = R(ind);
                maxP(k) = P(ind);
                thr1(k) = thr(ind);



                maxF_old(k) = max(F_old);
                ind = find(F_old==max(F_old));
                ind = ind(1);
                maxR_old(k) = R_old(ind);
                maxP_old(k) = P_old(ind);
                thr_old(k) = thr(ind);

                maxF_new_max(k) = max(F_new_max);
                ind = find(F_new_max==max(F_new_max));
                ind = ind(1);
                maxR_new_max(k) = R_new_max(ind);
                maxP_new_max(k) = P_new_max(ind);
                maxF_bestGT(k) = Fmax_whichGTs(ind);
                thr_new_max(k) = thr(ind);



                maxF_new_mean(k) = max(F_new_meanA);
                ind = find(F_new_meanA==max(F_new_meanA));
                ind = ind(1);
                maxR_new_mean(k) = R_new_mean(ind);
                maxP_new_mean(k) = P_new_mean(ind);
                maxF_new_std(k) = F_new_stdA(ind);   % std across groundtruths as threshold that gave max avgF value.
                thr_new_mean(k) = thr(ind);


                maxF_unionGT(k) = max(F_unionGT);
                ind = find(F_unionGT==max(F_unionGT));
                ind = ind(1);
                maxR_unionGT(k) = R_unionGT(ind);
                maxP_unionGT(k) = P(ind);
                thr_unionGT(k) = thr(ind);

                imgPtchName{k} = files(k).name(1:end-8);


                k

            end % looping through image patches
            
            
            
            % Take Mean & Std across image patches of Precision,Recall,F-measure computed in different ways.

            % (1). The way it was originally implemented in the benchmark code - Not Self Consistent.
            justMethod.maxF_meanAccImgs(i,j,A) = mean(maxF);
            justMethod.maxF_stdAccImgs(i,j,A) = std(maxF);
            justMethod.maxF_skewAccImgs(i,j,A) = skewness(maxF);
            %
            justMethod.maxR_meanAccImgs(i,j,A) = mean(maxR);
            justMethod.maxR_stdAccImgs(i,j,A) = std(maxR);
            %
            justMethod.maxP_meanAccImgs(i,j,A) = mean(maxP);
            justMethod.maxP_stdAccImgs(i,j,A) = std(maxP);



            % (2). Our reiteration of the original implementation (should be same as 1)
            justMethod.maxF_old_meanAccImgs(i,j,A) = mean(maxF_old);
            justMethod.maxF_old_stdAccImgs(i,j,A) = std(maxF_old);
            justMethod.maxF_old_skewAccImgs(i,j,A) = skewness(maxF_old);
            %
            justMethod.maxR_old_meanAccImgs(i,j,A) = mean(maxR_old);
            justMethod.maxR_old_stdAccImgs(i,j,A) = std(maxR_old);
            %
            justMethod.maxP_old_meanAccImgs(i,j,A) = mean(maxP_old);
            justMethod.maxP_old_stdAccImgs(i,j,A) = std(maxP_old);


            % (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
            justMethod.maxF_newMean_meanAccImgs(i,j,A) = mean(maxF_new_mean);
            justMethod.maxF_newMean_stdAccImgs(i,j,A) = std(maxF_new_mean);
            justMethod.maxF_newMean_skewAccImgs(i,j,A) = skewness(maxF_new_mean);
            justMethod.maxF_newStd_meanAccImgs(i,j,A) = mean(maxF_new_std);
            %
            justMethod.maxR_newMean_meanAccImgs(i,j,A) = mean(maxR_new_mean);
            justMethod.maxR_newMean_stdAccImgs(i,j,A) = std(maxR_new_mean);
            justMethod.maxR_newStd_meanAccImgs(i,j,A) = mean(maxR_new_std);
            %
            justMethod.maxP_newMean_meanAccImgs(i,j,A) = mean(maxP_new_mean);
            justMethod.maxP_newMean_stdAccImgs(i,j,A) = std(maxP_new_mean);
            justMethod.maxP_newStd_meanAccImgs(i,j,A) = mean(maxP_new_std);


            % (4). New way to compute P-R-F - (Max value using best single ground truther.). 
            justMethod.maxF_newMax_meanAccImgs(i,j,A) = mean(maxF_new_max);
            justMethod.maxF_newMax_stdAccImgs(i,j,A) = std(maxF_new_max);
            justMethod.maxF_newMax_skewAccImgs(i,j,A) = skewness(maxF_new_max);
            %
            justMethod.maxR_newMax_meanAccImgs(i,j,A) = mean(maxR_new_max);
            justMethod.maxR_newMax_stdAccImgs(i,j,A) = std(maxR_new_max);
            %
            justMethod.maxP_newMax_meanAccImgs(i,j,A) = mean(maxP_new_max);
            justMethod.maxP_newMax_stdAccImgs(i,j,A) = std(maxP_new_max);


            % (5). New way to compute R & F - Kinda matching what benchmark was
            % doing to compute Precision originally. Unioning all Ground Truthers
            % into one GroundTruth and computing P & R & F using that unioned GT.
            justMethod.maxF_unionGT_meanAccImgs(i,j,A) = mean(maxF_unionGT);
            justMethod.maxF_unionGT_stdAccImgs(i,j,A) = std(maxF_unionGT);
            justMethod.maxF_unionGT_skewAccImgs(i,j,A) = skewness(maxF_unionGT);
            %
            justMethod.maxR_unionGT_meanAccImgs(i,j,A) = mean(maxR_unionGT);
            justMethod.maxR_unionGT_stdAccImgs(i,j,A) = std(maxR_unionGT);
            %
            justMethod.maxP_unionGT_meanAccImgs(i,j,A) = mean(maxP_unionGT);
            justMethod.maxP_unionGT_stdAccImgs(i,j,A) = std(maxP_unionGT);
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            %% Compute different Precision/Recall values for Raw Image Pixels
            % Currently, doing this inside method & params for loops because I dont know how else to line up ImPix & ImBlur to have 
            % statistics relative to those.
            imPixDir = [dirPre,'images/BSDS_patch/101x101_ds1/benchmark_results/'];
            
            
            
            maxF_imPix = zeros(1,numel(files));            % original benchmark way of computing P/R/F - not self consistent
            maxR_imPix = zeros(1,numel(files)); 
            maxP_imPix = zeros(1,numel(files)); 
            thr1_imPix = zeros(1,numel(files));
            
            maxF_old_imPix = zeros(1,numel(files));        % Our reimplementation of benchmark P/R/F (should be same as above)
            maxR_old_imPix = zeros(1,numel(files));
            maxP_old_imPix = zeros(1,numel(files));
            thr_old_imPix  = zeros(1,numel(files));

            maxF_unionGT_imPix = zeros(1,numel(files));    % Take union of all K groundtruths as single groundtruth before computing P/R/F
            maxR_unionGT_imPix = zeros(1,numel(files));
            maxP_unionGT_imPix = zeros(1,numel(files));
            thr_unionGT_imPix  = zeros(1,numel(files));

            maxF_new_max_imPix = zeros(1,numel(files));    % For these next 3, we compute P/R/F for each groundtruth associated with an
            maxR_new_max_imPix = zeros(1,numel(files));    % image patch.  Here we just grab the max value.  What is the algorithm's
            maxP_new_max_imPix = zeros(1,numel(files));    % best match to any single human?
            maxF_bestGT_imPix  =  zeros(1,numel(files));
            thr_new_max_imPix  = zeros(1,numel(files));

            maxF_new_mean_imPix = zeros(1,numel(files));   % For a single image patch, the average P/R/F across the different groundtruths
            maxR_new_mean_imPix = zeros(1,numel(files));
            maxP_new_mean_imPix = zeros(1,numel(files));
            thr_new_mean_imPix  = zeros(1,numel(files));
            %
            maxF_new_std_imPix = zeros(1,numel(files));    % Standard Deviation of P/R/F across different groundtruths for single image ptch
            maxR_new_std_imPix = zeros(1,numel(files));
            maxP_new_std_imPix = zeros(1,numel(files));

            %numFiles(i,j,A) = numel(files);              % number of image patches processed for a given method or blur

            for k = 1:numel(files)

                filename = fullfile(imPixDir,files(k).name);
                AA  = dlmread(filename);
                thr = AA(:, 1); 
                %cntR = AA(:, 2);
                sumR = AA(:, 3);
                %cntP = AA(:, 4);
                sumP = AA(:, 5);
                %
                cntR0 = AA(:, 6);                      % CW: added these new ways to compute pixel correspondence in evaluation_bdry_image.
                cntP0 = AA(:, 7);

                sumR0 = AA(:, 8);

                R_old = AA(:, 9);
                P_old = AA(:, 10);

                R_new_mean = AA(:, 11);
                R_new_std = AA(:, 12);
                R_new_max = AA(:, 13);

                P_new_mean = AA(:, 14);
                P_new_std = AA(:, 15);
                P_new_max = AA(:, 16);
                
                
                F_new_meanA = AA(:, 17);
                F_new_stdA = AA(:, 18);
                F_new_maxA = AA(:, 19);

                R_unionGT = AA(:, 20);

                num_gTs = AA(:, 21);
                Fmax_whichGTs = AA(:, 22);
                

                R = cntR0 ./ (sumR + (sumR==0)); % Note: this (sumR==0) in the denominator is just to ensure we do not
                P = cntP0 ./ (sumP + (sumP==0)); %       divide by zero. If sumR/P == 0, then we divide by 1. (CW).

                F = 2*P.*R./(P+R+((P+R)==0));    % F-measure.

                F_old = 2*P_old.*R_old./(P_old+R_old+((P_old+R_old)==0));

                F_new_max = 2*P_new_max.*R_new_max./(P_new_max+R_new_max+((P_new_max+R_new_max)==0)); 

                F_new_mean = 2*P_new_mean.*R_new_mean./(P_new_mean+R_new_mean+((P_new_mean+R_new_mean)==0)); 

                F_new_std = 2*P_new_std.*R_new_std./(P_new_std+R_new_std+((P_new_std+R_new_std)==0)); 

                F_unionGT = 2*P.*R_unionGT./(P+R_unionGT+((P+R_unionGT)==0));


                % Collect max F (computed different ways) vs threshold into long vector.
                maxF_imPix(k) = max(F);
                ind = find(F==max(F));
                ind = ind(1);
                maxR_imPix(k) = R(ind);
                maxP_imPix(k) = P(ind);
                thr1_imPix(k) = thr(ind);



                maxF_old_imPix(k) = max(F_old);
                ind = find(F_old==max(F_old));
                ind = ind(1);
                maxR_old_imPix(k) = R_old(ind);
                maxP_old_imPix(k) = P_old(ind);
                thr_old_imPix(k) = thr(ind);


                maxF_new_max_imPix(k) = max(F_new_max);
                ind = find(F_new_max==max(F_new_max));
                ind = ind(1);
                maxR_new_max_imPix(k) = R_new_max(ind);
                maxP_new_max_imPix(k) = P_new_max(ind);
                maxF_bestGT_imPix(k) = Fmax_whichGTs(ind);
                thr_new_max_imPix(k) = thr(ind);
                


                maxF_new_mean_imPix(k) = max(F_new_meanA);
                ind = find(F_new_meanA==max(F_new_meanA));
                ind = ind(1);
                maxR_new_mean_imPix(k) = R_new_mean(ind);
                maxP_new_mean_imPix(k) = P_new_mean(ind);
                maxF_new_std_imPix(k) = F_new_stdA(ind);   % std across groundtruths as threshold that gave max avgF value.
                thr_new_mean_imPix(k) = thr(ind);


                maxF_unionGT_imPix(k) = max(F_unionGT);
                ind = find(F_unionGT==max(F_unionGT));
                ind = ind(1);
                maxR_unionGT_imPix(k) = R_unionGT(ind);
                maxP_unionGT_imPix(k) = P(ind);
                thr_unionGT_imPix(k) = thr(ind);




                k

            end % looping through image patches
            
            
            
            
            % Compute statistics of this method & parameters relative to imPix.
            
            % (1). The way it was originally implemented in the benchmark code - Not Self Consistent.
            relImPix.maxF_meanAccImgs(i,j,A) = mean(maxF-maxF_imPix);
            relImPix.maxF_stdAccImgs(i,j,A) = std(maxF-maxF_imPix);
            relImPix.maxF_skewAccImgs(i,j,A) = skewness(maxF-maxF_imPix);
            %
            relImPix.maxR_meanAccImgs(i,j,A) = mean(maxR-maxR_imPix);
            relImPix.maxR_stdAccImgs(i,j,A) = std(maxR-maxR_imPix);
            %
            relImPix.maxP_meanAccImgs(i,j,A) = mean(maxP-maxP_imPix);
            relImPix.maxP_stdAccImgs(i,j,A) = std(maxP-maxP_imPix);



            % (2). Our reiteration of the original implementation (should be same as 1)
            relImPix.maxF_old_meanAccImgs(i,j,A) = mean(maxF_old-maxF_old_imPix);
            relImPix.maxF_old_stdAccImgs(i,j,A) = std(maxF_old-maxF_old_imPix);
            relImPix.maxF_old_skewAccImgs(i,j,A) = skewness(maxF_old-maxF_old_imPix);
            %
            relImPix.maxR_old_meanAccImgs(i,j,A) = mean(maxR_old-maxR_old_imPix);
            relImPix.maxR_old_stdAccImgs(i,j,A) = std(maxR_old-maxR_old_imPix);
            %
            relImPix.maxP_old_meanAccImgs(i,j,A) = mean(maxP_old-maxP_old_imPix);
            relImPix.maxP_old_stdAccImgs(i,j,A) = std(maxP_old-maxP_old_imPix);


            % (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
            relImPix.maxF_newMean_meanAccImgs(i,j,A) = mean(maxF_new_mean-maxF_new_mean_imPix);
            relImPix.maxF_newMean_stdAccImgs(i,j,A) = std(maxF_new_mean-maxF_new_mean_imPix);
            relImPix.maxF_newMean_skewAccImgs(i,j,A) = skewness(maxF_new_mean-maxF_new_mean_imPix);
            relImPix.maxF_newStd_meanAccImgs(i,j,A) = mean(maxF_new_std-maxF_new_std_imPix);
            %
            relImPix.maxR_newMean_meanAccImgs(i,j,A) = mean(maxR_new_mean-maxF_new_mean_imPix);
            relImPix.maxR_newMean_stdAccImgs(i,j,A) = std(maxR_new_mean-maxR_new_mean_imPix);
            relImPix.maxR_newStd_meanAccImgs(i,j,A) = mean(maxR_new_std-maxR_new_std_imPix);
            %
            relImPix.maxP_newMean_meanAccImgs(i,j,A) = mean(maxP_new_mean-maxP_new_mean_imPix);
            relImPix.maxP_newMean_stdAccImgs(i,j,A) = std(maxP_new_mean-maxP_new_mean_imPix);
            relImPix.maxP_newStd_meanAccImgs(i,j,A) = mean(maxP_new_std-maxP_new_std_imPix);


            % (4). New way to compute P-R-F - (Max value using best single ground truther.). 
            relImPix.maxF_newMax_meanAccImgs(i,j,A) = mean(maxF_new_max-maxF_new_max_imPix);
            relImPix.maxF_newMax_stdAccImgs(i,j,A) = std(maxF_new_max-maxF_new_max_imPix);
            relImPix.maxF_newMax_skewAccImgs(i,j,A) = skewness(maxF_new_max-maxF_new_max_imPix);
            %
            relImPix.maxR_newMax_meanAccImgs(i,j,A) = mean(maxR_new_max-maxR_new_max_imPix);
            relImPix.maxR_newMax_stdAccImgs(i,j,A) = std(maxR_new_max-maxR_new_max_imPix);
            %
            relImPix.maxP_newMax_meanAccImgs(i,j,A) = mean(maxP_new_max-maxP_new_max_imPix);
            relImPix.maxP_newMax_stdAccImgs(i,j,A) = std(maxP_new_max-maxP_new_max_imPix);


            % (5). New way to compute R & F - Kinda matching what benchmark was
            % doing to compute Precision originally. Unioning all Ground Truthers
            % into one GroundTruth and computing P & R & F using that unioned GT.
            relImPix.maxF_unionGT_meanAccImgs(i,j,A) = mean(maxF_unionGT-maxF_unionGT_imPix);
            relImPix.maxF_unionGT_stdAccImgs(i,j,A) = std(maxF_unionGT-maxF_unionGT_imPix);
            relImPix.maxF_unionGT_skewAccImgs(i,j,A) = skewness(maxF_unionGT-maxF_unionGT_imPix);
            %
            relImPix.maxR_unionGT_meanAccImgs(i,j,A) = mean(maxR_unionGT-maxR_unionGT_imPix);
            relImPix.maxR_unionGT_stdAccImgs(i,j,A) = std(maxR_unionGT-maxR_unionGT_imPix);
            %
            relImPix.maxP_unionGT_meanAccImgs(i,j,A) = mean(maxP_unionGT-maxP_unionGT_imPix);
            relImPix.maxP_unionGT_stdAccImgs(i,j,A) = std(maxP_unionGT-maxP_unionGT_imPix);
            
            
            
            
            
            
            %% Compute different Precision/Recall values for Optimal Blurring
            imBlurDir = [dirPre,'images/BSDS_patch/101x101_ds1/',blur_tag_I,'benchmark_results/'];
            
            
            
            maxF_imBlur = zeros(1,numel(files));            % original benchmark way of computing P/R/F - not self consistent
            maxR_imBlur = zeros(1,numel(files)); 
            maxP_imBlur = zeros(1,numel(files)); 
            thr1_imBlur = zeros(1,numel(files)); 

            maxF_old_imBlur = zeros(1,numel(files));        % Our reimplementation of benchmark P/R/F (should be same as above)
            maxR_old_imBlur = zeros(1,numel(files));
            maxP_old_imBlur = zeros(1,numel(files));
            thr_old_imBlur  = zeros(1,numel(files)); 

            maxF_unionGT_imBlur = zeros(1,numel(files));    % Take union of all K groundtruths as single groundtruth before computing P/R/F
            maxR_unionGT_imBlur = zeros(1,numel(files));
            maxP_unionGT_imBlur = zeros(1,numel(files));
            thr_unionGT_imBlur  = zeros(1,numel(files)); 

            maxF_new_max_imBlur = zeros(1,numel(files));    % For these next 3, we compute P/R/F for each groundtruth associated with an
            maxR_new_max_imBlur = zeros(1,numel(files));    % image patch.  Here we just grab the max value.  What is the algorithm's
            maxP_new_max_imBlur = zeros(1,numel(files));    % best match to any single human?
            maxF_bestGT_imBlur  =  zeros(1,numel(files));
            thr_new_max_imBlur  = zeros(1,numel(files)); 

            maxF_new_mean_imBlur = zeros(1,numel(files));   % For a single image patch, the average P/R/F across the different groundtruths
            maxR_new_mean_imBlur = zeros(1,numel(files));
            maxP_new_mean_imBlur = zeros(1,numel(files));
            thr_new_mean_imBlur  = zeros(1,numel(files)); 

            maxF_new_std_imBlur = zeros(1,numel(files));    % Standard Deviation of P/R/F across different groundtruths for single image ptch
            maxR_new_std_imBlur = zeros(1,numel(files));
            maxP_new_std_imBlur = zeros(1,numel(files));

            %numFiles(i,j,A) = numel(files);              % number of image patches processed for a given method or blur

            for k = 1:numel(files)

                filename = fullfile(imBlurDir,files(k).name);
                AA  = dlmread(filename);
                thr = AA(:, 1); 
                %cntR = AA(:, 2);
                sumR = AA(:, 3);
                %cntP = AA(:, 4);
                sumP = AA(:, 5);
                %
                cntR0 = AA(:, 6);                      % CW: added these new ways to compute pixel correspondence in evaluation_bdry_image.
                cntP0 = AA(:, 7);

                sumR0 = AA(:, 8);

                R_old = AA(:, 9);
                P_old = AA(:, 10);

                R_new_mean = AA(:, 11);
                R_new_std = AA(:, 12);
                R_new_max = AA(:, 13);

                P_new_mean = AA(:, 14);
                P_new_std = AA(:, 15);
                P_new_max = AA(:, 16);
                
                
                F_new_meanA = AA(:, 17);
                F_new_stdA = AA(:, 18);
                F_new_maxA = AA(:, 19);

                R_unionGT = AA(:, 20);

                num_gTs = AA(:, 21);
                Fmax_whichGTs = AA(:, 22);

                R = cntR0 ./ (sumR + (sumR==0)); % Note: this (sumR==0) in the denominator is just to ensure we do not
                P = cntP0 ./ (sumP + (sumP==0)); %       divide by zero. If sumR/P == 0, then we divide by 1. (CW).

                F = 2*P.*R./(P+R+((P+R)==0));    % F-measure.

                F_old = 2*P_old.*R_old./(P_old+R_old+((P_old+R_old)==0));

                F_new_max = 2*P_new_max.*R_new_max./(P_new_max+R_new_max+((P_new_max+R_new_max)==0)); 

                F_new_mean = 2*P_new_mean.*R_new_mean./(P_new_mean+R_new_mean+((P_new_mean+R_new_mean)==0)); 

                F_new_std = 2*P_new_std.*R_new_std./(P_new_std+R_new_std+((P_new_std+R_new_std)==0)); 

                F_unionGT = 2*P.*R_unionGT./(P+R_unionGT+((P+R_unionGT)==0));


                % Collect max F (computed different ways) vs threshold into long vector.
                maxF_imBlur(k) = max(F);
                ind = find(F==max(F));
                ind = ind(1);
                maxR_imBlur(k) = R(ind);
                maxP_imBlur(k) = P(ind);
                thr1_imBlur(k) = thr(ind);



                maxF_old_imBlur(k) = max(F_old);
                ind = find(F_old==max(F_old));
                ind = ind(1);
                maxR_old_imBlur(k) = R_old(ind);
                maxP_old_imBlur(k) = P_old(ind);
                thr_old_imBlur(k) = thr(ind);


                maxF_new_max_imBlur(k) = max(F_new_max);
                ind = find(F_new_max==max(F_new_max));
                ind = ind(1);
                maxR_new_max_imBlur(k) = R_new_max(ind);
                maxP_new_max_imBlur(k) = P_new_max(ind);
                maxF_bestGT_imBlur(k) = Fmax_whichGTs(ind);
                thr_new_max_imBlur(k) = thr(ind);



                maxF_new_mean_imBlur(k) = max(F_new_meanA);
                ind = find(F_new_meanA==max(F_new_meanA));
                ind = ind(1);
                maxR_new_mean_imBlur(k) = R_new_mean(ind);
                maxP_new_mean_imBlur(k) = P_new_mean(ind);
                maxF_new_std_imBlur(k) = F_new_stdA(ind);   % std across groundtruths as threshold that gave max avgF value.
                thr_new_mean_imBlur(k) = thr(ind);


                maxF_unionGT_imBlur(k) = max(F_unionGT);
                ind = find(F_unionGT==max(F_unionGT));
                ind = ind(1);
                maxR_unionGT_imBlur(k) = R_unionGT(ind);
                maxP_unionGT_imBlur(k) = P(ind);
                thr_unionGT_imBlur(k) = thr(ind);




                k

            end % looping through image patches
            
            
            
            
            
            
            % Compute statistics of this method & parameters relative to imBlur.
            
            % (1). The way it was originally implemented in the benchmark code - Not Self Consistent.
            relImBlur.maxF_meanAccImgs(i,j,A) = mean(maxF-maxF_imBlur);
            relImBlur.maxF_stdAccImgs(i,j,A) = std(maxF-maxF_imBlur);
            relImBlur.maxF_skewAccImgs(i,j,A) = skewness(maxF-maxF_imBlur);
            %
            relImBlur.maxR_meanAccImgs(i,j,A) = mean(maxR-maxR_imBlur);
            relImBlur.maxR_stdAccImgs(i,j,A) = std(maxR-maxR_imBlur);
            %
            relImBlur.maxP_meanAccImgs(i,j,A) = mean(maxP-maxP_imBlur);
            relImBlur.maxP_stdAccImgs(i,j,A) = std(maxP-maxP_imBlur);



            % (2). Our reiteration of the original implementation (should be same as 1)
            relImBlur.maxF_old_meanAccImgs(i,j,A) = mean(maxF_old-maxF_old_imBlur);
            relImBlur.maxF_old_stdAccImgs(i,j,A) = std(maxF_old-maxF_old_imBlur);
            relImBlur.maxF_old_skewAccImgs(i,j,A) = skewness(maxF_old-maxF_old_imBlur);
            %
            relImBlur.maxR_old_meanAccImgs(i,j,A) = mean(maxR_old-maxR_old_imBlur);
            relImBlur.maxR_old_stdAccImgs(i,j,A) = std(maxR_old-maxR_old_imBlur);
            %
            relImBlur.maxP_old_meanAccImgs(i,j,A) = mean(maxP_old-maxP_old_imBlur);
            relImBlur.maxP_old_stdAccImgs(i,j,A) = std(maxP_old-maxP_old_imBlur);


            % (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
            relImBlur.maxF_newMean_meanAccImgs(i,j,A) = mean(maxF_new_mean-maxF_new_mean_imBlur);
            relImBlur.maxF_newMean_stdAccImgs(i,j,A) = std(maxF_new_mean-maxF_new_mean_imBlur);
            relImBlur.maxF_newMean_skewAccImgs(i,j,A) = skewness(maxF_new_mean-maxF_new_mean_imBlur);
            relImBlur.maxF_newStd_meanAccImgs(i,j,A) = mean(maxF_new_std-maxF_new_std_imBlur);
            %
            relImBlur.maxR_newMean_meanAccImgs(i,j,A) = mean(maxR_new_mean-maxR_new_mean_imBlur);
            relImBlur.maxR_newMean_stdAccImgs(i,j,A) = std(maxR_new_mean-maxR_new_mean_imBlur);
            relImBlur.maxR_newStd_meanAccImgs(i,j,A) = mean(maxR_new_std-maxR_new_std_imBlur);
            %
            relImBlur.maxP_newMean_meanAccImgs(i,j,A) = mean(maxP_new_mean-maxP_new_mean_imBlur);
            relImBlur.maxP_newMean_stdAccImgs(i,j,A) = std(maxP_new_mean-maxP_new_mean_imBlur);
            relImBlur.maxP_newStd_meanAccImgs(i,j,A) = mean(maxP_new_std-maxP_new_std_imBlur);


            % (4). New way to compute P-R-F - (Max value using best single ground truther.). 
            relImBlur.maxF_newMax_meanAccImgs(i,j,A) = mean(maxF_new_max-maxF_new_max_imBlur);
            relImBlur.maxF_newMax_stdAccImgs(i,j,A) = std(maxF_new_max-maxF_new_max_imBlur);
            relImBlur.maxF_newMax_skewAccImgs(i,j,A) = skewness(maxF_new_max-maxF_new_max_imBlur);
            %
            relImBlur.maxR_newMax_meanAccImgs(i,j,A) = mean(maxR_new_max-maxR_new_max_imBlur);
            relImBlur.maxR_newMax_stdAccImgs(i,j,A) = std(maxR_new_max-maxR_new_max_imBlur);
            %
            relImBlur.maxP_newMax_meanAccImgs(i,j,A) = mean(maxP_new_max-maxP_new_max_imBlur);
            relImBlur.maxP_newMax_stdAccImgs(i,j,A) = std(maxP_new_max-maxP_new_max_imBlur);


            % (5). New way to compute R & F - Kinda matching what benchmark was
            % doing to compute Precision originally. Unioning all Ground Truthers
            % into one GroundTruth and computing P & R & F using that unioned GT.
            relImBlur.maxF_unionGT_meanAccImgs(i,j,A) = mean(maxF_unionGT-maxF_unionGT_imBlur);
            relImBlur.maxF_unionGT_stdAccImgs(i,j,A) = std(maxF_unionGT-maxF_unionGT_imBlur);
            relImBlur.maxF_unionGT_skewAccImgs(i,j,A) = skewness(maxF_unionGT-maxF_unionGT_imBlur);
            %
            relImBlur.maxR_unionGT_meanAccImgs(i,j,A) = mean(maxR_unionGT-maxR_unionGT_imBlur);
            relImBlur.maxR_unionGT_stdAccImgs(i,j,A) = std(maxR_unionGT-maxR_unionGT_imBlur);
            %
            relImBlur.maxP_unionGT_meanAccImgs(i,j,A) = mean(maxP_unionGT-maxP_unionGT_imBlur);
            relImBlur.maxP_unionGT_stdAccImgs(i,j,A) = std(maxP_unionGT-maxP_unionGT_imBlur);
            
            
            
            
            
            
            
            
            
            
            
            
            
            % Create a data structure which holds all F-measure values (not just mean and std) so we can make Violin plots.
            maxF_old_struct_method_only{i,j,A} = maxF_old';
            maxF_old_struct_imBlur_only{i,j,A} = maxF_old_imBlur'; 
            maxF_old_struct_imPix_only{i,j,A}  = maxF_old_imPix';   
            %
            maxR_old_struct_method_only{i,j,A} = maxR_old';
            maxR_old_struct_imBlur_only{i,j,A} = maxR_old_imBlur';
            maxR_old_struct_imPix_only{i,j,A}  = maxR_old_imPix';
            %
            maxP_old_struct_method_only{i,j,A} = maxP_old';
            maxP_old_struct_imBlur_only{i,j,A} = maxP_old_imBlur';
            maxP_old_struct_imPix_only{i,j,A}  = maxP_old_imPix';
            %
            thr_old_struct_method_only{i,j,A} = thr_old';
            thr_old_struct_imBlur_only{i,j,A} = thr_old_imBlur';
            thr_old_struct_imPix_only{i,j,A}  = thr_old_imPix';
            %
            % %
            %
            maxF_new_mean_struct_method_only{i,j,A} = maxF_new_mean';
            maxF_new_mean_struct_imBlur_only{i,j,A} = maxF_new_mean_imBlur';
            maxF_new_mean_struct_imPix_only{i,j,A} = maxF_new_mean_imPix';
            %
            maxR_new_mean_struct_method_only{i,j,A} = maxR_new_mean';
            maxR_new_mean_struct_imBlur_only{i,j,A} = maxR_new_mean_imBlur';
            maxR_new_mean_struct_imPix_only{i,j,A} = maxR_new_mean_imPix';
            %
            maxP_new_mean_struct_method_only{i,j,A} = maxP_new_mean';
            maxP_new_mean_struct_imBlur_only{i,j,A} = maxP_new_mean_imBlur';
            maxP_new_mean_struct_imPix_only{i,j,A} = maxP_new_mean_imPix';
            %
            thr_new_mean_struct_method_only{i,j,A} = thr_new_mean';
            thr_new_mean_struct_imBlur_only{i,j,A} = thr_new_mean_imBlur';
            thr_new_mean_struct_imPix_only{i,j,A}  = thr_new_mean_imPix';
            %
            % %
            %
            maxF_new_max_struct_method_only{i,j,A} = maxF_new_max';
            maxF_new_max_struct_imBlur_only{i,j,A} = maxF_new_max_imBlur';
            maxF_new_max_struct_imPix_only{i,j,A} = maxF_new_max_imPix';
            %
            maxR_new_max_struct_method_only{i,j,A} = maxR_new_max';
            maxR_new_max_struct_imBlur_only{i,j,A} = maxR_new_max_imBlur';
            maxR_new_max_struct_imPix_only{i,j,A} = maxR_new_max_imPix';
            %
            maxP_new_max_struct_method_only{i,j,A} = maxP_new_max';
            maxP_new_max_struct_imBlur_only{i,j,A} = maxP_new_max_imBlur';
            maxP_new_max_struct_imPix_only{i,j,A} = maxP_new_max_imPix';
            %
            thr_new_max_struct_method_only{i,j,A} = thr_new_max';
            thr_new_max_struct_imBlur_only{i,j,A} = thr_new_max_imBlur';
            thr_new_max_struct_imPix_only{i,j,A}  = thr_new_max_imPix';
            
            bestGT_new_max_struct_method_only{i,j,A} = maxF_bestGT';
            bestGT_new_max_struct_imBlur_only{i,j,A} = maxF_bestGT_imBlur';
            bestGT_new_max_struct_imPix_only{i,j,A} = maxF_bestGT_imPix';
            
            %
            % %
            %
            maxF_unionGT_struct_method_only{i,j,A} = maxF_unionGT';
            maxF_unionGT_struct_imBlur_only{i,j,A} = maxF_unionGT_imBlur';
            maxF_unionGT_struct_imPix_only{i,j,A} = maxF_unionGT_imPix';
            %
            maxR_unionGT_struct_method_only{i,j,A} = maxR_unionGT';
            maxR_unionGT_struct_imBlur_only{i,j,A} = maxR_unionGT_imBlur';
            maxR_unionGT_struct_imPix_only{i,j,A} = maxR_unionGT_imPix';
            %
            maxP_unionGT_struct_method_only{i,j,A} = maxP_unionGT';
            maxP_unionGT_struct_imBlur_only{i,j,A} = maxP_unionGT_imBlur';
            maxP_unionGT_struct_imPix_only{i,j,A} = maxP_unionGT_imPix';
            %
            thr_unionGT_struct_method_only{i,j,A} = thr_unionGT';
            thr_unionGT_struct_imBlur_only{i,j,A} = thr_unionGT_imBlur';
            thr_unionGT_struct_imPix_only{i,j,A}  = thr_unionGT_imPix';
            

            % Also save name of image patch so I can go back and visualize them next to the performance.
            img_ptch_name_struct{i,j,A} = imgPtchName;
            

        end % loop over KS parameter (for single method)
    end % loop over rM parameter (for single method)
    
    
  
    
    
    
    
    % Package together the different P/R/F computations for ImPix & ImBlur alone - to plot for comparison later.
    %
    % Or actually, I think it is best to do this outside of the Methods loop.
    %
    % And this sorting and plotting stuff below can happen outside the
    % methods loop and inside another loop over methods...
    
    
    
    
    
    
    
    % Sort average F-measure results (for different parameter settings)
    
    % (1). Original benchmark P/R/F computation - (not self-consistent)
    x1a = justMethod.maxF_meanAccImgs(:,:,A);
    x1b = justMethod.maxF_stdAccImgs(:,:,A);
    x1c = justMethod.maxF_skewAccImgs(:,:,A);
    [px1,qx1] = sort(x1a(:),'descend');
    %
    x2a = relImPix.maxF_meanAccImgs(:,:,A);
    x2b = relImPix.maxF_stdAccImgs(:,:,A);
    x2c = relImPix.maxF_skewAccImgs(:,:,A);
    [px2,qx2] = sort(x2a(:),'descend');
    %
    x3a = relImBlur.maxF_meanAccImgs(:,:,A);
    x3b = relImBlur.maxF_stdAccImgs(:,:,A);
    x3c = relImBlur.maxF_skewAccImgs(:,:,A);
    [px3,qx3] = sort(x3a(:),'descend');
    %
    % for Recall
    x1aR = justMethod.maxR_meanAccImgs(:,:,A);
    x1bR = justMethod.maxR_stdAccImgs(:,:,A);
    %
    x2aR = relImPix.maxR_meanAccImgs(:,:,A);
    x2bR = relImPix.maxR_stdAccImgs(:,:,A);
    %
    x3aR = relImBlur.maxR_meanAccImgs(:,:,A);
    x3bR = relImBlur.maxR_stdAccImgs(:,:,A);
    %
    % for Precision
    x1aP = justMethod.maxP_meanAccImgs(:,:,A);
    x1bP = justMethod.maxP_stdAccImgs(:,:,A);
    %
    x2aP = relImPix.maxP_meanAccImgs(:,:,A);
    x2bP = relImPix.maxP_stdAccImgs(:,:,A);
    %
    x3aP = relImBlur.maxP_meanAccImgs(:,:,A);
    x3bP = relImBlur.maxP_stdAccImgs(:,:,A);
    
    
    
    % (2). Computing average (across GT's) P/R/F - our implementation.
    y1a = justMethod.maxF_newMean_meanAccImgs(:,:,A);
    y1b = justMethod.maxF_newMean_stdAccImgs(:,:,A);
    y1c = justMethod.maxF_newMean_skewAccImgs(:,:,A);
    [py1,qy1] = sort(y1a(:),'descend');
    %
    y2a = relImPix.maxF_newMean_meanAccImgs(:,:,A);
    y2b = relImPix.maxF_newMean_stdAccImgs(:,:,A);
    y2c = relImPix.maxF_newMean_skewAccImgs(:,:,A);
    [py2,qy2] = sort(y2a(:),'descend');
    %
    y3a = relImBlur.maxF_newMean_meanAccImgs(:,:,A);
    y3b = relImBlur.maxF_newMean_stdAccImgs(:,:,A);
    y3c = relImBlur.maxF_newMean_skewAccImgs(:,:,A);
    [py3,qy3] = sort(y3a(:),'descend');
    %
    % for Recall
    y1aR = justMethod.maxR_newMean_meanAccImgs(:,:,A);
    y1bR = justMethod.maxR_newMean_stdAccImgs(:,:,A);
    %
    y2aR = relImPix.maxR_newMean_meanAccImgs(:,:,A);
    y2bR = relImPix.maxR_newMean_stdAccImgs(:,:,A);
    %
    y3aR = relImBlur.maxR_newMean_meanAccImgs(:,:,A);
    y3bR = relImBlur.maxR_newMean_stdAccImgs(:,:,A);
    %
    % for Precision
    y1aP = justMethod.maxP_newMean_meanAccImgs(:,:,A);
    y1bP = justMethod.maxP_newMean_stdAccImgs(:,:,A);
    %
    y2aP = relImPix.maxP_newMean_meanAccImgs(:,:,A);
    y2bP = relImPix.maxP_newMean_stdAccImgs(:,:,A);
    %
    y3aP = relImBlur.maxP_newMean_meanAccImgs(:,:,A);
    y3bP = relImBlur.maxP_newMean_stdAccImgs(:,:,A);
    
    
    
    % (3). Computing best P/R/F by using only 1 (best) GT for each image patch.
    z1a = justMethod.maxF_newMax_meanAccImgs(:,:,A);
    z1b = justMethod.maxF_newMax_stdAccImgs(:,:,A);
    z1c = justMethod.maxF_newMax_skewAccImgs(:,:,A);
    [pz1,qz1] = sort(z1a(:),'descend');
    %
    z2a = relImPix.maxF_newMax_meanAccImgs(:,:,A);
    z2b = relImPix.maxF_newMax_stdAccImgs(:,:,A);
    z2c = relImPix.maxF_newMax_skewAccImgs(:,:,A);
    [pz2,qz2] = sort(z2a(:),'descend');
    %
    z3a = relImBlur.maxF_newMax_meanAccImgs(:,:,A);
    z3b = relImBlur.maxF_newMax_stdAccImgs(:,:,A);
    z3c = relImBlur.maxF_newMax_skewAccImgs(:,:,A);
    [pz3,qz3] = sort(z3a(:),'descend');
    %
    % for Recall
    z1aR = justMethod.maxR_newMax_meanAccImgs(:,:,A);
    z1bR = justMethod.maxR_newMax_stdAccImgs(:,:,A);
    %
    z2aR = relImPix.maxR_newMax_meanAccImgs(:,:,A);
    z2bR = relImPix.maxR_newMax_stdAccImgs(:,:,A);
    %
    z3aR = relImBlur.maxR_newMax_meanAccImgs(:,:,A);
    z3bR = relImBlur.maxR_newMax_stdAccImgs(:,:,A);
    %
    % for Precision
    z1aP = justMethod.maxP_newMax_meanAccImgs(:,:,A);
    z1bP = justMethod.maxP_newMax_stdAccImgs(:,:,A);
    %
    z2aP = relImPix.maxP_newMax_meanAccImgs(:,:,A);
    z2bP = relImPix.maxP_newMax_stdAccImgs(:,:,A);
    %
    z3aP = relImBlur.maxP_newMax_meanAccImgs(:,:,A);
    z3bP = relImBlur.maxP_newMax_stdAccImgs(:,:,A);
    
    
    % (4). Computing P/R/F by first unioning all GT's into one single GT.
    a1a = justMethod.maxF_unionGT_meanAccImgs(:,:,A);
    a1b = justMethod.maxF_unionGT_stdAccImgs(:,:,A);
    a1c = justMethod.maxF_unionGT_skewAccImgs(:,:,A);
    [pa1,qa1] = sort(a1a(:),'descend');
    %
    a2a = relImPix.maxF_unionGT_meanAccImgs(:,:,A);
    a2b = relImPix.maxF_unionGT_stdAccImgs(:,:,A);
    a2c = relImPix.maxF_unionGT_skewAccImgs(:,:,A);
    [pa2,qa2] = sort(a2a(:),'descend');
    %
    a3a = relImBlur.maxF_unionGT_meanAccImgs(:,:,A);
    a3b = relImBlur.maxF_unionGT_stdAccImgs(:,:,A);
    a3c = relImBlur.maxF_unionGT_skewAccImgs(:,:,A);
    [pa3,qa3] = sort(a3a(:),'descend');
    %
    % for Recall
    a1aR = justMethod.maxR_unionGT_meanAccImgs(:,:,A);
    a1bR = justMethod.maxR_unionGT_stdAccImgs(:,:,A);
    %
    a2aR = relImPix.maxR_unionGT_meanAccImgs(:,:,A);
    a2bR = relImPix.maxR_unionGT_stdAccImgs(:,:,A);
    %
    a3aR = relImBlur.maxR_unionGT_meanAccImgs(:,:,A);
    a3bR = relImBlur.maxR_unionGT_stdAccImgs(:,:,A);
    %
    % for Precision
    a1aP = justMethod.maxP_unionGT_meanAccImgs(:,:,A);
    a1bP = justMethod.maxP_unionGT_stdAccImgs(:,:,A);
    %
    a2aP = relImPix.maxP_unionGT_meanAccImgs(:,:,A);
    a2bP = relImPix.maxP_unionGT_stdAccImgs(:,:,A);
    %
    a3aP = relImBlur.maxP_unionGT_meanAccImgs(:,:,A);
    a3bP = relImBlur.maxP_unionGT_stdAccImgs(:,:,A);
    

    % Flags for plots below. Set to 1 to plot.
    plot_F_violins_all_params_each_method = 1;     
    plot_FRP_errorbars_all_params_each_method = 1;
    
        
    % Plot mean & std F-measure. Benchmark Implementation: Just Method.
    if(plot_F_violins_all_params_each_method)
        H1 = figure; hold on
        maxF_temp = maxF_old_struct_method_only(:,:,A);
        for xx = 1:numel(qx1)
            violin(maxF_temp(qx1(xx)),'x',[xx xx+1])
        end
        errorbar(x1a(qx1),x1b(qx1),'LineStyle','none','Marker','o','LineWidth',2);
        for xx = 1:numel(qx1)
            text(xx+0.1,x1a(qx1(xx))-2*x1b(qx1(xx)),['s=',num2str(x1c(qx1(xx)),'%+5.2f')]) % indicate the skew
        end
        axis([0.5 numel(qx1)+2.5 min(x1a(qx1)-3*x1b(qx1)) max(x1a(qx1)+3*x1b(qx1))])
        %
        set(gca,'XTick',1:numel(qx1)+2,'XTickLabel',[param_matrix(qx1);'ImPix';'ImBlur'],'FontSize',12,'FontWeight','Bold')
        title([method{A},' ',blur_tit,' : benchmark F :  just Method'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        %
        saveGoodImg(H1,[dirPre,'../Documentation/Cosyne_2016/F_benchmark_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_justMethod.jpg'],sizeGoodIm)
        close(H1)
    end
    
    
    % Plot Precision & Recall Separately (instead of combining them together as F-measure)
    % because maybe network methods effect one more than the other or something.
    if(plot_FRP_errorbars_all_params_each_method)
        H1_RP = figure; hold on
        errorbar([1:numel(qx1)]-0.3, x1aR(qx1), x1bR(qx1),'LineStyle','none','Marker','o','LineWidth',2,'Color','Red');  % Recall
        errorbar([1:numel(qx1)]+0.3, x1aP(qx1), x1bP(qx1),'LineStyle','none','Marker','o','LineWidth',2,'Color','Blue'); % Precision
        errorbar([1:numel(qx1)],     x1a(qx1),  x1b(qx1),'LineStyle','none','Marker','o','LineWidth',2,'Color','Green'); % F-measure
        %
        axis tight %([0.5 numel(qx1)+2.5 min([x1aR(qx1)-x1bR(qx1); x1aP(qx1)-x1bP(qx1); x1a(qx1)-x1b(qx1)]) max([x1aR(qx1)+x1bR(qx1); x1aP(qx1)+x1bP(qx1); x1a(qx1)+x1b(qx1)]) ])
        plot([0.5 numel(qx1)+2.5],[0 0],'k--')
        set(gca,'XTick',1:numel(qx1)+2,'XTickLabel',[param_matrix(qx1);'ImPix';'ImBlur'],'FontSize',12,'FontWeight','Bold')
        title([method{A},' ',blur_tit,' : benchmark F :  just Method'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        legend('Recall','Precision','F-measure')
        %
        saveGoodImg(H1_RP,[dirPre,'../Documentation/Cosyne_2016/FRP_benchmark_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_justMethod.jpg'],sizeGoodIm)
        close(H1_RP)
    end
    

    
    
    
    % Plot mean & std F-measure. Benchmark Implementation: Improvement over image pixels.
    if(plot_F_violins_all_params_each_method)
        H2 = figure; hold on
        % subtract all elements of vectors within these 2 cells (have same number of elements) to get relImPix.
        maxF_temp = cellfun(@minus,maxF_old_struct_method_only(:,:,A),maxF_old_struct_imPix_only(:,:,A),'UniformOutput',false);
        for xx = 1:numel(qx2)
            violin(maxF_temp(qx2(xx)),'x',[xx xx+1])
        end
        errorbar(x2a(qx2),x2b(qx2),'LineStyle','none','Marker','o','LineWidth',2);
        for xx = 1:numel(qx2)
            text(xx+0.1,x2a(qx2(xx))-2*x2b(qx2(xx)),['s=',num2str(x2c(qx2(xx)),'%+5.2f')]) % indicate the skew
        end
        axis([0.5 numel(qx2)+0.5 min(x2a(qx2)-3*x2b(qx2)) max(x2a(qx2)+3*x2b(qx2))])
        plot([1 numel(qx2)],[0 0],'k--')
        %
        set(gca,'XTick',1:numel(qx2),'XTickLabel',[param_matrix(qx2);],'FontSize',12,'FontWeight','Bold')
        title([method{A},' ',blur_tit,' : benchmark F :  rel Im Pix'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        %
        saveGoodImg(H2,[dirPre,'../Documentation/Cosyne_2016/F_benchmark_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_relImPix.jpg'],sizeGoodIm)
        close(H2)
    end
    
    
    
    
    % Plot Precision & Recall Separately (instead of combining them together as F-measure)
    % because maybe network methods effect one more than the other or something.
    if(plot_FRP_errorbars_all_params_each_method)
        H2_RP = figure; hold on
        errorbar([1:numel(qx2)]-0.3, x2aR(qx2), x2bR(qx2),'LineStyle','none','Marker','o','LineWidth',2,'Color','Blue');  % Recall
        errorbar([1:numel(qx2)]+0.3, x2aP(qx2), x2bP(qx2),'LineStyle','none','Marker','o','LineWidth',2,'Color','Red'); % Precision
        errorbar([1:numel(qx2)],     x2a(qx2),  x2b(qx2),'LineStyle','none','Marker','o','LineWidth',2,'Color','Green'); % F-measure
        %
        axis tight %([0.5 numel(qx2)+2.5 min([x2aR(qx2)-x2bR(qx2); x2aP(qx2)-x2bP(qx2); x2a(qx2)-x2b(qx2)]) max([x2aR(qx2)+x2bR(qx2); x2aP(qx2)+x2bP(qx2); x2a(qx2)+x2b(qx2)]) ])
        plot([0.5 numel(qx2)+2.5],[0 0],'k--')
        set(gca,'XTick',1:numel(qx2)+2,'XTickLabel',[param_matrix(qx2);'ImPix';'ImBlur'],'FontSize',12,'FontWeight','Bold')
        title([method{A},' ',blur_tit,' : benchmark F :  relImPix'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        legend('Recall','Precision','F-measure')
        %
        saveGoodImg(H2_RP,[dirPre,'../Documentation/Cosyne_2016/FRP_benchmark_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_relImPix.jpg'],sizeGoodIm)
        close(H2_RP)
    end
    
    

    
    % Plot mean & std F-measure. Benchmark Implementation: Improvement over image blur.
    if(plot_F_violins_all_params_each_method)
        H3 = figure; hold on
        % subtract all elements of vectors within these 2 cells (have same number of elements) to get relImBlur.
        maxF_temp = cellfun(@minus,maxF_old_struct_method_only(:,:,A),maxF_old_struct_imBlur_only(:,:,A),'UniformOutput',false);
        for xx = 1:numel(qx3)
            violin(maxF_temp(qx3(xx)),'x',[xx xx+1])
        end
        errorbar(x3a(qx3),x3b(qx3),'LineStyle','none','Marker','o','LineWidth',2);
        for xx = 1:numel(qx3)
            text(xx+0.1,x3a(qx3(xx))-2*x3b(qx3(xx)),['s=',num2str(x3c(qx3(xx)),'%+5.2f')]) % indicate the skew
        end
        axis([0.5 numel(qx3)+0.5 min(x3a(qx3)-3*x3b(qx3)) max(x3a(qx3)+3*x3b(qx3))])
        plot([1 numel(qx3)],[0 0],'k--')
        %
        set(gca,'XTick',1:numel(qx3),'XTickLabel',[param_matrix(qx3);],'FontSize',12,'FontWeight','Bold')
        title([method{A},' ',blur_tit,' : benchmark F :  rel Im Blur'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        saveGoodImg(H3,[dirPre,'../Documentation/Cosyne_2016/F_benchmark_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_relImBlur.jpg'],sizeGoodIm)
        close(H3)
    end
    
    
    
    
    
    % Plot Precision & Recall Separately (instead of combining them together as F-measure)
    % because maybe network methods effect one more than the other or something.
    if(plot_FRP_errorbars_all_params_each_method)
        H3_RP = figure; hold on
        errorbar([1:numel(qx3)]-0.3, x3aR(qx3), x3bR(qx3),'LineStyle','none','Marker','o','LineWidth',2,'Color','Blue');  % Recall
        errorbar([1:numel(qx3)]+0.3, x3aP(qx3), x3bP(qx3),'LineStyle','none','Marker','o','LineWidth',2,'Color','Red'); % Precision
        errorbar([1:numel(qx3)],     x3a(qx3),  x3b(qx3),'LineStyle','none','Marker','o','LineWidth',2,'Color','Green'); % F-measure
        %
        axis tight %([0.5 numel(qx3)+2.5 min([x3aR(qx3)-x3bR(qx3); x3aP(qx3)-x3bP(qx3); x3a(qx3)-x3b(qx3)]) max([x3aR(qx3)+x3bR(qx3); x3aP(qx3)+x3bP(qx3); x3a(qx3)+x3b(qx3)]) ])
        plot([0.5 numel(qx3)+2.5],[0 0],'k--')
        set(gca,'XTick',1:numel(qx3)+2,'XTickLabel',[param_matrix(qx3);'ImPix';'ImBlur'],'FontSize',12,'FontWeight','Bold')
        title([method{A},' ',blur_tit,' : benchmark F :  relImBlur'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        legend('Recall','Precision','F-measure')
        %
        saveGoodImg(H3_RP,[dirPre,'../Documentation/Cosyne_2016/FRP_benchmark_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_relImBlur.jpg'],sizeGoodIm)
        close(H3_RP)
    end
    
    
    
    
    % Plot mean & std F-measure. Mean (across GT) F: Just Method.
    if(plot_F_violins_all_params_each_method)
        H4 = figure; hold on
        maxF_temp = maxF_new_mean_struct_method_only(:,:,A);
        for xx = 1:numel(qy1)
            violin(maxF_temp(qy1(xx)),'x',[xx xx+1])
        end
        errorbar(y1a(qy1),y1b(qy1),'LineStyle','none','Marker','o','LineWidth',2);
        for xx = 1:numel(qy1)
            text(xx+0.1,y1a(qy1(xx))-2*y1b(qy1(xx)),['s=',num2str(y1c(qy1(xx)),'%+5.2f')]) % indicate the skew
        end
        axis([0.5 numel(qy1)+2.5 min(y1a(qy1)-3*y1b(qy1)) max(y1a(qy1)+3*y1b(qy1))])
        %
        set(gca,'XTick',1:numel(qy1)+2,'XTickLabel',[param_matrix(qy1);'ImPix';'ImBlur'],'FontSize',12,'FontWeight','Bold')
        title([method{A},' ',blur_tit,' : mean (across GT) F :  just Method'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        saveGoodImg(H4,[dirPre,'../Documentation/Cosyne_2016/F_meanAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_justMethod.jpg'],sizeGoodIm)
        close(H4)
    end

    
    
    % Plot Precision & Recall Separately (instead of combining them together as F-measure)
    % because maybe network methods effect one more than the other or something.
    if(plot_FRP_errorbars_all_params_each_method)
        H4_RP = figure; hold on
        errorbar([1:numel(qy1)]-0.3, y1aR(qy1), y1bR(qy1),'LineStyle','none','Marker','o','LineWidth',2,'Color','Blue');  % Recall
        errorbar([1:numel(qy1)]+0.3, y1aP(qy1), y1bP(qy1),'LineStyle','none','Marker','o','LineWidth',2,'Color','Red'); % Precision
        errorbar([1:numel(qy1)],     y1a(qy1),  y1b(qy1),'LineStyle','none','Marker','o','LineWidth',2,'Color','Green'); % F-measure
        %
        axis tight % ([0.5 numel(qy1)+2.5 0 1])
        plot([0.5 numel(qy1)+2.5],[0 0],'k--')
        set(gca,'XTick',1:numel(qy1)+2,'XTickLabel',[param_matrix(qy1);'ImPix';'ImBlur'],'FontSize',12,'FontWeight','Bold')
        title([method{A},' ',blur_tit,' : mean (across GT) F :  just Method'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        legend('Recall','Precision','F-measure')
        %
        saveGoodImg(H4_RP,[dirPre,'../Documentation/Cosyne_2016/FRP_meanAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_justMethod.jpg'],sizeGoodIm)
        close(H4_RP)
    end
    
    
    
    % Plot mean & std F-measure. Mean (across GT) F: Improvement over image pixels.
    if(plot_F_violins_all_params_each_method)
        H5 = figure; hold on
        % subtract all elements of vectors within these 2 cells (have same number of elements) to get relImPix.
        maxF_temp = cellfun(@minus,maxF_new_mean_struct_method_only(:,:,A),maxF_new_mean_struct_imPix_only(:,:,A),'UniformOutput',false);
        
        for xx = 1:numel(qy2)
            violin(maxF_temp(qy2(xx)),'x',[xx xx+1])
        end
        errorbar(y2a(qy2),y2b(qy2),'LineStyle','none','Marker','o','LineWidth',2);
        for xx = 1:numel(qy2)
            text(xx+0.1,y2a(qy2(xx))-2*y2b(qy2(xx)),['s=',num2str(y2c(qy2(xx)),'%+5.2f')]) % indicate the skew
        end
        axis([0.5 numel(qy2)+0.5 min(y2a(qy2)-3*y2b(qy2)) max(y2a(qy2)+3*y2b(qy2))])
        plot([1 numel(qy2)],[0 0],'k--')
        %
        set(gca,'XTick',1:numel(qy2),'XTickLabel',[param_matrix(qy2);],'FontSize',12,'FontWeight','Bold')
        title([method{A},' ',blur_tit,' : mean (across GT) F :  rel Im Pix'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        saveGoodImg(H5,[dirPre,'../Documentation/Cosyne_2016/F_meanAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_relImPix.jpg'],sizeGoodIm)
        close(H5)
    end
    
    
    
    
    % Plot Precision & Recall Separately (instead of combining them together as F-measure)
    % because maybe network methods effect one more than the other or something.
    if(plot_FRP_errorbars_all_params_each_method)
        H5_RP = figure; hold on
        errorbar([1:numel(qy2)]-0.3, y2aR(qy2), y2bR(qy2),'LineStyle','none','Marker','o','LineWidth',2,'Color','Blue');  % Recall
        errorbar([1:numel(qy2)]+0.3, y2aP(qy2), y2bP(qy2),'LineStyle','none','Marker','o','LineWidth',2,'Color','Red'); % Precision
        errorbar([1:numel(qy2)],     y2a(qy2),  y2b(qy2),'LineStyle','none','Marker','o','LineWidth',2,'Color','Green'); % F-measure
        %
        axis tight %([0.5 numel(qy2)+2.5 -0.6 0.6])
        plot([0.5 numel(qy2)+2.5],[0 0],'k--')
        set(gca,'XTick',1:numel(qy2)+2,'XTickLabel',[param_matrix(qy2);'ImPix';'ImBlur'],'FontSize',12,'FontWeight','Bold')
        title([method{A},' ',blur_tit,' : mean (across GT) F :  relImPix'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        legend('Recall','Precision','F-measure')
        %
        saveGoodImg(H5_RP,[dirPre,'../Documentation/Cosyne_2016/FRP_meanAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_relImPix.jpg'],sizeGoodIm)
        close(H5_RP)
    end
    

    
    % Plot mean & std F-measure. Mean (across GT) F: Improvement over image blur.
    if(plot_F_violins_all_params_each_method)
        H6 = figure; hold on
        % subtract all elements of vectors within these 2 cells (have same number of elements) to get relImBlur.
        maxF_temp = cellfun(@minus,maxF_new_mean_struct_method_only(:,:,A),maxF_new_mean_struct_imBlur_only(:,:,A),'UniformOutput',false);

        for xx = 1:numel(qy3)
            violin(maxF_temp(qy3(xx)),'x',[xx xx+1])
        end
        errorbar(y3a(qy3),x3b(qy3),'LineStyle','none','Marker','o','LineWidth',2);
        for xx = 1:numel(qy3)
            text(xx+0.1,y3a(qy3(xx))-2*y3b(qy3(xx)),['s=',num2str(y3c(qy3(xx)),'%+5.2f')]) % indicate the skew
        end
        axis([0.5 numel(qy3)+0.5 min(y3a(qy3)-3*y3b(qy3)) max(y3a(qy3)+3*y3b(qy3))])
        plot([1 numel(qy3)],[0 0],'k--')
        %
        set(gca,'XTick',1:numel(qy3),'XTickLabel',[param_matrix(qy3);],'FontSize',12,'FontWeight','Bold')
        title([method{A},' ',blur_tit,' : mean (across GT) F :  rel Im Blur'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        saveGoodImg(H6,[dirPre,'../Documentation/Cosyne_2016/F_meanAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_relImBlur.jpg'],sizeGoodIm)
        close(H6)
    end
    
    
    
    % Plot Precision & Recall Separately (instead of combining them together as F-measure)
    % because maybe network methods effect one more than the other or something.
    if(plot_FRP_errorbars_all_params_each_method)
        H6_RP = figure; hold on
        errorbar([1:numel(qy3)]-0.3, y3aR(qy3), y3bR(qy3),'LineStyle','none','Marker','o','LineWidth',2,'Color','Blue');  % Recall
        errorbar([1:numel(qy3)]+0.3, y3aP(qy3), y3bP(qy3),'LineStyle','none','Marker','o','LineWidth',2,'Color','Red'); % Precision
        errorbar([1:numel(qy3)],     y3a(qy3),  y3b(qy3),'LineStyle','none','Marker','o','LineWidth',2,'Color','Green'); % F-measure
        %
        axis tight %([0.5 numel(qy3)+2.5 -0.6 0.6])
        plot([0.5 numel(qy3)+2.5],[0 0],'k--')
        set(gca,'XTick',1:numel(qy3)+2,'XTickLabel',[param_matrix(qy3);'ImPix';'ImBlur'],'FontSize',12,'FontWeight','Bold')
        title([method{A},' ',blur_tit,' : mean (across GT) F :  relImBlur'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        legend('Recall','Precision','F-measure')
        %
        saveGoodImg(H6_RP,[dirPre,'../Documentation/Cosyne_2016/FRP_meanAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_relImBlur.jpg'],sizeGoodIm)
        close(H6_RP)
    end
    
    
    
    
    
    % Plot mean & std F-measure. Max (across GT) F: Just Method.
    if(plot_F_violins_all_params_each_method)
        H7 = figure; hold on
        maxF_temp = maxF_new_max_struct_method_only(:,:,A);
        for xx = 1:numel(qz1)
            violin(maxF_temp(qz1(xx)),'x',[xx xx+1])
        end
        errorbar(z1a(qz1),z1b(qz1),'LineStyle','none','Marker','o','LineWidth',2);
        for xx = 1:numel(qz1)
            text(xx+0.1,z1a(qz1(xx))-2*z1b(qz1(xx)),['s=',num2str(z1c(qz1(xx)),'%+5.2f')]) % indicate the skew
        end
        axis([0.5 numel(qz1)+2.5 min(z1a(qz1)-3*z1b(qz1)) max(z1a(qz1)+3*z1b(qz1))])
        %
        set(gca,'XTick',1:numel(qz1)+2,'XTickLabel',[param_matrix(qz1);'ImPix';'ImBlur'],'FontSize',12,'FontWeight','Bold')
        title([method{A},' ',blur_tit,' : max (across GT) F :  just Method'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        %
        saveGoodImg(H7,[dirPre,'../Documentation/Cosyne_2016/F_maxAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_justMethod.jpg'],sizeGoodIm)
        close(H7)
    end
    
    
    
    % Plot Precision & Recall Separately (instead of combining them together as F-measure)
    % because maybe network methods effect one more than the other or something.
    if(plot_FRP_errorbars_all_params_each_method)
        H7_RP = figure; hold on
        errorbar([1:numel(qz1)]-0.3, z1aR(qz1), z1bR(qz1),'LineStyle','none','Marker','o','LineWidth',2,'Color','Blue');  % Recall
        errorbar([1:numel(qz1)]+0.3, z1aP(qz1), z1bP(qz1),'LineStyle','none','Marker','o','LineWidth',2,'Color','Red'); % Precision
        errorbar([1:numel(qz1)],     z1a(qz1),  z1b(qz1),'LineStyle','none','Marker','o','LineWidth',2,'Color','Green'); % F-measure
        %
        axis tight %([0.5 numel(qz1)+2.5 -0.6 0.6])
        plot([0.5 numel(qz1)+2.5],[0 0],'k--')
        set(gca,'XTick',1:numel(qz1)+2,'XTickLabel',[param_matrix(qz1);'ImPix';'ImBlur'],'FontSize',12,'FontWeight','Bold')
        title([method{A},' ',blur_tit,' : max (across GT) F : just Method'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        legend('Recall','Precision','F-measure')
        %
        saveGoodImg(H7_RP,[dirPre,'../Documentation/Cosyne_2016/FRP_maxAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_justMethod.jpg'],sizeGoodIm)
        close(H7_RP)
    end
    

    
    % Plot mean & std F-measure. Max (across GT) F: Improvement over image pixels.
    if(plot_F_violins_all_params_each_method)
        H8 = figure; hold on
        % subtract all elements of vectors within these 2 cells (have same number of elements) to get relImPix.
        maxF_temp = cellfun(@minus,maxF_new_max_struct_method_only(:,:,A),maxF_new_max_struct_imPix_only(:,:,A),'UniformOutput',false);
        
        for xx = 1:numel(qz2)
            violin(maxF_temp(qz2(xx)),'x',[xx xx+1])
        end
        errorbar(z2a(qz2),z2b(qz2),'LineStyle','none','Marker','o','LineWidth',2);
        for xx = 1:numel(qz2)
            text(xx+0.1,z2a(qz2(xx))-2*z2b(qz2(xx)),['s=',num2str(z2c(qz2(xx)),'%+5.2f')]) % indicate the skew
        end
        axis([0.5 numel(qz2)+0.5 min(z2a(qz2)-3*z2b(qz2)) max(z2a(qz2)+3*z2b(qz2))])
        plot([1 numel(qz2)],[0 0],'k--')
        %
        set(gca,'XTick',1:numel(qz2),'XTickLabel',[param_matrix(qz2);],'FontSize',12,'FontWeight','Bold')
        title([method{A},' ',blur_tit,' : max (across GT) F :  rel Im Pix'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        %
        saveGoodImg(H8,[dirPre,'../Documentation/Cosyne_2016/F_maxAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_relImPix.jpg'],sizeGoodIm)
        close(H8)
    end
    
    
    
    % Plot Precision & Recall Separately (instead of combining them together as F-measure)
    % because maybe network methods effect one more than the other or something.
    if(plot_FRP_errorbars_all_params_each_method)
        H8_RP = figure; hold on
        errorbar([1:numel(qz2)]-0.3, z2aR(qz2), z2bR(qz2),'LineStyle','none','Marker','o','LineWidth',2,'Color','Blue');  % Recall
        errorbar([1:numel(qz2)]+0.3, z2aP(qz2), z2bP(qz2),'LineStyle','none','Marker','o','LineWidth',2,'Color','Red'); % Precision
        errorbar([1:numel(qz2)],     z2a(qz2),  z2b(qz2),'LineStyle','none','Marker','o','LineWidth',2,'Color','Green'); % F-measure
        %
        axis tight %([0.5 numel(qz2)+2.5 -0.6 0.6])
        plot([0.5 numel(qz2)+2.5],[0 0],'k--')
        set(gca,'XTick',1:numel(qz2)+2,'XTickLabel',[param_matrix(qz2);'ImPix';'ImBlur'],'FontSize',12,'FontWeight','Bold')
        title([method{A},' ',blur_tit,' : max (across GT) F : relImPix'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        legend('Recall','Precision','F-measure')
        %
        saveGoodImg(H8_RP,[dirPre,'../Documentation/Cosyne_2016/FRP_maxAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_relImPix.jpg'],sizeGoodIm)
        close(H8_RP)
    end
    
    

    
    % Plot mean & std F-measure. Max (across GT) F: Improvement over image blur.
    if(plot_F_violins_all_params_each_method)
        H9 = figure; hold on
        % subtract all elements of vectors within these 2 cells (have same number of elements) to get relImBlur.
        maxF_temp = cellfun(@minus,maxF_new_max_struct_method_only(:,:,A),maxF_new_max_struct_imBlur_only(:,:,A),'UniformOutput',false);
        for xx = 1:numel(qz3)
            violin(maxF_temp(qz3(xx)),'x',[xx xx+1])
        end
        errorbar(z3a(qz3),z3b(qz3),'LineStyle','none','Marker','o','LineWidth',2);
        for xx = 1:numel(qz3)
            text(xx+0.1,z3a(qz3(xx))-2*z3b(qz3(xx)),['s=',num2str(z3c(qz3(xx)),'%+5.2f')]) % indicate the skew
        end
        axis([0.5 numel(qz3)+0.5 min(z3a(qz3)-3*z3b(qz3)) max(z3a(qz3)+3*z3b(qz3))])
        plot([1 numel(qz3)],[0 0],'k--')
        %
        set(gca,'XTick',1:numel(qz3),'XTickLabel',[param_matrix(qz3);],'FontSize',12,'FontWeight','Bold')
        title([method{A},' ',blur_tit,' : max (across GT) F :  rel Im Blur'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        %
        saveGoodImg(H9,[dirPre,'../Documentation/Cosyne_2016/F_maxAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_relImBlur.jpg'],sizeGoodIm)
        close(H9)
    end
    
    
    
    % Plot Precision & Recall Separately (instead of combining them together as F-measure)
    % because maybe network methods effect one more than the other or something.
    if(plot_FRP_errorbars_all_params_each_method)
        H9_RP = figure; hold on
        errorbar([1:numel(qz3)]-0.3, z3aR(qz3), z3bR(qz3),'LineStyle','none','Marker','o','LineWidth',2,'Color','Blue');  % Recall
        errorbar([1:numel(qz3)]+0.3, z3aP(qz3), z3bP(qz3),'LineStyle','none','Marker','o','LineWidth',2,'Color','Red'); % Precision
        errorbar([1:numel(qz3)],     z3a(qz3),  z3b(qz3),'LineStyle','none','Marker','o','LineWidth',2,'Color','Green'); % F-measure
        %
        axis tight %([0.5 numel(qz3)+2.5 -0.6 0.6])
        plot([0.5 numel(qz3)+2.5],[0 0],'k--')
        set(gca,'XTick',1:numel(qz3)+2,'XTickLabel',[param_matrix(qz3);'ImPix';'ImBlur'],'FontSize',12,'FontWeight','Bold')
        title([method{A},' ',blur_tit,' : max (across GT) F : relImBlur'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        legend('Recall','Precision','F-measure')
        %
        saveGoodImg(H9_RP,[dirPre,'../Documentation/Cosyne_2016/FRP_maxAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_relImBlur.jpg'],sizeGoodIm)
        close(H9_RP)
    end
    
    
    
    % Plot mean & std F-measure. Union GT: Just Method.
    if(plot_F_violins_all_params_each_method)
        H10 = figure; hold on
        maxF_temp = maxF_unionGT_struct_method_only(:,:,A);
        for xx = 1:numel(qa1)
            violin(maxF_temp(qa1(xx)),'x',[xx xx+1])
        end
        errorbar(a1a(qa1),a1b(qa1),'LineStyle','none','Marker','o','LineWidth',2);
        for xx = 1:numel(qa1)
            text(xx+0.1,a1a(qa1(xx))-2*a1b(qa1(xx)),['s=',num2str(a1c(qa1(xx)),'%+5.2f')]) % indicate the skew
        end
        axis([0.5 numel(qa1)+2.5 min(a1a(qa1)-3*a1b(qa1)) max(a1a(qa1)+3*a1b(qa1))])
        %
        set(gca,'XTick',1:numel(qa1)+2,'XTickLabel',[param_matrix(qa1);'ImPix';'ImBlur'],'FontSize',12,'FontWeight','Bold')
        title([method{A},' ',blur_tit,' : Union GT :  just Method'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        %
        saveGoodImg(H10,[dirPre,'../Documentation/Cosyne_2016/F_unionGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_justMethod.jpg'],sizeGoodIm)
        close(H10)
    end

    
    
    
    % Plot Precision & Recall Separately (instead of combining them together as F-measure)
    % because maybe network methods effect one more than the other or something.
    if(plot_FRP_errorbars_all_params_each_method)
        H10_RP = figure; hold on
        errorbar([1:numel(qa1)]-0.3, a1aR(qa1), a1bR(qa1),'LineStyle','none','Marker','o','LineWidth',2,'Color','Blue');  % Recall
        errorbar([1:numel(qa1)]+0.3, a1aP(qa1), a1bP(qa1),'LineStyle','none','Marker','o','LineWidth',2,'Color','Red'); % Precision
        errorbar([1:numel(qa1)],     a1a(qa1),  a1b(qa1),'LineStyle','none','Marker','o','LineWidth',2,'Color','Green'); % F-measure
        %
        axis tight %([0.5 numel(qa1)+2.5 -0.6 0.6])
        plot([0.5 numel(qa1)+2.5],[0 0],'k--')
        set(gca,'XTick',1:numel(qa1)+2,'XTickLabel',[param_matrix(qa1);'ImPix';'ImBlur'],'FontSize',12,'FontWeight','Bold')
        title([method{A},' ',blur_tit,' : Union GT : justMethod'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        legend('Recall','Precision','F-measure')
        %
        saveGoodImg(H10_RP,[dirPre,'../Documentation/Cosyne_2016/FRP_unionGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_justMethod.jpg'],sizeGoodIm)
        close(H10_RP)
    end
    
    
    
    % Plot mean & std F-measure. Union GT: Improvement over image pixels.
    if(plot_F_violins_all_params_each_method)
        H11 = figure; hold on
        % subtract all elements of vectors within these 2 cells (have same number of elements) to get relImPix.
        maxF_temp = cellfun(@minus,maxF_unionGT_struct_method_only(:,:,A),maxF_unionGT_struct_imPix_only(:,:,A),'UniformOutput',false);
        
        for xx = 1:numel(qa2)
            violin(maxF_temp(qa2(xx)),'x',[xx xx+1])
        end
        errorbar(a2a(qa2),a2b(qa2),'LineStyle','none','Marker','o','LineWidth',2);
        for xx = 1:numel(qa2)
            text(xx+0.1,a2a(qa2(xx))-2*a2b(qa2(xx)),['s=',num2str(a2c(qa2(xx)),'%+5.2f')]) % indicate the skew
        end
        axis([0.5 numel(qa2)+0.5 min(a2a(qa2)-3*a2b(qa2)) max(a2a(qa2)+3*a2b(qa2))])
        plot([1 numel(qa2)],[0 0],'k--')
        %
        set(gca,'XTick',1:numel(qa2),'XTickLabel',[param_matrix(qa2);],'FontSize',12,'FontWeight','Bold')
        title([method{A},' ',blur_tit,' : Union GT :  rel Im Pix'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        %
        saveGoodImg(H11,[dirPre,'../Documentation/Cosyne_2016/F_unionGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_relImPix.jpg'],sizeGoodIm)
        close(H11)
    end
    
    
    
    
    % Plot Precision & Recall Separately (instead of combining them together as F-measure)
    % because maybe network methods effect one more than the other or something.
    if(plot_FRP_errorbars_all_params_each_method)
        H11_RP = figure; hold on
        errorbar([1:numel(qa2)]-0.3, a2aR(qa2), a2bR(qa2),'LineStyle','none','Marker','o','LineWidth',2,'Color','Blue');  % Recall
        errorbar([1:numel(qa2)]+0.3, a2aP(qa2), a2bP(qa2),'LineStyle','none','Marker','o','LineWidth',2,'Color','Red'); % Precision
        errorbar([1:numel(qa2)],     a2a(qa2),  a2b(qa2),'LineStyle','none','Marker','o','LineWidth',2,'Color','Green'); % F-measure
        %
        axis tight %([0.5 numel(qa2)+2.5 -0.6 0.6])
        plot([0.5 numel(qa2)+2.5],[0 0],'k--')
        set(gca,'XTick',1:numel(qa2)+2,'XTickLabel',[param_matrix(qa2);'ImPix';'ImBlur'],'FontSize',12,'FontWeight','Bold')
        title([method{A},' ',blur_tit,' : Union GT : relImPix'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        legend('Recall','Precision','F-measure')
        %
        saveGoodImg(H11_RP,[dirPre,'../Documentation/Cosyne_2016/FRP_unionGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_relImPix.jpg'],sizeGoodIm)
        close(H11_RP)
    end
    
    

    
    % Plot mean & std F-measure. Union GT: Improvement over image blur.
    if(plot_F_violins_all_params_each_method)
        H12 = figure; hold on
        % subtract all elements of vectors within these 2 cells (have same number of elements) to get relImBlur.
        maxF_temp = cellfun(@minus,maxF_unionGT_struct_method_only(:,:,A),maxF_unionGT_struct_imBlur_only(:,:,A),'UniformOutput',false);
        for xx = 1:numel(qa3)
            violin(maxF_temp(qa3(xx)),'x',[xx xx+1])
        end
        errorbar(a3a(qa3),a3b(qa3),'LineStyle','none','Marker','o','LineWidth',2);
        for xx = 1:numel(qa3)
            text(xx+0.1,a3a(qa3(xx))-2*a3b(qa3(xx)),['s=',num2str(a3c(qa3(xx)),'%+5.2f')]) % indicate the skew
        end
        axis([0.5 numel(qa3)+0.5 min(a3a(qa3)-3*a3b(qa3)) max(a3a(qa3)+3*a3b(qa3))])
        plot([1 numel(qa3)],[0 0],'k--')
        %
        set(gca,'XTick',1:numel(qa3),'XTickLabel',[param_matrix(qa3);],'FontSize',12,'FontWeight','Bold')
        title([method{A},' ',blur_tit,' : Union GT :  rel Im Blur'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        %
        saveGoodImg(H12,[dirPre,'../Documentation/Cosyne_2016/F_unionGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_relImBlur.jpg'],sizeGoodIm)
        close(H12)
    end
    
    
    
    % Plot Precision & Recall Separately (instead of combining them together as F-measure)
    % because maybe network methods effect one more than the other or something.
    if(plot_FRP_errorbars_all_params_each_method)
        H12_RP = figure; hold on
        errorbar([1:numel(qa3)]-0.3, a3aR(qa3), a3bR(qa3),'LineStyle','none','Marker','o','LineWidth',2,'Color','Blue');  % Recall
        errorbar([1:numel(qa3)]+0.3, a3aP(qa3), a3bP(qa3),'LineStyle','none','Marker','o','LineWidth',2,'Color','Red'); % Precision
        errorbar([1:numel(qa3)],     a3a(qa3),  a3b(qa2),'LineStyle','none','Marker','o','LineWidth',2,'Color','Green'); % F-measure
        %
        axis tight %([0.5 numel(qa3)+2.5 -0.6 0.6])
        plot([0.5 numel(qa3)+2.5],[0 0],'k--')
        set(gca,'XTick',1:numel(qa3)+2,'XTickLabel',[param_matrix(qa3);'ImPix';'ImBlur'],'FontSize',12,'FontWeight','Bold')
        title([method{A},' ',blur_tit,' : Union GT : relImBlur'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        legend('Recall','Precision','F-measure')
        %
        saveGoodImg(H12_RP,[dirPre,'../Documentation/Cosyne_2016/FRP_unionGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_relImBlur.jpg'],sizeGoodIm)
        close(H12_RP)
    end
    
    
    
    
    % Now get best parameter settings for given method A.
    %     for (justMethod, relImPix, relImBlur)
    % and for (benchmark, meanGT, maxGT, unionGT)
    %
    % Collect meanF, stdF, rM, ks.
    %
    % Should be lots of redundancy (repeated rM,ks param values) in this I hope.
    %
    % Benchmark P/R/F
    h = 1;
    maxF_meanAccImgs(A,h) = max(max(justMethod.maxF_meanAccImgs(:,:,A)));                                    
    [rM_max(A,h),ks_max(A,h)] = find( justMethod.maxF_meanAccImgs(:,:,A) == maxF_meanAccImgs(A,h) );
    maxR_meanAccImgs(A,h) = justMethod.maxR_meanAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxP_meanAccImgs(A,h) = justMethod.maxP_meanAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxF_stdAccImgs(A,h) = justMethod.maxF_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxR_stdAccImgs(A,h) = justMethod.maxR_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxP_stdAccImgs(A,h) = justMethod.maxP_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxF_skewAccImgs(A,h) = justMethod.maxF_skewAccImgs(rM_max(A,h),ks_max(A,h),A);
    %
    h = 2;
    maxF_meanAccImgs(A,h) = max(max(relImPix.maxF_meanAccImgs(:,:,A)));                                    
    [rM_max(A,h),ks_max(A,h)] = find( relImPix.maxF_meanAccImgs(:,:,A) == maxF_meanAccImgs(A,h) );
    maxR_meanAccImgs(A,h) = relImPix.maxR_meanAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxP_meanAccImgs(A,h) = relImPix.maxP_meanAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxF_stdAccImgs(A,h) = relImPix.maxF_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxR_stdAccImgs(A,h) = relImPix.maxR_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxP_stdAccImgs(A,h) = relImPix.maxP_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxF_skewAccImgs(A,h) = relImPix.maxF_skewAccImgs(rM_max(A,h),ks_max(A,h),A);
    %
    h = 3;
    maxF_meanAccImgs(A,h) = max(max(relImBlur.maxF_meanAccImgs(:,:,A)));                                    
    [rM_max(A,h),ks_max(A,h)] = find( relImBlur.maxF_meanAccImgs(:,:,A) == maxF_meanAccImgs(A,h) );
    maxR_meanAccImgs(A,h) = relImBlur.maxR_meanAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxP_meanAccImgs(A,h) = relImBlur.maxP_meanAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxF_stdAccImgs(A,h) = relImBlur.maxF_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxR_stdAccImgs(A,h) = relImBlur.maxR_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxP_stdAccImgs(A,h) = relImBlur.maxP_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxF_skewAccImgs(A,h) = relImBlur.maxF_skewAccImgs(rM_max(A,h),ks_max(A,h),A);
    %
    % mean (across GT's) of P/R/F
    h = 4;
    maxF_meanAccImgs(A,h) = max(max(justMethod.maxF_newMean_meanAccImgs(:,:,A)));                                    
    [rM_max(A,h),ks_max(A,h)] = find( justMethod.maxF_newMean_meanAccImgs(:,:,A) == maxF_meanAccImgs(A,h) );
    maxR_meanAccImgs(A,h) = justMethod.maxR_newMean_meanAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxP_meanAccImgs(A,h) = justMethod.maxP_newMean_meanAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxF_stdAccImgs(A,h) = justMethod.maxF_newMean_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxR_stdAccImgs(A,h) = justMethod.maxR_newMean_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxP_stdAccImgs(A,h) = justMethod.maxP_newMean_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxF_skewAccImgs(A,h) = justMethod.maxF_newMean_skewAccImgs(rM_max(A,h),ks_max(A,h),A);
    %
    h = 5;
    maxF_meanAccImgs(A,h) = max(max(relImPix.maxF_newMean_meanAccImgs(:,:,A)));                                    
    [rM_max(A,h),ks_max(A,h)] = find( relImPix.maxF_newMean_meanAccImgs(:,:,A) == maxF_meanAccImgs(A,h) );
    maxR_meanAccImgs(A,h) = relImPix.maxR_newMean_meanAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxP_meanAccImgs(A,h) = relImPix.maxP_newMean_meanAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxF_stdAccImgs(A,h) = relImPix.maxF_newMean_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxR_stdAccImgs(A,h) = relImPix.maxR_newMean_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxP_stdAccImgs(A,h) = relImPix.maxP_newMean_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxF_skewAccImgs(A,h) = relImPix.maxF_newMean_skewAccImgs(rM_max(A,h),ks_max(A,h),A);
    %
    h = 6;
    maxF_meanAccImgs(A,h) = max(max(relImBlur.maxF_newMean_meanAccImgs(:,:,A)));                                    
    [rM_max(A,h),ks_max(A,h)] = find( relImBlur.maxF_newMean_meanAccImgs(:,:,A) == maxF_meanAccImgs(A,h) );
    maxR_meanAccImgs(A,h) = relImBlur.maxR_newMean_meanAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxP_meanAccImgs(A,h) = relImBlur.maxP_newMean_meanAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxF_stdAccImgs(A,h) = relImBlur.maxF_newMean_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxR_stdAccImgs(A,h) = relImBlur.maxR_newMean_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxP_stdAccImgs(A,h) = relImBlur.maxP_newMean_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxF_skewAccImgs(A,h) = relImBlur.maxF_newMean_skewAccImgs(rM_max(A,h),ks_max(A,h),A);
    %
    % max (across GT's) of P/R/F
    h = 7;
    maxF_meanAccImgs(A,h) = max(max(justMethod.maxF_newMax_meanAccImgs(:,:,A)));                                    
    [rM_max(A,h),ks_max(A,h)] = find( justMethod.maxF_newMax_meanAccImgs(:,:,A) == maxF_meanAccImgs(A,h) );
    maxR_meanAccImgs(A,h) = justMethod.maxR_newMax_meanAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxP_meanAccImgs(A,h) = justMethod.maxP_newMax_meanAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxF_stdAccImgs(A,h) = justMethod.maxF_newMax_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxR_stdAccImgs(A,h) = justMethod.maxR_newMax_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxP_stdAccImgs(A,h) = justMethod.maxP_newMax_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxF_skewAccImgs(A,h) = justMethod.maxF_newMax_skewAccImgs(rM_max(A,h),ks_max(A,h),A);
    %
    h = 8;
    maxF_meanAccImgs(A,h) = max(max(relImPix.maxF_newMax_meanAccImgs(:,:,A)));                                    
    [rM_max(A,h),ks_max(A,h)] = find( relImPix.maxF_newMax_meanAccImgs(:,:,A) == maxF_meanAccImgs(A,h) );
    maxR_meanAccImgs(A,h) = relImPix.maxR_newMax_meanAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxP_meanAccImgs(A,h) = relImPix.maxP_newMax_meanAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxF_stdAccImgs(A,h) = relImPix.maxF_newMax_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxR_stdAccImgs(A,h) = relImPix.maxR_newMax_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxP_stdAccImgs(A,h) = relImPix.maxP_newMax_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxF_skewAccImgs(A,h) = relImPix.maxF_newMax_skewAccImgs(rM_max(A,h),ks_max(A,h),A);
    %
    h = 9;
    maxF_meanAccImgs(A,h) = max(max(relImBlur.maxF_newMax_meanAccImgs(:,:,A)));                                    
    [rM_max(A,h),ks_max(A,h)] = find( relImBlur.maxF_newMax_meanAccImgs(:,:,A) == maxF_meanAccImgs(A,h) );
    maxR_meanAccImgs(A,h) = relImBlur.maxR_newMax_meanAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxP_meanAccImgs(A,h) = relImBlur.maxP_newMax_meanAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxF_stdAccImgs(A,h) = relImBlur.maxF_newMax_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxR_stdAccImgs(A,h) = relImBlur.maxR_newMax_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxP_stdAccImgs(A,h) = relImBlur.maxP_newMax_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxF_skewAccImgs(A,h) = relImBlur.maxF_newMax_skewAccImgs(rM_max(A,h),ks_max(A,h),A);
    %
    % union of GTs to compute P/R/F
    h = 10;
    maxF_meanAccImgs(A,h) = max(max(justMethod.maxF_unionGT_meanAccImgs(:,:,A)));                                    
    [rM_max(A,h),ks_max(A,h)] = find( justMethod.maxF_unionGT_meanAccImgs(:,:,A) == maxF_meanAccImgs(A,h) );
    maxR_meanAccImgs(A,h) = justMethod.maxR_unionGT_meanAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxP_meanAccImgs(A,h) = justMethod.maxP_unionGT_meanAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxF_stdAccImgs(A,h) = justMethod.maxF_unionGT_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxR_stdAccImgs(A,h) = justMethod.maxR_unionGT_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxP_stdAccImgs(A,h) = justMethod.maxP_unionGT_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxF_skewAccImgs(A,h) = justMethod.maxF_unionGT_skewAccImgs(rM_max(A,h),ks_max(A,h),A);
    %
    h = 11;
    maxF_meanAccImgs(A,h) = max(max(relImPix.maxF_unionGT_meanAccImgs(:,:,A)));                                    
    [rM_max(A,h),ks_max(A,h)] = find( relImPix.maxF_unionGT_meanAccImgs(:,:,A) == maxF_meanAccImgs(A,h) );
    maxR_meanAccImgs(A,h) = relImPix.maxR_unionGT_meanAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxP_meanAccImgs(A,h) = relImPix.maxP_unionGT_meanAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxF_stdAccImgs(A,h) = relImPix.maxF_unionGT_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxR_stdAccImgs(A,h) = relImPix.maxR_unionGT_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxP_stdAccImgs(A,h) = relImPix.maxP_unionGT_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxF_skewAccImgs(A,h) = relImPix.maxF_unionGT_skewAccImgs(rM_max(A,h),ks_max(A,h),A);
    %
    h = 12;
    maxF_meanAccImgs(A,h) = max(max(relImBlur.maxF_unionGT_meanAccImgs(:,:,A)));                                    
    [rM_max(A,h),ks_max(A,h)] = find( relImBlur.maxF_unionGT_meanAccImgs(:,:,A) == maxF_meanAccImgs(A,h) );
    maxR_meanAccImgs(A,h) = relImBlur.maxR_unionGT_meanAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxP_meanAccImgs(A,h) = relImBlur.maxP_unionGT_meanAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxF_stdAccImgs(A,h) = relImBlur.maxF_unionGT_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxR_stdAccImgs(A,h) = relImBlur.maxR_unionGT_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxP_stdAccImgs(A,h) = relImBlur.maxP_unionGT_stdAccImgs(rM_max(A,h),ks_max(A,h),A);
    maxF_skewAccImgs(A,h) = relImBlur.maxF_unionGT_skewAccImgs(rM_max(A,h),ks_max(A,h),A);
    
    
    [rM_max; ks_max; maxF_meanAccImgs; maxF_stdAccImgs; maxF_skewAccImgs]
    
    

    
    
    
    
    
    
    % No Longer doing this because using our ("unweighted") average F-measure in figures Hmax & Hmax2 below...
    if(0)
        % Plot all against the Precision Recall curve of boundaries in the Raw Image Pixel Spatial Gradients
        evalDir = [dirPre,'images/BSDS_patch/101x301_ds1/benchmark_results/'];
        [H, evalRes] =  plot_eval(evalDir,'ko--',H,plot_max);

        disp(['Image Pixels - Strawman'])
        disp(['      Same TH -- R:',num2str(evalRes(2),2),' P:',num2str(evalRes(3),2),' F:',num2str(evalRes(4),2)])
        disp(['       Any TH -- R:',num2str(evalRes(5),2),' P:',num2str(evalRes(6),2),' F:',num2str(evalRes(7),2)])

        % Plot all against the Precision Recall curve of boundaries in the Blurred Image Pixel Spatial Gradients
        evalDir = [dirPre,'images/BSDS_patch/101x101_ds1/',blur_tag_I,'benchmark_results/'];
        [H, evalResB] =  plot_eval(evalDir,'go--',H,plot_max);

        disp(['Blurred Image Pixels - Strawman II'])
        disp(['      Same TH -- R:',num2str(evalResB(2),2),' P:',num2str(evalResB(3),2),' F:',num2str(evalResB(4),2)])
        disp(['       Any TH -- R:',num2str(evalResB(5),2),' P:',num2str(evalResB(6),2),' F:',num2str(evalResB(7),2)])

        F_max(numel(method)+1) = evalRes(7);
        F_max(numel(method)+2) = evalResB(7);
    end
    
    
    
    
    
    
    
    
    % For now, no longer plotting this.  Making Bar & ErrorBar plots.
    % Also, I am plotting mean & std using our non-weighted average.
    % Previously, the benchmark implements a weighted average where you add
    % up numerators & denominators of Precision & Recall before taking
    % ratio.
    if(plot_IsoF_flag)

        % Make up my own legend here...
        figure(H);
        hold on
        %
        text(0.67,1,{'\color{red}{rM1}','\color{blue}{rM3}','\color{cyan}{rM5}','\color{magenta}{rM10}'},'VerticalAlignment','top','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
        %
        scatter(0.65,0.85,150,'kx','filled')
        text(0.67,0.85,'= kssml','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
        %
        scatter(0.65,0.81,150,'ko','filled')
        text(0.67,0.81,'= ksmid','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
        %
        scatter(0.65,0.77,150,'k.','filled')
        text(0.67,0.77,'= kslrg','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
        %

        if(plot_max)
            % open vs. filled
            scatter(0.65,0.76,200,'ko','LineWidth',2)
            text(0.67,0.76,'= const THs','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
            %
            scatter(0.65,0.74,'ko','filled','LineWidth',2)
            text(0.67,0.74,'= vary THs','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
        else
            scatter(0.65,0.74,200,'ko','LineWidth',2)
            text(0.67,0.74,'= max F','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
        end

        %
        plot([0.63,0.67],[0.69,0.69],'g-','LineWidth',1.5)
        text(0.67,0.69,' \color{green}{= Iso F}','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
        %
        plot([0.63,0.67],[0.65,0.65],'ko--','LineWidth',1.5)
        text(0.67,0.65,'= ImPix','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
        %
        plot([0.63,0.67],[0.60,0.60],'go--','LineWidth',1.5)
        text(0.67,0.60,'= ImBlur','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')


        method_tag = method{A};
        method_tag(method_tag=='_')=' ';


        if(plot_max)
            max_tag = 'best&max';
        else
            max_tag = 'curve';
        end


        title([method_tag,' ',blur_tit,' : Kuramoto Coupled Oscilltor Sim.'],'FontSize',18,'FontWeight','Bold')

        hold off

        saveGoodImg(H,[dirPre,'../Documentation/Cosyne_2016/Overall_PR_curve_results_thinpbOFF/',method{A},blur_tag_M,'_',max_tag,'_Kur_allParams.jpg'],sizeGoodIm)
        close(H)
    end
    
    
    disp(method(A))
    
    
    % Loop thru GT files for image patches in each optimized method and build a matching vector that quantifies degree of gT agreement.
    for h = 1:numel(which_F_computation)

        if strcmp(method{A},'IsoDiff')
            evalDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/Kur_PIF_Fourier1/',method{A},'/benchmark_results/rM',rM{rM_max(A,h)},'/NF_60_0/ks',ks{ks_max(A,h)},'/'];
        else
            evalDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/Kur_PIF_Fourier1/',method{A},'/benchmark_results/rM',rM{rM_max(A,h)},'/sDInf/sP0p2/NF_60_0/ks',ks{ks_max(A,h)},'/'];
        end


        % Loop through each image patch and grab maxF. So later I can compute their mean and std.
        files = dir([evalDir,'*_ev1.txt']);


        for k = 1:numel(files)

            load([gTdir,files(k).name(1:end-8),'.mat'])

            gT_agreement{A,h}(k,1) = F_tot_stats(1);
            gT_agreement{A,h}(k,2) = F_tot_stats(2);
            %
            gT_agreement2{A,h}(k,1) = F_tot_stats2(1);
            gT_agreement2{A,h}(k,2) = F_tot_stats2(2);

        end

    end
    
    

    

end % loop over method A = 1:5


% F_all
% F_max
% rM(rM_max)
% ks(ks_max)


[rM_max; ks_max; maxF_meanAccImgs; maxF_stdAccImgs]













%% All this to get min and max for Precision & Recall for plotting 2D Histograms to be consistent across different methods.
RPFlims = zeros(6,3); % numel(which_F_computation)); % Just recording 3 different min & max values for the
% different relative_to_what values. Saving max & min lims on Recall, Precision, F-measure. 

for pp = 1:numel(which_F_computation) % loop thru combinations of F computation & what relative to (12)
    for xx = 1:numel(method) % loop thru different methods or methods_all.


        switch which_F_computation{pp}
            %
            case 'benchmarkF'
                %
                switch relative_to_what{pp}
                    case 'justMethod'
                        %disp('benchmark justMethod')
                        maxF_pre = maxF_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxR_pre = maxR_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxP_pre = maxP_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        r=1; % this indexes into max & min limits
                    case 'relImPix'
                        %disp('benchmark relPix')
                        maxF_pre = maxF_old_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxR_pre = maxR_old_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxP_pre = maxP_old_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        r=2; % this indexes into max & min limits
                    case 'relImBlur'
                        %disp('benchmark relBlur')
                        maxF_pre = maxF_old_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxR_pre = maxR_old_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxP_pre = maxP_old_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        r=3; % this indexes into max & min limits
                end 
                %
                maxF_post = maxF_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                maxR_post = maxR_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                maxP_post = maxP_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                %
                thr_method = thr_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                thr_imPix  = thr_old_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                thr_imBlur = thr_old_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
            %
            case 'meanGT'
                %
                switch relative_to_what{pp}
                    case 'justMethod'
                        %disp('meanGT justMethod')
                        maxF_pre = maxF_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxR_pre = maxR_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxP_pre = maxP_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        r=1; % this indexes into max & min limits
                    case 'relImPix'
                        %disp('meanGT relPix')
                        maxF_pre = maxF_new_mean_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxR_pre = maxR_new_mean_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxP_pre = maxP_new_mean_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        r=2; % this indexes into max & min limits
                    case 'relImBlur'
                        %disp('meanGT relBlur')
                        maxF_pre = maxF_new_mean_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxR_pre = maxR_new_mean_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxP_pre = maxP_new_mean_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        r=3; % this indexes into max & min limits
                end
                %
                maxF_post = maxF_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                maxR_post = maxR_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                maxP_post = maxP_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                %
                thr_method = thr_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                thr_imPix  = thr_new_mean_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                thr_imBlur = thr_new_mean_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
            %
            case 'maxGT'
                %
                switch relative_to_what{pp}
                    case 'justMethod'
                        %disp('maxGT justMethod')
                        maxF_pre = maxF_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxR_pre = maxR_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxP_pre = maxP_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        r=1; % this indexes into max & min limits
                    case 'relImPix'
                        %disp('maxGT relPix')
                        maxF_pre = maxF_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxR_pre = maxR_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxP_pre = maxP_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        r=2; % this indexes into max & min limits
                    case 'relImBlur'
                        %disp('maxGT relBlur')
                        maxF_pre = maxF_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxR_pre = maxR_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxP_pre = maxP_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        r=3; % this indexes into max & min limits
                end
                %
                maxF_post = maxF_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                maxR_post = maxR_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                maxP_post = maxP_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                %
                thr_method = thr_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                thr_imPix  = thr_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                thr_imBlur = thr_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                %
                bestGT_method = bestGT_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                bestGT_imPix = bestGT_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                bestGT_imBlur = bestGT_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};  
            %
            case 'unionGT'
                %
                switch relative_to_what{pp}
                    case 'justMethod'
                        %disp('unionGT justMethod')
                        maxF_pre = maxF_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxR_pre = maxR_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxP_pre = maxP_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        r=1; % this indexes into max & min limits
                    case 'relImPix'
                        %disp('unionGT relPix')
                        maxF_pre = maxF_unionGT_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxR_pre = maxR_unionGT_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxP_pre = maxP_unionGT_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        r=2; % this indexes into max & min limits
                    case 'relImBlur'
                        %disp('unionGT relBlur')
                        maxF_pre = maxF_unionGT_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxR_pre = maxR_unionGT_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxP_pre = maxP_unionGT_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        r=3; % this indexes into max & min limits
                end
                %
                maxF_post = maxF_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                maxR_post = maxR_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                maxP_post = maxP_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                %
                thr_method = thr_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                thr_imPix  = thr_unionGT_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                thr_imBlur = thr_unionGT_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
            %
        end
        %
        maxF_temp = maxF_post - maxF_pre; % this relImPix or relImBlur.
        maxR_temp = maxR_post - maxR_pre; % if justMethod, this will be all zeros
        maxP_temp = maxP_post - maxP_pre;



        % These are max & min values across different methods and F-measure computations.  
        % The r keeps track of different relative-to-what values.
        RPFlims(1,r) = max([ RPFlims(1,r); maxR_temp ]); % max R
        RPFlims(2,r) = min([ RPFlims(2,r); maxR_temp ]); % min R
        RPFlims(3,r) = max([ RPFlims(3,r); maxP_temp ]); % max P
        RPFlims(4,r) = min([ RPFlims(4,r); maxP_temp ]); % min P
        RPFlims(5,r) = max([ RPFlims(5,r); maxF_temp ]); % max F
        RPFlims(6,r) = min([ RPFlims(6,r); maxF_temp ]); % min F

    end % loop over xx = 1:5 method

end % loop over pp = 1:12 which_F_computation & relative_to_what





%% Here, we plot Violin Plots of F-measure comparing different optimized methods as well as
%  A 2D scatter plot of delta R vs delta P for all individual image patches.
%  Compute Statistical Significance too.
if(1)
    
    
    methods_all = {method{:}, 'ImPix', 'ImBlur'};
    
    % Ok, now that we have solved for consistent limits for plotting, we plot F-measure with 
    % error bars and a 2D Histogram in R-P space.
    
    for pp = 1:numel(which_F_computation) % loop thru combinations of F computation & what relative to (12)
                   % relative_to_what
    
        H=figure; hold on
        plot([0.5 numel(methods_all)+0.5],[0 0],'k--','LineWidth',1.5)
        
        Pr = ones(1,numel(method));
        Hr = zeros(1,numel(method));
        
        for xx = 1:numel(method) % loop thru different methods or methods_all.

            switch which_F_computation{pp}
                %
                case 'benchmarkF'
                    %
                    switch relative_to_what{pp}
                        case 'justMethod'
                            disp('benchmark justMethod')
                            maxF_pre = maxF_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxP_pre = maxP_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=1; % this indexes into max & min limits
                        case 'relImPix'
                            disp('benchmark relPix')
                            maxF_pre = maxF_old_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_old_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxP_pre = maxP_old_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=2; % this indexes into max & min limits
                        case 'relImBlur'
                            disp('benchmark relBlur')
                            maxF_pre = maxF_old_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_old_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxP_pre = maxP_old_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=3; % this indexes into max & min limits
                    end 
                    %
                    maxF_post = maxF_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    maxR_post = maxR_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    maxP_post = maxP_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    %
                    thr_method = thr_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    thr_imPix  = thr_old_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    thr_imBlur = thr_old_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                %
                case 'meanGT'
                    %
                    switch relative_to_what{pp}
                        case 'justMethod'
                            disp('meanGT justMethod')
                            maxF_pre = maxF_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxP_pre = maxP_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=1; % this indexes into max & min limits
                        case 'relImPix'
                            disp('meanGT relPix')
                            maxF_pre = maxF_new_mean_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_new_mean_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxP_pre = maxP_new_mean_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=2; % this indexes into max & min limits
                        case 'relImBlur'
                            disp('meanGT relBlur')
                            maxF_pre = maxF_new_mean_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_new_mean_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxP_pre = maxP_new_mean_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=3; % this indexes into max & min limits
                    end
                    %
                    maxF_post = maxF_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    maxR_post = maxR_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    maxP_post = maxP_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    %
                    thr_method = thr_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    thr_imPix  = thr_new_mean_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    thr_imBlur = thr_new_mean_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                %
                case 'maxGT'
                    %
                    switch relative_to_what{pp}
                        case 'justMethod'
                            disp('maxGT justMethod')
                            maxF_pre = maxF_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxP_pre = maxP_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=1; % this indexes into max & min limits
                        case 'relImPix'
                            disp('maxGT relPix')
                            maxF_pre = maxF_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxP_pre = maxP_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=2; % this indexes into max & min limits
                        case 'relImBlur'
                            disp('maxGT relBlur')
                            maxF_pre = maxF_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxP_pre = maxP_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=3; % this indexes into max & min limits
                    end
                    %
                    maxF_post = maxF_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    maxR_post = maxR_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    maxP_post = maxP_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    %
                    thr_method = thr_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    thr_imPix  = thr_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    thr_imBlur = thr_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    %
                    bestGT_method = bestGT_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    bestGT_imPix  = bestGT_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    bestGT_imBlur = bestGT_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx}; 
                %
                case 'unionGT'
                    %
                    switch relative_to_what{pp}
                        case 'justMethod'
                            disp('unionGT justMethod')
                            maxF_pre = maxF_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxP_pre = maxP_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=1; % this indexes into max & min limits
                        case 'relImPix'
                            disp('unionGT relPix')
                            maxF_pre = maxF_unionGT_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_unionGT_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxP_pre = maxP_unionGT_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=2; % this indexes into max & min limits
                        case 'relImBlur'
                            disp('unionGT relBlur')
                            maxF_pre = maxF_unionGT_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_unionGT_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxP_pre = maxP_unionGT_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=3; % this indexes into max & min limits
                    end
                    %
                    maxF_post = maxF_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    maxR_post = maxR_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    maxP_post = maxP_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    %
                    thr_method = thr_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    thr_imPix  = thr_unionGT_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    thr_imBlur = thr_unionGT_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                %
            end
            %
            maxF_temp = maxF_post - maxF_pre; % this relImPix or relImBlur.
            maxR_temp = maxR_post - maxR_pre; % if justMethod, this will be all zeros
            maxP_temp = maxP_post - maxP_pre;


            violin(maxF_temp,'x',[xx xx+1])

            
            % Determine Statistical Significance of difference in Network Method F-measure distribution across images vs ImBlur or ImPix
            % distribution.  Use Mann-Whitney U (aka ranksum) test.
            [Pr(xx),Hr(xx)] = ranksum(maxF_pre,maxF_post);
            % best_vec_for_ranksum{1,i}, best_vec_for_ranksum{2,i}
            
            
        end % loop over method xx=1:5 
        

        
        errorbar( [1:numel(method)], maxF_meanAccImgs(:,pp),maxF_stdAccImgs(:,pp), 'bd', 'Linewidth', 2, 'Color','green' )
%         errorbar( [1:numel(method)]-0.3, maxR_meanAccImgs(:,pp),maxR_stdAccImgs(:,pp), 'bd', 'Linewidth', 2, 'Color','blue')
%         errorbar( [1:numel(method)]+0.3, maxP_meanAccImgs(:,pp),maxP_stdAccImgs(:,pp), 'bd', 'Linewidth', 2, 'Color','red' )
        
        % Labels for Precision, Recall & F-measure error bars.
        text(1+0.0-0.1, maxF_meanAccImgs(1,pp)+0.02, ['F'], 'Color','green','FontSize',16,'FontWeight','Bold')
%         text(1-0.3-0.1, maxR_meanAccImgs(1,pp)+0.02, ['R'], 'Color','blue','FontSize',16,'FontWeight','Bold')
%         text(1+0.3-0.1, maxP_meanAccImgs(1,pp)+0.02, ['P'], 'Color','red','FontSize',16,'FontWeight','Bold')

        
        for xx = 1:numel(method)
            text(xx+0.1, maxF_meanAccImgs(xx,pp)+2*maxF_stdAccImgs(xx,pp) , {['\mu=',num2str( maxF_meanAccImgs(xx,pp) ,'%+5.2f')], ...
                ['\sigma=',num2str( maxF_stdAccImgs(xx,pp) ,'%+5.2f')],['s=',num2str( maxF_skewAccImgs(xx,pp) ,'%+5.2f')]}, ...
                'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold') % indicate the mean, std & skew
            
            if ( Pr(xx) < 0.001)
                text(xx+0.1, maxF_meanAccImgs(xx,pp)-2*maxF_stdAccImgs(xx,pp) , {['{rM',rM{rM_max(xx,pp)},',ks',ks{ks_max(xx,pp)},'}***'],['p=',num2str(Pr(xx))]}, ...
                'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold') % optimized parameters
            elseif ( Pr(xx) < 0.01)
                text(xx+0.1, maxF_meanAccImgs(xx,pp)-2*maxF_stdAccImgs(xx,pp) , {['{rM',rM{rM_max(xx,pp)},',ks',ks{ks_max(xx,pp)},'}**'],['p=',num2str(Pr(xx))]}, ...
                'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold') % optimized parameters
            elseif ( Pr(xx) < 0.05)
                text(xx+0.1, maxF_meanAccImgs(xx,pp)-2*maxF_stdAccImgs(xx,pp) , {['{rM',rM{rM_max(xx,pp)},',ks',ks{ks_max(xx,pp)},'}*'],['p=',num2str(Pr(xx))]}, ...
                'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold') % optimized parameters
            else
                text(xx+0.1, maxF_meanAccImgs(xx,pp)-2*maxF_stdAccImgs(xx,pp) , {['{rM',rM{rM_max(xx,pp)},',ks',ks{ks_max(xx,pp)},'}'],['p=',num2str(Pr(xx))]}, ...
                'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold') % optimized parameters
            end
        
        end
        grid on
        set(gca,'XTick',1:numel(methods_all),'XTickLabel',methods_all,'FontSize',16,'FontWeight','Bold')
        title([which_F_computation{pp},' : ',relative_to_what{pp}, ' : ',blur_tit],'FontSize',20,'FontWeight','Bold')
        xlabel('Method with Optimized Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['F-measure (across ~',num2str(round(mean(numFiles(:)))),' imgs)'],'FontSize',18,'FontWeight','Bold')
        axis([0.5 numel(methods_all)+0.5 min(maxF_meanAccImgs(:,pp)-3*maxF_stdAccImgs(:,pp)) max(maxF_meanAccImgs(:,pp)+3*maxF_stdAccImgs(:,pp))])
        %
        saveGoodImg(H,[dirPre,'../Documentation/Cosyne_2016/Fmax_compareMethodsKurBestParams',blur_tag_M,'_',which_F_computation{pp},'_',relative_to_what{pp},'.jpg'],sizeGoodIm)
        close(H)
        
        
        

        % Scatter plot for Recall vs Precsion & Marginal histograms for each : change for each method and each F-computation
        nBins = 100;
        
        for xx = 1:numel(method)
            
            switch which_F_computation{pp}
                %
                case 'benchmarkF'
                    %
                    switch relative_to_what{pp}
                        case 'justMethod'
                            disp('benchmark justMethod')
                            maxF_pre = maxF_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxP_pre = maxP_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=1; % this indexes into max & min limits
                        case 'relImPix'
                            disp('benchmark relPix')
                            maxF_pre = maxF_old_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_old_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxP_pre = maxP_old_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=2; % this indexes into max & min limits
                        case 'relImBlur'
                            disp('benchmark relBlur')
                            maxF_pre = maxF_old_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_old_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxP_pre = maxP_old_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=3; % this indexes into max & min limits
                    end 
                    %
                    maxF_post = maxF_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    maxR_post = maxR_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    maxP_post = maxP_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    %
                    thr_method = thr_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    thr_imPix  = thr_old_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    thr_imBlur = thr_old_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                %
                case 'meanGT'
                    %
                    switch relative_to_what{pp}
                        case 'justMethod'
                            disp('meanGT justMethod')
                            maxF_pre = maxF_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxP_pre = maxP_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=1; % this indexes into max & min limits
                        case 'relImPix'
                            disp('meanGT relPix')
                            maxF_pre = maxF_new_mean_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_new_mean_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxP_pre = maxP_new_mean_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=2; % this indexes into max & min limits
                        case 'relImBlur'
                            disp('meanGT relBlur')
                            maxF_pre = maxF_new_mean_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_new_mean_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxP_pre = maxP_new_mean_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=3; % this indexes into max & min limits
                    end
                    %
                    maxF_post = maxF_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    maxR_post = maxR_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    maxP_post = maxP_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    %
                    thr_method = thr_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    thr_imPix  = thr_new_mean_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    thr_imBlur = thr_new_mean_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                %
                case 'maxGT'
                    %
                    switch relative_to_what{pp}
                        case 'justMethod'
                            disp('maxGT justMethod')
                            maxF_pre = maxF_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxP_pre = maxP_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=1; % this indexes into max & min limits
                        case 'relImPix'
                            disp('maxGT relPix')
                            maxF_pre = maxF_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxP_pre = maxP_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=2; % this indexes into max & min limits
                        case 'relImBlur'
                            disp('maxGT relBlur')
                            maxF_pre = maxF_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxP_pre = maxP_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=3; % this indexes into max & min limits
                    end
                    %
                    maxF_post = maxF_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    maxR_post = maxR_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    maxP_post = maxP_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    %
                    thr_method = thr_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    thr_imPix  = thr_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    thr_imBlur = thr_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    %
                    bestGT_method = bestGT_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    bestGT_imPix  = bestGT_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    bestGT_imBlur = bestGT_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx}; 
                %
                case 'unionGT'
                    %
                    switch relative_to_what{pp}
                        case 'justMethod'
                            disp('unionGT justMethod')
                            maxF_pre = maxF_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxP_pre = maxP_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=1; % this indexes into max & min limits
                        case 'relImPix'
                            disp('unionGT relPix')
                            maxF_pre = maxF_unionGT_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_unionGT_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxP_pre = maxP_unionGT_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=2; % this indexes into max & min limits
                        case 'relImBlur'
                            disp('unionGT relBlur')
                            maxF_pre = maxF_unionGT_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_unionGT_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxP_pre = maxP_unionGT_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=3; % this indexes into max & min limits
                    end
                    %
                    maxF_post = maxF_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    maxR_post = maxR_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    maxP_post = maxP_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    %
                    thr_method = thr_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    thr_imPix  = thr_unionGT_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    thr_imBlur = thr_unionGT_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                %
            end
            %
            maxF_temp = maxF_post - maxF_pre; % this relImPix or relImBlur.
            maxR_temp = maxR_post - maxR_pre; % if justMethod, this will be all zeros
            maxP_temp = maxP_post - maxP_pre;
            
            
            maxX = RPFlims(1,r); % max(maxR_temp{rM_max(xx,pp),ks_max(xx,pp),xx});
            minX = RPFlims(2,r); % min( [min(maxR_temp{rM_max(xx,pp),ks_max(xx,pp),xx}), 0]);
            maxY = RPFlims(3,r); % max(maxP_temp{rM_max(xx,pp),ks_max(xx,pp),xx});
            minY = RPFlims(4,r); % min( [min(maxP_temp{rM_max(xx,pp),ks_max(xx,pp),xx}), 0]);
            maxF = RPFlims(5,r);
            minF = RPFlims(6,r);
            %
            
            r_lims = linspace(minX,maxX,nBins); % Recall
            p_lims = linspace(minY,maxY,nBins)'; % Precision
            f_lims = linspace(minF,maxF,nBins)'; % F-measure

            [nFcnts,f_lims] = hist(maxF_temp,f_lims);
            mf = mean(maxF_temp);
            sf = std(maxF_temp);
            skf = skewness(maxF_temp);
            %
            % find entries with F (or delta F) values above & below mean.
            ind1 = find(maxF_temp>0);  % improvement over imBlur or imPix
            ind2 = find(maxF_temp<=0); % no improvement
            n_better = numel(ind1);
            n_worse = numel(ind2);
            n_total = numel(maxF_temp);
            %
            nPcnts = hist(maxP_temp,p_lims);
            nPcnts1 = hist(maxP_temp(ind1),p_lims);
            nPcnts2 = hist(maxP_temp(ind2),p_lims);
            mp = mean(maxP_temp);
            sp = std(maxP_temp);
            skp = skewness(maxP_temp);
            %
            mp1 = mean(maxP_temp(ind1));
            mp2 = mean(maxP_temp(ind2));
            %
            sp1 = std(maxP_temp(ind1));
            sp2 = std(maxP_temp(ind2));
            
            %
            nRcnts = hist(maxR_temp,r_lims);
            nRcnts1 = hist(maxR_temp(ind1),r_lims);
            nRcnts2 = hist(maxR_temp(ind2),r_lims);
            mr = mean(maxR_temp);
            sr = std(maxR_temp);
            skr = skewness(maxR_temp);
            %
            mr1 = mean(maxR_temp(ind1));
            mr2 = mean(maxR_temp(ind2));
            %
            sr1 = std(maxR_temp(ind1));
            sr2 = std(maxR_temp(ind2));
            
            %
            Hout = hist2d([maxR_temp,maxP_temp],nBins,nBins,[minX maxX],[minY maxY]);
            close
            
            
            % Hout0 = hist2d([maxR_temp,maxP_temp],nBins,nBins,[minX maxX],[minY maxY]);
            % close
            Hout1 = hist2d([maxR_temp(ind1),maxP_temp(ind1)],nBins,nBins,[minX maxX],[minY maxY]);
            close
            Hout2 = hist2d([maxR_temp(ind2),maxP_temp(ind2)],nBins,nBins,[minX maxX],[minY maxY]);
            close
            
            
            
            

            
            % If I wanna scatter plot [ Note: I dont wanna do this for just method. (only for relImPix & relImBlur) ]
            if(1 & r~=1) 
                
                
                H=figure;

                % Plot the Precision-Recall 2D Scatter here.
                subplot(5,5,[1:4,6:9,11:14,16:19]);
                hold on, 
            
                scatter_patches(maxR_temp(ind1), maxP_temp(ind1), 3, 'b','o', 'FaceAlpha',0.3, 'EdgeColor','none');
                scatter_patches(maxR_temp(ind2), maxP_temp(ind2), 3, 'r','o', 'FaceAlpha',0.3, 'EdgeColor','none');
                %
                scatter(mr1,mp1,150,'b','x','LineWidth',1.5)
                scatter(mr2,mp2,150,'r','x','LineWidth',1.5)
                scatter(mr,mp,150,'k','+','LineWidth',2)
                %
                plot([0 0],[minY maxY],'k--')
                plot([minX maxX],[0 0],'k--')
                plot([min([minX,minY]) max([maxX,maxY])], [min([minX,minY]) max([maxX,maxY])],'k--')
                axis([minX maxX minY maxY])
                title([method{xx},' : ',which_F_computation{pp},' : ',relative_to_what{pp}, ' : ',blur_tit],'FontSize',20,'FontWeight','Bold')

                text( minX, maxY, {method{xx},['{rM',rM{rM_max(xx,pp)},',ks',ks{ks_max(xx,pp)},'}']}, ...
                    'Color','black','FontSize',18,'FontWeight','Bold','VerticalAlignment','top','HorizontalAlignment','left')
                ylabel('(-)  \DeltaP  (+)','FontSize',18,'FontWeight','Bold')
                set(gca,'FontSize',16,'FontWeight','Bold')

                text( maxX, maxY, {['When \DeltaF > \mu_F (N = ',num2str(n_better),')'],...
                    ['\mu_{(\DeltaR,\DeltaP)}=(',num2str(mr1,'%+5.2f'),',',num2str(mp1,'%+5.2f'),')'] ...
                    ['\sigma_{(\DeltaR,\DeltaP)}=',num2str(sr1,'%+5.2f'),',',num2str(sp1,'%+5.2f'),')']},...
                    'VerticalAlignment','top','HorizontalAlignment','right','FontSize',16,'FontWeight','Bold','Color','Blue','LineStyle','-')

                text( maxX, minY,{['When \DeltaF \leq \mu_F (N = ',num2str(n_worse),')'],...
                    ['\mu_{(\DeltaR,\DeltaP)}=(',num2str(mr2,'%+5.2f'),',',num2str(mp2,'%+5.2f'),')'] ...
                    ['\sigma_{(\DeltaR,\DeltaP)}=',num2str(sr2,'%+5.2f'),',',num2str(sp2,'%+5.2f'),')']},...
                    'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16,'FontWeight','Bold','Color','Red','LineStyle','-')
                


                % Plot Precision Marginalization on vertical pane.
                subplot(5,5,[5,10,15,20]), 
                hold on
                plot(nPcnts./sum(nPcnts(:)),p_lims,'k','LineWidth',2)
                plot(nPcnts1./sum(nPcnts(:)),p_lims,'b','LineWidth',2)
                plot(nPcnts2./sum(nPcnts(:)),p_lims,'r','LineWidth',2)
                %
                plot([0.7*max(nPcnts)./sum(nPcnts(:)), 0.7*max(nPcnts)./sum(nPcnts(:))], [mp2-sp2, mp2+sp2], 'Color','red','LineWidth',1.5)
                scatter(0.7*max(nPcnts)./sum(nPcnts(:)), mp2, 300, 'rx','LineWidth',1.5)
                plot([0.8*max(nPcnts)./sum(nPcnts(:)), 0.8*max(nPcnts)./sum(nPcnts(:))], [mp-sp, mp+sp], 'Color','black','LineWidth',1.5)
                scatter(0.8*max(nPcnts)./sum(nPcnts(:)), mp, 300, 'kx','LineWidth',1.5)
                plot([0.9*max(nPcnts)./sum(nPcnts(:)), 0.9*max(nPcnts)./sum(nPcnts(:))], [mp1-sp1, mp1+sp1], 'Color','blue','LineWidth',1.5)
                scatter(0.9*max(nPcnts)./sum(nPcnts(:)), mp1, 300, 'bx','LineWidth',1.5)
                %
                plot([0 max(nPcnts)./sum(nPcnts(:))], [0 0], 'k--')
                axis([0 max(nPcnts)./sum(nPcnts(:)) minY maxY])
                set(gca,'FontSize',16,'FontWeight','Bold')
                text( 0.5*max(nPcnts)./sum(nPcnts(:)), p_lims(round(0.25.*numel(p_lims))), ...
                    {['\mu_P=',num2str(mp,'%+5.2f')],['\sigma_P=',num2str(sp,'%+5.2f')],['s_P=',num2str(skp,'%+5.2f')]},...
                    'VerticalAlignment','top','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
                text( 0.8*max(nPcnts)./sum(nPcnts(:)), p_lims(round(0.90.*numel(p_lims))), ...
                    ['\DeltaP'],'VerticalAlignment','top','HorizontalAlignment','left','FontSize',20,'FontWeight','Bold')



                % Plot Recall Marginalization on horizontal pane.
                subplot(5,5,[21:24])
                hold on
                plot(r_lims,nRcnts./sum(nRcnts(:)),'k','LineWidth',2)
                plot(r_lims,nRcnts1./sum(nRcnts(:)),'b','LineWidth',2)
                plot(r_lims,nRcnts2./sum(nRcnts(:)),'r','LineWidth',2)
                %
                plot([mr2-sr2, mr2+sr2], [0.7*max(nRcnts)./sum(nRcnts(:)), 0.7*max(nRcnts)./sum(nRcnts(:))],'Color','red','LineWidth',1.5)
                scatter(mr2, 0.7*max(nRcnts)./sum(nRcnts(:)), 300, 'rx','LineWidth',1.5)
                plot([mr-sr, mr+sr], [0.8*max(nRcnts)./sum(nRcnts(:)), 0.8*max(nRcnts)./sum(nRcnts(:))],'Color','black','LineWidth',1.5)
                scatter(mr, 0.8*max(nRcnts)./sum(nRcnts(:)), 300, 'kx','LineWidth',1.5)
                plot([mr1-sr1, mr1+sr1], [0.9*max(nRcnts)./sum(nRcnts(:)), 0.9*max(nRcnts)./sum(nRcnts(:))],'Color','blue','LineWidth',1.5)
                scatter(mr1, 0.9*max(nRcnts)./sum(nRcnts(:)), 300, 'bx','LineWidth',1.5)
                %
                plot([0 0], [0 max(nRcnts)./sum(nRcnts(:))], 'k--')
                xlabel('(-)  \DeltaR  (+)','FontSize',18,'FontWeight','Bold')
                axis([minX maxX 0 max(nRcnts)./sum(nRcnts(:))])
                set(gca,'FontSize',16,'FontWeight','Bold')
                text( r_lims(round(0.75.*numel(r_lims))), 0.5*max(nRcnts)./sum(nRcnts(:)), ...
                    {['\mu_R=',num2str(mr,'%+5.2f')],['\sigma_R=',num2str(sr,'%+5.2f')],['s_R=',num2str(skr,'%+5.2f')]},...
                    'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')



                text( r_lims(round(0.10.*numel(r_lims))), 0.8*max(nRcnts)./sum(nRcnts(:)), ...
                    ['\DeltaR'],'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',20,'FontWeight','Bold')


                % Plot F-measure Marginalization
                subplot(5,5,[25])
                hold on

                ind1 = find(f_lims > 0); % was mf
                ind2 = find(f_lims <= 0); % was mf

                plot(f_lims,nFcnts./sum(nFcnts(:)),'g','LineWidth',2)

                plot(f_lims(ind1),nFcnts(ind1)./sum(nFcnts(:)),'b--','LineWidth',2)
                plot(f_lims(ind2),nFcnts(ind2)./sum(nFcnts(:)),'r--','LineWidth',2)

                plot([mf-sf, mf+sf], [0.5*max(nFcnts)./sum(nFcnts(:)), 0.5*max(nFcnts)./sum(nFcnts(:))],'Color','black','LineWidth',1.5)
                scatter(mf, 0.5*max(nFcnts)./sum(nFcnts(:)), 300, 'kx','LineWidth',1.5)
                plot([0 0], [0 max(nFcnts)./sum(nFcnts(:))], 'k--')
                xlabel('(-)  \DeltaF  (+)','FontSize',18,'FontWeight','Bold')
                axis([minF maxF 0 max(nFcnts)./sum(nFcnts(:))])
                set(gca,'FontSize',16,'FontWeight','Bold')
                text( f_lims(round(0.75.*numel(f_lims))), 0.5*max(nFcnts)./sum(nFcnts(:)), ...
                    {['\mu_F=',num2str(mf,'%+5.2f')],['\sigma_F=',num2str(sf,'%+5.2f')],['s_F=',num2str(skf,'%+5.2f')]},...
                    'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')


    %             % to investigate full top bin in histogram of SKH relImBlur maxGT & meanGT 
    %             if( xx==4 && ( pp==6 | pp==9) ) % Mod_SKH relImBlur meanGT, maxGT
    %                 keyboard
    %             end



                saveGoodImg(H,[dirPre,'../Documentation/Cosyne_2016/RPmax_',which_F_computation{pp},'_',relative_to_what{pp},'_',method{xx},'_rM',rM{rM_max(xx,pp)},'_ks',ks{ks_max(xx,pp)},blur_tag_M,'.jpg'],sizeGoodIm)
                close(H)
            
            end
            
            
        end % loop over xx = 1:5 for method.
        
        
    end % loop over pp = 1:12 for which_F_computation & relative_to_what
    
    
    
end



%% Scatter plot in R-P Space -  plot the vector from (R,P)_blur to (R,P)_method. 
%  Similar as above, but plot this instead of plotting delta R & delta P in a 2D plane from -1 to 1.
%  I want to  with similar color coding for points with delta F > or <= 0.
if(0)
    
    
    for xx = 1:numel(method) % loop thru different methods or methods_all.
    
        for pp = 1:numel(which_F_computation) % loop thru combinations of F computation & what relative to (12)
                   % relative_to_what

                    
                   
            switch which_F_computation{pp}
                %
                case 'benchmarkF'
                    maxF_post = maxF_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    maxR_post = maxR_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    maxP_post = maxP_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    %
                    thr_method = thr_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    thr_imPix  = thr_old_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    thr_imBlur = thr_old_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    %
                    switch relative_to_what{pp}
                        case 'justMethod'
                            continue % dont want to do this justMethod analysis (move on to next pp loop hopefully).
                        case 'relImPix'
                            maxF_pre = maxF_old_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_old_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx}; % should be pix.
                            maxP_pre = maxP_old_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=1; % this indexes into max & min limits
                        case 'relImBlur'
                            maxF_pre = maxF_old_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_old_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx}; % should be blur
                            maxP_pre = maxP_old_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=1; % this indexes into max & min limits
                    end
                %    
                case 'meanGT'
                    maxF_post = maxF_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    maxR_post = maxR_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    maxP_post = maxP_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    %
                    thr_method = thr_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    thr_imPix  = thr_new_mean_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    thr_imBlur = thr_new_mean_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    %
                    switch relative_to_what{pp}
                        case 'justMethod'
                            continue % dont want to do this justMethod analysis (move on to next pp loop hopefully).
                        case 'relImPix'
                            maxF_pre = maxF_new_mean_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_new_mean_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx}; % should be pix.
                            maxP_pre = maxP_new_mean_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=1; % this indexes into max & min limits
                        case 'relImBlur'
                            maxF_pre = maxF_new_mean_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_new_mean_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx}; % should be blur
                            maxP_pre = maxP_new_mean_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=1; % this indexes into max & min limits
                    end
                %    
                case 'maxGT'
                    maxF_post = maxF_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    maxR_post = maxR_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    maxP_post = maxP_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    %
                    thr_method = thr_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    thr_imPix  = thr_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    thr_imBlur = thr_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    %
                    bestGT_method = bestGT_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    bestGT_imPix = bestGT_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    bestGT_imBlur = bestGT_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx}; 
                    %
                    switch relative_to_what{pp}
                        case 'justMethod'
                            continue % dont want to do this justMethod analysis (move on to next pp loop hopefully).
                        case 'relImPix'
                            maxF_pre = maxF_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx}; % should be pix.
                            maxP_pre = maxP_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=1; % this indexes into max & min limits
                        case 'relImBlur'
                            maxF_pre = maxF_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx}; % should be blur
                            maxP_pre = maxP_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=1; % this indexes into max & min limits
                    end
                %    
                case 'unionGT'
                    maxF_post = maxF_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    maxR_post = maxR_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    maxP_post = maxP_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    %
                    thr_method = thr_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    thr_imPix  = thr_unionGT_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    thr_imBlur = thr_unionGT_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                    %
                    switch relative_to_what{pp}
                        case 'justMethod'
                            continue % dont want to do this justMethod analysis (move on to next pp loop hopefully).
                        case 'relImPix'
                            maxF_pre = maxF_unionGT_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_unionGT_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx}; % should be pix.
                            maxP_pre = maxP_unionGT_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=1; % this indexes into max & min limits
                        case 'relImBlur'
                            maxF_pre = maxF_unionGT_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            maxR_pre = maxR_unionGT_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx}; % should be blur
                            maxP_pre = maxP_unionGT_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                            r=1; % this indexes into max & min limits
                    end
                    
            end % switch which_f_computation
            

            
            maxX = RPFlims(1,r); 
            minX = RPFlims(2,r); 
            maxY = RPFlims(3,r); 
            minY = RPFlims(4,r); 
            maxF = RPFlims(5,r);
            minF = RPFlims(6,r);
            %
            
            


            % find entries with F (or delta F) values above & below mean.
            ind1 = find(maxF_post-maxF_pre>0); % was mf
            ind2 = find(maxF_post-maxF_pre<=0); % was mf
            n_better = numel(ind1);
            n_worse = numel(ind2);
            n_total = numel(maxF_temp);



            nBins = 50;
            r_lims = linspace(minX,maxX,nBins);  % Recall
            p_lims = linspace(minY,maxY,nBins)'; % Precisionfigur,e 
            f_lims = linspace(minF,maxF,nBins)'; % F-measure

            [nFcnts_pre0,f_lims] = hist(maxF_pre,f_lims);
            [nFcnts_pre1,f_lims] = hist(maxF_pre(ind1),f_lims);
            [nFcnts_pre2,f_lims] = hist(maxF_pre(ind2),f_lims);
            %
            [nFcnts_post0,f_lims] = hist(maxF_post,f_lims);
            [nFcnts_post1,f_lims] = hist(maxF_post(ind1),f_lims);
            [nFcnts_post2,f_lims] = hist(maxF_post(ind2),f_lims);
            
      
%             mf = mean(maxF_temp);
%             sf = std(maxF_temp);
%             skf = skewness(maxF_temp);
%             %


            [nPcnts_pre0,p_lims] = hist(maxP_pre,p_lims);
            [nPcnts_pre1,p_lims] = hist(maxP_pre(ind1),p_lims);
            [nPcnts_pre2,p_lims] = hist(maxP_pre(ind2),p_lims);
            %
            [nPcnts_post0,p_lims] = hist(maxP_post,p_lims);
            [nPcnts_post1,p_lims] = hist(maxP_post(ind1),p_lims);
            [nPcnts_post2,p_lims] = hist(maxP_post(ind2),p_lims);
            
            
%             mp = mean(maxP_temp);
%             sp = std(maxP_temp);
%             skp = skewness(maxP_temp);
%             %
%             mp1 = mean(maxP_temp(ind1));
%             mp2 = mean(maxP_temp(ind2));
%             %
%             sp1 = std(maxP_temp(ind1));
%             sp2 = std(maxP_temp(ind2));
            

            [nRcnts_pre0,r_lims] = hist(maxR_pre,r_lims);
            [nRcnts_pre1,r_lims] = hist(maxR_pre(ind1),r_lims);
            [nRcnts_pre2,r_lims] = hist(maxR_pre(ind2),r_lims);
            %
            [nRcnts_post0,r_lims] = hist(maxR_post,r_lims);
            [nRcnts_post1,r_lims] = hist(maxR_post(ind1),r_lims);
            [nRcnts_post2,r_lims] = hist(maxR_post(ind2),r_lims);
            

            
            
%             mr = mean(maxR_temp);
%             sr = std(maxR_temp);
%             skr = skewness(maxR_temp);
%             %
%             mr1 = mean(maxR_temp(ind1));
%             mr2 = mean(maxR_temp(ind2));
%             %
%             sr1 = std(maxR_temp(ind1));
%             sr2 = std(maxR_temp(ind2));
            
            
            H=figure;
            
            
            % Plot the Precision-Recall 2D Histogram here.
            s = subplot(5,5,[1:4,6:9,11:14,16:19]);
            hold on, 

            
            h1 = openfig('isoF.fig','reuse'); % open figure
            ax1 = gca; % get handle to axes of figure
            fig1 = get(ax1,'children'); %get handle to all the children in the figure
            copyobj(fig1,s); %copy children to new parent axes i.e. the subplot axes
            close(h1)
            figure(H)
            subplot(s)
            hold on
            scatter(0.70,0.90,1500,'w.','Linewidth',2) % cover the green dot in Iso-F contours
            colormap('bone')
            %
            ylabel('Recall','FontSize',18,'FontWeight','Bold')
            xlabel('Precision','FontSize',18,'FontWeight','Bold')
            % axis square

            axis([minX maxX minY maxY])




            % 1st draw lines connecting pre & post
            line([maxP_pre(ind1)';maxR_post(ind1)'],[maxR_pre(ind1)';maxP_post(ind1)'],'Color','cyan');
            line([maxP_pre(ind2)';maxR_post(ind2)'],[maxR_pre(ind2)';maxP_post(ind2)'],'Color','yellow');
            %
            % draw pre points
            scatter_patches(maxP_pre(ind1), maxR_pre(ind1), 3, 'c','s', 'FaceAlpha',1, 'EdgeColor','none');
            scatter_patches(maxP_pre(ind2), maxR_pre(ind2), 3, 'm','s', 'FaceAlpha',1, 'EdgeColor','Black');

            %
            % draw post points
            scatter_patches( maxR_post(ind1), maxP_post(ind1),4, 'b','o', 'FaceAlpha',1, 'EdgeColor','none');
            scatter_patches(maxR_post(ind2),maxP_post(ind2),  4, 'r','o', 'FaceAlpha',1, 'EdgeColor','none');



            title([method{xx},' : ',which_F_computation{pp},' : ',relative_to_what{pp}, ' : ',blur_tit],'FontSize',20,'FontWeight','Bold')

            

            % Plot Marginalized Precision Distributions splitting them up by images with improved F and degraded F.
            subplot(5,5,[21:24])
            hold on
            
            plot(p_lims,nPcnts_pre0./n_total,'k--','LineWidth',2)
            plot(p_lims,nPcnts_pre1./n_total,'c--','LineWidth',1)
            plot(p_lims,nPcnts_post1./n_total,'b-','LineWidth',2)
            plot(p_lims,nPcnts_pre2./n_total,'m--','LineWidth',1)
            plot(p_lims,nPcnts_post2./n_total,'r-','LineWidth',2)
            %
            legend({'before all','before better','after better','before worse','after worse',}) %  
            text(f_lims(round(0.9*numel(f_lims))), 0.3*max(nPcnts_pre0./n_total), ...
                {['Ntot = ',num2str(n_total)],['\color{blue}Nbetter = ',num2str(n_better)],['\color{red}Nworse = ',num2str(n_worse)]})
            ylabel('% counts')
            xlabel('Precision')
            
            
            
            
         
            
            % Plot Marginalized Recall Distributions splitting them up by images with improved F and degraded F.
            subplot(5,5,[5,10,15,20]), 
            hold on
            plot(nRcnts_pre0./n_total,r_lims,'k--','LineWidth',2)
            plot(nRcnts_pre1./n_total,r_lims,'c--','LineWidth',1)
            plot(nRcnts_post1./n_total,r_lims,'b-','LineWidth',2)
            plot(nRcnts_pre2./n_total,r_lims,'m--','LineWidth',1)
            plot(nRcnts_post2./n_total,r_lims,'r-','LineWidth',2)
            %
            ylabel('% counts')
            xlabel('Recall')
            
            
            
            
            
            % Plot F-measure "Marginalization"  Distributions splitting them up by images with improved F and degraded F.
            subplot(5,5,[25]), hold on
            plot(f_lims,nFcnts_pre0./n_total,'k--','LineWidth',2.5)
            plot(f_lims,nFcnts_pre1./n_total,'c--','LineWidth',1)
            plot(f_lims,nFcnts_post1./n_total,'b-','LineWidth',2.5)
            plot(f_lims,nFcnts_pre2./n_total,'m--','LineWidth',1)
            plot(f_lims,nFcnts_post2./n_total,'r-','LineWidth',2.5)
            %
            ylabel('% counts')
            xlabel('F-measure')
            
           

%             % to investigate full top bin in histogram of SKH relImBlur maxGT & meanGT 
%             if( xx==4 && ( pp==6 | pp==9) ) % Mod_SKH relImBlur meanGT, maxGT
%                 keyboard
%             end
            
            
            
            saveGoodImg(H,[dirPre,'../Documentation/Cosyne_2016/RPmax_lines_',which_F_computation{pp},'_',relative_to_what{pp},'_',method{xx},'_rM',rM{rM_max(xx,pp)},'_ks',ks{ks_max(xx,pp)},blur_tag_M,'.jpg'],sizeGoodIm)
            close(H)
            
            
        end % loop over pp = 1:12 for which_F_computation & relative_to_what
    
    end % loop over xx = 1:5 for method.
    
    
    
end





%% Here, I plot <Delta> F-measure before and after method vs. GroundTruth (Dis-)Agreement.
if(1)
    
    maxDelF = -1.*ones(size(which_F_computation));
    minDelF = 1.*ones(size(which_F_computation));
    
    for pp = 9 % 1:numel(which_F_computation) % loop thru combinations of F computation & what relative to (12)
                   % relative_to_what
                   
       for xx = 1:numel(method)
           
           switch which_F_computation{pp}

                    case 'benchmarkF'
                        maxF_post = maxF_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxR_post = maxR_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxP_post = maxP_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        %
                        thr_method = thr_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        thr_imPix  = thr_old_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        thr_imBlur = thr_old_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        %
                        switch relative_to_what{pp}
                            case 'justMethod'
                                continue % dont want to do this justMethod analysis (move on to next pp loop hopefully).
                            case 'relImPix'
                                maxF_pre = maxF_old_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                maxR_pre = maxR_old_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx}; % should be pix.
                                maxP_pre = maxP_old_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                r=1; % this indexes into max & min limits
                            case 'relImBlur'
                                maxF_pre = maxF_old_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                maxR_pre = maxR_old_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx}; % should be blur
                                maxP_pre = maxP_old_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                r=1; % this indexes into max & min limits
                        end

                    case 'meanGT'
                        maxF_post = maxF_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxR_post = maxR_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxP_post = maxP_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        %
                        thr_method = thr_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        thr_imPix  = thr_new_mean_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        thr_imBlur = thr_new_mean_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        %
                        switch relative_to_what{pp}
                            case 'justMethod'
                                continue % dont want to do this justMethod analysis (move on to next pp loop hopefully).
                            case 'relImPix'
                                maxF_pre = maxF_new_mean_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                maxR_pre = maxR_new_mean_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx}; % should be pix.
                                maxP_pre = maxP_new_mean_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                r=1; % this indexes into max & min limits
                            case 'relImBlur'
                                maxF_pre = maxF_new_mean_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                maxR_pre = maxR_new_mean_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx}; % should be blur
                                maxP_pre = maxP_new_mean_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                r=1; % this indexes into max & min limits
                        end

                    case 'maxGT'
                        maxF_post = maxF_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxR_post = maxR_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxP_post = maxP_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        %
                        thr_method = thr_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        thr_imPix  = thr_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        thr_imBlur = thr_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        %
                        bestGT_method = bestGT_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        bestGT_imPix = bestGT_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        bestGT_imBlur = bestGT_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx}; 
                        %
                        switch relative_to_what{pp}
                            case 'justMethod'
                                continue % dont want to do this justMethod analysis (move on to next pp loop hopefully).
                            case 'relImPix'
                                maxF_pre = maxF_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                maxR_pre = maxR_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx}; % should be pix.
                                maxP_pre = maxP_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                r=1; % this indexes into max & min limits
                            case 'relImBlur'
                                maxF_pre = maxF_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                maxR_pre = maxR_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx}; % should be blur
                                maxP_pre = maxP_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                r=1; % this indexes into max & min limits
                        end

                    case 'unionGT'
                        maxF_post = maxF_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxR_post = maxR_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxP_post = maxP_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        %
                        thr_method = thr_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        thr_imPix  = thr_unionGT_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        thr_imBlur = thr_unionGT_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        %
                        switch relative_to_what{pp}
                            case 'justMethod'
                                continue % dont want to do this justMethod analysis (move on to next pp loop hopefully).
                            case 'relImPix'
                                maxF_pre = maxF_unionGT_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                maxR_pre = maxR_unionGT_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx}; % should be pix.
                                maxP_pre = maxP_unionGT_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                r=1; % this indexes into max & min limits
                            case 'relImBlur'
                                maxF_pre = maxF_unionGT_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                maxR_pre = maxR_unionGT_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx}; % should be blur
                                maxP_pre = maxP_unionGT_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                r=1; % this indexes into max & min limits
                        end

            end % switch which_f_computation


            % Want to set the Y-axis of these plots to be same across Optimized Methods
            maxDelF(pp) = max( [ maxDelF(pp); maxF_post-maxF_pre ] );
            minDelF(pp) = min( [ minDelF(pp); maxF_post-maxF_pre ] );

       end % end loop over xx = 1:numel(method)
       %      
            
          
    
        for xx = 4 %1:numel(method) % loop thru different methods or methods_all.
            
            % Directory to probabalistic boundary png image. 
            pbDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/Kur_PIF_Fourier1/',method{xx},'/pb_png/rM',rM{rM_max(xx,pp)},'/sDInf/sP0p2/NF_60_0/ks',ks{ks_max(xx,pp)},'/'];  
        
            switch which_F_computation{pp}

                    case 'benchmarkF'
                        maxF_post = maxF_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxR_post = maxR_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxP_post = maxP_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        %
                        thr_method = thr_old_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        thr_imPix  = thr_old_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        thr_imBlur = thr_old_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        %
                        switch relative_to_what{pp}
                            case 'justMethod'
                                continue % dont want to do this justMethod analysis (move on to next pp loop hopefully).
                            case 'relImPix'
                                maxF_pre = maxF_old_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                maxR_pre = maxR_old_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx}; % should be pix.
                                maxP_pre = maxP_old_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                r=1; % this indexes into max & min limits
                            case 'relImBlur'
                                maxF_pre = maxF_old_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                maxR_pre = maxR_old_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx}; % should be blur
                                maxP_pre = maxP_old_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                r=1; % this indexes into max & min limits
                        end

                    case 'meanGT'
                        maxF_post = maxF_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxR_post = maxR_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxP_post = maxP_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        %
                        thr_method = thr_new_mean_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        thr_imPix  = thr_new_mean_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        thr_imBlur = thr_new_mean_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        %
                        switch relative_to_what{pp}
                            case 'justMethod'
                                continue % dont want to do this justMethod analysis (move on to next pp loop hopefully).
                            case 'relImPix'
                                maxF_pre = maxF_new_mean_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                maxR_pre = maxR_new_mean_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx}; % should be pix.
                                maxP_pre = maxP_new_mean_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                r=1; % this indexes into max & min limits
                            case 'relImBlur'
                                maxF_pre = maxF_new_mean_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                maxR_pre = maxR_new_mean_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx}; % should be blur
                                maxP_pre = maxP_new_mean_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                r=1; % this indexes into max & min limits
                        end

                    case 'maxGT'
                        maxF_post = maxF_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxR_post = maxR_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxP_post = maxP_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        %
                        thr_method = thr_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        thr_imPix  = thr_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        thr_imBlur = thr_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        %
                        bestGT_method = bestGT_new_max_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        bestGT_imPix = bestGT_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        bestGT_imBlur = bestGT_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx}; 
                        %
                        switch relative_to_what{pp}
                            case 'justMethod'
                                continue % dont want to do this justMethod analysis (move on to next pp loop hopefully).
                            case 'relImPix'
                                maxF_pre = maxF_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                maxR_pre = maxR_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx}; % should be pix.
                                maxP_pre = maxP_new_max_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                r=1; % this indexes into max & min limits
                            case 'relImBlur'
                                maxF_pre = maxF_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                maxR_pre = maxR_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx}; % should be blur
                                maxP_pre = maxP_new_max_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                r=1; % this indexes into max & min limits
                        end

                    case 'unionGT'
                        maxF_post = maxF_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxR_post = maxR_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        maxP_post = maxP_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        %
                        thr_method = thr_unionGT_struct_method_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        thr_imPix  = thr_unionGT_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        thr_imBlur = thr_unionGT_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                        %
                        switch relative_to_what{pp}
                            case 'justMethod'
                                continue % dont want to do this justMethod analysis (move on to next pp loop hopefully).
                            case 'relImPix'
                                maxF_pre = maxF_unionGT_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                maxR_pre = maxR_unionGT_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx}; % should be pix.
                                maxP_pre = maxP_unionGT_struct_imPix_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                r=1; % this indexes into max & min limits
                            case 'relImBlur'
                                maxF_pre = maxF_unionGT_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                maxR_pre = maxR_unionGT_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx}; % should be blur
                                maxP_pre = maxP_unionGT_struct_imBlur_only{rM_max(xx,pp),ks_max(xx,pp),xx};
                                r=1; % this indexes into max & min limits
                        end

            end % switch which_f_computation



            % Plot F_blur vs F_method with points colorcoded by gT agreement (boundary coincidence)
            H=figure;
            hold on
            %scatter(maxF_pre(ind1), maxF_post(ind1),'b.')
            %scatter(maxF_pre(ind2), maxF_post(ind2),'r.')
            scatter(maxF_pre, maxF_post,20,gT_agreement{xx,pp}(:,1),'filled')
            cb=colorbar;
            set(get(cb,'ylabel'),'string','gT co-similarity (direct overlap)');
            axis([0 1 0 1])
            plot([0 1],[0 1],'k--')
            title([method{xx},' : ',which_F_computation{pp},' : ',relative_to_what{pp}, ' : ',blur_tit],'FontSize',20,'FontWeight','Bold')
            xlabel('F-measure before method','FontSize',18,'FontWeight','Bold')
            ylabel('F-measure after method','FontSize',18,'FontWeight','Bold')
            set(gca,'FontSize',16,'FontWeight','Bold')
            text(1,0,{['Ntot = ',num2str(n_total)],['\color{red}Nbetter = ',num2str(n_better)],...
                ['\color{blue}Nworse = ',num2str(n_worse)]},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16,'FontWeight','Bold')
            axis square
            grid on
            %
            [Pr,Hr] = ranksum(maxF_pre, maxF_post);
            text(0.1, 0.9, ['p = ',num2str(Pr)], 'VerticalAlignment','top','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
            %
            saveGoodImg(H,[dirPre,'../Documentation/Cosyne_2016/FvsF_gTagree1_',which_F_computation{pp},'_',relative_to_what{pp},'_',method{xx},'_rM',rM{rM_max(xx,pp)},'_ks',ks{ks_max(xx,pp)},blur_tag_M,'.jpg'],sizeGoodIm)
            close(H)
            
            
            % Plot <Delta> F (pre-blur to post-method) vs gT agreement (boundary coincidence)
            [p,S,mu] = polyfit(gT_agreement{xx,pp}(:,1),maxF_post-maxF_pre,1);
            x = linspace( min(gT_agreement{xx,pp}(:,1)), max(gT_agreement{xx,pp}(:,1)));
            %
            H=figure;
            hold on
            scatter(gT_agreement{xx,pp}(:,1),maxF_post - maxF_pre,'filled')
            title([method{xx},' : ',which_F_computation{pp},' : ',relative_to_what{pp}, ' : ',blur_tit],'FontSize',20,'FontWeight','Bold')
            ylabel('\Delta F-measure before vs after method','FontSize',18,'FontWeight','Bold')
            xlabel('Ground Truth Similarity Measure (Direct Overlap)','FontSize',18,'FontWeight','Bold')
            set(gca,'FontSize',16,'FontWeight','Bold')
            plot([0 max(gT_agreement{xx,pp}(:,1))],[0 0],'k--')
            plot(x, p(1) + p(2).*x,'r--','LineWidth',2 )
            text(min(gT_agreement2{xx,pp}(:,1)), max(maxF_post-maxF_pre), ...
                { ['\color{red}y0 = ',num2str(p(1),2)], ['\color{red}m = ',num2str(p(2),2)] },...
                'VerticalAlignment','top','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
            grid on
            ylim([minDelF(pp),maxDelF(pp)]);
            saveGoodImg(H,[dirPre,'../Documentation/Cosyne_2016/delF_gTagree1_',which_F_computation{pp},'_',relative_to_what{pp},'_',method{xx},'_rM',rM{rM_max(xx,pp)},'_ks',ks{ks_max(xx,pp)},blur_tag_M,'.jpg'],sizeGoodIm)
            close(H)
            
            
            % Plot F-blur vs gT agreement (boundary coincidence)
            [p,S,mu] = polyfit(gT_agreement{xx,pp}(:,1),maxF_pre,1);
            x = linspace( min(gT_agreement{xx,pp}(:,1)), max(gT_agreement{xx,pp}(:,1)));
            %
            H=figure;
            hold on
            scatter(gT_agreement{xx,pp}(:,1),maxF_pre,'filled')
            title([method{xx},' : ',which_F_computation{pp},' : ',relative_to_what{pp}, ' : ',blur_tit],'FontSize',20,'FontWeight','Bold')
            ylabel('F-measure after blur / before method','FontSize',18,'FontWeight','Bold')
            xlabel('Ground Truth Similarity Measure (Direct Overlap)','FontSize',18,'FontWeight','Bold')
            set(gca,'FontSize',16,'FontWeight','Bold')
            plot([0 max(gT_agreement{xx,pp}(:,1))],[0 0],'k--')
            plot(x, p(1) + p(2).*x,'r--','LineWidth',2 )
            text(min(gT_agreement{xx,pp}(:,1)), max(maxF_pre), ...
                { ['\color{red}y0 = ',num2str(p(1),2)], ['\color{red}m = ',num2str(p(2),2)] },...
                'VerticalAlignment','top','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
            grid on
            ylim([min(maxF_pre),max(maxF_pre)]);
            saveGoodImg(H,[dirPre,'../Documentation/Cosyne_2016/Fblur_v_gTagree1_',which_F_computation{pp},'_',relative_to_what{pp},'_',method{xx},'_rM',rM{rM_max(xx,pp)},'_ks',ks{ks_max(xx,pp)},blur_tag_M,'.jpg'],sizeGoodIm)
            close(H)
            
            
            
            
            % % % % % % %
            
            
            % Plot F_blur vs F_method with points colorcoded by gT agreement (d=0.0075)
            H=figure;
            hold on
            scatter(maxF_pre, maxF_post,20,gT_agreement2{xx,pp}(:,1),'filled')
            cb=colorbar;
            set(get(cb,'ylabel'),'string','gT co-similarity (boundary closeness d=0.0075)');
            axis([0 1 0 1])
            plot([0 1],[0 1],'k--')
            title([method{xx},' : ',which_F_computation{pp},' : ',relative_to_what{pp}, ' : ',blur_tit],'FontSize',20,'FontWeight','Bold')
            xlabel('F-measure before method','FontSize',18,'FontWeight','Bold')
            ylabel('F-measure after method','FontSize',18,'FontWeight','Bold')
            set(gca,'FontSize',16,'FontWeight','Bold')
            text(1,0,{['Ntot = ',num2str(n_total)],['\color{red}Nbetter = ',num2str(n_better)],...
                ['\color{blue}Nworse = ',num2str(n_worse)]},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16,'FontWeight','Bold')
            axis square
            grid on
            %
            [Pr,Hr] = ranksum(maxF_pre, maxF_post);
            text(0.1, 0.9, ['p = ',num2str(Pr)], 'VerticalAlignment','top','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
            %
            saveGoodImg(H,[dirPre,'../Documentation/Cosyne_2016/FvsF_gTagree2_',which_F_computation{pp},'_',relative_to_what{pp},'_',method{xx},'_rM',rM{rM_max(xx,pp)},'_ks',ks{ks_max(xx,pp)},blur_tag_M,'.jpg'],sizeGoodIm)
            close(H)
            
            
            
            % Plot <Delta> F (pre-blur to post-method) vs gT agreement (d=0.0075)
            [p,S,mu] = polyfit(gT_agreement2{xx,pp}(:,1),maxF_post-maxF_pre,1);
            x = linspace( min(gT_agreement2{xx,pp}(:,1)), max(gT_agreement2{xx,pp}(:,1)));
            %
            H=figure;
            hold on
            scatter(gT_agreement2{xx,pp}(:,1),maxF_post - maxF_pre,'filled')
            title([method{xx},' : ',which_F_computation{pp},' : ',relative_to_what{pp}, ' : ',blur_tit],'FontSize',20,'FontWeight','Bold')
            ylabel('\Delta F-measure before vs after method','FontSize',18,'FontWeight','Bold')
            xlabel('Ground Truth Similarity Measure (boundary closeness d=0.0075)','FontSize',18,'FontWeight','Bold')
            set(gca,'FontSize',16,'FontWeight','Bold')
            plot([0 max(gT_agreement2{xx,pp}(:,1))],[0 0],'k--')
            plot(x, p(1) + p(2).*x,'r--','LineWidth',2 )
            text(min(gT_agreement2{xx,pp}(:,1)), max(maxF_post-maxF_pre), ...
                { ['\color{red}y0 = ',num2str(p(1),2)], ['\color{red}m = ',num2str(p(2),2)] },...
                'VerticalAlignment','top','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
            grid on
            ylim([minDelF(pp),maxDelF(pp)]);
            saveGoodImg(H,[dirPre,'../Documentation/Cosyne_2016/delF_gTagree2_',which_F_computation{pp},'_',relative_to_what{pp},'_',method{xx},'_rM',rM{rM_max(xx,pp)},'_ks',ks{ks_max(xx,pp)},blur_tag_M,'.jpg'],sizeGoodIm)
            close(H)
            
            
            
            % Plot F-blur vs gT agreement (d=0.0075)
            [p,S,mu] = polyfit(gT_agreement2{xx,pp}(:,1),maxF_pre,1);
            x = linspace( min(gT_agreement2{xx,pp}(:,1)), max(gT_agreement2{xx,pp}(:,1)));
            %
            H=figure;
            hold on
            scatter(gT_agreement2{xx,pp}(:,1),maxF_pre,'filled')
            title([method{xx},' : ',which_F_computation{pp},' : ',relative_to_what{pp}, ' : ',blur_tit],'FontSize',20,'FontWeight','Bold')
            ylabel('F-measure after blur / before method','FontSize',18,'FontWeight','Bold')
            xlabel('Ground Truth Similarity Measure (boundary closeness d=0.0075)','FontSize',18,'FontWeight','Bold')
            set(gca,'FontSize',16,'FontWeight','Bold')
            plot([0 max(gT_agreement2{xx,pp}(:,1))],[0 0],'k--')
            plot(x, p(1) + p(2).*x,'r--','LineWidth',2 )
            text(min(gT_agreement2{xx,pp}(:,1)), max(maxF_pre), ...
                { ['\color{red}y0 = ',num2str(p(1),2)], ['\color{red}m = ',num2str(p(2),2)] },...
                'VerticalAlignment','top','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
            grid on
            ylim([min(maxF_pre),max(maxF_pre)]);
            saveGoodImg(H,[dirPre,'../Documentation/Cosyne_2016/Fblur_v_gTagree2_',which_F_computation{pp},'_',relative_to_what{pp},'_',method{xx},'_rM',rM{rM_max(xx,pp)},'_ks',ks{ks_max(xx,pp)},blur_tag_M,'.jpg'],sizeGoodIm)
            close(H)
            
           
            % % % % % %
            
            
            % Plot F-blur vs <Delta> F 
            [p,S,mu] = polyfit(maxF_pre,maxF_post-maxF_pre,1);
            x = linspace( min(maxF_pre), max(maxF_pre) );
            %
            H=figure;
            hold on
            F_before = maxF_pre;
            delF_aft = maxF_post-maxF_pre;
            scatter(F_before,delF_aft,20,gT_agreement2{xx,pp}(:,1),'filled')
            cb=colorbar;
            set(get(cb,'ylabel'),'string','gT co-similarity (boundary closeness d=0.0075)');
            title([method{xx},' : ',which_F_computation{pp},' : ',relative_to_what{pp}, ' : ',blur_tit],'FontSize',20,'FontWeight','Bold')
            ylabel('\Delta F after method','FontSize',18,'FontWeight','Bold')
            xlabel('F-measure after blur / before method','FontSize',18,'FontWeight','Bold')
            set(gca,'FontSize',16,'FontWeight','Bold')
            plot([0 max(maxF_pre)],[0 0],'k--')
            plot(x, p(1) + p(2).*x,'r--','LineWidth',2 )
            text(min(maxF_pre), max(maxF_post-maxF_pre), ...
                { ['\color{red}y0 = ',num2str(p(1),2)], ['\color{red}m = ',num2str(p(2),2)] },...
                'VerticalAlignment','top','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
            grid on
            ylim([minDelF(pp),maxDelF(pp)]);

            
            % Within this plot, number some points of certain improvement quality (good, mid or bad) spread across different F_blur quality
            % and then visualize the image patches for them in another figure. (NOTE: Only wanna do for [ Mod_SK : MaxGT : RelImBlur ]
            if(xx==4 && pp==9)  % if Mod_SK, MaxGT, relImBlur
                
                % break up F_pre into N bins and then grab N points with certain <Delta> F properties within each bin.
                which_delF = 'good'; % options are: 'good', 'bad', 'mid', 'assort'
                N=8;
                F_blur_bins = linspace( min(maxF_pre), max(maxF_pre), N+1 );
                ind_tot = [];
                ind_tot_loc = [];
                
                for i = 1:N % Looping over F-blur bins
                    
                    ind_blur = find(F_before > F_blur_bins(i) & F_before <= F_blur_bins(i+1));
                    
                    [ppp,ind_delF] = sort(delF_aft(ind_blur));
                    
                    try
                        switch which_delF
                            case 'good'
                                ind = ind_blur(ind_delF(end-N+1:end));
                            case 'bad'
                                ind = ind_blur(ind_delF(1:N));
                            case 'mid'
                                ind = ind_blur(ind_delF(round(numel(ind_blur)./2)-N/2+1:round(numel(ind_blur)./2)+N/2));
                            case 'assort'
                                ind = randsample(ind_blur(ind_delF),N);
                                  % sort ind according to <Delta> F after randomly choosing.
                                  % ... Do it later.  THis doesnt make sense yet.
                        end
                    catch
                        % I think the problem is there is not enough points in this bin.  So just take em all.
                        ind = ind_blur;
                    end
                    
                    ind_loc = N*(i-1) + [1:numel(ind)];
                    ind_tot_loc = [ind_tot_loc,ind_loc];
                    ind_tot = [ind_tot;ind];
                    
                    % Label the image patches in the scatter plot with numbers
                    for j = 1:numel(ind)
                        text( F_before(ind(j)), delF_aft(ind(j)), num2str(ind_loc(j)) );
                    end
                    
                    
                    scatter(F_blur_bins,zeros(size(F_blur_bins)),300,'k+')
                    
%                     for j = 1:numel(F_blur_bins)
%                         
%                     end
                    
                    
                    
                end % i = 1:N % Looping over F-blur bins
            
            
            end % if Mod_SK, MaxGT, relImBlur
            
            
            
            saveGoodImg(H,[dirPre,'../Documentation/Cosyne_2016/Fblur_v_delF_',which_F_computation{pp},'_',relative_to_what{pp},'_',method{xx},'_rM',rM{rM_max(xx,pp)},'_ks',ks{ks_max(xx,pp)},blur_tag_M,'_',which_delF,'.jpg'],sizeGoodIm)
            close(H)
            
            
            
            
            
            % Make a figure with image patches arranged according to F-blur and <Delta> F after applying method (using ind_tot from above).
            if(xx==4 && pp==9)  % if Mod_SK, MaxGT, relImBlur
            
                % (1). Plot image patch, boundary found by method & best matching Ground Truth
                H=figure; 
                ha = tight_subplot(N, N, [0.01 0], [0.05 0.05], [0.05 0]);
                
                im_Fnames = img_ptch_name_struct{rM_max(xx,pp),ks_max(xx,pp),xx};
                
                Fblur = maxF_pre;
                delF = maxF_post-maxF_pre;
                
                for i = 1:numel(ind_tot)

                    load([imDir, im_Fnames{ind_tot(i)}, '.mat']);
                    load([gTdir, im_Fnames{ind_tot(i)}, '.mat']);

%                         % which gT ?
%                         bestGT_method(ind_tot(i))
%                         bestGT_imPix(ind_tot(i))
%                         bestGT_imBlur(ind_tot(i))
%                         
%                         % threshold pb and plot.
%                         thr_method(ind_tot(i))
%                         thr_imPix(ind_tot(i))
%                         thr_imBlur(ind_tot(i))        

                    % load pb png image and threshold it to show what boundaries let to max F-measure with maxGT criteria.
                    pbFile = [pbDir,im_Fnames{ind_tot(i)}, '.png'];
                    pb = double(imread(pbFile))/255;

%                     % load pb png image and threshold it to show what boundaries let to max F-measure with maxGT criteria.
%                     pbBlurFile = [pbDirBlur,im_Fnames{ind_tot(i)}, '.png'];
%                     pbBlur = double(imread(pbBlurFile))/255;

                    subplot(ha(ind_tot_loc(i))), hold on
                    %
%                     Igt1 = ( pbBlur >=thr_imBlur(ind_tot(i)) ); % Boundaries determined from thresholded probabalistic boundary from blur (that gave maxF)
%                     im0 = imoverlay(im,Igt1,[0 1 0]);
                    %
                    Igt2 = ( pb >=thr_method(ind_tot(i)) ); % Boundaries determined from thresholded probabalistic boundary from method (that gave maxF)
                    im1 = imoverlay(im,Igt2,[0 0 1]);
                    %
                    Igt3 = groundTruth{bestGT_method(ind_tot(i))}.Boundaries; % GT that best matched threshold pb and gave max F-measure using maxGT.
                    im2 = imoverlay(im1,Igt3,[1 0 0]);
                    imagesc(im2), 
                    axis square tight ij
                    set(gca,'XTick',[],'YTick',[])
                    if( bestGT_method(ind_tot(i)) == bestGT_imBlur(ind_tot(i)) )
                        text(size(im,1),size(im,2),num2str(i),'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',16,'FontWeight','Bold')
                    else
                        text(size(im,1),size(im,2),{'  *',num2str(i)},'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',16,'FontWeight','Bold')
                    end
                    ylabel(['\color{blue}{',num2str(Fblur(ind_tot(i)),2),'}','\color{red}{',num2str(delF(ind_tot(i)),'%+5.2f'),'}'])
                        
                end

                subplot(ha(N^2));
                xlabel(['(-) <---   \color{red}{\Delta F}  \color{black}{ ---> (+)}'])
                ylabel(['(+) <---   \color{blue}{F blur}  \color{black}{ ---> (-)}'])
                
                annotation('textbox', [0 0.9 1 0.1],'String',[method{xx},' : ',which_F_computation{pp},' : ',relative_to_what{pp}, ' : ',blur_tit, ' : ',which_delF, ' : (METHOD)'], ...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold')
                
%                 annotation('textbox', [0 0.5 1 0.1],'String',['(-) <---   \color{blue}{F blur}  \color{black}{ ---> (+)}'], ...
%                 'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold','Rotation',90);
% 
%                 annotation('textbox', [0 0 1 0.1],'String',['(-) <---   \color{red}{\Delta F}  \color{black}{ ---> (+)}'], ...
%                 'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold')

                saveGoodImg(H,[dirPre,'../Documentation/Cosyne_2016/Fblur_v_delF_',which_F_computation{pp},'_',relative_to_what{pp},'_',method{xx},'_rM',rM{rM_max(xx,pp)},'_ks',ks{ks_max(xx,pp)},blur_tag_M,'_',which_delF,'_visPatchesMethod.jpg'],sizeGoodIm)
                close(H)
                
                
                
                % (2). Plot image patch, boundary found by imBlur & best matching Ground Truth
                H=figure; 
                ha = tight_subplot(N, N, [0.01 0], [0.05 0.05], [0.05 0]);
                
                im_Fnames = img_ptch_name_struct{rM_max(xx,pp),ks_max(xx,pp),xx};
                
                Fblur = maxF_pre;
                delF = maxF_post-maxF_pre;
                
                for i = 1:numel(ind_tot)

                    load([imDir, im_Fnames{ind_tot(i)}, '.mat']);
                    load([gTdir, im_Fnames{ind_tot(i)}, '.mat']);
     

%                     % load pb png image and threshold it to show what boundaries let to max F-measure with maxGT criteria.
%                     pbFile = [pbDir,im_Fnames{ind_tot(i)}, '.png'];
%                     pb = double(imread(pbFile))/255;

                    % load pb png image and threshold it to show what boundaries let to max F-measure with maxGT criteria.
                    pbBlurFile = [pbDirBlur,im_Fnames{ind_tot(i)}, '.png'];
                    pbBlur = double(imread(pbBlurFile))/255;

                    subplot(ha(ind_tot_loc(i))), hold on
                    %
                    Igt1 = ( pbBlur >=thr_imBlur(ind_tot(i)) ); % Boundaries determined from thresholded probabalistic boundary from blur (that gave maxF)
                    im0 = imoverlay(im,Igt1,[0 1 0]);
                    %
%                     Igt2 = ( pb >=thr_method(ind_tot(i)) ); % Boundaries determined from thresholded probabalistic boundary from method (that gave maxF)
%                     im1 = imoverlay(im0,Igt2,[0 0 1]);
                    %
                    Igt3 = groundTruth{bestGT_imBlur(ind_tot(i))}.Boundaries; % GT that best matched threshold pb and gave max F-measure using maxGT.
                    im2 = imoverlay(im0,Igt3,[1 0 0]);
                    imagesc(im2), 
                    axis square tight ij
                    set(gca,'XTick',[],'YTick',[])
                    if( bestGT_method(ind_tot(i)) == bestGT_imBlur(ind_tot(i)) )
                        text(size(im,1),size(im,2),num2str(i),'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',16,'FontWeight','Bold')
                    else
                        text(size(im,1),size(im,2),{'  *',num2str(i)},'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',16,'FontWeight','Bold')
                    end
                    ylabel(['\color{blue}{',num2str(Fblur(ind_tot(i)),2),'}','\color{red}{',num2str(delF(ind_tot(i)),'%+5.2f'),'}'])
                        
                end

                subplot(ha(N^2));
                xlabel(['(-) <---   \color{red}{\Delta F}  \color{black}{ ---> (+)}'])
                ylabel(['(+) <---   \color{blue}{F blur}  \color{black}{ ---> (-)}'])
                
                annotation('textbox', [0 0.9 1 0.1],'String',[method{xx},' : ',which_F_computation{pp},' : ',relative_to_what{pp}, ' : ',blur_tit, ' : ',which_delF, ' : (BLUR)'], ...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold')
                
%                 annotation('textbox', [0 0.5 1 0.1],'String',['(-) <---   \color{blue}{F blur}  \color{black}{ ---> (+)}'], ...
%                 'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold','Rotation',90);
% 
%                 annotation('textbox', [0 0 1 0.1],'String',['(-) <---   \color{red}{\Delta F}  \color{black}{ ---> (+)}'], ...
%                 'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold')

                saveGoodImg(H,[dirPre,'../Documentation/Cosyne_2016/Fblur_v_delF_',which_F_computation{pp},'_',relative_to_what{pp},'_',method{xx},'_rM',rM{rM_max(xx,pp)},'_ks',ks{ks_max(xx,pp)},blur_tag_M,'_',which_delF,'_visPatchesBlur.jpg'],sizeGoodIm)
                close(H)
                
            
            end  % if Mod_SK, MaxGT, relImBlur
            
            
            
            
            
            
        end % loop over xx = 4 (Mod_SKH).

    end % loop over pp = 9 maxGT & relImBlur
    
    
    
    
    
end





keyboard
