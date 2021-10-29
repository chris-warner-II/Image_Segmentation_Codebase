% This script/function written by CW 11/15 to plot multiple precision
% recall curves on the same isoF plot.  We can use this to compare
% different methods or different parameter settings for the same method.


[dirPre,sizeGoodIm] = onCluster;
addpath([dirPre,'images/BSDS_images/BSR/bench/benchmarks/'])




sig = {'0p5','1','1p5','2','2p5','5','10'};
colors =  {'r','b','y','c','k','m','g'};

sz =  {'7','13','17','21','27','49','95'}; 
shapes = {'v-','d-','v-','^-','v-','s-','o-'};


H=open('isoF.fig');

plot_max = 1; % set to 1 to plot average of Precision & Recall leading to largest F, regardless of threshold
              % set to 0 to plot average of Precision & Recall at each threshold value for each image



%% Preallocate memory for arrays to hold statistics (across image patches) for each different blur.

% justMethod - Absolute value of performance. Not relative to any strawman or null model.

% (1). The way it was originally implemented in the benchmark code - Not Self Consistent.
justMethod.maxF_meanAccImgs = zeros(1,numel(sig)+1);
justMethod.maxF_stdAccImgs = zeros(1,numel(sig)+1);
%
justMethod.maxR_meanAccImgs = zeros(1,numel(sig)+1);
justMethod.maxR_stdAccImgs = zeros(1,numel(sig)+1);
%
justMethod.maxP_meanAccImgs = zeros(1,numel(sig)+1);
justMethod.maxP_stdAccImgs = zeros(1,numel(sig)+1);


% (2). Our reiteration of the original implementation (should be same as 1)
justMethod.maxF_old_meanAccImgs = zeros(1,numel(sig)+1);
justMethod.maxF_old_stdAccImgs = zeros(1,numel(sig)+1);
%
justMethod.maxR_old_meanAccImgs = zeros(1,numel(sig)+1);
justMethod.maxR_old_stdAccImgs = zeros(1,numel(sig)+1);
%
justMethod.maxP_old_meanAccImgs = zeros(1,numel(sig)+1);
justMethod.maxP_old_stdAccImgs = zeros(1,numel(sig)+1);


% (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
justMethod.maxF_newMean_meanAccImgs = zeros(1,numel(sig)+1);
justMethod.maxF_newMean_stdAccImgs = zeros(1,numel(sig)+1);
justMethod.maxF_newStd_meanAccImgs = zeros(1,numel(sig)+1);
%
justMethod.maxR_newMean_meanAccImgs = zeros(1,numel(sig)+1);
justMethod.maxR_newMean_stdAccImgs = zeros(1,numel(sig)+1);
justMethod.maxR_newStd_meanAccImgs = zeros(1,numel(sig)+1);
%
justMethod.maxP_newMean_meanAccImgs = zeros(1,numel(sig)+1);
justMethod.maxP_newMean_stdAccImgs = zeros(1,numel(sig)+1);
justMethod.maxP_newStd_meanAccImgs = zeros(1,numel(sig)+1);


% (4). New way to compute P-R-F - (Max value using best single ground truther.). 
justMethod.maxF_newMax_meanAccImgs = zeros(1,numel(sig)+1);
justMethod.maxF_newMax_stdAccImgs = zeros(1,numel(sig)+1);
%
justMethod.maxR_newMax_meanAccImgs = zeros(1,numel(sig)+1);
justMethod.maxR_newMax_stdAccImgs = zeros(1,numel(sig)+1);
%
justMethod.maxP_newMax_meanAccImgs = zeros(1,numel(sig)+1);
justMethod.maxP_newMax_stdAccImgs = zeros(1,numel(sig)+1);


% (5). New way to compute R & F - Kinda matching what benchmark was
% doing to compute Precision originally. Unioning all Ground Truthers
% into one GroundTruth and computing P & R & F using that unioned GT.
justMethod.maxF_unionGT_meanAccImgs = zeros(1,numel(sig)+1);
justMethod.maxF_unionGT_stdAccImgs = zeros(1,numel(sig)+1);
%
justMethod.maxR_unionGT_meanAccImgs = zeros(1,numel(sig)+1);
justMethod.maxR_unionGT_stdAccImgs = zeros(1,numel(sig)+1);
%
justMethod.maxP_unionGT_meanAccImgs = zeros(1,numel(sig)+1);
justMethod.maxP_unionGT_stdAccImgs = zeros(1,numel(sig)+1);


% relImPix - Performance measures relative to Raw image Pixels.

% (1). The way it was originally implemented in the benchmark code - Not Self Consistent.
relImPix.maxF_meanAccImgs = zeros(1,numel(sig)+1);
relImPix.maxF_stdAccImgs = zeros(1,numel(sig)+1);
%
relImPix.maxR_meanAccImgs = zeros(1,numel(sig)+1);
relImPix.maxR_stdAccImgs = zeros(1,numel(sig)+1);
%
relImPix.maxP_meanAccImgs = zeros(1,numel(sig)+1);
relImPix.maxP_stdAccImgs = zeros(1,numel(sig)+1);


% (2). Our reiteration of the original implementation (should be same as 1)
relImPix.maxF_old_meanAccImgs = zeros(1,numel(sig)+1);
relImPix.maxF_old_stdAccImgs = zeros(1,numel(sig)+1);
%
relImPix.maxR_old_meanAccImgs = zeros(1,numel(sig)+1);
relImPix.maxR_old_stdAccImgs = zeros(1,numel(sig)+1);
%
relImPix.maxP_old_meanAccImgs = zeros(1,numel(sig)+1);
relImPix.maxP_old_stdAccImgs = zeros(1,numel(sig)+1);


% (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
relImPix.maxF_newMean_meanAccImgs = zeros(1,numel(sig)+1);
relImPix.maxF_newMean_stdAccImgs = zeros(1,numel(sig)+1);
relImPix.maxF_newStd_meanAccImgs = zeros(1,numel(sig)+1);
%
relImPix.maxR_newMean_meanAccImgs = zeros(1,numel(sig)+1);
relImPix.maxR_newMean_stdAccImgs = zeros(1,numel(sig)+1);
relImPix.maxR_newStd_meanAccImgs = zeros(1,numel(sig)+1);
%
relImPix.maxP_newMean_meanAccImgs = zeros(1,numel(sig)+1);
relImPix.maxP_newMean_stdAccImgs = zeros(1,numel(sig)+1);
relImPix.maxP_newStd_meanAccImgs = zeros(1,numel(sig)+1);


% (4). New way to compute P-R-F - (Max value using best single ground truther.). 
relImPix.maxF_newMax_meanAccImgs = zeros(1,numel(sig)+1);
relImPix.maxF_newMax_stdAccImgs = zeros(1,numel(sig)+1);
%
relImPix.maxR_newMax_meanAccImgs = zeros(1,numel(sig)+1);
relImPix.maxR_newMax_stdAccImgs = zeros(1,numel(sig)+1);
%
relImPix.maxP_newMax_meanAccImgs = zeros(1,numel(sig)+1);
relImPix.maxP_newMax_stdAccImgs = zeros(1,numel(sig)+1);


% (5). New way to compute R & F - Kinda matching what benchmark was
% doing to compute Precision originally. Unioning all Ground Truthers
% into one GroundTruth and computing P & R & F using that unioned GT.
relImPix.maxF_unionGT_meanAccImgs = zeros(1,numel(sig)+1);
relImPix.maxF_unionGT_stdAccImgs = zeros(1,numel(sig)+1);
%
relImPix.maxR_unionGT_meanAccImgs = zeros(1,numel(sig)+1);
relImPix.maxR_unionGT_stdAccImgs = zeros(1,numel(sig)+1);
%
relImPix.maxP_unionGT_meanAccImgs = zeros(1,numel(sig)+1);
relImPix.maxP_unionGT_stdAccImgs = zeros(1,numel(sig)+1);



numFiles = zeros(1,numel(sig)+1);



%% Do analysis on raw image pixels first (before doing blurring by different sigma)
%  NOTE: To be able to find measures relative to ImPix, you have to do this
%  inside the blur loop, I think.  No other way... Is there?
%
% Plot all against the Precision Recall curve of boundaries in the Unblurred Raw Image Pixel Spatial Gradients
evalDir = [dirPre,'images/BSDS_patch/101x101_ds1/benchmark_results/'];

i=1;

[H, evalRes] = plot_eval(evalDir,[colors{i},shapes{i}],H,plot_max); % 


% Loop through each image patch and grab maxF. So later I can compute their mean and std.
files = dir([evalDir,'*_ev1.txt']);

maxF_imPix = zeros(1,numel(files));            % original benchmark way of computing P/R/F - not self consistent
maxR_imPix = zeros(1,numel(files)); 
maxP_imPix = zeros(1,numel(files)); 

maxF_old_imPix = zeros(1,numel(files));        % Our reimplementation of benchmark P/R/F (should be same as above)
maxR_old_imPix = zeros(1,numel(files));
maxP_old_imPix = zeros(1,numel(files));

maxF_unionGT_imPix = zeros(1,numel(files));    % Take union of all K groundtruths as single groundtruth before computing P/R/F
maxR_unionGT_imPix = zeros(1,numel(files));
maxP_unionGT_imPix = zeros(1,numel(files));

maxF_new_max_imPix = zeros(1,numel(files));    % For these next 3, we compute P/R/F for each groundtruth associated with an
maxR_new_max_imPix = zeros(1,numel(files));    % image patch.  Here we just grab the max value.  What is the algorithm's
maxP_new_max_imPix = zeros(1,numel(files));    % best match to any single human?

maxF_new_mean_imPix = zeros(1,numel(files));   % For a single image patch, the average P/R/F across the different groundtruths
maxR_new_mean_imPix = zeros(1,numel(files));
maxP_new_mean_imPix = zeros(1,numel(files));

maxF_new_std_imPix = zeros(1,numel(files));    % Standard Deviation of P/R/F across different groundtruths for single image ptch
maxR_new_std_imPix = zeros(1,numel(files));
maxP_new_std_imPix = zeros(1,numel(files));

numFiles(i) = numel(files);              % number of image patches processed for a given method or blur

for k = 1:numel(files)

    filename = fullfile(evalDir,files(k).name);
    AA  = dlmread(filename);
    %cntR = AA(:, 2);
    sumR = AA(:, 3);
    %cntP = AA(:, 4);
    sumP = AA(:, 5);
    %
    cntR0 = AA(:, 6);                      % CW: added these new ways to compute pixel correspondence in evaluation_bdry_image.
    cntP0 = AA(:, 7);

    sumR0 = AA(:, 9);

    R_old = AA(:, 9);
    P_old = AA(:, 10);

    R_new_mean = AA(:, 11);
    R_new_std = AA(:, 12);
    R_new_max = AA(:, 13);

    P_new_mean = AA(:, 14);
    P_new_std = AA(:, 15);
    P_new_max = AA(:, 16);

    R_unionGT = AA(:, 17);

    num_gTs = AA(:, 18);

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



    maxF_old_imPix(k) = max(F_old);
    ind = find(F_old==max(F_old));
    ind = ind(1);
    maxR_old_imPix(k) = R_old(ind);
    maxP_old_imPix(k) = P_old(ind);


    maxF_new_max_imPix(k) = max(F_new_max);
    ind = find(F_new_max==max(F_new_max));
    ind = ind(1);
    maxR_new_max_imPix(k) = R_new_max(ind);
    maxP_new_max_imPix(k) = P_new_max(ind);



    maxF_new_mean_imPix(k) = max(F_new_mean);
    ind = find(F_new_mean==max(F_new_mean));
    ind = ind(1);
    maxR_new_mean_imPix(k) = R_new_mean(ind);
    maxP_new_mean_imPix(k) = P_new_mean(ind);
    maxF_new_std_imPix(k) = F_new_std(ind);   % std across groundtruths as threshold that gave max avgF value.



    maxF_unionGT_imPix(k) = max(F_unionGT);
    ind = find(F_unionGT==max(F_unionGT));
    ind = ind(1);
    maxR_unionGT_imPix(k) = R_unionGT(ind);
    maxP_unionGT_imPix(k) = P(ind);

    k

end % looping through image patches


% justMethod - Take Mean & Std across image patches of Precision,Recall,F-measure computed in different ways.

% (1). The way it was originally implemented in the benchmark code - Not Self Consistent.
justMethod.maxF_meanAccImgs(i) = mean(maxF_imPix);
justMethod.maxF_stdAccImgs(i) = std(maxF_imPix);
%
justMethod.maxR_meanAccImgs(i) = mean(maxR_imPix);
justMethod.maxR_stdAccImgs(i) = std(maxR_imPix);
%
justMethod.maxP_meanAccImgs(i) = mean(maxP_imPix);
justMethod.maxP_stdAccImgs(i) = std(maxP_imPix);



% (2). Our reiteration of the original implementation (should be same as 1)
justMethod.maxF_old_meanAccImgs(i) = mean(maxF_old_imPix);
justMethod.maxF_old_stdAccImgs(i) = std(maxF_old_imPix);
%
justMethod.maxR_old_meanAccImgs(i) = mean(maxR_old_imPix);
justMethod.maxR_old_stdAccImgs(i) = std(maxR_old_imPix);
%
justMethod.maxP_old_meanAccImgs(i) = mean(maxP_old_imPix);
justMethod.maxP_old_stdAccImgs(i) = std(maxP_old_imPix);


% (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
justMethod.maxF_newMean_meanAccImgs(i) = mean(maxF_new_mean_imPix);
justMethod.maxF_newMean_stdAccImgs(i) = std(maxF_new_mean_imPix);
justMethod.maxF_newStd_meanAccImgs(i) = mean(maxF_new_std_imPix);
%
justMethod.maxR_newMean_meanAccImgs(i) = mean(maxR_new_mean_imPix);
justMethod.maxR_newMean_stdAccImgs(i) = std(maxR_new_mean_imPix);
justMethod.maxR_newStd_meanAccImgs(i) = mean(maxR_new_std_imPix);
%
justMethod.maxP_newMean_meanAccImgs(i) = mean(maxP_new_mean_imPix);
justMethod.maxP_newMean_stdAccImgs(i) = std(maxP_new_mean_imPix);
justMethod.maxP_newStd_meanAccImgs(i) = mean(maxP_new_std_imPix);


% (4). New way to compute P-R-F - (Max value using best single ground truther.). 
justMethod.maxF_newMax_meanAccImgs(i) = mean(maxF_new_max_imPix);
justMethod.maxF_newMax_stdAccImgs(i) = std(maxF_new_max_imPix);
%
justMethod.maxR_newMax_meanAccImgs(i) = mean(maxR_new_max_imPix);
justMethod.maxR_newMax_stdAccImgs(i) = std(maxR_new_max_imPix);
%
justMethod.maxP_newMax_meanAccImgs(i) = mean(maxP_new_max_imPix);
justMethod.maxP_newMax_stdAccImgs(i) = std(maxP_new_max_imPix);


% (5). New way to compute R & F - Kinda matching what benchmark was
% doing to compute Precision originally. Unioning all Ground Truthers
% into one GroundTruth and computing P & R & F using that unioned GT.
justMethod.maxF_unionGT_meanAccImgs(i) = mean(maxF_unionGT_imPix);
justMethod.maxF_unionGT_stdAccImgs(i) = std(maxF_unionGT_imPix);
%
justMethod.maxR_unionGT_meanAccImgs(i) = mean(maxR_unionGT_imPix);
justMethod.maxR_unionGT_stdAccImgs(i) = std(maxR_unionGT_imPix);
%
justMethod.maxP_unionGT_meanAccImgs(i) = mean(maxP_unionGT_imPix);
justMethod.maxP_unionGT_stdAccImgs(i) = std(maxP_unionGT_imPix);






% relImPix - calculate P/R/F different ways of method relative to raw image pixels.

% (1). The way it was originally implemented in the benchmark code - Not Self Consistent.
relImPix.maxF_meanAccImgs(i) = 0;
relImPix.maxF_stdAccImgs(i) = 0;
%
relImPix.maxR_meanAccImgs(i) = 0;
relImPix.maxR_stdAccImgs(i) = 0;
%
relImPix.maxP_meanAccImgs(i) = 0;
relImPix.maxP_stdAccImgs(i) = 0;



% (2). Our reiteration of the original implementation (should be same as 1)
relImPix.maxF_old_meanAccImgs(i) = 0;
relImPix.maxF_old_stdAccImgs(i) = 0;
%
relImPix.maxR_old_meanAccImgs(i) = 0;
relImPix.maxR_old_stdAccImgs(i) = 0;
%
relImPix.maxP_old_meanAccImgs(i) = 0;
relImPix.maxP_old_stdAccImgs(i) = 0;


% (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
relImPix.maxF_newMean_meanAccImgs(i) = 0;
relImPix.maxF_newMean_stdAccImgs(i) = 0;
relImPix.maxF_newStd_meanAccImgs(i) = 0;
%
relImPix.maxR_newMean_meanAccImgs(i) = 0;
relImPix.maxR_newMean_stdAccImgs(i) = 0;
relImPix.maxR_newStd_meanAccImgs(i) = 0;
%
relImPix.maxP_newMean_meanAccImgs(i) = 0;
relImPix.maxP_newMean_stdAccImgs(i) = 0;
relImPix.maxP_newStd_meanAccImgs(i) = 0;


% (4). New way to compute P-R-F - (Max value using best single ground truther.). 
relImPix.maxF_newMax_meanAccImgs(i) = 0;
relImPix.maxF_newMax_stdAccImgs(i) = 0;
%
relImPix.maxR_newMax_meanAccImgs(i) = 0;
relImPix.maxR_newMax_stdAccImgs(i) = 0;
%
relImPix.maxP_newMax_meanAccImgs(i) = 0;
relImPix.maxP_newMax_stdAccImgs(i) = 0;


% (5). New way to compute R & F - Kinda matching what benchmark was
% doing to compute Precision originally. Unioning all Ground Truthers
% into one GroundTruth and computing P & R & F using that unioned GT.
relImPix.maxF_unionGT_meanAccImgs(i) = 0;
relImPix.maxF_unionGT_stdAccImgs(i) = 0;
%
relImPix.maxR_unionGT_meanAccImgs(i) = 0;
relImPix.maxR_unionGT_stdAccImgs(i) = 0;
%
relImPix.maxP_unionGT_meanAccImgs(i) = 0;
relImPix.maxP_unionGT_stdAccImgs(i) = 0;










%% Repeat analysis in For Loop for each different Blurring amount

for i = 2:numel(sig)+1

    evalDir = [dirPre,'images/BSDS_patch/101x101_ds1/blur_sz',sz{i-1},'_sig',sig{i-1},'/benchmark_results/'];

    [H, evalRes] = plot_eval(evalDir,[colors{i-1},shapes{i-1}],H,plot_max); % 
    
    
    
    
    % Loop through each image patch and grab maxF. So later I can compute their mean and std.
    files = dir([evalDir,'*_ev1.txt']);

    maxF = zeros(1,numel(files));            % original benchmark way of computing P/R/F - not self consistent
    maxR = zeros(1,numel(files)); 
    maxP = zeros(1,numel(files)); 
    
    maxF_old = zeros(1,numel(files));        % Our reimplementation of benchmark P/R/F (should be same as above)
    maxR_old = zeros(1,numel(files));
    maxP_old = zeros(1,numel(files));
    
    maxF_unionGT = zeros(1,numel(files));    % Take union of all K groundtruths as single groundtruth before computing P/R/F
    maxR_unionGT = zeros(1,numel(files));
    maxP_unionGT = zeros(1,numel(files));
    
    maxF_new_max = zeros(1,numel(files));    % For these next 3, we compute P/R/F for each groundtruth associated with an
    maxR_new_max = zeros(1,numel(files));    % image patch.  Here we just grab the max value.  What is the algorithm's
    maxP_new_max = zeros(1,numel(files));    % best match to any single human?
    
    maxF_new_mean = zeros(1,numel(files));   % For a single image patch, the average P/R/F across the different groundtruths
    maxR_new_mean = zeros(1,numel(files));
    maxP_new_mean = zeros(1,numel(files));
    
    maxF_new_std = zeros(1,numel(files));    % Standard Deviation of P/R/F across different groundtruths for single image ptch
    maxR_new_std = zeros(1,numel(files));
    maxP_new_std = zeros(1,numel(files));

    numFiles(i) = numel(files);              % number of image patches processed for a given method or blur

    for k = 1:numel(files)

        filename = fullfile(evalDir,files(k).name);
        AA  = dlmread(filename);
        %cntR = AA(:, 2);
        sumR = AA(:, 3);
        %cntP = AA(:, 4);
        sumP = AA(:, 5);
        %
        cntR0 = AA(:, 6);                      % CW: added these new ways to compute pixel correspondence in evaluation_bdry_image.
        cntP0 = AA(:, 7);
        
        sumR0 = AA(:, 9);
        
        R_old = AA(:, 9);
        P_old = AA(:, 10);
        
        R_new_mean = AA(:, 11);
        R_new_std = AA(:, 12);
        R_new_max = AA(:, 13);
        
        P_new_mean = AA(:, 14);
        P_new_std = AA(:, 15);
        P_new_max = AA(:, 16);
        
        R_unionGT = AA(:, 17);
        
        num_gTs = AA(:, 18);

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
        
        
        
        maxF_old(k) = max(F_old);
        ind = find(F_old==max(F_old));
        ind = ind(1);
        maxR_old(k) = R_old(ind);
        maxP_old(k) = P_old(ind);
        
        
        maxF_new_max(k) = max(F_new_max);
        ind = find(F_new_max==max(F_new_max));
        ind = ind(1);
        maxR_new_max(k) = R_new_max(ind);
        maxP_new_max(k) = P_new_max(ind);
        
        
        
        maxF_new_mean(k) = max(F_new_mean);
        ind = find(F_new_mean==max(F_new_mean));
        ind = ind(1);
        maxR_new_mean(k) = R_new_mean(ind);
        maxP_new_mean(k) = P_new_mean(ind);
        maxF_new_std(k) = F_new_std(ind);   % std across groundtruths as threshold that gave max avgF value.
        
        
        
        maxF_unionGT(k) = max(F_unionGT);
        ind = find(F_unionGT==max(F_unionGT));
        ind = ind(1);
        maxR_unionGT(k) = R_unionGT(ind);
        maxP_unionGT(k) = P(ind);
        
        
        

        k

    end % looping through image patches
    
    
    
    
    
    
    
    
    
    % justMethod - Take Mean & Std across image patches of Precision,Recall,F-measure computed in different ways.
    
    % (1). The way it was originally implemented in the benchmark code - Not Self Consistent.
    justMethod.maxF_meanAccImgs(i) = mean(maxF);
    justMethod.maxF_stdAccImgs(i) = std(maxF);
    %
    justMethod.maxR_meanAccImgs(i) = mean(maxR);
    justMethod.maxR_stdAccImgs(i) = std(maxR);
    %
    justMethod.maxP_meanAccImgs(i) = mean(maxP);
    justMethod.maxP_stdAccImgs(i) = std(maxP);
    
    % (2). Our reiteration of the original implementation (should be same as 1)
    justMethod.maxF_old_meanAccImgs(i) = mean(maxF_old);
    justMethod.maxF_old_stdAccImgs(i) = std(maxF_old);
    %
    justMethod.maxR_old_meanAccImgs(i) = mean(maxR_old);
    justMethod.maxR_old_stdAccImgs(i) = std(maxR_old);
    %
    justMethod.maxP_old_meanAccImgs(i) = mean(maxP_old);
    justMethod.maxP_old_stdAccImgs(i) = std(maxP_old);
    
    
    % (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
    justMethod.maxF_newMean_meanAccImgs(i) = mean(maxF_new_mean);
    justMethod.maxF_newMean_stdAccImgs(i) = std(maxF_new_mean);
    justMethod.maxF_newStd_meanAccImgs(i) = mean(maxF_new_std);
    %
    justMethod.maxR_newMean_meanAccImgs(i) = mean(maxR_new_mean);
    justMethod.maxR_newMean_stdAccImgs(i) = std(maxR_new_mean);
    justMethod.maxR_newStd_meanAccImgs(i) = mean(maxR_new_std);
    %
    justMethod.maxP_newMean_meanAccImgs(i) = mean(maxP_new_mean);
    justMethod.maxP_newMean_stdAccImgs(i) = std(maxP_new_mean);
    justMethod.maxP_newStd_meanAccImgs(i) = mean(maxP_new_std);
    
    
    % (4). New way to compute P-R-F - (Max value using best single ground truther.). 
    justMethod.maxF_newMax_meanAccImgs(i) = mean(maxF_new_max);
    justMethod.maxF_newMax_stdAccImgs(i) = std(maxF_new_max);
    %
    justMethod.maxR_newMax_meanAccImgs(i) = mean(maxR_new_max);
    justMethod.maxR_newMax_stdAccImgs(i) = std(maxR_new_max);
    %
    justMethod.maxP_newMax_meanAccImgs(i) = mean(maxP_new_max);
    justMethod.maxP_newMax_stdAccImgs(i) = std(maxP_new_max);
    
    
    % (5). New way to compute R & F - Kinda matching what benchmark was
    % doing to compute Precision originally. Unioning all Ground Truthers
    % into one GroundTruth and computing P & R & F using that unioned GT.
    justMethod.maxF_unionGT_meanAccImgs(i) = mean(maxF_unionGT);
    justMethod.maxF_unionGT_stdAccImgs(i) = std(maxF_unionGT);
    %
    justMethod.maxR_unionGT_meanAccImgs(i) = mean(maxR_unionGT);
    justMethod.maxR_unionGT_stdAccImgs(i) = std(maxR_unionGT);
    %
    justMethod.maxP_unionGT_meanAccImgs(i) = mean(maxP_unionGT);
    justMethod.maxP_unionGT_stdAccImgs(i) = std(maxP_unionGT);
    
    
    
    
    % relImPix - Different ways of computing P/R/F all relative to raw image pixels performance.
    
    % (1). The way it was originally implemented in the benchmark code - Not Self Consistent.
    relImPix.maxF_meanAccImgs(i) = mean(maxF-maxF_imPix);
    relImPix.maxF_stdAccImgs(i) = std(maxF-maxF_imPix);
    %
    relImPix.maxR_meanAccImgs(i) = mean(maxR-maxR_imPix);
    relImPix.maxR_stdAccImgs(i) = std(maxR-maxR_imPix);
    %
    relImPix.maxP_meanAccImgs(i) = mean(maxP-maxP_imPix);
    relImPix.maxP_stdAccImgs(i) = std(maxP-maxP_imPix);
    
    % (2). Our reiteration of the original implementation (should be same as 1)
    relImPix.maxF_old_meanAccImgs(i) = mean(maxF_old-maxF_old_imPix);
    relImPix.maxF_old_stdAccImgs(i) = std(maxF_old-maxF_old_imPix);
    %
    relImPix.maxR_old_meanAccImgs(i) = mean(maxR_old-maxR_old_imPix);
    relImPix.maxR_old_stdAccImgs(i) = std(maxR_old-maxR_old_imPix);
    %
    relImPix.maxP_old_meanAccImgs(i) = mean(maxP_old-maxP_old_imPix);
    relImPix.maxP_old_stdAccImgs(i) = std(maxP_old-maxP_old_imPix);
    
    
    % (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
    relImPix.maxF_newMean_meanAccImgs(i) = mean(maxF_new_mean-maxF_new_mean_imPix);
    relImPix.maxF_newMean_stdAccImgs(i) = std(maxF_new_mean-maxF_new_mean_imPix);
    relImPix.maxF_newStd_meanAccImgs(i) = mean(maxF_new_std-maxF_new_std_imPix);
    %
    relImPix.maxR_newMean_meanAccImgs(i) = mean(maxR_new_mean-maxR_new_mean_imPix);
    relImPix.maxR_newMean_stdAccImgs(i) = std(maxR_new_mean-maxR_new_mean_imPix);
    relImPix.maxR_newStd_meanAccImgs(i) = mean(maxR_new_std-maxR_new_std_imPix);
    %
    relImPix.maxP_newMean_meanAccImgs(i) = mean(maxP_new_mean-maxP_new_mean_imPix);
    relImPix.maxP_newMean_stdAccImgs(i) = std(maxP_new_mean-maxP_new_mean_imPix);
    relImPix.maxP_newStd_meanAccImgs(i) = mean(maxP_new_std-maxP_new_std_imPix);
    
    
    % (4). New way to compute P-R-F - (Max value using best single ground truther.). 
    relImPix.maxF_newMax_meanAccImgs(i) = mean(maxF_new_max-maxF_new_max_imPix);
    relImPix.maxF_newMax_stdAccImgs(i) = std(maxF_new_max-maxF_new_max_imPix);
    %
    relImPix.maxR_newMax_meanAccImgs(i) = mean(maxR_new_max-maxR_new_max_imPix);
    relImPix.maxR_newMax_stdAccImgs(i) = std(maxR_new_max-maxR_new_max_imPix);
    %
    relImPix.maxP_newMax_meanAccImgs(i) = mean(maxP_new_max-maxP_new_max_imPix);
    relImPix.maxP_newMax_stdAccImgs(i) = std(maxP_new_max-maxP_new_max_imPix);
    
    
    % (5). New way to compute R & F - Kinda matching what benchmark was
    % doing to compute Precision originally. Unioning all Ground Truthers
    % into one GroundTruth and computing P & R & F using that unioned GT.
    relImPix.maxF_unionGT_meanAccImgs(i) = mean(maxF_unionGT-maxF_unionGT_imPix);
    relImPix.maxF_unionGT_stdAccImgs(i) = std(maxF_unionGT-maxF_unionGT_imPix);
    %
    relImPix.maxR_unionGT_meanAccImgs(i) = mean(maxR_unionGT-maxR_unionGT_imPix);
    relImPix.maxR_unionGT_stdAccImgs(i) = std(maxR_unionGT-maxR_unionGT_imPix);
    %
    relImPix.maxP_unionGT_meanAccImgs(i) = mean(maxP_unionGT-maxP_unionGT_imPix);
    relImPix.maxP_unionGT_stdAccImgs(i) = std(maxP_unionGT-maxP_unionGT_imPix);
    
    
    
    
    
    
    
    
    
    
    
    
    % Plotting histograms of all different P/R/F calculations for each blur.
    
    if(0)
        F_hist_bins = linspace(0,1,20);

        N_maxF = hist(maxF,F_hist_bins);
        N_maxF_old = hist(maxF_old,F_hist_bins);

        N_maxF_new_max =  hist(maxF_new_max,F_hist_bins);
        N_maxF_new_mean =  hist(maxF_new_mean,F_hist_bins);

        N_maxF_unionGT =  hist(maxF_unionGT,F_hist_bins);



        figure, hold on,
        plot(F_hist_bins, N_maxF,'k')
        plot(F_hist_bins, N_maxF_old,'r--')

        plot(F_hist_bins, N_maxF_new_max,'b')
        plot(F_hist_bins, N_maxF_new_mean,'g')

        plot(F_hist_bins, N_maxF_unionGT,'c')

        xlabel('F-measure')
        ylabel('counts')
        title(['Blurring \sigma = ',sig{i-1}])

        legend('bench','bench2','best gT','mean across gT','union gT')
    end
    

    

    % Note: evalRes = [bestT,bestR,bestP,bestF,R_max,P_max,F_max,Area_PR]
    disp(['sig:',sig{i-1},' sz:',sz{i-1}])
    disp(['      Same TH -- R:',num2str(evalRes(2),2),' P:',num2str(evalRes(3),2),' F:',num2str(evalRes(4),2)])
    disp(['       Any TH -- R:',num2str(evalRes(5),2),' P:',num2str(evalRes(6),2),' F:',num2str(evalRes(7),2)])


end % Looping over different Gaussian Blur sigma values (i)






%% Plot mean & std of different P/R/F measures for justMethod.
H1=figure; hold on
errorbar([1:numel(sig)+1]+0.0, justMethod.maxF_meanAccImgs, justMethod.maxF_stdAccImgs,'k.-','LineWidth',2)
errorbar([1:numel(sig)+1]+0.05, justMethod.maxF_old_meanAccImgs, justMethod.maxF_old_stdAccImgs,'r^--','LineWidth',2)
errorbar([1:numel(sig)+1]+0.1, justMethod.maxF_newMax_meanAccImgs, justMethod.maxF_newMax_stdAccImgs,'bs--','LineWidth',2)
errorbar([1:numel(sig)+1]+0.15, justMethod.maxF_newMean_meanAccImgs, justMethod.maxF_newMean_stdAccImgs,'gs--','LineWidth',2)
errorbar([1:numel(sig)+1]+0.2, justMethod.maxF_unionGT_meanAccImgs, justMethod.maxF_unionGT_stdAccImgs,'co--','LineWidth',2)
legend('bench','bench2','best gT','mean across gT','union gT')
set(gca,'XTick',1:numel(sig)+1,'XTickLabel',{'0',sig{:}},'FontSize',16,'FontWeight','Bold')
grid on
xlabel('\sigma Value for Gaussian Blurring','FontSize',18,'FontWeight','Bold')
ylabel('F-measure (across ~500 images)','FontSize',18,'FontWeight','Bold')
title('Comparing Different Precision - Recall - Fmeasure computations vs. Image Blurring','FontSize',20,'FontWeight','Bold')
saveGoodImg(H1,[dirPre,'../Documentation/Cosyne_2016/GaussianBlur_1.jpg'],sizeGoodIm)




%% Plot mean & std of different P/R/F measures for relImPix.
H2=figure; hold on
errorbar([1:numel(sig)+1]+0.0, relImPix.maxF_meanAccImgs, relImPix.maxF_stdAccImgs,'k.-','LineWidth',2)
errorbar([1:numel(sig)+1]+0.05, relImPix.maxF_old_meanAccImgs, relImPix.maxF_old_stdAccImgs,'r^--','LineWidth',2)
errorbar([1:numel(sig)+1]+0.1, relImPix.maxF_newMax_meanAccImgs, relImPix.maxF_newMax_stdAccImgs,'bs--','LineWidth',2)
errorbar([1:numel(sig)+1]+0.15, relImPix.maxF_newMean_meanAccImgs, relImPix.maxF_newMean_stdAccImgs,'gs--','LineWidth',2)
errorbar([1:numel(sig)+1]+0.2, relImPix.maxF_unionGT_meanAccImgs, relImPix.maxF_unionGT_stdAccImgs,'co--','LineWidth',2)
legend('bench','bench2','best gT','mean across gT','union gT')
set(gca,'XTick',1:numel(sig)+1,'XTickLabel',{'0',sig{:}},'FontSize',16,'FontWeight','Bold')
grid on
xlabel('\sigma Value for Gaussian Blurring','FontSize',18,'FontWeight','Bold')
ylabel('F-measure (across ~500 images)','FontSize',18,'FontWeight','Bold')
title('Comparing Different Precision - Recall - Fmeasure computations vs. Image Blurring','FontSize',20,'FontWeight','Bold')
saveGoodImg(H2,[dirPre,'../Documentation/Cosyne_2016/GaussianBlur_2.jpg'],sizeGoodIm)




%% Plot mean (using justMethod) & std (using relImPix) of different P/R/F measures
%  NOTE: THIS MAY BE SLIGHTLY DISINGENUOUS...
H3=figure; hold on
errorbar([1:numel(sig)+1]+0.0, justMethod.maxF_meanAccImgs, relImPix.maxF_stdAccImgs,'k.-','LineWidth',2)
errorbar([1:numel(sig)+1]+0.05, justMethod.maxF_old_meanAccImgs, relImPix.maxF_old_stdAccImgs,'r^--','LineWidth',2)
errorbar([1:numel(sig)+1]+0.1, justMethod.maxF_newMax_meanAccImgs, relImPix.maxF_newMax_stdAccImgs,'bs--','LineWidth',2)
errorbar([1:numel(sig)+1]+0.15, justMethod.maxF_newMean_meanAccImgs, relImPix.maxF_newMean_stdAccImgs,'gs--','LineWidth',2)
errorbar([1:numel(sig)+1]+0.2, justMethod.maxF_unionGT_meanAccImgs, relImPix.maxF_unionGT_stdAccImgs,'co--','LineWidth',2)
legend('bench','bench2','best gT','mean across gT','union gT')
set(gca,'XTick',1:numel(sig)+1,'XTickLabel',{'0',sig{:}},'FontSize',16,'FontWeight','Bold')
grid on
xlabel('\sigma Value for Gaussian Blurring','FontSize',18,'FontWeight','Bold')
ylabel('F-measure (across ~500 images)','FontSize',18,'FontWeight','Bold')
title('Comparing Different Precision - Recall - Fmeasure computations vs. Image Blurring','FontSize',20,'FontWeight','Bold')
saveGoodImg(H3,[dirPre,'../Documentation/Cosyne_2016/GaussianBlur_3.jpg'],sizeGoodIm)


keyboard






% Make up my own legend here...

figure(H);
hold on
%
text(0.87,1,{'\color{red}{\sigma=0.5}','\color{blue}{\sigma=1}','\color{yellow}{\sigma=1.5}','\color{cyan}{\sigma=2}','\color{black}{\sigma=2.5}','\color{magenta}{\sigma=5}','\color{green}{\sigma=10}'},'VerticalAlignment','top','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')

if(plot_max)
    % open vs. filled
    scatter(0.65,0.76,200,'ko','LineWidth',2)
    text(0.67,0.76,'= const THs','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
    %
    scatter(0.65,0.74,200,'ko','filled','LineWidth',2)
    text(0.67,0.74,'= vary THs','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
else
    scatter(0.65,0.70,200,'ko','LineWidth',2)
    text(0.67,0.70,'= max F','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
end

%
plot([0.63,0.67],[0.65,0.65],'g-','LineWidth',1.5)
text(0.67,0.65,' \color{green}{= Iso F}','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
%
plot([0.63,0.87],[0.60,0.60],'ko--','LineWidth',1.5)
text(0.67,0.60,'= ImPix','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')




if(plot_max)
    max_tag = 'best&max';
else
    max_tag = 'curve';
end


title(['Gaussian Blurred Image. ',max_tag],'FontSize',18,'FontWeight','Bold')

hold off


saveGoodImg(H,[dirPre,'../Documentation/Cosyne_2016/Overall_PR_curve_results_thinpbOFF/GaussianBlur_',max_tag,'_Kur_allParams.jpg'],sizeGoodIm)