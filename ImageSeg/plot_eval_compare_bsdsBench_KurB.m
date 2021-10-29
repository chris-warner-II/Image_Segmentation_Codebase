% This script/function written by CW 11/15 to plot multiple precision
% recall curves on the same isoF plot.  We can use this to compare
% different methods or different parameter settings for the same method.


[dirPre,sizeGoodIm] = onCluster;
addpath([dirPre,'images/BSDS_images/BSR/bench/benchmarks/'])


method = {'IsoDiff','AAnrm','GLnrm','Mod_SKHEuc','Mod_SKHAdj','Mod_N&G'}; % 'IsoDiff','AAnrm','GLnrm','Mod_SKHAdj','Mod_N&G';
methodStr =   {'ISO','AA','GL','TM 2D','TM 1D','M'};
legendMethods = {'Isotropic Diffusion', 'Normalized Average Association [Sarker]', 'Normalized Graph Laplacian [Shi]',...
       'Topographic Modularity [2D]','Topographic Modularity [1D]', 'Nontopographic Modularity [Newman]'};

   
   
method_colors = {'k','y','g','r','m','b'};
method_shapes = {'o','^','s','*','d','d'};


rM = {'1','3','5','10'};
rMcolors =  {'r','b','c','m'};

ks = {'sml','mid','lrg'}; 
shapes = {'v-','d-','^-'};

plot_max = 1; % set to 1 to plot average of Precision & Recall leading to largest F, regardless of threshold
              % set to 0 to plot average of Precision & Recall at each threshold value for each image
              
which_delF = {'better', 'average', 'worse'}; 
M = 4;	 % number of f-blur bins for image patches to mark and display when showing delF vs. F
N = 1; 	 % number of image patches within an f-blur bin to mark and display when showing delF vs. F
              
which_errbars = 'sem'; % either 'std' for standard deviation or 'sem' for standard error.    

switch which_errbars
    case 'sem'
        plot_name = 'SE';
    case'std'
        plot_name = '\sigma';
end


% Flags for plots below. Set to 1 to plot.
plot_F_errbars_params_optimize_each_method = 0;     %     
plot_FRP_errorbars_params_optimize_each_method = 0; % Ignore these for now. Come back to it.



              

blur_flg = 1; % 2 if DoG filter, 1 if blur gaussian filter.
%
if(blur_flg==2)
    blur_tag_M = '_blur_sigC1_S8_Kr0p01';
    blur_tit = 'w/ DoG filtering (\sigma_C=1, \sigma_S=8, K_r=0.01)';
    which_F_computation = {'meanGT','meanGT','meanGT','maxGT','maxGT','maxGT'};
    relative_to_what = {'justMethod','relRawPix','relGaussRF','justMethod','relRawPix','relGaussRF'};
elseif(blur_flg==1)
    blur_tag_M = '_blur_sig1';
    blur_tit = 'w/ GaussRF features (\sigma=1)';
    which_F_computation = {'meanGT','meanGT','meanGT','maxGT','maxGT','maxGT'};
    relative_to_what = {'justMethod','relRawPix','relGaussRF','justMethod','relRawPix','relGaussRF'};
else
    blur_tag_M = ''; % if we are not blurring.
    blur_tit = 'w/ RawPix features.';
    which_F_computation = {'meanGT','meanGT','maxGT','maxGT'};
    relative_to_what = {'justMethod','relRawPix','justMethod','relRawPix'};
end


% Hardcoding to do fewer options. Comment out to do all depending on blur option.
which_F_computation = {'meanGT','maxGT'};
relative_to_what = {'relGaussRF','relGaussRF'};



blur_tag_I = 'blur_sz13_sig1/';        % note: only using these to access pb images for visualization and benchmark results for justBlur and justDoG.
DoG_tag_I = 'blur_sigC1_S8_Kr0p01/';


% directory to ground truth mat files.
gTdir = [dirPre,'images/BSDS_patch/101x101_ds1/groundTruth/'];

% directory to image patches mat files.
imDir = [dirPre,'images/BSDS_patch/101x101_ds1/'];

% directory to probabalistic boundary png images for RawPix & Blurring
pbDirDoG = [dirPre,'images/BSDS_patch/101x101_ds1/',DoG_tag_I,'pb_png/'];
pbDirBlur = [dirPre,'images/BSDS_patch/101x101_ds1/',blur_tag_I,'pb_png/'];
pbDirIm = [dirPre,'images/BSDS_patch/101x101_ds1/pb_png/'];


num_cPdist = 4; % number of distances allowed in correspondPixelsB function
cPd_colors = {'r','g','b','k'};



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% load(['Kur_plot_results',blur_tag_M,'_',which_errbars])

% % matrix of parameter values for nice labeling later.
% param_matrix = cell(numel(rM),numel(ks));
% param_matrix_wd = cell(numel(rM),numel(ks),num_cPdist);
% for i = 1:numel(rM)
%     for j = 1:numel(ks)
%         param_matrix{i,j} = ['rM',rM{i},':ks',ks{j}];
%         for d = 1:num_cPdist
%             param_matrix_wd{i,j,d} = ['rM',rM{i},':ks',ks{j},':d',num2str(d)];
%         end
%     end   
% end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

try
    %noooo
     load(['Kur_plot_results',blur_tag_M])
     disp(['Loading in saved file called: Kur_plot_results',blur_tag_M])
catch

    %% Preallocate memory to hold all the data. Statistics for all methods of P, R and F
    %  absolute, relative to RawPix (RawPix) and GaussRF (GaussRF) Independent Sensors models.
    %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % (1). For just the method - the absolute value. Not comparing to any baseline or null model.
    %
    % Mean & STD across different ground truthers
    justMethod.maxF_meanGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    justMethod.maxF_meanGT_semAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    justMethod.maxF_meanGT_stdAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    justMethod.maxF_stdGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    %
    justMethod.maxR_meanGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    justMethod.maxR_meanGT_semAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    justMethod.maxR_meanGT_stdAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    justMethod.maxR_stdGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    %
    justMethod.maxP_meanGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    justMethod.maxP_meanGT_semAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    justMethod.maxP_meanGT_stdAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    justMethod.maxP_stdGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    %
    % Max value using best single ground truther
    justMethod.maxF_maxGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    justMethod.maxF_maxGT_semAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    justMethod.maxF_maxGT_stdAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    %
    justMethod.maxR_maxGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    justMethod.maxR_maxGT_semAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    justMethod.maxR_maxGT_stdAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    %
    justMethod.maxP_maxGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    justMethod.maxP_maxGT_semAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    justMethod.maxP_maxGT_stdAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % (2). justRawPix - mean & sem performance across 500 image patches using raw 
    %                   image pixel spatial gradients.
    %
    % Mean & STD across different ground truthers
    justRawPix.maxF_meanGT_meanAccImgs = zeros(1,num_cPdist);
    justRawPix.maxF_meanGT_semAccImgs = zeros(1,num_cPdist);
    justRawPix.maxF_meanGT_stdAccImgs = zeros(1,num_cPdist);
    justRawPix.maxF_stdGT_meanAccImgs = zeros(1,num_cPdist);
    %
    justRawPix.maxR_meanGT_meanAccImgs = zeros(1,num_cPdist);
    justRawPix.maxR_meanGT_semAccImgs = zeros(1,num_cPdist);
    justRawPix.maxR_meanGT_stdAccImgs = zeros(1,num_cPdist);
    justRawPix.maxR_stdGT_meanAccImgs = zeros(1,num_cPdist);
    %
    justRawPix.maxP_meanGT_meanAccImgs = zeros(1,num_cPdist);
    justRawPix.maxP_meanGT_semAccImgs = zeros(1,num_cPdist);
    justRawPix.maxP_meanGT_stdAccImgs = zeros(1,num_cPdist);
    justRawPix.maxP_stdGT_meanAccImgs = zeros(1,num_cPdist);
    %
    % Max value using best single ground truther.
    justRawPix.maxF_maxGT_meanAccImgs = zeros(1,num_cPdist);
    justRawPix.maxF_maxGT_semAccImgs = zeros(1,num_cPdist);
    justRawPix.maxF_maxGT_stdAccImgs = zeros(1,num_cPdist);
    %
    justRawPix.maxR_maxGT_meanAccImgs = zeros(1,num_cPdist);
    justRawPix.maxR_maxGT_semAccImgs = zeros(1,num_cPdist);
    justRawPix.maxR_maxGT_stdAccImgs = zeros(1,num_cPdist);
    %
    justRawPix.maxP_maxGT_meanAccImgs = zeros(1,num_cPdist);
    justRawPix.maxP_maxGT_semAccImgs = zeros(1,num_cPdist);
    justRawPix.maxP_maxGT_stdAccImgs = zeros(1,num_cPdist);
    %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % (3). relRawPix - values of P/R/F of method relative to just using Raw Image Pixels.
    %
    % Mean & STD across different ground truthers
    relRawPix.maxF_meanGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relRawPix.maxF_meanGT_semAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relRawPix.maxF_meanGT_stdAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relRawPix.maxF_stdGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    %
    relRawPix.maxR_meanGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relRawPix.maxR_meanGT_semAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relRawPix.maxR_meanGT_stdAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relRawPix.maxR_stdGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    %
    relRawPix.maxP_meanGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relRawPix.maxP_meanGT_semAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relRawPix.maxP_meanGT_stdAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relRawPix.maxP_stdGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    %
    % Max value using best single ground truther
    relRawPix.maxF_maxGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relRawPix.maxF_maxGT_semAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relRawPix.maxF_maxGT_stdAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    %
    relRawPix.maxR_maxGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relRawPix.maxR_maxGT_semAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relRawPix.maxR_maxGT_stdAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    %
    relRawPix.maxP_maxGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relRawPix.maxP_maxGT_semAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relRawPix.maxP_maxGT_stdAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % (4). justGaussRF - mean & sem performance across 500 image patches using blurred 
    %                   image pixel spatial gradients.
    %
    % Mean & STD across different ground truthers 
    justGaussRF.maxF_meanGT_meanAccImgs = zeros(1,num_cPdist);
    justGaussRF.maxF_meanGT_semAccImgs = zeros(1,num_cPdist);
    justGaussRF.maxF_meanGT_stdAccImgs = zeros(1,num_cPdist);
    justGaussRF.maxF_stdGT_meanAccImgs = zeros(1,num_cPdist);
    %
    justGaussRF.maxR_meanGT_meanAccImgs = zeros(1,num_cPdist);
    justGaussRF.maxR_meanGT_semAccImgs = zeros(1,num_cPdist);
    justGaussRF.maxR_meanGT_stdAccImgs = zeros(1,num_cPdist);
    justGaussRF.maxR_stdGT_meanAccImgs = zeros(1,num_cPdist);
    %
    justGaussRF.maxP_meanGT_meanAccImgs = zeros(1,num_cPdist);
    justGaussRF.maxP_meanGT_semAccImgs = zeros(1,num_cPdist);
    justGaussRF.maxP_meanGT_stdAccImgs = zeros(1,num_cPdist);
    justGaussRF.maxP_stdGT_meanAccImgs = zeros(1,num_cPdist);
    %
    % Max value using best single ground truther 
    justGaussRF.maxF_maxGT_meanAccImgs = zeros(1,num_cPdist);
    justGaussRF.maxF_maxGT_semAccImgs = zeros(1,num_cPdist);
    justGaussRF.maxF_maxGT_stdAccImgs = zeros(1,num_cPdist);
    %
    justGaussRF.maxR_maxGT_meanAccImgs = zeros(1,num_cPdist);
    justGaussRF.maxR_maxGT_semAccImgs = zeros(1,num_cPdist);
    justGaussRF.maxR_maxGT_stdAccImgs = zeros(1,num_cPdist);
    %
    justGaussRF.maxP_maxGT_meanAccImgs = zeros(1,num_cPdist);
    justGaussRF.maxP_maxGT_semAccImgs = zeros(1,num_cPdist);
    justGaussRF.maxP_maxGT_stdAccImgs = zeros(1,num_cPdist);
    %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % (5). relGaussRF - values of P/R/F of method relative to optimal Gaussian Blurring (sig=1).
    %
    % Mean & STD across different ground truthers
    relGaussRF.maxF_meanGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relGaussRF.maxF_meanGT_semAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relGaussRF.maxF_meanGT_stdAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relGaussRF.maxF_stdGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    %
    relGaussRF.maxR_meanGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relGaussRF.maxR_meanGT_semAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relGaussRF.maxR_meanGT_stdAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relGaussRF.maxR_stdGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    %
    relGaussRF.maxP_meanGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relGaussRF.maxP_meanGT_semAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relGaussRF.maxP_meanGT_stdAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relGaussRF.maxP_stdGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    %
    % Max value using best single ground truther
    relGaussRF.maxF_maxGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relGaussRF.maxF_maxGT_semAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relGaussRF.maxF_maxGT_stdAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    %
    relGaussRF.maxR_maxGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relGaussRF.maxR_maxGT_semAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relGaussRF.maxR_maxGT_stdAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    %
    relGaussRF.maxP_maxGT_meanAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relGaussRF.maxP_maxGT_semAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    relGaussRF.maxP_maxGT_stdAccImgs = zeros(numel(rM),numel(ks),num_cPdist,numel(method));
    %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % (6). RawPixVsGaussRF - mean & sem performance across 500 image patches directly 
    %                       comparing RawPix & GaussRF
    %
    % Mean & STD across different ground truthers
    RawPixVsGaussRF.maxF_meanGT_meanAccImgs = zeros(1,num_cPdist);
    RawPixVsGaussRF.maxF_meanGT_semAccImgs = zeros(1,num_cPdist);
    RawPixVsGaussRF.maxF_meanGT_stdAccImgs = zeros(1,num_cPdist);
    RawPixVsGaussRF.maxF_stdGT_meanAccImgs = zeros(1,num_cPdist);
    %
    RawPixVsGaussRF.maxR_meanGT_meanAccImgs = zeros(1,num_cPdist);
    RawPixVsGaussRF.maxR_meanGT_semAccImgs = zeros(1,num_cPdist);
    RawPixVsGaussRF.maxR_meanGT_stdAccImgs = zeros(1,num_cPdist);
    RawPixVsGaussRF.maxR_stdGT_meanAccImgs = zeros(1,num_cPdist);
    %
    RawPixVsGaussRF.maxP_meanGT_meanAccImgs = zeros(1,num_cPdist);
    RawPixVsGaussRF.maxP_meanGT_semAccImgs = zeros(1,num_cPdist);
    RawPixVsGaussRF.maxP_meanGT_stdAccImgs = zeros(1,num_cPdist);
    RawPixVsGaussRF.maxP_stdGT_meanAccImgs = zeros(1,num_cPdist);
    %
    % Max value using best single ground truther
    RawPixVsGaussRF.maxF_maxGT_meanAccImgs = zeros(1,num_cPdist);
    RawPixVsGaussRF.maxF_maxGT_semAccImgs = zeros(1,num_cPdist);
    RawPixVsGaussRF.maxF_maxGT_stdAccImgs = zeros(1,num_cPdist);
    %
    RawPixVsGaussRF.maxR_maxGT_meanAccImgs = zeros(1,num_cPdist);
    RawPixVsGaussRF.maxR_maxGT_semAccImgs = zeros(1,num_cPdist);
    RawPixVsGaussRF.maxR_maxGT_stdAccImgs = zeros(1,num_cPdist);
    %
    RawPixVsGaussRF.maxP_maxGT_meanAccImgs = zeros(1,num_cPdist);
    RawPixVsGaussRF.maxP_maxGT_semAccImgs = zeros(1,num_cPdist);
    RawPixVsGaussRF.maxP_maxGT_stdAccImgs = zeros(1,num_cPdist);
    %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % (7). Preallocate memory to store a mish-mash of other things.
    %
    % to store mean & std of each method with optimal
    % parameter settings for different ways of computing F-measure.
    maxF_meanAccImgs = zeros(num_cPdist,numel(which_F_computation),numel(method)); 
    maxR_meanAccImgs = zeros(num_cPdist,numel(which_F_computation),numel(method)); 
    maxP_meanAccImgs = zeros(num_cPdist,numel(which_F_computation),numel(method)); 
    %
    maxF_semAccImgs = zeros(num_cPdist,numel(which_F_computation),numel(method)); 
    maxR_semAccImgs = zeros(num_cPdist,numel(which_F_computation),numel(method)); 
    maxP_semAccImgs = zeros(num_cPdist,numel(which_F_computation),numel(method)); 
    %
    maxF_stdAccImgs = zeros(num_cPdist,numel(which_F_computation),numel(method)); 
    maxR_stdAccImgs = zeros(num_cPdist,numel(which_F_computation),numel(method)); 
    maxP_stdAccImgs = zeros(num_cPdist,numel(which_F_computation),numel(method)); 
    %
    % to store overlap between groundtruthers as a measure of how easy
    % segmentation was. Based on human consensus of what boundaries were.
    gT_agreement = cell(numel(method),numel(which_F_computation));
    gT_agreement2 = cell(numel(method),numel(which_F_computation));


    numFiles = zeros(numel(rM),numel(ks),numel(method));

    %% Loop through Different Methods & Parameter settings to calculate statistics across all available image patches on
    % Precision/Recall/F-measure calculated in a couple different ways and on those things relative to RawPix and Optimal GaussRF.


    % matrix of parameter values for nice labeling later.
    param_matrix = cell(numel(rM),numel(ks));
    param_matrix_wd = cell(numel(rM),numel(ks),num_cPdist);
    for i = 1:numel(rM)
        for j = 1:numel(ks)
            param_matrix{i,j} = ['rM',rM{i},':ks',ks{j}];
            for d = 1:num_cPdist
                param_matrix_wd{i,j,d} = ['rM',rM{i},':ks',ks{j},':d',num2str(d)];
            end
        end   
    end



    plot_IsoF_flag = 0; % no longer plotting this for now.



    for A = 1:numel(method)

        cntr = 0;

        % H=open('isoF.fig');

        method{A}

        for i = 1:numel(rM)      % loop over rM parameter
            for j = 1:numel(ks)     % loop over ks parameter
                for d = 1:num_cPdist   % loop over correspondPixels distance allowance.

                    cntr=cntr+1;

                    disp( ['rM: ',num2str(i),' / ',num2str(numel(rM)),', ks: ',num2str(j),' / ',num2str(numel(ks)),' , cP: ',num2str(d),' / ',num2str(num_cPdist)] )


                    if strcmp(method{A},'IsoDiff')
                        evalDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/Kur_PIF_Fourier1/',method{A},'/benchmark_results/rM',rM{i},'/NF_60_0/ks',ks{j},'/'];
                    else
                        evalDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/Kur_PIF_Fourier1/',method{A},'/benchmark_results/rM',rM{i},'/sDInf/sP0p2/NF_60_0/ks',ks{j},'/'];
                    end



                    % Loop through each image patch and grab maxF. So later I can compute their mean and std.
                    files = dir([evalDir,'*_d',num2str(d),'_ev2.txt']);

                    maxF_maxGT = zeros(1,numel(files));    % For these next 3, we compute P/R/F for each groundtruth associated with an
                    maxR_maxGT = zeros(1,numel(files));    % image patch.  Here we just grab the max value.  What is the algorithm's
                    maxP_maxGT = zeros(1,numel(files));    % best match to any single human?
                    maxF_bestGT  = zeros(1,numel(files));
                    thr_maxGT  = zeros(1,numel(files)); 

                    maxF_meanGT = zeros(1,numel(files));   % For a single image patch, the average P/R/F across the different groundtruths
                    maxR_meanGT = zeros(1,numel(files));
                    maxP_meanGT = zeros(1,numel(files));
                    thr_meanGT  = zeros(1,numel(files)); 
                    %
                    maxF_stdGT = zeros(1,numel(files));    % Standard Deviation of P/R/F across different groundtruths for single image ptch
                    maxR_stdGT = zeros(1,numel(files));
                    maxP_stdGT = zeros(1,numel(files));

                    numFiles(i,j,A) = numel(files);              % number of image patches processed for a given method or blur

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


                        maxF_maxGT(k) = max(Fmax);
                        ind = find(Fmax==max(Fmax));
                        ind = ind(1);
                        maxR_maxGT(k) = Rmax(ind);
                        maxP_maxGT(k) = Pmax(ind);
                        maxF_bestGT(k) = Fmax_whichGTs(ind);
                        thr_maxGT(k) = thr(ind);

                        maxF_meanGT(k) = max(Fmean);
                        ind = find(Fmean==max(Fmean));
                        ind = ind(1);
                        maxR_meanGT(k) = Rmean(ind);
                        maxP_meanGT(k) = Pmean(ind);
                        maxF_stdGT(k) = Fstd(ind);   % std across groundtruths as threshold that gave max avgF value.
                        thr_meanGT(k) = thr(ind);


                        imgPtchName{k} = files(k).name(1:end-8);


                        %k

                    end % looping through image patches



                    % Take Mean & Std across image patches of Precision,Recall,F-measure computed in different ways.
                    %
                    % (3). New way to compute P-R-F - (Mean & SEM across different ground truthers). 
                    justMethod.maxF_meanGT_meanAccImgs(i,j,d,A) = nanmean(maxF_meanGT);
                    justMethod.maxF_meanGT_semAccImgs(i,j,d,A) = sem(maxF_meanGT);
                    justMethod.maxF_meanGT_stdAccImgs(i,j,d,A) = nanstd(maxF_meanGT);
                    justMethod.maxF_stdGT_meanAccImgs(i,j,d,A) = nanmean(maxF_stdGT);
                    %
                    justMethod.maxR_meanGT_meanAccImgs(i,j,d,A) = nanmean(maxR_meanGT);
                    justMethod.maxR_meanGT_semAccImgs(i,j,d,A) = sem(maxR_meanGT);
                    justMethod.maxR_meanGT_stdAccImgs(i,j,d,A) = nanstd(maxR_meanGT);
                    justMethod.maxR_stdGT_meanAccImgs(i,j,d,A) = nanmean(maxR_stdGT);
                    %
                    justMethod.maxP_meanGT_meanAccImgs(i,j,d,A) = nanmean(maxP_meanGT);
                    justMethod.maxP_meanGT_semAccImgs(i,j,d,A) = sem(maxP_meanGT);
                    justMethod.maxP_meanGT_stdAccImgs(i,j,d,A) = nanstd(maxP_meanGT);
                    justMethod.maxP_stdGT_meanAccImgs(i,j,d,A) = nanmean(maxP_stdGT);
                    %
                    % (4). New way to compute P-R-F - (Max value using best single ground truther.). 
                    justMethod.maxF_maxGT_meanAccImgs(i,j,d,A) = nanmean(maxF_maxGT);
                    justMethod.maxF_maxGT_semAccImgs(i,j,d,A) =sem(maxF_maxGT);
                    justMethod.maxF_maxGT_stdAccImgs(i,j,d,A) = nanstd(maxF_maxGT);
                    %
                    justMethod.maxR_maxGT_meanAccImgs(i,j,d,A) = nanmean(maxR_maxGT);
                    justMethod.maxR_maxGT_semAccImgs(i,j,d,A) = sem(maxR_maxGT);
                    justMethod.maxR_maxGT_stdAccImgs(i,j,d,A) = nanstd(maxR_maxGT);
                    %
                    justMethod.maxP_maxGT_meanAccImgs(i,j,d,A) = nanmean(maxP_maxGT);
                    justMethod.maxP_maxGT_semAccImgs(i,j,d,A) = sem(maxP_maxGT);
                    justMethod.maxP_maxGT_stdAccImgs(i,j,d,A) = nanstd(maxP_maxGT);



                    %% Compute different Precision/Recall values for Raw Image Pixels
                    % Currently, doing this inside method & params for loops because I dont know how else to line up RawPix & GaussRF to have 
                    % statistics relative to those.
                    RawPixDir = [dirPre,'images/BSDS_patch/101x101_ds1/benchmark_results/'];

                    maxF_maxGT_RawPix = zeros(1,numel(files));    % For these next 3, we compute P/R/F for each groundtruth associated with an
                    maxR_maxGT_RawPix = zeros(1,numel(files));    % image patch.  Here we just grab the max value.  What is the algorithm's
                    maxP_maxGT_RawPix = zeros(1,numel(files));    % best match to any single human?
                    maxF_bestGT_RawPix  =  zeros(1,numel(files));
                    thr_maxGT_RawPix  = zeros(1,numel(files));

                    maxF_meanGT_RawPix = zeros(1,numel(files));   % For a single image patch, the average P/R/F across the different groundtruths
                    maxR_meanGT_RawPix = zeros(1,numel(files));
                    maxP_meanGT_RawPix = zeros(1,numel(files));
                    thr_meanGT_RawPix  = zeros(1,numel(files));
                    %
                    maxF_stdGT_RawPix = zeros(1,numel(files));    % Standard Deviation of P/R/F across different groundtruths for single image ptch
                    maxR_stdGT_RawPix = zeros(1,numel(files));
                    maxP_stdGT_RawPix = zeros(1,numel(files));

                    %numFiles(i,j,A) = numel(files);              % number of image patches processed for a given method or blur

                    for k = 1:numel(files)

                        filename = fullfile(RawPixDir,files(k).name);
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


                        maxF_maxGT_RawPix(k) = max(Fmax);
                        ind = find(Fmax==max(Fmax));
                        ind = ind(1);
                        maxR_maxGT_RawPix(k) = Rmax(ind);
                        maxP_maxGT_RawPix(k) = Pmax(ind);
                        maxF_bestGT_RawPix(k) = Fmax_whichGTs(ind);
                        thr_maxGT_RawPix(k) = thr(ind);

                        maxF_meanGT_RawPix(k) = max(Fmean);
                        ind = find(Fmean==max(Fmean));
                        ind = ind(1);
                        maxR_meanGT_RawPix(k) = Rmean(ind);
                        maxP_meanGT_RawPix(k) = Pmean(ind);
                        maxF_stdGT_RawPix(k) = Fstd(ind);   % std across groundtruths as threshold that gave max avgF value.
                        thr_meanGT_RawPix(k) = thr(ind);

                        %k

                    end % looping through image patches



                    % justRawPix - mean & sem performance across 500 image patches using just image pixel spatial gradients.
                    %
                    % (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
                    if(A==1 && i==i && j==1)
                        justRawPix.maxF_meanGT_meanAccImgs(d) = nanmean(maxF_meanGT_RawPix);
                        justRawPix.maxF_meanGT_semAccImgs(d) = sem(maxF_meanGT_RawPix);
                        justRawPix.maxF_meanGT_stdAccImgs(d) = nanstd(maxF_meanGT_RawPix);
                        justRawPix.maxF_stdGT_meanAccImgs(d) = nanmean(maxF_stdGT_RawPix);
                        %
                        justRawPix.maxR_meanGT_meanAccImgs(d) = nanmean(maxR_meanGT_RawPix);
                        justRawPix.maxR_meanGT_semAccImgs(d) = sem(maxR_meanGT_RawPix);
                        justRawPix.maxR_meanGT_stdAccImgs(d) = nanstd(maxR_meanGT_RawPix);
                        justRawPix.maxR_stdGT_meanAccImgs(d) = nanmean(maxR_stdGT_RawPix);
                        %
                        justRawPix.maxP_meanGT_meanAccImgs(d) = nanmean(maxP_meanGT_RawPix);
                        justRawPix.maxP_meanGT_semAccImgs(d) = sem(maxP_meanGT_RawPix);
                        justRawPix.maxP_meanGT_stdAccImgs(d) = nanstd(maxP_meanGT_RawPix);
                        justRawPix.maxP_stdGT_meanAccImgs(d) = nanmean(maxP_stdGT_RawPix);
                        %
                        % (4). New way to compute P-R-F - (Max value using best single ground truther.). 
                        justRawPix.maxF_maxGT_meanAccImgs(d) = nanmean(maxF_maxGT_RawPix);
                        justRawPix.maxF_maxGT_semAccImgs(d) = sem(maxF_maxGT_RawPix);
                        justRawPix.maxF_maxGT_stdAccImgs(d) = nanstd(maxF_maxGT_RawPix);
                        %
                        justRawPix.maxR_maxGT_meanAccImgs(d) = nanmean(maxR_maxGT_RawPix);
                        justRawPix.maxR_maxGT_semAccImgs(d) = sem(maxR_maxGT_RawPix);
                        justRawPix.maxR_maxGT_stdAccImgs(d) = nanstd(maxR_maxGT_RawPix);
                        %
                        justRawPix.maxP_maxGT_meanAccImgs(d) = nanmean(maxP_maxGT_RawPix);
                        justRawPix.maxP_maxGT_semAccImgs(d) = sem(maxP_maxGT_RawPix);
                        justRawPix.maxP_maxGT_stdAccImgs(d) = nanstd(maxP_maxGT_RawPix);
                    end




                    % Compute statistics of this method & parameters relative to RawPix.
                    %
                    % (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
                    relRawPix.maxF_meanGT_meanAccImgs(i,j,d,A) = nanmean(maxF_meanGT-maxF_meanGT_RawPix);
                    relRawPix.maxF_meanGT_semAccImgs(i,j,d,A) = sem(maxF_meanGT-maxF_meanGT_RawPix);
                    relRawPix.maxF_meanGT_semAccImgs(i,j,d,A) = sem(maxF_meanGT-maxF_meanGT_RawPix);
                    relRawPix.maxF_stdGT_meanAccImgs(i,j,d,A) = nanmean(maxF_stdGT-maxF_stdGT_RawPix);
                    %
                    relRawPix.maxR_meanGT_meanAccImgs(i,j,d,A) = nanmean(maxR_meanGT-maxF_meanGT_RawPix);
                    relRawPix.maxR_meanGT_semAccImgs(i,j,d,A) = sem(maxR_meanGT-maxR_meanGT_RawPix);
                    relRawPix.maxR_meanGT_stdAccImgs(i,j,d,A) = nanstd(maxR_meanGT-maxR_meanGT_RawPix);
                    relRawPix.maxR_stdGT_meanAccImgs(i,j,d,A) = nanmean(maxR_stdGT-maxR_stdGT_RawPix);
                    %
                    relRawPix.maxP_meanGT_meanAccImgs(i,j,d,A) = nanmean(maxP_meanGT-maxP_meanGT_RawPix);
                    relRawPix.maxP_meanGT_semAccImgs(i,j,d,A) = sem(maxP_meanGT-maxP_meanGT_RawPix);
                    relRawPix.maxP_meanGT_stdAccImgs(i,j,d,A) = nanstd(maxP_meanGT-maxP_meanGT_RawPix);
                    relRawPix.maxP_stdGT_meanAccImgs(i,j,d,A) = nanmean(maxP_stdGT-maxP_stdGT_RawPix);
                    %
                    % (4). New way to compute P-R-F - (Max value using best single ground truther.). 
                    relRawPix.maxF_maxGT_meanAccImgs(i,j,d,A) = nanmean(maxF_maxGT-maxF_maxGT_RawPix);
                    relRawPix.maxF_maxGT_semAccImgs(i,j,d,A) = sem(maxF_maxGT-maxF_maxGT_RawPix);
                    relRawPix.maxF_maxGT_stdAccImgs(i,j,d,A) = nanstd(maxF_maxGT-maxF_maxGT_RawPix);
                    %
                    relRawPix.maxR_maxGT_meanAccImgs(i,j,d,A) = nanmean(maxR_maxGT-maxR_maxGT_RawPix);
                    relRawPix.maxR_maxGT_semAccImgs(i,j,d,A) = sem(maxR_maxGT-maxR_maxGT_RawPix);
                    relRawPix.maxR_maxGT_stdAccImgs(i,j,d,A) = nanstd(maxR_maxGT-maxR_maxGT_RawPix);
                    %
                    relRawPix.maxP_maxGT_meanAccImgs(i,j,d,A) = mean(maxP_maxGT-maxP_maxGT_RawPix);
                    relRawPix.maxP_maxGT_semAccImgs(i,j,d,A) = sem(maxP_maxGT-maxP_maxGT_RawPix);
                    relRawPix.maxP_maxGT_stdAccImgs(i,j,d,A) = nanstd(maxP_maxGT-maxP_maxGT_RawPix);



                    %% Compute different Precision/Recall values for Optimal Blurring
                    %if(blur_flg)

                    GaussRFDir = [dirPre,'images/BSDS_patch/101x101_ds1/',blur_tag_I,'benchmark_results/'];

                    maxF_maxGT_GaussRF = zeros(1,numel(files));    % For these next 3, we compute P/R/F for each groundtruth associated with an
                    maxR_maxGT_GaussRF = zeros(1,numel(files));    % image patch.  Here we just grab the max value.  What is the algorithm's
                    maxP_maxGT_GaussRF = zeros(1,numel(files));    % best match to any single human?
                    maxF_bestGT_GaussRF  =  zeros(1,numel(files));
                    thr_maxGT_GaussRF  = zeros(1,numel(files)); 

                    maxF_meanGT_GaussRF = zeros(1,numel(files));   % For a single image patch, the average P/R/F across the different groundtruths
                    maxR_meanGT_GaussRF = zeros(1,numel(files));
                    maxP_meanGT_GaussRF = zeros(1,numel(files));
                    thr_meanGT_GaussRF  = zeros(1,numel(files)); 

                    maxF_stdGT_GaussRF = zeros(1,numel(files));    % Standard Deviation of P/R/F across different groundtruths for single image ptch
                    maxR_stdGT_GaussRF = zeros(1,numel(files));
                    maxP_stdGT_GaussRF = zeros(1,numel(files));

                    %numFiles(i,j,A) = numel(files);              % number of image patches processed for a given method or blur

                    for k = 1:numel(files)

                        filename = fullfile(GaussRFDir,files(k).name);
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

                        maxF_maxGT_GaussRF(k) = max(Fmax);
                        ind = find(Fmax==max(Fmax));
                        ind = ind(1);
                        maxR_maxGT_GaussRF(k) = Rmax(ind);
                        maxP_maxGT_GaussRF(k) = Pmax(ind);
                        maxF_bestGT_GaussRF(k) = Fmax_whichGTs(ind);
                        thr_maxGT_GaussRF(k) = thr(ind);

                        maxF_meanGT_GaussRF(k) = max(Fmean);
                        ind = find(Fmean==max(Fmean));
                        ind = ind(1);
                        maxR_meanGT_GaussRF(k) = Rmean(ind);
                        maxP_meanGT_GaussRF(k) = Pmean(ind);
                        maxF_stdGT_GaussRF(k) = Fstd(ind);   % std across groundtruths as threshold that gave max avgF value.
                        thr_meanGT_GaussRF(k) = thr(ind);

                        %k

                    end % looping through image patches


                    if(A==1 && i==i && j==1)
                        % justGaussRF - mean & sem performance across 500 image patches using blurred image pixel spatial gradients.
                        %
                        % (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
                        justGaussRF.maxF_meanGT_meanAccImgs(d) = nanmean(maxF_meanGT_GaussRF);
                        justGaussRF.maxF_meanGT_semAccImgs(d) = sem(maxF_meanGT_GaussRF);
                        justGaussRF.maxF_meanGT_stdAccImgs(d) = nanstd(maxF_meanGT_GaussRF);
                        justGaussRF.maxF_stdGT_meanAccImgs(d) = nanmean(maxF_stdGT_GaussRF);
                        %
                        justGaussRF.maxR_meanGT_meanAccImgs(d) = nanmean(maxR_meanGT_GaussRF);
                        justGaussRF.maxR_meanGT_semAccImgs(d) = sem(maxR_meanGT_GaussRF);
                        justGaussRF.maxR_meanGT_stdAccImgs(d) = nanstd(maxR_meanGT_GaussRF);
                        justGaussRF.maxR_stdGT_meanAccImgs(d) = nanmean(maxR_stdGT_GaussRF);
                        %
                        justGaussRF.maxP_meanGT_meanAccImgs(d) = nanmean(maxP_meanGT_GaussRF);
                        justGaussRF.maxP_meanGT_semAccImgs(d) = sem(maxP_meanGT_GaussRF);
                        justGaussRF.maxP_meanGT_stdAccImgs(d) = nanstd(maxP_meanGT_GaussRF);
                        justGaussRF.maxP_stdGT_meanAccImgs(d) = nanmean(maxP_stdGT_GaussRF);
                        %
                        % (4). New way to compute P-R-F - (Max value using best single ground truther.). 
                        justGaussRF.maxF_maxGT_meanAccImgs(d) = nanmean(maxF_maxGT_GaussRF);
                        justGaussRF.maxF_maxGT_semAccImgs(d) = sem(maxF_maxGT_GaussRF);
                        justGaussRF.maxF_maxGT_stdAccImgs(d) = nanstd(maxF_maxGT_GaussRF);
                        %
                        justGaussRF.maxR_maxGT_meanAccImgs(d) = nanmean(maxR_maxGT_GaussRF);
                        justGaussRF.maxR_maxGT_semAccImgs(d) = sem(maxR_maxGT_GaussRF);
                        justGaussRF.maxR_maxGT_stdAccImgs(d) = nanstd(maxR_maxGT_GaussRF);
                        %
                        justGaussRF.maxP_maxGT_meanAccImgs(d) = nanmean(maxP_maxGT_GaussRF);
                        justGaussRF.maxP_maxGT_semAccImgs(d) = sem(maxP_maxGT_GaussRF);
                        justGaussRF.maxP_maxGT_stdAccImgs(d) = nanstd(maxP_maxGT_GaussRF);
                    end




                    %RawPixVsGaussRF - keeps track of the direct comparisonbetween RawPix & GaussRF.  
                    % This is GaussRF relRawPix and the negative of RawPix relGaussRF.
                    if(A==1 && i==i && j==1)
                        % (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
                        RawPixVsGaussRF.maxF_meanGT_meanAccImgs(d) = nanmean(maxF_meanGT_GaussRF-maxF_meanGT_RawPix);
                        RawPixVsGaussRF.maxF_meanGT_semAccImgs(d) = sem(maxF_meanGT_GaussRF-maxF_meanGT_RawPix);
                        RawPixVsGaussRF.maxF_meanGT_stdAccImgs(d) = nanstd(maxF_meanGT_GaussRF-maxF_meanGT_RawPix);
                        RawPixVsGaussRF.maxF_stdGT_meanAccImgs(d) = nanmean(maxF_stdGT_GaussRF-maxF_stdGT_RawPix);
                        %
                        RawPixVsGaussRF.maxR_meanGT_meanAccImgs(d) = nanmean(maxR_meanGT_GaussRF-maxR_meanGT_RawPix);
                        RawPixVsGaussRF.maxR_meanGT_semAccImgs(d) = sem(maxR_meanGT_GaussRF-maxR_meanGT_RawPix);
                        RawPixVsGaussRF.maxR_meanGT_stdAccImgs(d) = nanstd(maxR_meanGT_GaussRF-maxR_meanGT_RawPix);
                        RawPixVsGaussRF.maxR_stdGT_meanAccImgs(d) = nanmean(maxR_stdGT_GaussRF-maxR_stdGT_RawPix);
                        %
                        RawPixVsGaussRF.maxP_meanGT_meanAccImgs(d) = nanmean(maxP_meanGT_GaussRF-maxP_meanGT_RawPix);
                        RawPixVsGaussRF.maxP_meanGT_semAccImgs(d) = sem(maxP_meanGT_GaussRF-maxP_meanGT_RawPix);
                        RawPixVsGaussRF.maxP_meanGT_stdAccImgs(d) = nanstd(maxP_meanGT_GaussRF-maxP_meanGT_RawPix);
                        RawPixVsGaussRF.maxP_stdGT_meanAccImgs(d) = nanmean(maxP_stdGT_GaussRF-maxP_stdGT_RawPix);
                        %
                        % (4). New way to compute P-R-F - (Max value using best single ground truther.). 
                        RawPixVsGaussRF.maxF_maxGT_meanAccImgs(d) = nanmean(maxF_maxGT_GaussRF-maxF_maxGT_RawPix);
                        RawPixVsGaussRF.maxF_maxGT_semAccImgs(d) = sem(maxF_maxGT_GaussRF-maxF_maxGT_RawPix);
                        RawPixVsGaussRF.maxF_maxGT_stdAccImgs(d) = nanstd(maxF_maxGT_GaussRF-maxF_maxGT_RawPix);
                        %
                        RawPixVsGaussRF.maxR_maxGT_meanAccImgs(d) = nanmean(maxR_maxGT_GaussRF-maxR_maxGT_RawPix);
                        RawPixVsGaussRF.maxR_maxGT_semAccImgs(d) = sem(maxR_maxGT_GaussRF-maxR_maxGT_RawPix);
                        RawPixVsGaussRF.maxR_maxGT_stdAccImgs(d) = nanstd(maxR_maxGT_GaussRF-maxR_maxGT_RawPix);
                        %
                        RawPixVsGaussRF.maxP_maxGT_meanAccImgs(d) = nanmean(maxP_maxGT_GaussRF-maxP_maxGT_RawPix);
                        RawPixVsGaussRF.maxP_maxGT_semAccImgs(d) = sem(maxP_maxGT_GaussRF-maxP_maxGT_RawPix);
                        RawPixVsGaussRF.maxP_maxGT_stdAccImgs(d) = nanstd(maxP_maxGT_GaussRF-maxP_maxGT_RawPix);
                    end




                    % Compute statistics of this method & parameters relative to GaussRF.
                    %
                    % (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
                    relGaussRF.maxF_meanGT_meanAccImgs(i,j,d,A) = nanmean(maxF_meanGT-maxF_meanGT_GaussRF);
                    relGaussRF.maxF_meanGT_semAccImgs(i,j,d,A) = sem(maxF_meanGT-maxF_meanGT_GaussRF);
                    relGaussRF.maxF_meanGT_stdAccImgs(i,j,d,A) = nanstd(maxF_meanGT-maxF_meanGT_GaussRF);
                    relGaussRF.maxF_stdGT_meanAccImgs(i,j,d,A) = nanmean(maxF_stdGT-maxF_stdGT_GaussRF);
                    %
                    relGaussRF.maxR_meanGT_meanAccImgs(i,j,d,A) = nanmean(maxR_meanGT-maxR_meanGT_GaussRF);
                    relGaussRF.maxR_meanGT_semAccImgs(i,j,d,A) = sem(maxR_meanGT-maxR_meanGT_GaussRF);
                    relGaussRF.maxR_meanGT_stdAccImgs(i,j,d,A) = nanstd(maxR_meanGT-maxR_meanGT_GaussRF);
                    relGaussRF.maxR_stdGT_meanAccImgs(i,j,d,A) = nanmean(maxR_stdGT-maxR_stdGT_GaussRF);
                    %
                    relGaussRF.maxP_meanGT_meanAccImgs(i,j,d,A) = nanmean(maxP_meanGT-maxP_meanGT_GaussRF);
                    relGaussRF.maxP_meanGT_semAccImgs(i,j,d,A) = sem(maxP_meanGT-maxP_meanGT_GaussRF);
                    relGaussRF.maxP_meanGT_stdAccImgs(i,j,d,A) = nanstd(maxP_meanGT-maxP_meanGT_GaussRF);
                    relGaussRF.maxP_stdGT_meanAccImgs(i,j,d,A) = nanmean(maxP_stdGT-maxP_stdGT_GaussRF);
                    %
                    % (4). New way to compute P-R-F - (Max value using best single ground truther.). 
                    relGaussRF.maxF_maxGT_meanAccImgs(i,j,d,A) = nanmean(maxF_maxGT-maxF_maxGT_GaussRF);
                    relGaussRF.maxF_maxGT_semAccImgs(i,j,d,A) = sem(maxF_maxGT-maxF_maxGT_GaussRF);
                    relGaussRF.maxF_maxGT_stdAccImgs(i,j,d,A) = nanstd(maxF_maxGT-maxF_maxGT_GaussRF);
                    %
                    relGaussRF.maxR_maxGT_meanAccImgs(i,j,d,A) = nanmean(maxR_maxGT-maxR_maxGT_GaussRF);
                    relGaussRF.maxR_maxGT_semAccImgs(i,j,d,A) = sem(maxR_maxGT-maxR_maxGT_GaussRF);
                    relGaussRF.maxR_maxGT_stdAccImgs(i,j,d,A) = nanstd(maxR_maxGT-maxR_maxGT_GaussRF);
                    %
                    relGaussRF.maxP_maxGT_meanAccImgs(i,j,d,A) = nanmean(maxP_maxGT-maxP_maxGT_GaussRF);
                    relGaussRF.maxP_maxGT_semAccImgs(i,j,d,A) = sem(maxP_maxGT-maxP_maxGT_GaussRF);
                    relGaussRF.maxP_maxGT_stdAccImgs(i,j,d,A) = nanstd(maxP_maxGT-maxP_maxGT_GaussRF);


                    % Create a data structure which holds F-measure values for all images
                    % (not just mean and std) so we can make Violin plots.
                    maxF_meanGT_struct_method_only{i,j,d,A} = maxF_meanGT';
                    maxR_meanGT_struct_method_only{i,j,d,A} = maxR_meanGT';
                    maxP_meanGT_struct_method_only{i,j,d,A} = maxP_meanGT';
                    thr_meanGT_struct_method_only{i,j,d,A}      = thr_meanGT';
                    %
                    maxF_maxGT_struct_method_only{i,j,d,A}  = maxF_maxGT';
                    maxR_maxGT_struct_method_only{i,j,d,A}  = maxR_maxGT';
                    maxP_maxGT_struct_method_only{i,j,d,A}  = maxP_maxGT';
                    thr_maxGT_struct_method_only{i,j,d,A}       = thr_maxGT';
                    bestGT_maxGT_struct_method_only{i,j,d,A} = maxF_bestGT';
                    %
                    % % %
                    %
                    maxF_meanGT_struct_GaussRF_only{i,j,d,A} = maxF_meanGT_GaussRF';
                    maxR_meanGT_struct_GaussRF_only{i,j,d,A} = maxR_meanGT_GaussRF';
                    maxP_meanGT_struct_GaussRF_only{i,j,d,A} = maxP_meanGT_GaussRF';
                    thr_meanGT_struct_GaussRF_only{i,j,d,A}     = thr_meanGT_GaussRF';
                    %
                    maxF_maxGT_struct_GaussRF_only{i,j,d,A} = maxF_maxGT_GaussRF';
                    maxR_maxGT_struct_GaussRF_only{i,j,d,A} = maxR_maxGT_GaussRF';
                    maxP_maxGT_struct_GaussRF_only{i,j,d,A} = maxP_maxGT_GaussRF';
                    thr_maxGT_struct_GaussRF_only{i,j,d,A}      = thr_maxGT_GaussRF';
                    bestGT_maxGT_struct_GaussRF_only{i,j,d,A} = maxF_bestGT_GaussRF';
                    %
                    % % %
                    %
                    maxF_meanGT_struct_RawPix_only{i,j,d,A} = maxF_meanGT_RawPix';
                    maxR_meanGT_struct_RawPix_only{i,j,d,A} = maxR_meanGT_RawPix';
                    maxP_meanGT_struct_RawPix_only{i,j,d,A} = maxP_meanGT_RawPix';
                    thr_meanGT_struct_RawPix_only{i,j,d,A}  = thr_meanGT_RawPix';
                    %
                    maxF_maxGT_struct_RawPix_only{i,j,d,A} = maxF_maxGT_RawPix';
                    maxR_maxGT_struct_RawPix_only{i,j,d,A} = maxR_maxGT_RawPix';
                    maxP_maxGT_struct_RawPix_only{i,j,d,A} = maxP_maxGT_RawPix';
                    thr_maxGT_struct_RawPix_only{i,j,d,A}  = thr_maxGT_RawPix';
                    bestGT_maxGT_struct_RawPix_only{i,j,d,A} = maxF_bestGT_RawPix';

                    % Also save name of image patch so I can go back and visualize them next to the performance.
                    img_ptch_name_struct{i,j,A} = imgPtchName;

                    % Do same for DoG as we did with blur.
                    imDoGDir = [dirPre,'images/BSDS_patch/101x101_ds1/',DoG_tag_I,'benchmark_results/'];


                end % end looping over correspondPixels distance
            end % loop over KS parameter (for single method)
        end % loop over rM parameter (for single method)

        % Free up some memory space
        clear maxF_bestGT  maxF_bestGT_GaussRF maxF_bestGT_RawPix ...
              maxF_maxGT maxF_maxGT_GaussRF maxF_maxGT_RawPix ...
              maxF_meanGT maxF_meanGT_GaussRF maxF_meanGT_RawPix ...
              maxF_stdGT maxF_stdGT_GaussRF maxF_stdGT_RawPix ...
              maxP_maxGT maxP_maxGT_GaussRF maxP_maxGT_RawPix ...
              maxP_meanGT maxP_meanGT_GaussRF maxP_meanGT_RawPix ...
              maxP_stdGT maxP_stdGT_GaussRF maxP_stdGT_RawPix ...
              maxR_maxGT maxR_maxGT_GaussRF maxR_maxGT_RawPix ... 
              maxR_meanGT maxR_meanGT_GaussRF maxR_meanGT_RawPix ...
              maxR_stdGT maxR_stdGT_GaussRF maxR_stdGT_RawPix

    end   





    % Now get best parameter settings for each of the methods A.
    %
    % Params are {rM, ks, d}. 
    % The output structures are of dimension size: num_cPdist x 6 x numel(method).
    %
    % The 8 different output structures are max{P,R,F}_{mean,sem}AccImgs and 
    % rM_max, ks_max.
    %
    % The 6 values in the 2nd dim of the output structures go thru.
    % 1 = justMethod,	meanGT
    % 2 = relRawPix,     meanGT 
    % 3 = relGaussRF,    meanGT
    % 4 = justMethod,   maxGT
    % 5 = relRawPix,     maxGT
    % 6 = relGaussRF,    maxGT
    %
    % Should be lots of redundancy (repeated rM,ks param values) in this I hope.

    for A = 1:numel(method)
        for d = 1:num_cPdist
            %
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            % mean (across GT's) of P/R/F : Just Method
            h = 1;                                                            %(rM,ks,cPd,Meth)
            maxF_meanAccImgs(d,h,A) = max(max(justMethod.maxF_meanGT_meanAccImgs(:,:,d,A)));                                    
            [rM_max(d,h,A),ks_max(d,h,A)] = find( justMethod.maxF_meanGT_meanAccImgs(:,:,d,A) == maxF_meanAccImgs(d,h,A) );
            maxR_meanAccImgs(d,h,A) = justMethod.maxR_meanGT_meanAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxP_meanAccImgs(d,h,A) = justMethod.maxP_meanGT_meanAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            %
            maxF_semAccImgs(d,h,A) = justMethod.maxF_meanGT_semAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxR_semAccImgs(d,h,A) = justMethod.maxR_meanGT_semAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxP_semAccImgs(d,h,A) = justMethod.maxP_meanGT_semAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            %
            maxF_stdAccImgs(d,h,A) = justMethod.maxF_meanGT_stdAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxR_stdAccImgs(d,h,A) = justMethod.maxR_meanGT_stdAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxP_stdAccImgs(d,h,A) = justMethod.maxP_meanGT_stdAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            %
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            % mean (across GT's) of P/R/F : Relative to RawPix
            h = 2;
            maxF_meanAccImgs(d,h,A) = max(max(relRawPix.maxF_meanGT_meanAccImgs(:,:,d,A)));                                    
            [rM_max(d,h,A),ks_max(d,h,A)] = find( relRawPix.maxF_meanGT_meanAccImgs(:,:,d,A) == maxF_meanAccImgs(d,h,A) );
            maxR_meanAccImgs(d,h,A) = relRawPix.maxR_meanGT_meanAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxP_meanAccImgs(d,h,A) = relRawPix.maxP_meanGT_meanAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            %
            maxF_semAccImgs(d,h,A) = relRawPix.maxF_meanGT_semAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxR_semAccImgs(d,h,A) = relRawPix.maxR_meanGT_semAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxP_semAccImgs(d,h,A) = relRawPix.maxP_meanGT_semAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            %
            maxF_stdAccImgs(d,h,A) = relRawPix.maxF_meanGT_stdAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxR_stdAccImgs(d,h,A) = relRawPix.maxR_meanGT_stdAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxP_stdAccImgs(d,h,A) = relRawPix.maxP_meanGT_stdAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            %
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            % mean (across GT's) of P/R/F : Relative to GaussRF
            h = 3;
            maxF_meanAccImgs(d,h,A) = max(max(relGaussRF.maxF_meanGT_meanAccImgs(:,:,d,A)));                                    
            [rM_max(d,h,A),ks_max(d,h,A)] = find( relGaussRF.maxF_meanGT_meanAccImgs(:,:,d,A) == maxF_meanAccImgs(d,h,A) );
            maxR_meanAccImgs(d,h,A) = relGaussRF.maxR_meanGT_meanAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxP_meanAccImgs(d,h,A) = relGaussRF.maxP_meanGT_meanAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            %
            maxF_semAccImgs(d,h,A) = relGaussRF.maxF_meanGT_semAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxR_semAccImgs(d,h,A) = relGaussRF.maxR_meanGT_semAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxP_semAccImgs(d,h,A) = relGaussRF.maxP_meanGT_semAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            %
            maxF_stdAccImgs(d,h,A) = relGaussRF.maxF_meanGT_stdAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxR_stdAccImgs(d,h,A) = relGaussRF.maxR_meanGT_stdAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxP_stdAccImgs(d,h,A) = relGaussRF.maxP_meanGT_stdAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            %
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            % max (across GT's) of P/R/F : Just Method.
            h = 4;
            maxF_meanAccImgs(d,h,A) = max(max(justMethod.maxF_maxGT_meanAccImgs(:,:,d,A)));                                    
            [rM_max(d,h,A),ks_max(d,h,A)] = find( justMethod.maxF_maxGT_meanAccImgs(:,:,d,A) == maxF_meanAccImgs(d,h,A) );
            maxR_meanAccImgs(d,h,A) = justMethod.maxR_maxGT_meanAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxP_meanAccImgs(d,h,A) = justMethod.maxP_maxGT_meanAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            %
            maxF_semAccImgs(d,h,A) = justMethod.maxF_maxGT_semAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxR_semAccImgs(d,h,A) = justMethod.maxR_maxGT_semAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxP_semAccImgs(d,h,A) = justMethod.maxP_maxGT_semAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            %
            maxF_stdAccImgs(d,h,A) = justMethod.maxF_maxGT_stdAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxR_stdAccImgs(d,h,A) = justMethod.maxR_maxGT_stdAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxP_stdAccImgs(d,h,A) = justMethod.maxP_maxGT_stdAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            %
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            % max (across GT's) of P/R/F : Relative to RawPix.
            h = 5;
            maxF_meanAccImgs(d,h,A) = max(max(relRawPix.maxF_maxGT_meanAccImgs(:,:,d,A)));                                    
            [rM_max(d,h,A),ks_max(d,h,A)] = find( relRawPix.maxF_maxGT_meanAccImgs(:,:,d,A) == maxF_meanAccImgs(d,h,A) );
            maxR_meanAccImgs(d,h,A) = relRawPix.maxR_maxGT_meanAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxP_meanAccImgs(d,h,A) = relRawPix.maxP_maxGT_meanAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            %
            maxF_semAccImgs(d,h,A) = relRawPix.maxF_maxGT_semAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxR_semAccImgs(d,h,A) = relRawPix.maxR_maxGT_semAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxP_semAccImgs(d,h,A) = relRawPix.maxP_maxGT_semAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            %
            maxF_stdAccImgs(d,h,A) = relRawPix.maxF_maxGT_stdAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxR_stdAccImgs(d,h,A) = relRawPix.maxR_maxGT_stdAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxP_stdAccImgs(d,h,A) = relRawPix.maxP_maxGT_stdAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            %
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            % max (across GT's) of P/R/F : Relative to GaussRF.
            h = 6;
            maxF_meanAccImgs(d,h,A) = max(max(relGaussRF.maxF_maxGT_meanAccImgs(:,:,d,A)));                                    
            [rM_max(d,h,A),ks_max(d,h,A)] = find( relGaussRF.maxF_maxGT_meanAccImgs(:,:,d,A) == maxF_meanAccImgs(d,h,A) );
            maxR_meanAccImgs(d,h,A) = relGaussRF.maxR_maxGT_meanAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxP_meanAccImgs(d,h,A) = relGaussRF.maxP_maxGT_meanAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            %
            maxF_semAccImgs(d,h,A) = relGaussRF.maxF_maxGT_semAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxR_semAccImgs(d,h,A) = relGaussRF.maxR_maxGT_semAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxP_semAccImgs(d,h,A) = relGaussRF.maxP_maxGT_semAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            %
            maxF_stdAccImgs(d,h,A) = relGaussRF.maxF_maxGT_stdAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxR_stdAccImgs(d,h,A) = relGaussRF.maxR_maxGT_stdAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
            maxP_stdAccImgs(d,h,A) = relGaussRF.maxP_maxGT_stdAccImgs(rM_max(d,h,A),ks_max(d,h,A),d,A);
        end

        disp(method(A))


        % REINSTATE LATER: THIS IS PROBABLY STILL IMPORTANT !
        if(0)
            % Loop thru GT files for image patches in each optimized method and
            % build a matching vector that quantifies degree of gT agreement. ks_max(d,h,A)
            % NOTE: I do not have to do this across all different methods
            % because the ground truths will be the same regardless of method.
            for h = 1:numel(which_F_computation)

                if strcmp(method{A},'IsoDiff')
                    evalDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/Kur_PIF_Fourier1/',method{A},'/benchmark_results/rM',rM{rM_max(d,h,A)},'/NF_60_0/ks',ks{ks_max(d,h,A)},'/'];
                else
                    evalDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/Kur_PIF_Fourier1/',method{A},'/benchmark_results/rM',rM{rM_max(d,h,A)},'/sDInf/sP0p2/NF_60_0/ks',ks{ks_max(d,h,A)},'/'];
                end
                %
                % Loop through each image patch and grab maxF. So later I can compute their mean and std.
                files = dir([evalDir,'*_d1_ev2.txt']);
                %
                for k = 1:numel(files)
                    load([gTdir,files(k).name(1:end-11),'.mat'])
                    %
                    gT_agreement{A,h}(k,1) = F_tot_stats(1);
                    gT_agreement{A,h}(k,2) = F_tot_stats(2);
                    %
                    gT_agreement2{A,h}(k,1) = F_tot_stats2(1);
                    gT_agreement2{A,h}(k,2) = F_tot_stats2(2);
                end
            end
        end

    end % loop over method A = 1:6


    % F_all
    % F_max
    % rM(rM_max)
    % ks(ks_max)


    % rM_max
    % ks_max
    % maxF_meanAccImgs
    % maxF_semAccImgs


    % Save into a mat file all the things I plot in 4 & 5 above so I can plot them with Kuramoto stuff.
    save(['Kur_plot_results',blur_tag_M,'_',which_errbars], 'justRawPix','justGaussRF','justMethod',...
    'relRawPix','relGaussRF','RawPixVsGaussRF','method','which_errbars','blur_tag_M','blur_tit',  ...
    'relative_to_what','which_F_computation','num_cPdist', ...
    'maxF_meanAccImgs','maxF_semAccImgs', 'maxF_stdAccImgs', ...
    'maxR_meanAccImgs','maxR_semAccImgs','maxR_stdAccImgs', ...
    'maxP_meanAccImgs','maxP_semAccImgs','maxP_stdAccImgs', ...
    'rM_max', 'ks_max', 'gT_agreement', 'gT_agreement2', ...
    'maxF_meanGT_struct_method_only', 'maxR_meanGT_struct_method_only', ...
    'maxP_meanGT_struct_method_only', 'thr_meanGT_struct_method_only', ...
    'maxF_maxGT_struct_method_only', 'maxR_maxGT_struct_method_only', ...
    'maxP_maxGT_struct_method_only', 'thr_maxGT_struct_method_only', ...
    'bestGT_maxGT_struct_method_only', ...
    'maxF_meanGT_struct_GaussRF_only', 'maxR_meanGT_struct_GaussRF_only', ...
    'maxP_meanGT_struct_GaussRF_only', 'thr_meanGT_struct_GaussRF_only', ...
    'maxF_maxGT_struct_GaussRF_only', 'maxR_maxGT_struct_GaussRF_only', ...
    'maxP_maxGT_struct_GaussRF_only', 'thr_maxGT_struct_GaussRF_only', ...
    'bestGT_maxGT_struct_GaussRF_only', ...
    'maxF_meanGT_struct_RawPix_only', 'maxR_meanGT_struct_RawPix_only', ...
    'maxP_meanGT_struct_RawPix_only', 'thr_meanGT_struct_RawPix_only', ...
    'maxF_maxGT_struct_RawPix_only', 'maxR_maxGT_struct_RawPix_only', ...
    'maxP_maxGT_struct_RawPix_only', 'thr_maxGT_struct_RawPix_only', ...
    'bestGT_maxGT_struct_RawPix_only', ...
    'img_ptch_name_struct', ...
    'param_matrix', 'param_matrix_wd')

    % idea: average of F_maxGT across cPd's and of F_meanGT across cPd's.

end % try to load in already saved file.

















% Compare the 5 methods (SK wins).

% (6). meanMethod meanRawPix meanBlur maxMethod maxRawPix maxBlur
%  iso
%     d=1
%     d=2
%     d=3
%     d=4
%   aa
%     d
%   sk
%     d

if(0)
    disp('idea: average of F_maxGT across cPds and of F_meanGT across cPds')
    squeeze(mean(maxF_meanAccImgs,1))'
    disp('idea: maybe take max too of F_maxGT across cPds and of F_meanGT across cPds')
    squeeze(max(maxF_meanAccImgs,1))'
end








%% NOW, DO SOME PLOTTING FOR DIFFERENT METHOD, RM, KS COMBINATION COMPARISONS.
%
%
% Package together the different P/R/F computations for Network methods
% absolute performance as well as relRawPix and relGaussRF
for A =1:numel(method)    
    
    for d = 1:num_cPdist % loop thru different d values

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % For best mean value across all ground truthers (GTmean)
        %
        % For F-measure (sort by F-measure)
        y1aF = justMethod.maxF_meanGT_meanAccImgs(:,:,d,A);
        y1bF = justMethod.maxF_meanGT_semAccImgs(:,:,d,A);
        [py1(:,d),qy1(:,d)] = sort(y1aF(:),'descend');   % sort by mean (accImgs) F measure
        py1(:,d) = y1aF(:);                              % mean (unsorted)
        ry1(:,d) = y1bF(:);                              % sem or std (unsorted)
        %
        y2aF = relRawPix.maxF_meanGT_meanAccImgs(:,:,d,A);
        y2bF = relRawPix.maxF_meanGT_semAccImgs(:,:,d,A);
        [py2(:,d),qy2(:,d)] = sort(y2aF(:),'descend');   % mean (sort by mean F measure)
        py2(:,d) = y2aF(:);                              % mean (unsorted)
        ry2(:,d) = y2bF(:);                              % sem or std (unsorted)
        %
        y3aF = relGaussRF.maxF_meanGT_meanAccImgs(:,:,d,A);
        y3bF = relGaussRF.maxF_meanGT_semAccImgs(:,:,d,A);
        [py3(:,d),qy3(:,d)] = sort(y3aF(:),'descend');   % mean (sort by mean F measure)
        py3(:,d) = y3aF(:);                              % mean (unsorted)
        ry3(:,d) = y3bF(:);                              % sem or std (unsorted)
        %
        % for Recall
        y1aR = justMethod.maxR_meanGT_meanAccImgs(:,:,d,A);
        y1bR = justMethod.maxR_meanGT_semAccImgs(:,:,d,A);
        %
        y2aR = relRawPix.maxR_meanGT_meanAccImgs(:,:,d,A);
        y2bR = relRawPix.maxR_meanGT_semAccImgs(:,:,d,A);
        %
        y3aR = relGaussRF.maxR_meanGT_meanAccImgs(:,:,d,A);
        y3bR = relGaussRF.maxR_meanGT_semAccImgs(:,:,d,A);
        %
        % for Precision
        y1aP = justMethod.maxP_meanGT_meanAccImgs(:,:,d,A);
        y1bP = justMethod.maxP_meanGT_semAccImgs(:,:,d,A);
        %
        y2aP = relRawPix.maxP_meanGT_meanAccImgs(:,:,d,A);
        y2bP = relRawPix.maxP_meanGT_semAccImgs(:,:,d,A);
        %
        y3aP = relGaussRF.maxP_meanGT_meanAccImgs(:,:,d,A);
        y3bP = relGaussRF.maxP_meanGT_semAccImgs(:,:,d,A);
        %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % For best value for single best matching ground truth (GTmax)
        %
        % For F-measure (sort by F-measure)
        z1aF = justMethod.maxF_maxGT_meanAccImgs(:,:,d,A);
        z1bF = justMethod.maxF_maxGT_semAccImgs(:,:,d,A);
        [pz1(:,d),qz1(:,d)] = sort(z1aF(:),'descend');   % mean (sort by mean F measure)
        pz1(:,d) = z1aF(:);                              % mean (unsorted)
        rz1(:,d) = z1bF(:);                              % sem or std (unsorted)
        %
        z2aF = relRawPix.maxF_maxGT_meanAccImgs(:,:,d,A);
        z2bF = relRawPix.maxF_maxGT_semAccImgs(:,:,d,A);
        [pz2(:,d),qz2(:,d)] = sort(z2aF(:),'descend');   % mean (sort by mean F measure)
        pz2(:,d) = z2aF(:);                              % mean (unsorted)
        rz2(:,d) = z2bF(:);                              % sem or std (unsorted)
        %
        z3aF = relGaussRF.maxF_maxGT_meanAccImgs(:,:,d,A);
        z3bF = relGaussRF.maxF_maxGT_semAccImgs(:,:,d,A);
        [pz3(:,d),qz3(:,d)] = sort(z3aF(:),'descend'); 	% mean (sort by mean F measure)
        pz3(:,d) = z3aF(:);                              % mean (unsorted)
        rz3(:,d) = z3bF(:);                              % sem or std (unsorted)
        %
        % for Recall
        z1aR = justMethod.maxR_maxGT_meanAccImgs(:,:,d,A);
        z1bR = justMethod.maxR_maxGT_semAccImgs(:,:,d,A);
        %
        z2aR = relRawPix.maxR_maxGT_meanAccImgs(:,:,d,A);
        z2bR = relRawPix.maxR_maxGT_semAccImgs(:,:,d,A);
        %
        z3aR = relGaussRF.maxR_maxGT_meanAccImgs(:,:,d,A);
        z3bR = relGaussRF.maxR_maxGT_semAccImgs(:,:,d,A);
        %
        % for Precision
        z1aP = justMethod.maxP_maxGT_meanAccImgs(:,:,d,A);
        z1bP = justMethod.maxP_maxGT_semAccImgs(:,:,d,A);
        %
        z2aP = relRawPix.maxP_maxGT_meanAccImgs(:,:,d,A);
        z2bP = relRawPix.maxP_maxGT_semAccImgs(:,:,d,A);
        %
        z3aP = relGaussRF.maxP_maxGT_meanAccImgs(:,:,d,A);
        z3bP = relGaussRF.maxP_maxGT_semAccImgs(:,:,d,A);

    end
    %
    %
    % Compute mode across all cp_dist d values. 
    % To see which params scored best across different cp_dist d values. 
    qy1m = mode(qy1,2); % meanGT just method
    qz1m = mode(qz1,2); % maxGT just method
    %
    qy2m = mode(qy2,2); % meanGT relRawPix
    qz2m = mode(qz2,2); % maxGT relRawPix
    %
    qy3m = mode(qy3,2); % meanGT relGaussRF
    qz3m = mode(qz3,2); % maxGT relGaussRF
    

    if(plot_F_errbars_params_optimize_each_method)
        %
        % Plot mean & std F-measure. Mean (across GT) F: Just Method.
        H4 = figure; hold on
        for d = 1:num_cPdist
            errorbar([1:numel(qy1m)]+0.1*d,py1(qy1m,d),ry1(qy1m,d),'LineStyle','none','Marker','o','LineWidth',2,'Color',cPd_colors{d});            % errorbar mean & std
        end
        axis([0.5 numel(qy1m)+2.5 min(min(py1-1.2*ry1)) max(max(py1+1.2*ry1))])
        set(gca,'XTick',1:numel(qy1m)+2,'XTickLabel',[param_matrix(qy1m);'RawPix';'GaussRF'],'FontSize',12,'FontWeight','Bold')
        title([methodStr{A},' ',blur_tit,' : mean (across GT) F :  just Method'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['F-measure (avg across image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        legend('d=0','d=1','d=1.4','d=2')
        saveGoodImg(H4,[dirPre,'../Documentation/Cosyne_2016/F_meanAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_justMethod.jpg'],sizeGoodIm)
        close(H4)
        %
        %
        % Plot mean & std F-measure. Mean (across GT) F: Improvement over image pixels.
        H5 = figure; hold on
        for d = 1:num_cPdist
            errorbar([1:numel(qy2m)]+0.1*d,py2(qy2m,d),ry2(qy2m,d),'LineStyle','none','Marker','o','LineWidth',2,'Color',cPd_colors{d});            % errorbar mean & std
        end
        plot( [0 numel(qy2m)+3], [0 0],'k--')
        axis([0.5 numel(qy2m)+2.5 min(min(py2-1.2*ry2)) max(max(py2+1.2*ry2))])
        set(gca,'XTick',1:numel(qy2m)+2,'XTickLabel',[param_matrix(qy2m);'RawPix';'GaussRF'],'FontSize',12,'FontWeight','Bold')
        title([methodStr{A},' ',blur_tit,' : mean (across GT) F :  relRawPix'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['F-measure (avg across image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        legend('d=0','d=1','d=1.4','d=2')
        saveGoodImg(H5,[dirPre,'../Documentation/Cosyne_2016/F_meanAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_relRawPix.jpg'],sizeGoodIm)
        close(H5)
        %
        %
        % Plot mean & std F-measure. Mean (across GT) F: Improvement over image blur.
        H6 = figure; hold on
        for d = 1:num_cPdist
            errorbar([1:numel(qy3m)]+0.1*d,py3(qy3m,d),ry3(qy3m,d),'LineStyle','none','Marker','o','LineWidth',2,'Color',cPd_colors{d});            % errorbar mean & std
        end
        plot( [0 numel(qy3m)+3], [0 0],'k--')
        axis([0.5 numel(qy3m)+2.5 min(min(py3-1.2*ry3)) max(max(py3+1.2*ry3))])
        set(gca,'XTick',1:numel(qy3m)+2,'XTickLabel',[param_matrix(qy3m);'RawPix';'GaussRF'],'FontSize',12,'FontWeight','Bold')
        title([methodStr{A},' ',blur_tit,' : mean (across GT) F :  rel Im Blur'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['F-measure (avg across image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        legend('d=0','d=1','d=1.4','d=2')
        saveGoodImg(H6,[dirPre,'../Documentation/Cosyne_2016/F_meanAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_relGaussRF.jpg'],sizeGoodIm)
        close(H6)
        %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        %
        % Plot mean & std F-measure. Max (across GT) F: Just Method.
        H7 = figure; hold on
        for d = 1:num_cPdist
            errorbar([1:numel(qz1m)]+0.1*d,pz1(qz1m,d),rz1(qz1m,d),'LineStyle','none','Marker','o','LineWidth',2,'Color',cPd_colors{d});            % errorbar mean & std
        end
        plot( [0 numel(qz1m)+3], [0 0],'k--')
        axis([0.5 numel(qz1m)+2.5 min(min(pz1-1.2*rz1)) max(max(pz1+1.2*rz1))])
        set(gca,'XTick',1:numel(qz1m)+2,'XTickLabel',[param_matrix(qz1m);'RawPix';'GaussRF'],'FontSize',12,'FontWeight','Bold')
        title([methodStr{A},' ',blur_tit,' : max (across GT) F :  just Method'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['F-measure (avg across image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        legend('d=0','d=1','d=1.4','d=2')
        saveGoodImg(H7,[dirPre,'../Documentation/Cosyne_2016/F_maxAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_justMethod.jpg'],sizeGoodIm)
        close(H7)
        %
        %
        % Plot mean & std F-measure. Max (across GT) F: Improvement over image pixels.
        H8 = figure; hold on
        for d = 1:num_cPdist
            errorbar([1:numel(qz2m)]+0.1*d,pz2(qz2m,d),rz2(qz2m,d),'LineStyle','none','Marker','o','LineWidth',2,'Color',cPd_colors{d});            % errorbar mean & std
        end
        plot( [0 numel(qz2m)+3], [0 0],'k--')
        axis([0.5 numel(qz2m)+2.5 min(min(pz2-1.2*rz2)) max(max(pz2+1.2*rz2))])
        set(gca,'XTick',1:numel(qz2m)+2,'XTickLabel',[param_matrix(qz2m);'RawPix';'GaussRF'],'FontSize',12,'FontWeight','Bold')
        title([methodStr{A},' ',blur_tit,' : max (across GT) F :  relRawPix'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['F-measure (avg across image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        legend('d=0','d=1','d=1.4','d=2')
        saveGoodImg(H8,[dirPre,'../Documentation/Cosyne_2016/F_maxAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_relRawPix.jpg'],sizeGoodIm)
        close(H8)
        %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % Plot mean & std F-measure. Max (across GT) F: Improvement over image blur.
        H9 = figure; hold on
        for d = 1:num_cPdist
            errorbar([1:numel(qz3m)]+0.1*d,pz3(qz3m,d),rz3(qz3m,d),'LineStyle','none','Marker','o','LineWidth',2,'Color',cPd_colors{d});            % errorbar mean & std
        end
        plot( [0 numel(qy3m)+3], [0 0],'k--')
        axis([0.5 numel(qz3m)+2.5 min(min(pz3-1.2*rz3)) max(max(pz3+1.2*rz3))])
        set(gca,'XTick',1:numel(qz3m)+2,'XTickLabel',[param_matrix(qz3m);'RawPix';'GaussRF'],'FontSize',12,'FontWeight','Bold')
        title([methodStr{A},' ',blur_tit,' : max (across GT) F :  rel Im Blur'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['F-measure (avg across image patches)'],'FontSize',18,'FontWeight','Bold')
        grid on
        legend('d=0','d=1','d=1.4','d=2')
        saveGoodImg(H9,[dirPre,'../Documentation/Cosyne_2016/F_maxAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_relGaussRF.jpg'],sizeGoodIm)
        close(H9)
    end
    
    
    
    
    % Plot Precision & Recall Separately (instead of combining them together as F-measure)
    % because maybe network methods effect one more than the other or something.
    if(plot_FRP_errorbars_params_optimize_each_method)
        %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % Plot mean & std F-measure,Precision & Recall. Mean (across GT) F: Just Method.
        H4_RP = figure; hold on
        errorbar([1:numel(qy1)]-0.3, reshape(y1aR(qy1),1,[]), reshape(y1bR(qy1),1,[]),'LineStyle','none','Marker','o','LineWidth',2,'Color','Blue');  % Recall
        errorbar([1:numel(qy1)]+0.3, reshape(y1aP(qy1),1,[]), reshape(y1bP(qy1),1,[]),'LineStyle','none','Marker','o','LineWidth',2,'Color','Red'); % Precision
        errorbar([1:numel(qy1)],     reshape(y1aF(qy1),1,[]), reshape(y1bF(qy1),1,[]),'LineStyle','none','Marker','o','LineWidth',2,'Color','Green'); % F-measure
        %
        axis tight % ([0.5 numel(qy1)+2.5 0 1])
        plot([0.5 numel(qy1)+2.5],[0 0],'k--')
        set(gca,'XTick',1:numel(qy1)+2,'XTickLabel',[reshape(param_matrix_wd(qy1),1,[]),'RawPix','GaussRF'],'FontSize',12,'FontWeight','Bold')
        title([methodStr{A},' ',blur_tit,' : mean (across GT) F :  just Method'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['F-measure (avg across image patches)'],'FontSize',18,'FontWeight','Bold')
        xtickangle(90)
        grid on
        legend('Recall','Precision','F-measure')
        %
        saveGoodImg(H4_RP,[dirPre,'../Documentation/Cosyne_2016/FRP_meanAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_justMethod.jpg'],sizeGoodIm)
        close(H4_RP)
        %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % Plot mean & std F-measure,Precision & Recall. Mean (across GT) F: Improvement over RawPix.
        H5_RP = figure; hold on
        errorbar([1:numel(qy2)]-0.3, reshape(y2aR(qy2),1,[]), reshape(y2bR(qy2),1,[]),'LineStyle','none','Marker','o','LineWidth',2,'Color','Blue');  % Recall
        errorbar([1:numel(qy2)]+0.3, reshape(y2aP(qy2),1,[]), reshape(y2bP(qy2),1,[]),'LineStyle','none','Marker','o','LineWidth',2,'Color','Red'); % Precision
        errorbar([1:numel(qy2)],     reshape(y2aF(qy2),1,[]), reshape(y2bF(qy2),1,[]),'LineStyle','none','Marker','o','LineWidth',2,'Color','Green'); % F-measure
        %
        axis tight %([0.5 numel(qy2)+2.5 -0.6 0.6])
        plot( [0 numel(qy2)+3], [0 0],'k--')
        plot([0.5 numel(qy2)+2.5],[0 0],'k--')
        set(gca,'XTick',1:numel(qy2)+2,'XTickLabel',[reshape(param_matrix_wd(qy2),1,[]),'RawPix','GaussRF'],'FontSize',12,'FontWeight','Bold')
        title([methodStr{A},' ',blur_tit,' : mean (across GT) F :  relRawPix'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['F-measure (avg across image patches)'],'FontSize',18,'FontWeight','Bold')
        xtickangle(90)
        grid on
        legend('Recall','Precision','F-measure')
        %
        saveGoodImg(H5_RP,[dirPre,'../Documentation/Cosyne_2016/FRP_meanAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_relRawPix.jpg'],sizeGoodIm)
        close(H5_RP)
        %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % Plot mean & std F-measure,Precision & Recall. Mean (across GT) F: Improvement over GaussRF.
        H6_RP = figure; hold on
        errorbar([1:numel(qy3)]-0.3, reshape(y3aR(qy3),1,[]), reshape(y3bR(qy3),1,[]),'LineStyle','none','Marker','o','LineWidth',2,'Color','Blue');  % Recall
        errorbar([1:numel(qy3)]+0.3, reshape(y3aP(qy3),1,[]), reshape(y3bP(qy3),1,[]),'LineStyle','none','Marker','o','LineWidth',2,'Color','Red'); % Precision
        errorbar([1:numel(qy3)],     reshape(y3aF(qy3),1,[]), reshape(y3bF(qy3),1,[]),'LineStyle','none','Marker','o','LineWidth',2,'Color','Green'); % F-measure
        %
        axis tight %([0.5 numel(qy3)+2.5 -0.6 0.6])
        plot([0 numel(qy3)+3],[0 0],'k--')
        set(gca,'XTick',1:numel(qy3)+2,'XTickLabel',[reshape(param_matrix_wd(qy3),1,[]),'RawPix','GaussRF'],'FontSize',12,'FontWeight','Bold')
        title([methodStr{A},' ',blur_tit,' : mean (across GT) F :  relGaussRF'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['F-measure (avg across image patches)'],'FontSize',18,'FontWeight','Bold')
        xtickangle(90)
        grid on
        legend('Recall','Precision','F-measure')
        %
        saveGoodImg(H6_RP,[dirPre,'../Documentation/Cosyne_2016/FRP_meanAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_relGaussRF.jpg'],sizeGoodIm)
        close(H6_RP)
        %       
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % Plot mean & std F-measure,Precision & Recall. Max (across GT) F: Just Method.
        H7_RP = figure; hold on
        errorbar([1:numel(qz1)]-0.3, reshape(z1aR(qz1),1,[]), reshape(z1bR(qz1),1,[]),'LineStyle','none','Marker','o','LineWidth',2,'Color','Blue');  % Recall
        errorbar([1:numel(qz1)]+0.3, reshape(z1aP(qz1),1,[]), reshape(z1bP(qz1),1,[]),'LineStyle','none','Marker','o','LineWidth',2,'Color','Red'); % Precision
        errorbar([1:numel(qz1)],     reshape(z1aF(qz1),1,[]), reshape(z1bF(qz1),1,[]),'LineStyle','none','Marker','o','LineWidth',2,'Color','Green'); % F-measure
        %
        axis tight %([0.5 numel(qz1)+2.5 -0.6 0.6])
        plot([0 numel(qz1)+3],[0 0],'k--')
        set(gca,'XTick',1:numel(qz1)+2,'XTickLabel',[reshape(param_matrix_wd(qz1),1,[]),'RawPix','GaussRF'],'FontSize',12,'FontWeight','Bold')
        title([methodStr{A},' ',blur_tit,' : max (across GT) F : just Method'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['F-measure (avg across image patches)'],'FontSize',18,'FontWeight','Bold')
        xtickangle(90)
        grid on
        legend('Recall','Precision','F-measure')
        %
        saveGoodImg(H7_RP,[dirPre,'../Documentation/Cosyne_2016/FRP_maxAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_justMethod.jpg'],sizeGoodIm)
        close(H7_RP)
        %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % Plot mean & std F-measure,Precision & Recall. Max (across GT) F: Improvement over RawPix.
        H8_RP = figure; hold on
        errorbar([1:numel(qz2)]-0.3, reshape(z2aR(qz2),1,[]), reshape(z2bR(qz2),1,[]),'LineStyle','none','Marker','o','LineWidth',2,'Color','Blue');  % Recall
        errorbar([1:numel(qz2)]+0.3, reshape(z2aP(qz2),1,[]), reshape(z2bP(qz2),1,[]),'LineStyle','none','Marker','o','LineWidth',2,'Color','Red'); % Precision
        errorbar([1:numel(qz2)],     reshape(z2aF(qz2),1,[]), reshape(z2bF(qz2),1,[]),'LineStyle','none','Marker','o','LineWidth',2,'Color','Green'); % F-measure
        %
        axis tight %([0.5 numel(qz2)+2.5 -0.6 0.6])
        plot([0 numel(qz2)+3],[0 0],'k--')
        set(gca,'XTick',1:numel(qz2)+2,'XTickLabel',[reshape(param_matrix_wd(qz2),1,[]),'RawPix','GaussRF'],'FontSize',12,'FontWeight','Bold')
        title([methodStr{A},' ',blur_tit,' : max (across GT) F : relRawPix'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['F-measure (avg across image patches)'],'FontSize',18,'FontWeight','Bold')
        xtickangle(90)
        grid on
        legend('Recall','Precision','F-measure')
        %
        saveGoodImg(H8_RP,[dirPre,'../Documentation/Cosyne_2016/FRP_maxAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_relRawPix.jpg'],sizeGoodIm)
        close(H8_RP)
        %
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % Plot mean & std F-measure,Precision & Recall. Max (across GT) F: Improvment over GaussRF.
        H9_RP = figure; hold on
        errorbar([1:numel(qz3)]-0.3, reshape(z3aR(qz3),1,[]), reshape(z3bR(qz3),1,[]),'LineStyle','none','Marker','o','LineWidth',2,'Color','Blue');  % Recall
        errorbar([1:numel(qz3)]+0.3, reshape(z3aP(qz3),1,[]), reshape(z3bP(qz3),1,[]),'LineStyle','none','Marker','o','LineWidth',2,'Color','Red'); % Precision
        errorbar([1:numel(qz3)],     reshape(z3aF(qz3),1,[]), reshape(z3bF(qz3),1,[]),'LineStyle','none','Marker','o','LineWidth',2,'Color','Green'); % F-measure
        %
        axis tight %([0.5 numel(qz3)+2.5 -0.6 0.6])
        plot([0 numel(qz3)+3],[0 0],'k--')
        set(gca,'XTick',1:numel(qz3)+2,'XTickLabel',[reshape(param_matrix_wd(qz3),1,[]),'RawPix','GaussRF'],'FontSize',12,'FontWeight','Bold')
        title([methodStr{A},' ',blur_tit,' : max (across GT) F : relGaussRF'],'FontSize',18,'FontWeight','Bold')
        xlabel('Parameters','FontSize',18,'FontWeight','Bold')
        ylabel(['F-measure (avg across image patches)'],'FontSize',18,'FontWeight','Bold')
        xtickangle(90)
        grid on
        legend('Recall','Precision','F-measure')
        %
        saveGoodImg(H9_RP,[dirPre,'../Documentation/Cosyne_2016/FRP_maxAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_allParams_relGaussRF.jpg'],sizeGoodIm)
        close(H9_RP)
    end
    
end
    











%% Get min and max for Precision & Recall for plotting 2D Histograms to be consistent across different methods.
%
% RPFlims is max and min values for Precision, Recall & F across all 6 Methods.
%         size is (maxminRPF=6, relative_to_what=3, num_cPdist=4)
%
%           maxminRPF = {maxR, minR, maxP, minP, maxF, minF}
%           relative_to_what = {JustMethod, relRawPix, relGaussRF}
%           cPdists = {0, 1, 1.4, 2}
% 
%
RPFlims = zeros(6,3,num_cPdist); % numel(which_F_computation)); % Just recording 3 different min & max values for the
% different relative_to_what values. Saving max & min lims on Recall, Precision, F-measure. 

for d = 1:num_cPdist
    for h = 1:numel(which_F_computation) % loop thru combinations of F computation & what relative to (12)
        for A = 1:numel(method) % loop thru different methods or methods_all.

            switch which_F_computation{h}
                %
                case 'meanGT'
                    %
                    switch relative_to_what{h}
                        case 'justMethod'
                            %disp('meanGT justMethod')
                            maxF_pre = 0; %maxF_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                            maxR_pre = 0; %maxR_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                            maxP_pre = 0; %maxP_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                            r=1; % this indexes into max & min limits
                        case 'relRawPix'
                            %disp('meanGT relPix')
                            maxF_pre = maxF_meanGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                            maxR_pre = maxR_meanGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                            maxP_pre = maxP_meanGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                            r=2; % this indexes into max & min limits
                        case 'relGaussRF'
                            %disp('meanGT relBlur')
                            maxF_pre = maxF_meanGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                            maxR_pre = maxR_meanGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                            maxP_pre = maxP_meanGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                            r=3; % this indexes into max & min limits
                    end
                    %
                    maxF_post = maxF_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                    maxR_post = maxR_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                    maxP_post = maxP_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                    %
                    thr_method = thr_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                    thr_RawPix  = thr_meanGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                    %if(blur_flg)
                    thr_GaussRF = thr_meanGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                    %end
                %
                case 'maxGT'
                    %
                    switch relative_to_what{h}
                        case 'justMethod'
                            %disp('maxGT justMethod')
                            maxF_pre = 0; %maxF_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                            maxR_pre = 0; %maxR_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                            maxP_pre = 0; %maxP_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                            r=1; % this indexes into max & min limits
                        case 'relRawPix'
                            %disp('maxGT relPix')
                            maxF_pre = maxF_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                            maxR_pre = maxR_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                            maxP_pre = maxP_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                            r=2; % this indexes into max & min limits
                        case 'relGaussRF'
                            %disp('maxGT relBlur')
                            maxF_pre = maxF_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                            maxR_pre = maxR_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                            maxP_pre = maxP_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                            r=3; % this indexes into max & min limits
                    end
                    %
                    maxF_post = maxF_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                    maxR_post = maxR_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                    maxP_post = maxP_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                    %
                    thr_method = thr_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                    thr_RawPix  = thr_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                    %if(blur_flg)
                    thr_GaussRF = thr_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                    %end
                    %
                    bestGT_method = bestGT_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                    bestGT_RawPix = bestGT_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                    %if(blur_flg)
                    bestGT_GaussRF = bestGT_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A}; 
                    %end
                %
            end
            % 
            % maxF_pre is for relRawPix or relGaussRF. If justMethod, it will be all zeros.
            maxF_temp = maxF_post - maxF_pre; % 
            maxR_temp = maxR_post - maxR_pre; % 
            maxP_temp = maxP_post - maxP_pre;



            % These are max & min values across different methods and F-measure computations.  
            % The r keeps track of different relative-to-what values.
            RPFlims(1,r,d) = max([ RPFlims(1,r,d); maxR_temp ]); % max R
            RPFlims(2,r,d) = min([ RPFlims(2,r,d); maxR_temp ]); % min R
            RPFlims(3,r,d) = max([ RPFlims(3,r,d); maxP_temp ]); % max P
            RPFlims(4,r,d) = min([ RPFlims(4,r,d); maxP_temp ]); % min P
            RPFlims(5,r,d) = max([ RPFlims(5,r,d); maxF_temp ]); % max F
            RPFlims(6,r,d) = min([ RPFlims(6,r,d); maxF_temp ]); % min F

        end % loop over A = 1:5 method

    end % loop over h = 1:12 which_F_computation & relative_to_what
    
end % loop over d = 1:num_cPdists 














%% Here, we plot Violin Plots of F-measure comparing different optimized methods and  Statistical Significance 
if(0)

    methods_all = {method{:}}; % , 'RawPix', 'GaussRF'
    
    % Ok, now that we have solved for consistent limits for plotting, we plot F-measure with 
    % error bars and a 2D Histogram in R-P space.
    
    for d = 1:num_cPdist % different correspondPixel distances
    
        for h = 1:numel(which_F_computation) % loop thru combinations of F computation & what relative to (12)
                       % relative_to_what

            H=figure; hold on
            plot([0.5 numel(methods_all)+0.5],[0 0],'k--','LineWidth',1.5)

            Pr = ones(numel(method),num_cPdist);  % significance value (from ranksum)
            Hr = zeros(numel(method),num_cPdist); % whether significant (from ranksum)

            for A = 1:numel(method) % loop thru different methods or methods_all.

                switch which_F_computation{h}
                    %
                    case 'meanGT'
                        %
                        switch relative_to_what{h}
                            case 'justMethod'
                                disp('meanGT justMethod')
                                maxF_pre = 0; %maxF_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                maxR_pre = 0; %maxR_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                maxP_pre = 0; %maxP_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                r=1; % this indexes into max & min limits
                            case 'relRawPix'
                                disp('meanGT relPix')
                                maxF_pre = maxF_meanGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                maxR_pre = maxR_meanGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                maxP_pre = maxP_meanGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                r=2; % this indexes into max & min limits
                            case 'relGaussRF'
                                disp('meanGT relBlur')
                                maxF_pre = maxF_meanGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                maxR_pre = maxR_meanGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                maxP_pre = maxP_meanGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                r=3; % this indexes into max & min limits
                        end
                        %
                        maxF_post = maxF_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        maxR_post = maxR_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        maxP_post = maxP_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        %
                        thr_method = thr_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        thr_RawPix  = thr_meanGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        thr_GaussRF = thr_meanGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                    %
                    case 'maxGT'
                        %
                        switch relative_to_what{h}
                            case 'justMethod'
                                disp('maxGT justMethod')
                                maxF_pre = 0; %maxF_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                maxR_pre = 0; %maxR_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                maxP_pre = 0; %maxP_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                r=1; % this indexes into max & min limits
                            case 'relRawPix'
                                disp('maxGT relPix')
                                maxF_pre = maxF_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                maxR_pre = maxR_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                maxP_pre = maxP_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                r=2; % this indexes into max & min limits
                            case 'relGaussRF'
                                disp('maxGT relBlur')
                                maxF_pre = maxF_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                maxR_pre = maxR_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                maxP_pre = maxP_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                r=3; % this indexes into max & min limits
                        end
                        %
                        maxF_post = maxF_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        maxR_post = maxR_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        maxP_post = maxP_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        %
                        thr_method = thr_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        thr_RawPix  = thr_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        thr_GaussRF = thr_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        %
                        bestGT_method = bestGT_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        bestGT_RawPix  = bestGT_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        bestGT_GaussRF = bestGT_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A}; 
                    %
                end
                %
                maxF_temp = maxF_post - maxF_pre; % this relRawPix or relGaussRF.
                maxR_temp = maxR_post - maxR_pre; % if justMethod, this will be all zeros
                maxP_temp = maxP_post - maxP_pre;

                violin(maxF_temp,'x',[A A+1])

                % Determine Statistical Significance of difference in Network Method F-measure distribution across images vs GaussRF or RawPix
                % distribution.  Use Mann-Whitney U (aka ranksum) test.
                [Pr(A,d),Hr(A,d)] = ranksum(maxF_pre,maxF_post);
                % best_vec_for_ranksum{1,i}, best_vec_for_ranksum{2,i}

            end % loop over method A=1:numel(method) 


            errorbar( [1:numel(method)],     squeeze(maxF_meanAccImgs(d,h,:)), squeeze(maxF_semAccImgs(d,h,:)), 'bd', 'Linewidth', 3, 'Color','green' )
            errorbar( [1:numel(method)]-0.3, squeeze(maxR_meanAccImgs(d,h,:)), squeeze(maxR_semAccImgs(d,h,:)), 'bd', 'Linewidth', 3, 'Color','blue')
            errorbar( [1:numel(method)]+0.3, squeeze(maxP_meanAccImgs(d,h,:)), squeeze(maxP_semAccImgs(d,h,:)), 'bd', 'Linewidth', 3, 'Color','red' )
            %
            errorbar( [1:numel(method)],     squeeze(maxF_meanAccImgs(d,h,:)), squeeze(maxF_stdAccImgs(d,h,:)), 'bd', 'Linewidth', 3, 'Color','green' )
            errorbar( [1:numel(method)]-0.3, squeeze(maxR_meanAccImgs(d,h,:)), squeeze(maxR_stdAccImgs(d,h,:)), 'bd', 'Linewidth', 3, 'Color','blue')
            errorbar( [1:numel(method)]+0.3, squeeze(maxP_meanAccImgs(d,h,:)), squeeze(maxP_stdAccImgs(d,h,:)), 'bd', 'Linewidth', 3, 'Color','red' )
            %
            % Labels for Precision, Recall & F-measure error bars.
            text(1+0.0, maxF_meanAccImgs(d,h,1) - maxF_stdAccImgs(d,h,1) - 0.02, ['F'], 'Color','green','VerticalAlignment','top','HorizontalAlignment','center','FontSize',18,'FontWeight','Bold')
            text(1-0.3, maxR_meanAccImgs(d,h,1) - maxR_stdAccImgs(d,h,1) - 0.02, ['R'], 'Color','blue','VerticalAlignment','top','HorizontalAlignment','center','FontSize',18,'FontWeight','Bold')
            text(1+0.3, maxP_meanAccImgs(d,h,1) - maxP_stdAccImgs(d,h,1) - 0.02, ['P'], 'Color','red','VerticalAlignment','top','HorizontalAlignment','center','FontSize',18,'FontWeight','Bold')


            % Fill in whether distributions are Statistically Significant on plot
            for A = 1:numel(method)
                pTxtBeg = {['\mu=',num2str( maxF_meanAccImgs(d,h,A) ,'%+5.2f')],['\sigma=',num2str( maxF_semAccImgs(d,h,A) ,'%+5.2f')],['SE=',num2str( maxF_semAccImgs(d,h,A) ,'%+5.2f')]}
                %text(A+0.1, maxF_meanAccImgs(d,h,A)+2*maxF_semAccImgs(d,h,A) , {['\mu=',num2str( maxF_meanAccImgs(d,h,A) ,'%+5.2f')],[plot_name,'=',num2str( maxF_semAccImgs(d,h,A) ,'%+5.2f')]}, ...
                %'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold') % indicate the mean, std & skew
                if ( Pr(A,d) < 0.00001)
                    pTxtFin = {['{rM',rM{rM_max(d,h,A)},',ks',ks{ks_max(d,h,A)},'}***'],['p=',num2str(Pr(A,d),'%5.2e')]}
                    
                elseif ( Pr(A,d) < 0.001)
                    pTxtFin = {['{rM',rM{rM_max(d,h,A)},',ks',ks{ks_max(d,h,A)},'}***'],['p=',num2str(Pr(A,d),'%5.4f')]}
                    %text(A+0.1, maxF_meanAccImgs(d,h,A)-2*maxF_semAccImgs(d,h,A) , {['{rM',rM{rM_max(d,h,A)},',ks',ks{ks_max(d,h,A)},'}***'],['p=',num2str(Pr(A,d),'%5.4f')]}, ...
                    %'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold') % optimized parameters
                elseif ( Pr(A,d) < 0.01)
                    pTxtFin = {['{rM',rM{rM_max(d,h,A)},',ks',ks{ks_max(d,h,A)},'}**'],['p=',num2str(Pr(A,d),'%5.3f')]}
                    %text(A+0.1, maxF_meanAccImgs(d,h,A)-2*maxF_semAccImgs(d,h,A) , {['{rM',rM{rM_max(d,h,A)},',ks',ks{ks_max(d,h,A)},'}**'],['p=',num2str(Pr(A,d),'%5.3f')]}, ...
                    %'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold') % optimized parameters
                elseif ( Pr(A,d) < 0.05)
                    pTxtFin = {['{rM',rM{rM_max(d,h,A)},',ks',ks{ks_max(d,h,A)},'}*'],['p=',num2str(Pr(A,d),'%5.2f')]}
                    %text(A+0.1, maxF_meanAccImgs(d,h,A)-2*maxF_semAccImgs(d,h,A) , {['{rM',rM{rM_max(d,h,A)},',ks',ks{ks_max(d,h,A)},'}*'],['p=',num2str(Pr(A,d),'%5.2f')]}, ...
                    %'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold') % optimized parameters
                else
                    pTxtFin = {['{rM',rM{rM_max(d,h,A)},',ks',ks{ks_max(d,h,A)},'}'],['p=',num2str(Pr(A,d),'%5.2f')]}
                    %text(A+0.1, maxF_meanAccImgs(d,h,A)-2*maxF_semAccImgs(d,h,A) , {['{rM',rM{rM_max(d,h,A)},',ks',ks{ks_max(d,h,A)},'}'],['p=',num2str(Pr(A,d),'%5.2f')]}, ...
                    %'VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold') % optimized parameters
                end
                %
                text( A+0.1, max( maxP_meanAccImgs(d,h,A)+maxP_stdAccImgs(d,h,A), maxR_meanAccImgs(d,h,A)+maxR_stdAccImgs(d,h,A) ) + 0.02, ...
                    { pTxtFin{1}, pTxtFin{2}, pTxtBeg{1}, pTxtBeg{2} }, 'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold') % optimized parameters
            end

            set(gca,'XTick',1:numel(method),'XTickLabel',methodStr,'FontSize',18,'FontWeight','Bold')
            title([which_F_computation{h},' : ',relative_to_what{h},' : ',blur_tit,' : (d=',num2str(d),')'],'FontSize',24,'FontWeight','Bold')
            xlabel('Method with Optimized Parameters','FontSize',20,'FontWeight','Bold')
            ylabel(['F-measure (avg across image patches)'],'FontSize',20,'FontWeight','Bold')
            % axis([0.5 numel(methods_all)+0.5 min(maxF_meanAccImgs(d,h,A)-3*maxF_stdAccImgs(d,h,A)) max(maxF_meanAccImgs(d,h,A)+3*maxF_stdAccImgs(d,h,A))])
            axis([ 0.5 numel(method)+0.7 min(RPFlims(6,r,:)) max(RPFlims(5,r,:)) ]) % this will keep plot axes consistent across different d values
            grid on
            % NOTE: CAN MAYBE JUST DO 'AXIS TIGHT' INSTEAD OF USING RPFLIMS THING... (?)

            
            %plot2svg([dirPre,'../Documentation/Cosyne_2016/Fmax_compareMethodsKurBestParams',blur_tag_M,'_',which_F_computation{h},'_',relative_to_what{h},'_d',num2str(d),'_errbars',which_errbars,'.svg'],H)
            saveGoodImg(H,[dirPre,'../Documentation/Cosyne_2016/Fmax_compareMethodsKurBestParams',blur_tag_M,'_',which_F_computation{h},'_',relative_to_what{h},'_d',num2str(d),'_errbars',which_errbars,'.jpg'],sizeGoodIm)
            close(H)

        end
    end
end
        



%% Scatter plot in R-P Space -  plot the vector from (R,P)_blur to (R,P)_method. 
%  Similar as above, but plot this instead of plotting delta R & delta P in a 2D plane from -1 to 1.
%  I want to  with similar color coding for points with delta F > or <= 0.
if(1)
    
    
    n_better    = zeros(num_cPdist, numel(which_F_computation), numel(method));
    n_worse     = zeros(num_cPdist, numel(which_F_computation), numel(method));
    n_total     = zeros(num_cPdist, numel(which_F_computation), numel(method));
    
    
    
    % Ok, now that we have solved for consistent limits for plotting, we plot F-measure with 
    % error bars and a 2D Histogram in R-P space.
    
    for d = 4 %1:num_cPdist % different correspondPixel distances
    
        for h = 1:numel(which_F_computation) % loop thru combinations of F computation & what relative to (12)
                       % relative_to_what
                       
                       
            % preallocate memory to hold distributions for Precision, Recall, F-measure before & after method.
            nBins = 10;
            nFcnts_post0 = zeros(nBins,numel(method));
            mf = zeros(3,numel(method)+1);
            nPcnts_post0 = zeros(nBins,numel(method));
            mp = zeros(3,numel(method)+1);
            nRcnts_post0 = zeros(nBins,numel(method));
            mr = zeros(3,numel(method)+1);
            

            for A = 4 % 1:numel(method) % loop thru different methods or methods_all.

                switch which_F_computation{h}
                    %    
                    case 'meanGT'
                        maxF_post = maxF_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        maxR_post = maxR_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        maxP_post = maxP_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        %
                        thr_method = thr_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        thr_RawPix  = thr_meanGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        thr_GaussRF = thr_meanGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        %
                        switch relative_to_what{h}
                            case 'justMethod'
                                maxF_pre = zeros(size(maxF_post));
                                maxR_pre = zeros(size(maxR_post));
                                maxP_pre = zeros(size(maxP_post));
                                r=1; % this indexes into max & min limits
                            case 'relRawPix'
                                maxF_pre = maxF_meanGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                maxR_pre = maxR_meanGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A}; % should be pix.
                                maxP_pre = maxP_meanGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                r=1; % this indexes into max & min limits
                            case 'relGaussRF'
                                maxF_pre = maxF_meanGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                maxR_pre = maxR_meanGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A}; % should be blur
                                maxP_pre = maxP_meanGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                r=1; % this indexes into max & min limits
                                
                        end
                    %    
                    case 'maxGT'
                        maxF_post = maxF_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        maxR_post = maxR_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        maxP_post = maxP_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        %
                        thr_method = thr_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        thr_RawPix  = thr_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        thr_GaussRF = thr_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        %
                        bestGT_method = bestGT_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        bestGT_RawPix = bestGT_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        bestGT_GaussRF = bestGT_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A}; 
                        %
                        switch relative_to_what{h}
                            case 'justMethod'
                                maxF_pre = zeros(size(maxF_post));
                                maxR_pre = zeros(size(maxR_post));
                                maxP_pre = zeros(size(maxP_post));
                                r=1; % this indexes into max & min limits
                            case 'relRawPix'
                                maxF_pre = maxF_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                maxR_pre = maxR_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A}; % should be pix.
                                maxP_pre = maxP_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                r=1; % this indexes into max & min limits
                            case 'relGaussRF'
                                maxF_pre = maxF_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                maxR_pre = maxR_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A}; % should be blur
                                maxP_pre = maxP_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                r=1; % this indexes into max & min limits
                        end
                    %    
                end % switch which_f_computation


               

                maxX = 1%RPFlims(1,r); 
                minX = 0%RPFlims(2,r); 
                maxY = 1%RPFlims(3,r); 
                minY = 0%RPFlims(4,r); 
                maxF = 1%RPFlims(5,r);
                minF = 0%RPFlims(6,r);
                %




                % find entries with F (or delta F) values above & below mean.
                ind1 = find(maxF_post-maxF_pre>0); % was mf
                ind2 = find(maxF_post-maxF_pre<=0); % was mf
                n_better(d,h,A) = numel(ind1);
                n_worse(d,h,A) = numel(ind2);
                n_total(d,h,A) = numel(maxF_temp);



                
                r_lims = linspace(minX,maxX,nBins);  % Recall
                p_lims = linspace(minY,maxY,nBins)'; % Precision
                f_lims = linspace(minF,maxF,nBins)'; % F-measure

                
                
                
                
                % F-measure
                [nFcnts_pre0,f_lims] = hist(maxF_pre,f_lims);
                [nFcnts_pre1,f_lims] = hist(maxF_pre(ind1),f_lims);
                [nFcnts_pre2,f_lims] = hist(maxF_pre(ind2),f_lims);
                %
                [nFcnts_post0(:,A),f_lims] = hist(maxF_post,f_lims);
                [nFcnts_post1,f_lims] = hist(maxF_post(ind1),f_lims);
                [nFcnts_post2,f_lims] = hist(maxF_post(ind2),f_lims);
                %
                mf(1,1) = mean(maxF_pre);
                mf(2,1) = std(maxF_pre);
                mf(3,1) = sem(maxF_pre);
                mf(1,A+1) = mean(maxF_post);
                mf(2,A+1) = std(maxF_post);
                mf(3,A+1) = sem(maxF_post);


                % Precision
                [nPcnts_pre0,p_lims] = hist(maxP_pre,p_lims);
                [nPcnts_pre1,p_lims] = hist(maxP_pre(ind1),p_lims);
                [nPcnts_pre2,p_lims] = hist(maxP_pre(ind2),p_lims);
                %
                [nPcnts_post0(:,A),p_lims] = hist(maxP_post,p_lims);
                [nPcnts_post1,p_lims] = hist(maxP_post(ind1),p_lims);
                [nPcnts_post2,p_lims] = hist(maxP_post(ind2),p_lims);
                %
                mp(1,1) = mean(maxP_pre);
                mp(2,1) = std(maxP_pre);
                mp(3,1) = sem(maxP_pre);
                mp(1,A+1) = mean(maxP_post);
                mp(2,A+1) = std(maxP_post);
                mp(3,A+1) = sem(maxP_post);

                % Recall
                [nRcnts_pre0,r_lims] = hist(maxR_pre,r_lims);
                [nRcnts_pre1,r_lims] = hist(maxR_pre(ind1),r_lims);
                [nRcnts_pre2,r_lims] = hist(maxR_pre(ind2),r_lims);
                %
                [nRcnts_post0(:,A),r_lims] = hist(maxR_post,r_lims);
                [nRcnts_post1,r_lims] = hist(maxR_post(ind1),r_lims);
                [nRcnts_post2,r_lims] = hist(maxR_post(ind2),r_lims);
                %
                mr(1,1) = mean(maxR_pre);
                mr(2,1) = std(maxR_pre);
                mr(3,1) = sem(maxR_pre);
                mr(1,A+1) = mean(maxR_post);
                mr(2,A+1) = std(maxR_post);
                mr(3,A+1) = sem(maxR_post);


                
                
                % NOTE: I HAVE COME IN HERE AND COMMENTED OUT THINGS AND
                % MADE SOME HARDCODED CHANGES TO PLOT IMPROVEMENT AND
                % DEGRADING IMAGE PATCHES SEPARATELY TO MAKE THINGS MORE
                % CLEAR FOR THE PAPER. JUST FYI...
                if(1)
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
                    scatter(0.70,0.90,1700,'w.','Linewidth',2) % cover the green dot in Iso-F contours
                    colormap('bone')
                    %
                    ylabel('Recall','FontSize',18,'FontWeight','Bold')
                    xlabel('Precision')
                    % axis square

                    %axis([minX maxX minY maxY])
                    axis([0 1 0 1])



                    % 1st draw lines connecting pre & post
                    %quiver(maxP_pre(ind1)',    maxR_pre(ind1)',   maxP_post(ind1)'-maxP_pre(ind1)',   maxR_post(ind1)'-maxR_pre(ind1)', 0 ,'Color','blue', 'MaxHeadSize',0.05);
                    quiver(maxP_pre(ind2)',    maxR_pre(ind2)',   maxP_post(ind2)'-maxP_pre(ind2)',   maxR_post(ind2)'-maxR_pre(ind2)', 0 ,'Color','red', 'MaxHeadSize',0.05);

                    alpha 1
                   
                    %
                    % draw pre points
                    %scatter_patches(maxP_pre(ind1), maxR_pre(ind1), 2, 'c','o', 'EdgeColor','Black'); % 'FaceAlpha',abs(maxF_post(ind1)-maxF_pre(ind1))/max(abs(maxF_post-maxF_pre)),
                    %scatter_patches(maxP_pre(ind2), maxR_pre(ind2), 2, 'm','o', 'EdgeColor','Black'); % 'FaceAlpha',abs(maxF_post(ind2)-maxF_pre(ind2))/max(abs(maxF_post-maxF_pre)), 

                    %
                    % draw post points
                    %scatter_patches(maxP_post(ind1), maxR_post(ind1), 6, 'b','^', 'EdgeColor','Black'); % 'FaceAlpha',abs(maxF_post(ind1)-maxF_pre(ind1))/max(abs(maxF_post-maxF_pre)),
                    %scatter_patches(maxP_post(ind2), maxR_post(ind2), 16, 'r','v',  'EdgeColor','Black'); % 'FaceAlpha',abs(maxF_post(ind2)-maxF_pre(ind2))/max(abs(maxF_post-maxF_pre)),
                    
                    grid on
                    % Ntot is actually... ,num2str(n_total(d,h,A)).  wtf ever.
                    text(.97, 0.03 , {['Ntot = 500'],['\color{blue}Nbetter = ',num2str(n_better(d,h,A))], ['\color{red}Nworse = ',num2str(n_worse(d,h,A))]}, ...
                         'verticalalignment','bottom' ,'horizontalalignment','right' ,'BackgroundColor','w')
                    
                    title([method{A},' : ',which_F_computation{h},' : ',relative_to_what{h}, ' : ',blur_tit],'FontSize',20,'FontWeight','Bold')



                    % Plot Marginalized Precision Distributions splitting them up by images with improved F and degraded F.
                    subplot(5,5,[21:24])
                    hold on

%                     plot(p_lims,nPcnts_pre0./n_total(d,h,A),'k--','LineWidth',2)
%                     plot(p_lims,nPcnts_pre1./n_total(d,h,A),'c--','LineWidth',1)
%                     plot(p_lims,nPcnts_post1./n_total(d,h,A),'b-','LineWidth',2)
%                     plot(p_lims,nPcnts_pre2./n_total(d,h,A),'m--','LineWidth',1)
%                     plot(p_lims,nPcnts_post2./n_total(d,h,A),'r-','LineWidth',2)
                    %
                     
%                     plot(p_lims,nPcnts_pre1,'c--','LineWidth',1)
%                     plot(p_lims,nPcnts_post1,'b-','LineWidth',2)
                    plot(p_lims,nPcnts_pre2,'m--','LineWidth',1)
                    plot(p_lims,nPcnts_post2,'r-','LineWidth',2)
                    grid on
                    %
                   % legend({'before all','before improve','after improve','before degrade','after degrade'}) %  
                    legend({'before network','after network'}) %  
                    ylabel('counts')
                    xlabel('Precision','FontSize',18,'FontWeight','Bold')


                    % Plot Marginalized Recall Distributions splitting them up by images with improved F and degraded F.
                    subplot(5,5,[5,10,15,20]), 
                    hold on
%                     plot(nRcnts_pre0./n_total(d,h,A),r_lims,'k--','LineWidth',2)
%                     plot(nRcnts_pre1./n_total(d,h,A),r_lims,'c--','LineWidth',1)
%                     plot(nRcnts_post1./n_total(d,h,A),r_lims,'b-','LineWidth',2)
%                     plot(nRcnts_pre2./n_total(d,h,A),r_lims,'m--','LineWidth',1)
%                     plot(nRcnts_post2./n_total(d,h,A),r_lims,'r-','LineWidth',2)
                    %
%                     plot(nRcnts_pre1,r_lims,'c--','LineWidth',1)
%                     plot(nRcnts_post1,r_lims,'b-','LineWidth',2)
                    plot(nRcnts_pre2,r_lims,'m--','LineWidth',1)
                    plot(nRcnts_post2,r_lims,'r-','LineWidth',2)
%                     grid on
                    %
                    xlabel('counts')
                    ylabel('Recall')


                    % Plot F-measure "Marginalization"  Distributions splitting them up by images with improved F and degraded F.
                    subplot(5,5,[25]), hold on
%                     plot(f_lims,nFcnts_pre0./n_total(d,h,A),'k--','LineWidth',2.5)
%                     plot(f_lims,nFcnts_pre1./n_total(d,h,A),'c--','LineWidth',1)
%                     plot(f_lims,nFcnts_post1./n_total(d,h,A),'b-','LineWidth',2.5)
%                     plot(f_lims,nFcnts_pre2./n_total(d,h,A),'m--','LineWidth',1)
%                     plot(f_lims,nFcnts_post2./n_total(d,h,A),'r-','LineWidth',2.5)
                    %
%                     plot(f_lims,nFcnts_pre1,'c--','LineWidth',1)
%                     plot(f_lims,nFcnts_post1,'b-','LineWidth',2.5)
                    plot(f_lims,nFcnts_pre2,'m--','LineWidth',1)
                    plot(f_lims,nFcnts_post2,'r-','LineWidth',2.5)
                    grid on
                    %
                    ylabel('counts')
                    xlabel('F-measure','FontSize',18,'FontWeight','Bold')

                    saveGoodImg(H,[dirPre,'../Documentation/Cosyne_2016/RPmax_lines_',which_F_computation{h},'_',relative_to_what{h},'_',method{A},'_d',num2str(d),'_rM',rM{rM_max(d,h,A)},'_ks',ks{ks_max(d,h,A)},blur_tag_M,'.png'], [0 0 0.6 1]) 
                    close(H)
                end


            end % loop over A = 1:numel(method)
            
            
            
            %% # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
            %  Plot distributions of F,P,R for GaussRF and all 6 network models. 
            %  Also plot mean, std & sem.
            
            if(0)
                % F-measure
                Hf=figure; hold on
                plot(f_lims,nFcnts_pre0,'c--','LineWidth',10)
                %
                heb=herrorbar(mf(1,1), max(nFcnts_post0(:))+4, mf(2,1), 'c^'); % mean + std. NOTE: mf(3,:) I bet is Standard Error.
                plot( [mf(1,1), mf(1,1)], [0, max(nFcnts_post0(:))+4]   ,'c--','LineWidth',2)
                set(heb,'LineWidth',6);
                %
                for A = 1:numel(method)
                    plot(f_lims,nFcnts_post0(:,A),method_colors{A},'LineWidth',6)
                    heb=herrorbar(mf(1,A+1), max(nFcnts_post0(:))+4*(A+1), mf(2,A+1), [method_colors{A},'^']);
                    set(heb,'LineWidth',6);
                end
                %
                xlabel('F-measure','FontSize',18,'FontWeight','Bold')
                ylabel('counts','FontSize',18,'FontWeight','Bold')
                set(gca,'FontSize',16,'FontWeight','Bold')
                grid on
                %
                %plot2svg([dirPre,'../Documentation/Cosyne_2016/CompareDistributions_F_',relative_to_what{h},blur_tag_M,'_d',num2str(d),'.svg'],Hf)
                saveGoodImg(Hf,[dirPre,'../Documentation/Cosyne_2016/CompareDistributions_F_',which_F_computation{h},'_',relative_to_what{h},blur_tag_M,'_d',num2str(d),'.png'],sizeGoodIm)
                close(Hf)
                %
                % %
                %
                % Precision
                Hp=figure; hold on
                plot(p_lims,nPcnts_pre0,'c--','LineWidth',10)
                %
                heb=herrorbar(mp(1,1), max(nPcnts_post0(:))+4, mp(2,1), 'c^');
                plot( [mp(1,1), mp(1,1)], [0, max(nPcnts_post0(:))+4]   ,'c--','LineWidth',2)
                set(heb,'LineWidth',6);
                %
                for A = 1:numel(method)
                    plot(p_lims,nPcnts_post0(:,A),method_colors{A},'LineWidth',6)
                    heb=herrorbar(mp(1,A+1), max(nPcnts_post0(:))+4*(A+1), mp(2,A+1), [method_colors{A},'^']);
                    set(heb,'LineWidth',6);
                end
                %
                xlabel('Precision','FontSize',18,'FontWeight','Bold')
                ylabel('counts','FontSize',18,'FontWeight','Bold')
                set(gca,'FontSize',16,'FontWeight','Bold')
                grid on
                %             
                %plot2svg([dirPre,'../Documentation/Cosyne_2016/CompareDistributions_P_',relative_to_what{h},blur_tag_M,'_d',num2str(d),'.svg'],Hp)
                saveGoodImg(Hp,[dirPre,'../Documentation/Cosyne_2016/CompareDistributions_P_',which_F_computation{h},'_',relative_to_what{h},blur_tag_M,'_d',num2str(d),'.png'],sizeGoodIm)
                close(Hp)
                %
                % %
                %
                % Recall
                Hr=figure; hold on
                plot(r_lims,nRcnts_pre0,'c--','LineWidth',10)
                %
                heb=herrorbar(mr(1,1), max(nRcnts_post0(:))+4, mr(2,1), 'c^');
                plot( [mr(1,1), mr(1,1)], [0, max(nRcnts_post0(:))+4]   ,'c--','LineWidth',2)
                set(heb,'LineWidth',6);
                %
                for A = 1:numel(method)
                    plot(r_lims,nRcnts_post0(:,A),method_colors{A},'LineWidth',6)
                    heb=herrorbar(mr(1,A+1), max(nRcnts_post0(:))+4*(A+1), mr(2,A+1), [method_colors{A},'^']);
                    set(heb,'LineWidth',6);
                end
                %
                xlabel('Recall','FontSize',18,'FontWeight','Bold')
                ylabel('counts','FontSize',18,'FontWeight','Bold')
                set(gca,'FontSize',16,'FontWeight','Bold')
                grid on
                %
                %plot2svg([dirPre,'../Documentation/Cosyne_2016/CompareDistributions_P_',relative_to_what{h},blur_tag_M,'_d',num2str(d),'.svg'],Hp)
                saveGoodImg(Hr,[dirPre,'../Documentation/Cosyne_2016/CompareDistributions_R_',which_F_computation{h},'_',relative_to_what{h},blur_tag_M,'_d',num2str(d),'.png'],sizeGoodIm)
                close(Hr)
            end
            
       
        end % loop over h = 1:6 for which_F_computation & relative_to_what
    
    end % loop over d = 1:numel(cpD)
    
end % if(1)





%% # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
% 
% <Delta>F vs cP_dist for all methods on single plot.  Show errorbars for mean and std
% of distribution across 500 image patches. Note: statistical significance.
% Do this for justMethod, relRawPix & relGaussRF.
if(0)
    
    % TO DO:
    
    % TO FRITZ: Black line IsoDiff kinda same in Blur v NotBlur imgs. Can
    % put that in plot somehow? TRY IT
    
    % put violins back in? with standard error?
    
%     which_F_computation = {'meanGT','meanGT','meanGT','maxGT','maxGT','maxGT'};
%     relative_to_what = {'justMethod','relRawPix','relGaussRF','justMethod','relRawPix','relGaussRF'};
   
    
    for h = 1:3 % half the size. trying to plot maxGT & meanGT on same plot.
        
        H=figure; hold on
        
        % NOTE: TO DO: IF IT IS JUST METHOD, PLOT RawPix (MAGENTA) & GaussRF (CYAN) TOO FOR COMPARISON.
        switch relative_to_what{h}
            case 'justMethod'
                errorbar( [1:num_cPdist] +0.01*(-1),   justRawPix.maxF_maxGT_meanAccImgs, justRawPix.maxF_maxGT_semAccImgs,['m-'], 'LineWidth', 4)
                errorbar( [1:num_cPdist] +0.01*(-1),   justRawPix.maxF_meanGT_meanAccImgs, justRawPix.maxF_meanGT_semAccImgs,['m--'], 'LineWidth', 4)
                %
                errorbar( [1:num_cPdist] +0.01*(-2),   justGaussRF.maxF_maxGT_meanAccImgs, justGaussRF.maxF_maxGT_semAccImgs,['c-'], 'LineWidth', 4)
                errorbar( [1:num_cPdist] +0.01*(-2),   justGaussRF.maxF_meanGT_meanAccImgs, justGaussRF.maxF_meanGT_semAccImgs,['c--'], 'LineWidth', 4)
                % %
                scatter( [1:num_cPdist] +0.01*(-1),   justRawPix.maxF_maxGT_meanAccImgs, 100, 'k','o', 'filled')
                scatter( [1:num_cPdist] +0.01*(-1),   justRawPix.maxF_meanGT_meanAccImgs, 100, 'k','o', 'filled')
                %
                scatter( [1:num_cPdist] +0.01*(-2),   justGaussRF.maxF_maxGT_meanAccImgs, 100, 'k','o', 'filled')
                scatter( [1:num_cPdist] +0.01*(-2),   justGaussRF.maxF_meanGT_meanAccImgs, 100, 'k','o', 'filled')
                
            case 'relRawPix'
                errorbar( [1:num_cPdist] +0.01*(-1), [0 0 0 0], [0 0 0 0], ['m-'], 'LineWidth', 4)
                errorbar( [1:num_cPdist] +0.01*(-1), [0 0 0 0], [0 0 0 0], ['m--'], 'LineWidth', 4)
                %
                errorbar( [1:num_cPdist] +0.01*(-2), RawPixVsGaussRF.maxF_maxGT_meanAccImgs, RawPixVsGaussRF.maxF_maxGT_semAccImgs, ['c-'], 'LineWidth', 4)
                errorbar( [1:num_cPdist] +0.01*(-2), RawPixVsGaussRF.maxF_meanGT_meanAccImgs, RawPixVsGaussRF.maxF_meanGT_semAccImgs, ['c--'], 'LineWidth', 4)
                % %
                scatter( [1:num_cPdist] +0.01*(-1),   [0 0 0 0], 100, 'k','o', 'filled')
                scatter( [1:num_cPdist] +0.01*(-1),   [0 0 0 0], 100, 'k','o', 'filled')
                %
                scatter( [1:num_cPdist] +0.01*(-2),   RawPixVsGaussRF.maxF_maxGT_meanAccImgs, 100, 'k','o', 'filled')
                scatter( [1:num_cPdist] +0.01*(-2),   RawPixVsGaussRF.maxF_meanGT_meanAccImgs, 100, 'k','o', 'filled')
                
                
            case 'relGaussRF'
                errorbar( [1:num_cPdist] +0.01*(-1), -RawPixVsGaussRF.maxF_maxGT_meanAccImgs, RawPixVsGaussRF.maxF_maxGT_semAccImgs, ['m-'], 'LineWidth', 4)
                errorbar( [1:num_cPdist] +0.01*(-1), -RawPixVsGaussRF.maxF_meanGT_meanAccImgs, RawPixVsGaussRF.maxF_meanGT_semAccImgs, ['m--'], 'LineWidth', 4)
                %
                errorbar( [1:num_cPdist] +0.01*(-2), [0 0 0 0], [0 0 0 0], ['c-'], 'LineWidth', 4)
                errorbar( [1:num_cPdist] +0.01*(-2), [0 0 0 0], [0 0 0 0], ['c--'], 'LineWidth', 4)
                % %
                scatter( [1:num_cPdist] +0.01*(-1),   -RawPixVsGaussRF.maxF_maxGT_meanAccImgs, 100, 'k','o', 'filled')
                scatter( [1:num_cPdist] +0.01*(-1),   -RawPixVsGaussRF.maxF_meanGT_meanAccImgs, 100, 'k','o', 'filled')
                %
                scatter( [1:num_cPdist] +0.01*(-2),   [0 0 0 0], 100, 'k','o', 'filled')
                scatter( [1:num_cPdist] +0.01*(-2),   [0 0 0 0], 100, 'k','o', 'filled')
                
        end
        
        
        
        
        for A = 1:numel(method)
            errorbar( [1:num_cPdist] +0.01*A,   squeeze(maxF_meanAccImgs(:,3+h,A)), squeeze(maxF_semAccImgs(:,3+h,A)), ... % 'maxGT'
                [method_colors{A},'-'], 'LineWidth', 3) % ,method_shapes{A}
            
            errorbar( [1:num_cPdist] +0.01*A,   squeeze(maxF_meanAccImgs(:,h,A)), squeeze(maxF_semAccImgs(:,h,A)), ...     % 'meanGT'
                [method_colors{A},'--'], 'LineWidth', 3) % method_shapes{A},
        end
        
        
        %legend(legendMethods{1},[],legendMethods{2},[],legendMethods{3},[],legendMethods{4},[],legendMethods{5},[],...
        %'Location','Best','FontSize',16,'FontWeight','Bold') % Put in good legend.
        % title([relative_to_what{i},' ',blur_tit],'FontSize',20,'FontWeight','Bold')
        xlabel(['d_t'],'FontSize',18,'FontWeight','Bold')
        ylabel(['\Delta F  '],'FontSize',18,'FontWeight','Bold')
        set(gca,'XTick',[1:num_cPdist],'XTickLabel',[0 1 1.4 2],'FontSize',16,'FontWeight','Bold')

        plot( [0.5 num_cPdist+0.5], [0 0],'k.-.','LineWidth',2)
        grid on
        
        % Hardcode y-axis limits.
        switch i
            case 1
                if ~isempty(strmatch(which_errbars,'sem'))
                    ylim([0 0.45]); % justMethod w/ sem.
                else
                    ylim([0 0.60]); % justMethod w/ std.
                end
            case 2
                if ~isempty(strmatch(which_errbars,'sem'))
                    ylim([-0.01 0.09]); % relRawPix w/ sem. 
                else
                    ylim([-0.10 0.18]); % relRawPix w/ std.
                end
                
            case 3
                if ~isempty(strmatch(which_errbars,'sem'))
                    ylim([-0.05 0.06]); % relGaussRF w/ sem. 
                else
                    ylim([-0.11 0.11]); % relGaussRF w/ std.
                end
        end
        
%        % Would use this but it doesnt compare y-axis for Pre-Blur vs Not.
%         y_vals = [maxF_meanAccImgs(:,i,:),maxF_meanAccImgs(:,3+i,:)]; 
%         ylim( [ min(y_vals(:)) , max(y_vals(:)) ] ); 

        disp([dirPre,'../Documentation/Cosyne_2016/F_vs_cpD_optimized_methods_',relative_to_what{h},blur_tag_M,'_errbars_',which_errbars,'.svg'])
        %plot2svg([dirPre,'../Documentation/Cosyne_2016/F_vs_cpD_optimized_methods_',relative_to_what{h},blur_tag_M,'_errbars_',which_errbars,'.svg'],H)

        saveGoodImg(H,[dirPre,'../Documentation/Cosyne_2016/F_vs_cpD_optimized_methods_',relative_to_what{h},blur_tag_M,'_errbars_',which_errbars,'.jpg'],sizeGoodIm)
        close(H)
        
        
        % TODO:  MAYBE ALSO MAKE A SIMILAR PLOT FOR R AND P SO WE CAN SAY
        % SOMETHING ABOUT WHAT THE MODELS ARE CHANGING MORE.
        
    end
    
end




%% # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
% Here, I plot <Delta> F-measure before and after method 
if(0)
    
    d = 4;
    
    maxDelF = -1.*ones(size(which_F_computation));
    minDelF = 1.*ones(size(which_F_computation));
    
    for h = 1:numel(which_F_computation) % loop thru combinations of F computation & what relative to 
        % which_F_computation   =  {'meanGT',     'meanGT',   'meanGT',    'maxGT',      'maxGT',    'maxGT'    }         
        % relative_to_what      =  {'justMethod', 'relRawPix', 'relGaussRF', 'justMethod', 'relRawPix', 'relGaussRF'}
                   
       for A = 1:numel(method)
           
           switch which_F_computation{h}
                % % % % % % % % % % % % % %
                %
                case 'meanGT'
                maxF_post = maxF_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                maxR_post = maxR_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                maxP_post = maxP_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                %
                thr_method = thr_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                thr_RawPix  = thr_meanGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                thr_GaussRF = thr_meanGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                %
                switch relative_to_what{h}
                    case 'justMethod'
                        continue % dont want to do this justMethod analysis (move on to next pp loop hopefully).
                    case 'relRawPix'
                        maxF_pre = maxF_meanGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        maxR_pre = maxR_meanGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A}; % should be pix.
                        maxP_pre = maxP_meanGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        r=1; % this indexes into max & min limits
                    case 'relGaussRF'
                        maxF_pre = maxF_meanGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        maxR_pre = maxR_meanGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A}; % should be blur
                        maxP_pre = maxP_meanGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        r=1; % this indexes into max & min limits
                end
                % % % % % % % % % % % % % %
                %
                case 'maxGT'
                maxF_post = maxF_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                maxR_post = maxR_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                maxP_post = maxP_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                %
                thr_method = thr_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                thr_RawPix  = thr_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                thr_GaussRF = thr_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                %
                bestGT_method = bestGT_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                bestGT_RawPix = bestGT_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                bestGT_GaussRF = bestGT_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A}; 
                %
                switch relative_to_what{h}
                    case 'justMethod'
                        continue % dont want to do this justMethod analysis (move on to next h loop hopefully).
                    case 'relRawPix'
                        maxF_pre = maxF_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        maxR_pre = maxR_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A}; % should be RawPix.
                        maxP_pre = maxP_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        r=1; % this indexes into max & min limits
                    case 'relGaussRF'
                        maxF_pre = maxF_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        maxR_pre = maxR_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A}; % should be GaussRF
                        maxP_pre = maxP_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        r=1; % this indexes into max & min limits
                end
                %
            end % switch which_f_computation
            %
            % Want to set the Y-axis of these plots to be same across Optimized Methods
            maxDelF(h) = max( [ maxDelF(h); maxF_post-maxF_pre ] );
            minDelF(h) = min( [ minDelF(h); maxF_post-maxF_pre ] );
            %
       end % end loop over A = 1:numel(method)
       %      
    end % end loop over h = 1:numel(which_F_computation)
            
   
    
    
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %for f = 1:numel(which_delF)
        %
        for h = 1:numel(which_F_computation) % loop thru combinations of F computation & what relative to 
            % which_F_computation   =  {'meanGT',     'meanGT',   'meanGT',    'maxGT',      'maxGT',    'maxGT'    }         
            % relative_to_what      =  {'justMethod', 'relRawPix', 'relGaussRF', 'justMethod', 'relRawPix', 'relGaussRF'}
            %
            for A = 1:numel(method) % loop thru different methods or methods_all.
                %
                % Directory to probabalistic boundary png image. 
                pbDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/Kur_PIF_Fourier1/',...
                            method{A},'/pb_png/rM',rM{rM_max(d,h,A)},'/sDInf/sP0p2/NF_60_0/ks',ks{ks_max(d,h,A)},'/'];  

                switch which_F_computation{h}
                    %
                    case 'meanGT'
                        maxF_post = maxF_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        maxR_post = maxR_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        maxP_post = maxP_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        %
                        thr_method = thr_meanGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        thr_RawPix  = thr_meanGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        thr_GaussRF = thr_meanGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        %
                        switch relative_to_what{h}
                            case 'justMethod'
                                continue % dont want to do this justMethod analysis (move on to next pp loop hopefully).
                            case 'relRawPix'
                                maxF_pre = maxF_meanGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                maxR_pre = maxR_meanGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A}; % should be pix.
                                maxP_pre = maxP_meanGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                r=1; % this indexes into max & min limits
                            case 'relGaussRF'
                                maxF_pre = maxF_meanGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                maxR_pre = maxR_meanGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A}; % should be blur
                                maxP_pre = maxP_meanGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                r=1; % this indexes into max & min limits
                        end

                    case 'maxGT'
                        maxF_post = maxF_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        maxR_post = maxR_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        maxP_post = maxP_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        %
                        thr_method = thr_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        thr_RawPix  = thr_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        thr_GaussRF = thr_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        %
                        bestGT_method = bestGT_maxGT_struct_method_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        bestGT_RawPix = bestGT_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                        bestGT_GaussRF = bestGT_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A}; 
                        %
                        switch relative_to_what{h}
                            case 'justMethod'
                                continue % dont want to do this justMethod analysis (move on to next pp loop hopefully).
                            case 'relRawPix'
                                maxF_pre = maxF_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                maxR_pre = maxR_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A}; % should be pix.
                                maxP_pre = maxP_maxGT_struct_RawPix_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                r=1; % this indexes into max & min limits
                            case 'relGaussRF'
                                maxF_pre = maxF_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                maxR_pre = maxR_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A}; % should be blur
                                maxP_pre = maxP_maxGT_struct_GaussRF_only{rM_max(d,h,A),ks_max(d,h,A),d,A};
                                r=1; % this indexes into max & min limits
                        end

                end %  switch which_F_computation{h}



                % Plot F_blur vs F_method
                if(0)
                    H=figure;
                    hold on
                    %scatter(maxF_pre(ind1), maxF_post(ind1),'b.')
                    %scatter(maxF_pre(ind2), maxF_post(ind2),'r.')
                    s = scatter(maxF_pre, maxF_post, 20, 'k', 'filled');
                    s.MarkerFaceAlpha = 0.5;
                    axis([0 1 0 1])
                    plot([0 1],[0 1],'k--')
                    title([methodStr{A},' : ',which_F_computation{h},' : ',relative_to_what{h}, ' : ',blur_tit],'FontSize',20,'FontWeight','Bold') % ,' : ',which_delF{f},' performance'
                    xlabel('F-measure before method','FontSize',18,'FontWeight','Bold')
                    ylabel('F-measure after method','FontSize',18,'FontWeight','Bold')
                    set(gca,'FontSize',16,'FontWeight','Bold')
                    text(1,0,{['Ntot = ',num2str(n_total(d,h,A))],['\color{red}Nbetter = ',num2str(n_better(d,h,A))],...
                        ['\color{blue}Nworse = ',num2str(n_worse(d,h,A))]},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16,'FontWeight','Bold')
                    axis square
                    grid on
                    %
                    [Pr,Hr] = ranksum(maxF_pre, maxF_post);
                    text(0.1, 0.9, ['p = ',num2str(Pr)], 'VerticalAlignment','top','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
                    %
                    saveGoodImg(H,[dirPre,'../Documentation/Cosyne_2016/FvsF_',which_F_computation{h},'_',relative_to_what{h},'_',method{A},'_rM',rM{rM_max(d,h,A)},'_ks',ks{ks_max(d,h,A)},blur_tag_M,'.jpg'],sizeGoodIm)
                    close(H)
                end

                % % % % % %
                
                
                
                
               if(0) 
               % Here, just plot scatter of Gaussian RF indep. sensors F vs
               % F for each of the 3 modularity methods on one scatter plot
               % with p-values. Can do for h==6 :: maxGT and h==3 meanGT.
               %
               % A==4 :: TMod 1D
               % A==5 :: TMod 2D
               % A==6 :: Modularity
               % 
               % PLOT TM1D and TM2D scatter of F_{IS} vs. F_{TM} on same
               % axis with scatter points colored. Hardcoding it.
               % if( (A==4 | A==5 | A==6) && (h==6) ) 
               
               
                   for hh = [3,6]
                       
                       H=figure;
                       hold on
                       cols ={'red','green','blue'};
                       for AA = [4,5,6]


                               maxF_pre = maxF_maxGT_struct_GaussRF_only{rM_max(d,hh,AA),ks_max(d,hh,AA),d,AA};
                               maxF_post = maxF_maxGT_struct_method_only{rM_max(d,hh,AA),ks_max(d,hh,AA),d,AA};
                               s = scatter(maxF_pre, maxF_post, 20, cols{AA-3}, 'filled');
                               s.MarkerFaceAlpha = 0.8;
                               %
                               [Pr,Hr] = ranksum(maxF_pre, maxF_post);
                               text(0.1, 1 - (AA-3)*0.05, ['\color{',cols{AA-3},'}', methodStr{AA}, ' p = ',num2str(Pr)], 'VerticalAlignment','top','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
                               % 
                       end
                       
                        axis([0 1 0 1])
                        plot([0 1],[0 1],'k--')
                        %title([methodStr{A},' : ',which_F_computation{h},' : ',relative_to_what{h}, ' : ',blur_tit],'FontSize',20,'FontWeight','Bold') % ,' : ',which_delF{f},' performance'
                        xlabel('F-measure before method','FontSize',18,'FontWeight','Bold')
                        ylabel('F-measure after method','FontSize',18,'FontWeight','Bold')
                        set(gca,'FontSize',16,'FontWeight','Bold')
                        %text(1,0,{['Ntot = ',num2str(n_total(d,h,A))],['\color{red}Nbetter = ',num2str(n_better(d,h,A))],...
                        %['\color{blue}Nworse = ',num2str(n_worse(d,h,A))]},'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16,'FontWeight','Bold')
                        axis square
                        grid on

                        saveGoodImg(H,[dirPre,'../Documentation/Cosyne_2016/FvsF_',which_F_computation{hh},'_',relative_to_what{hh},'_3Modularities.jpg'],sizeGoodIm)
                        close(H)

                   end
                   
               end
                
                
                
                
                
                
                
                


                % Plot F-gaussRF vs <Delta>F 
                if(0)
                    [p,S,mu] = polyfit(maxF_pre,maxF_post-maxF_pre,1);
                    x = linspace( min(maxF_pre), max(maxF_pre) );
                    %
                    H1=figure;
                    hold on
                    F_before = maxF_pre;
                    delF_aft = maxF_post-maxF_pre;
                    s = scatter(F_before, delF_aft, 30, 'k',  'filled');
                    s.MarkerFaceAlpha = 0.3;
                    ylabel(['\Delta F ', legendMethods{A}],'FontSize',20,'FontWeight','Bold')
                    xlabel('F_G Gaussian RF Independent Sensors','FontSize',20,'FontWeight','Bold')
                    set(gca,'FontSize',18,'FontWeight','Bold')
                    plot([0 max(maxF_pre)],[0 0],'k--')
                    plot(x, p(1) + p(2).*x,'m--','LineWidth',1.5 )
                    text( max(maxF_pre), 0, ... %min(maxF_post-maxF_pre), ...
                        ['\color{magenta}Line Fit: \color{black} \DeltaF = ',num2str(p(2),2),' F_G + ',num2str(p(1),2)],...
                        'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',16,'FontWeight','Bold')
                    %
                    text( max(maxF_pre), max(maxF_post-maxF_pre) ,{['Ntot = ',num2str(n_total(d,h,A))],['\color{red}Nbetter = ',num2str(n_better(d,h,A))],...
                        ['\color{blue}Nworse = ',num2str(n_worse(d,h,A))]},'VerticalAlignment','top','HorizontalAlignment','right','FontSize',16,'FontWeight','Bold')
                    grid on
                    %ylim([minDelF(h),maxDelF(h)]);
                    %title([ legendMethods{A}],'FontSize',24,'FontWeight','Bold')
                    axis('tight')    


                    
                    
                    % Within this plot, number some points of certain improvement quality (better, average or worse) spread across different F_blur quality
                    % and then visualize the image patches for them in another figure. (NOTE: Only wanna do for [ Mod_SK : MaxGT : relGaussRF ]
                    if( (A==4 | A==5) && (h==3 | h==6) )  % if Mod_SK, MaxGT, relGaussRF

                        % break up F_pre into N bins and then grab N points with certain <Delta> F properties within each bin.
                        % N=3;
                        F_blur_bins = linspace( min(maxF_pre), max(maxF_pre), M+1 );
                        ind_tot = [];
                        ind_tot_loc = [];

                        fm=0;
                        for f = 1:numel(which_delF)
                            for i = 1:M % Looping over F-blur bins

                                ind_blur = find(F_before > F_blur_bins(i) & F_before <= F_blur_bins(i+1));
                                [hp,ind_delF] = sort(delF_aft(ind_blur));

                            
                                %
                                try
                                    switch which_delF{f}
                                        case 'better'
                                            ind = ind_blur(ind_delF(end-N+1:end));
                                        case 'worse'
                                            ind = ind_blur(ind_delF(1:N));
                                        case 'average'
                                            ind = ind_blur(ind_delF(round(numel(ind_blur)./2)-N/2+1:round(numel(ind_blur)./2)+N/2));
                                    end
                                catch
                                    % I think the problem is there is not enough points in this bin.  So just take em all.
                                    ind = ind_blur;
                                end

                                ind_loc = N*(i-1) + [1:numel(ind)];
                                ind_tot_loc = [ind_tot_loc,ind_loc];
                                ind_tot = [ind_tot;ind];
                                
                                delF_colors = {'red','green','blue'};
                                
                                
                                
                                % Label the image patches in the scatter plot with numbers
                                for j = 1:numel(ind)
                                    fm = fm+1;
                                    text( F_before(ind(j)), delF_aft(ind(j)), ['\color{',delF_colors{f},'}', num2str(fm) ], 'fontsize',20,'fontweight','bold');
                                end

                                %scatter(F_blur_bins,zeros(size(F_blur_bins)),300,'k+','linewidth',3)

                            end % f = numel(which_delF) = better, worse, average

                        end % i = 1:M % Looping over F-blur bins and number of individual images showing pbs and boundaries 

                        saveGoodImg(H1,[dirPre,'../Documentation/Cosyne_2016/Fblur_v_delF_',which_F_computation{h},'_',relative_to_what{h},'_',method{A},'_rM',rM{rM_max(d,h,A)},'_ks',ks{ks_max(d,h,A)},blur_tag_M,'.jpg'],[0 0 0.8 1])
                        


                        % Make a figure with image patches arranged according to F-blur and <Delta> F after applying method (using ind_tot from above).
                        % TO DO: LOWER CONTRAST OF IMAGE PATCHES SO LINES STAND OUT MORE.
                        % (1). Plot image patch, boundary found by method & best matching Ground Truth
                        for f = 1:numel(which_delF)
                            
%                             ind_tot_loc = [];
%                             ind_tot  = [];

                            H2=figure; 
                            ha = tight_subplot( 3, M*N, [0.01 0], [0.05 0.05], [0.05 0]);
                            im_Fnames = img_ptch_name_struct{rM_max(d,h,A),ks_max(d,h,A),A};

                            Fblur = maxF_pre;
                            delF = maxF_post-maxF_pre;

                            for i = 1:numel(ind_tot)

                                load([imDir, im_Fnames{ind_tot(i)}(1:end-3), '.mat']); % method
                                load([gTdir, im_Fnames{ind_tot(i)}(1:end-3), '.mat']); % ground truth
                                %    
                                % load pb png image and threshold it to show what boundaries let to max F-measure with maxGT criteria.
                                pbFile = [pbDir,im_Fnames{ind_tot(i)}(1:end-3), '.png'];
                                pb = double(imread(pbFile))/255;
                                %
                                % load pb png image and threshold it to show what boundaries let to max F-measure with maxGT criteria.
                                pbBlurFile = [pbDirBlur,im_Fnames{ind_tot(i)}(1:end-3), '.png'];
                                pbBlur = double(imread(pbBlurFile))/255;


                                subplot( ha(i) ), hold on
                                %disp( [min(im(:)), max(im(:))]  )   % NOTE: I THINK IM GOES FROM 0 to 1 BY CONSTUCTION.
                                imagesc(im, [-0.5,1.5]), colormap(bone), freezeColors % Here set clim to lower contrast of image patch
                                axis square tight ij
                                set(gca,'XTick',[],'YTick',[])
                               
                                %
                                Igt1 = groundTruth{bestGT_method(ind_tot(i))}.Boundaries; % GT that best matched threshold pb and gave max F-measure using maxGT.
                                [x1,y1] = find(Igt1);
                                %
                                Igt2 = ( pbBlur >=thr_GaussRF(ind_tot(i)) ); % Boundaries determined from thresholded probabalistic boundary from blur (that gave maxF)
                                [x2,y2] = find(Igt2);
                                %
                                Igt3 = ( pb >=thr_method(ind_tot(i)) ); % Boundaries determined from thresholded probabalistic boundary from method (that gave maxF)
                                [x3,y3] = find(Igt3);
                                %
                                Igt4 = ( Igt2 & Igt3 ); % Boundaries determined from thresholded probabalistic boundary from method (that gave maxF)
                                [x4,y4] = find(Igt4);
                                %
                                subplot( ha(i) ), hold on
                                scatter(y1,x1,5,'k','.')
                                scatter(y2,x2,25,'c','.')
                                scatter(y3,x3,25,'r','.')
                                scatter(y4,x4,23,'y','.')
                                axis square ij
                                axis([0 size(im,1) 0 size(im,2)])
                                set(gca,'XTick',[],'YTick',[])
%                                 text(1,size(im,2),[num2str(i),'). \color{cyan}{',num2str(Fblur(ind_tot(i)),2),'}','\color{red}{',num2str(delF(ind_tot(i)),'%+5.2f'),'}'], ...
%                                     'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',16,'FontWeight','Bold')
                                %
                                ylabel([num2str(i),'). \color{cyan}{',num2str(Fblur(ind_tot(i)),2),'}','\color{red}{',num2str(delF(ind_tot(i)),'%+5.2f'),'}'],'FontSize',16,'FontWeight','Bold')
                                %
                            end                    
                            %
                            annotation('textbox', [0 0.9 1 0.1],'String',[legendMethods{A},'  ::  \color{cyan}{F_G} \color{red}{ + \Delta F}'], ...
                            'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold')

                            saveGoodImg(H2,[dirPre,'../Documentation/Cosyne_2016/Fblur_v_delF_',which_F_computation{h},'_',relative_to_what{h},'_',method{A},'_rM',rM{rM_max(d,h,A)},'_ks',ks{ks_max(d,h,A)},blur_tag_M,'_visPatchesMethod.jpg'],[0 0 0.7 1])
                            close(H2)

                        end % for f = 1:which_delF = better, worse, average performance.

                    end % if( (A==4 | A==5) && (h==3 | h==6) ) ie. (if TM1D or TM2D) and (MaxGT or meanGT) relGaussRF

                end % if(1) - Plot these image patches.

                try
                    close(H1)
                catch
                    disp('meh')
                end
                
            end % loop over A = 1:numel(method)

        end % loop over h = 1:numel(which_F_computation) maxGT & relGaussRF

    end % if(1) - plot image patches and scatter plots.



%% TO DO: Visualize statistically how things change for Mod_SKHEuc as we vary Rmax and Ks parameters.


% (1). Fix rM to best value for TM 2D. Vary ks. Plot F-measure stats.
% Plot mean & std F-measure. Max (across GT) F: Just relGaussRF.
if(0)
    h = 6; % maxGT
    A = 4; % TM 2D - best params were {rM=10,ks=mid}
    r = 4; % index into rM value that was best for TM 2D
    %
    pz3 = zeros(num_cPdist,numel(ks));
    qz3 = zeros(num_cPdist,numel(ks));
    rz3 = zeros(num_cPdist,numel(ks));
    sz3 = zeros(num_cPdist,numel(ks));
    
    H7 = figure; hold on

    for d = 1:num_cPdist % loop thru different d values

%         % For best mean value of F-measure across all Imgs (meanGT)
%         y3aF = relGaussRF.maxF_meanGT_meanAccImgs(r,:,d,A);
%         y3bF = relGaussRF.maxF_meanGT_semAccImgs(r,:,d,A);
%         [py3(:,d),qy3(:,d)] = sort(y3aF(:),'descend');   % mean (sort by mean F measure)
%         py3(:,d) = y3aF(:);                              % mean (unsorted)
%         ry3(:,d) = y3bF(:);                              % sem or std (unsorted)
        %
        % For best mean value of F-measure across all Imgs (maxGT)
        z3aF = relGaussRF.maxF_maxGT_meanAccImgs(r,:,d,A);
        z3bF = relGaussRF.maxF_maxGT_semAccImgs(r,:,d,A);
        z3cF = relGaussRF.maxF_maxGT_stdAccImgs(r,:,d,A);
        [pz3(d,:),qz3(d,:)] = sort(z3aF,'descend'); 	 % mean (sort by mean F measure)    
        %
        pz3(d,:) = z3aF;                              % mean (unsorted)
        rz3(d,:) = z3bF;                              % sem(unsorted)
        sz3(d,:) = z3cF;                              % std (unsorted)
    end
    %
    % Plot mean & std F-measure. Max (across GT) F: Just Method.
    for d = 1:num_cPdist % loop thru different d values
        %errorbar([1:numel(ks)]+0.1*d-0.25,pz3(d,:),rz3(d,:),'LineStyle','none','Marker','+','MarkerSize',12,'LineWidth',4,'Color',cPd_colors{d});
        errorbar([1:numel(ks)]+0.1*d-0.25,pz3(d,:),sz3(d,:),'LineStyle','none','Marker','+','MarkerSize',12,'LineWidth',6,'Color',cPd_colors{d});
    end
    %
    plot( [0 numel(ks)+1], [0 0],'k--')
    %axis([0.5 numel(ks)+0.5 1.1*min(min(pz3-1.2*rz3)) 1.1*max(max(pz3+1.2*rz3))])
    axis([0.5 numel(ks)+0.5 1.1*min(min(pz3-1.2*sz3)) 1.1*max(max(pz3+1.2*sz3))])
    set(gca,'XTick',1:numel(ks),'XTickLabel',ks,'FontSize',16,'FontWeight','Bold')
    title([methodStr{A},' ',blur_tit,' : max (across GT)'],'FontSize',18,'FontWeight','Bold')
    xlabel(['Spring Constant scaling (ks) w/ optimized rM = ',rM{r}],'FontSize',18,'FontWeight','Bold')
    ylabel(['\Delta F_G (\mu , \sigma)'],'FontSize',18,'FontWeight','Bold')
    grid on
    legend('d_t=0','d_t=1','d_t=1.4','d_t=2','Location','Best')
    
    saveGoodImg(H7,[dirPre,'../Documentation/Cosyne_2016/F_maxAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_rMopt_relGaussRF.jpg'],[0 0 0.5 1])
    close(H7)
    
end



% (2). Fix ks to best value for TM 2D. Vary rM. Plot F-measure stats.
% Plot mean & std F-measure. Max (across GT) F: relGaussRF.
if(0)
    h = 6; % maxGT
    A = 4; % TM 2D - best params were {rM=10,ks=mid}
    k = 2; % index into ks value that was best for TM 2D
    %
    pz3 = zeros(num_cPdist,numel(rM));
    qz3 = zeros(num_cPdist,numel(rM));
    rz3 = zeros(num_cPdist,numel(rM));
    sz3 = zeros(num_cPdist,numel(rM));
    
    H7 = figure; hold on
    
    for d = 1:num_cPdist % loop thru different d values

%         % For best mean value of F-measure across all Imgs (meanGT)
%         y3aF = relGaussRF.maxF_meanGT_meanAccImgs(r,:,d,A);
%         y3bF = relGaussRF.maxF_meanGT_semAccImgs(r,:,d,A);
%         [py3(:,d),qy3(:,d)] = sort(y3aF(:),'descend');   % mean (sort by mean F measure)
%         py3(:,d) = y3aF(:);                              % mean (unsorted)
%         ry3(:,d) = y3bF(:);                              % sem or std (unsorted)
        %
        % For best mean value of F-measure across all Imgs (maxGT)
        z3aF = relGaussRF.maxF_maxGT_meanAccImgs(:,k,d,A);
        z3bF = relGaussRF.maxF_maxGT_semAccImgs(:,k,d,A);
        z3cF = relGaussRF.maxF_maxGT_stdAccImgs(:,k,d,A);
        [pz3(d,:),qz3(d,:)] = sort(z3aF,'descend'); 	 % mean (sort by mean F measure)    
        %
        pz3(d,:) = z3aF;                              % mean (unsorted)
        rz3(d,:) = z3bF;                              % sem(unsorted)
        sz3(d,:) = z3cF;                              % std (unsorted)
    end
    %
    % Plot mean & std F-measure. Max (across GT) F: Just Method.
    for d = 1:num_cPdist % loop thru different d values
        %errorbar([1:numel(rM)]+0.1*d-0.25,pz3(d,:),rz3(d,:),'LineStyle','none','Marker','+','MarkerSize',12,'LineWidth',4,'Color',cPd_colors{d});
        errorbar([1:numel(rM)]+0.1*d-0.25,pz3(d,:),sz3(d,:),'LineStyle','none','Marker','+','MarkerSize',12,'LineWidth',6,'Color',cPd_colors{d});
    end
    
    plot( [0 numel(rM)+1], [0 0],'k--')
    %axis([0.5 numel(rM)+0.5 1.1*min(min(pz3-1.2*rz3)) 1.1*max(max(pz3+1.2*rz3))])
    axis([0.5 numel(rM)+0.5 1.1*min(min(pz3-1.2*sz3)) 1.1*max(max(pz3+1.2*sz3))])
    set(gca,'XTick',1:numel(rM),'XTickLabel',rM,'FontSize',16,'FontWeight','Bold')
    title([methodStr{A},' ',blur_tit,' : max (across GT) F : relGaussRF'],'FontSize',18,'FontWeight','Bold')
    xlabel(['Network Neighborhood Radius (rM) w/ optimized ks = ',ks{k}],'FontSize',18,'FontWeight','Bold')
    ylabel(['\Delta F_G (\mu , \sigma)'],'FontSize',18,'FontWeight','Bold')
    grid on
    legend({'d_t=0','d_t=1','d_t=1.4','d_t=2'},'Location','Best')
    
    saveGoodImg(H7,[dirPre,'../Documentation/Cosyne_2016/F_maxAccGT_ErrBar_',method{A},blur_tag_M,'_Kur_ksopt_relGaussRF.jpg'],[0 0 0.5 1])
    close(H7)
    
end


%% Visualize for single image examples how things change for Mod_SKHEuc as we vary Rmax and Ks parameters.



if(0)
    % (1). Plot image patch, boundary found by method & best matching Ground Truth
    A = 4; % TM 2D
    h = 6; % relGaussRF, maxGT
    d = 4; % d_t = 2 pix

    for i = 1:500
        
        H2=figure; 
        ha = tight_subplot(numel(ks), numel(rM), [0.01 0.01], [0.01 0.01], [0.01 0.01]);

        hind = 0;
        for k = 1:numel(ks)
            for r = 1:numel(rM)
                hind = hind+1;
                im_Fname = img_ptch_name_struct{r,k,A}(i);

                try
                    % load pb png image and threshold it to show what boundaries let to max F-measure with maxGT criteria.
                    pbDir = ['./output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1_blur_sig1/data/Kur_PIF_Fourier1/',method{A},'/pb_png/rM',rM{r},'/sDInf/sP0p2/NF_60_0/ks',ks{k},'/'];
                    pbFile = [pbDir,im_Fname{:}(1:end-3), '.png'];
                    pb = double(imread(pbFile))/255;
                catch
                    disp('Cant find..')
                    break
                end

                ax = subplot( ha(hind) ); hold on
                imagesc(pb)
                colormap(bone)
                axis square tight ij
                set(gca,'XTick',[],'YTick',[])
                %
                if r==1
                    text(0,size(pb,2),['\color{red} ks = ',ks{k}],'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',18,'FontWeight','Bold')
                end
                if k==1
                    text(0,0,['\color{red} rM = ',rM{r}],'HorizontalAlignment','left','VerticalAlignment','top','FontSize',18,'FontWeight','Bold')
                end

                box(ax,'on');
                    
                if r==rM_max(d,h,A) & k==ks_max(d,h,A)
                    box(ax,'on');
                    ax.XColor='red';
                    ax.YColor='red';
                    ax.LineWidth=10; 
                end
                
            end
        end
        %
        saveGoodImg(H2,[dirPre,'../Documentation/Cosyne_2016/single_pb_ksVrM',num2str(i),'.jpg'],[0,0,0.7,1])
        close(H2)
        
    end
end
                    
                    
