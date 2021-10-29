
edit edit_all_Kuramoto
%
%edit write_sbatch_loop_Kuramoto_Handcooked % 
%edit loop_main_Kuramoto_HandCookedNetwork  % a wrapper around main_Kuramoto to loop thru different parameter values
edit main_Kuramoto                         % wrapper function that constructs K and w, runs Kuramoto simulation and plots results
%edit main_EvecSeg                          % a function that takes coupling matrix K and computes a segmentation using eigenvector method.

%edit phase_init                            % Initialize oscillator phase based on pixel intensity values

%edit setup_Kuramoto_HandCookedNetwork      % I moved all the defining of runflags & runparams into this script for cleanness sake.
% edit setup_Kuramoto_ImgSeg                 % runflags & runparams set when calling Kuramoto simulation from SegmentMethod function



edit kurflags_setup
edit kurParams_setup
%edit netParams_setup_HCN
edit netParams_setup_IMG
edit netflags_setup_IMG

%edit kuramoto_output_dirs_setup_HCN
edit kuramoto_output_dirs_setup_IMG

%
if(0)
    edit KurOsc_addingFeatures.m               % a script where I am playing around with adding features to oscillators instead of specifying Win/Wout explicitly
    edit explore_transient % script (for now) to look at initial transient dynamics of system to pull out coupling by that
    %                      % Not Finished.  Need to Explore More !!!            
end

if(1)
    edit reduce_footprint_of_ModNG_mat_files
    edit write_sbatch_loop_reduce_footprint_of_ModNG_mat_files
end



%edit Kuramoto      % Kuramoto simulation loop.
edit KuramotoB
% edit Kuramoto2_test % trying out RK4 numerical integration using ODE45 function.

%% From Image Segmentation project (reusing the code).
if(1)
    
    % edit graph_construct_test
    
    % Different Methods to calculate coupling matrix K from image
%     edit calc_weights
%     edit compute_AvgAssociation
%     edit compute_Laplacian
%     edit compute_Modularity
%     edit compute_null_model
    %
    edit calc_weightsB
    edit compute_AvgAssociationB
    edit compute_LaplacianB
    edit compute_ModularityB
    edit compute_null_modelB
    %
%     % Generate or extract image to use to calculate coupling matrix
%     edit GaussianBoxLoopGen         % Box inside field with pixels drawn from 2 different probability distributions.
%         edit GaussianBox            %
%     edit GradientBoxLoopGen         % Box inside field each containing gradient in pixel values going opposite directions.
%         edit GradientBox            %
%     edit HalfSplitLoopGen           % Simple image with half pixels light and half pixels dark.
%         edit HalfSplit              %
%     edit BSDSimsLoopGen             % Grabs patches from Berkeley Segmentation Data Set images
end


%% For Hand Cooked Network, generate K. And generate w's randomly
if(0)
    edit build_K2 % similar to build_K1, except you dont have full connectivity within clusters (nearest-neighbor-like connections)
    edit build_w1 % build up natural frequency distribution randomly.
end

%% Different Phase Interaction Functions
if(0)
    edit PIF_Fourier1 % Original Kuramoto Sine interaction Function
    edit PIF_Fourier2 % Fourier (Sine) basis, adding a 2nd term with desynchronizing effect
    %
    edit PIF_dGauss1  % 1st derivative of Gaussian (sine wave like but with tunable concentration)
    edit PIF_dGauss3  % 3rd derivative of Gaussian - similar to 1st but with 4 humps
    edit PIF_dGauss5  % 5th derivative of Gaussian - similar to 1st but with 6 humps
    edit PIF_dGauss7  % 7th derivative of Gaussian - similar to 1st but with 8 humps
    %
    edit PIF_VonMises % A Gaussian-like curve that is tunable by shape parameter
    edit PIF_MexHat   % Center-Surround PIF.  Not sure it make sense actually.
    %
    edit pick_PIF     % function to run correct PIF function given a plag

    edit vizPhaseLock % function to visualize phase of different oscillators
    edit viz_PIF      % function to plot the Phase Interaction Function
end

%% Different Clustering-Quality Metrics
if(0)
    
    %edit metaClusterAnalysis  % something like relative margin to do instead of clustering!
    %edit metaClusterAnalysisB % advancement.  Obselete.
    
    edit metaClusterAnalysisC_circ % splitting circular for phase variables.
    edit metaClusterAnalysisC_linr % splitting linear for eigen & image pix.
    
    edit metaClusterVisualize % takes data file output from main_Kuramoto and plots phase and metric vs. T for each run.
    edit parseSeparationStruct % function to take in Separation Struct in MetaCluster dir that was produced by metaClusterVisualize
                                % and make some analysis plots that show how Separation metric changes along axes of parameter space
end    

%% Different plottings to visualize Separation measure for different runs
if(0)
    
    edit plot_Separation_vs_EvecVizNonLin
    % edit plot_Separation_vs_netParams % old version
    
%     edit write_sbatch_loop_explore_SepVPar
%     edit explore_Separation_vs_Parameters
%     edit write_sbatch_loop_SepVPar_avgOvrImgPatches
%     edit exp_SepVPar_avgOvrImgPatches
%     edit cleanup_Kur_metaSummary_files % obselete. Replaced by explore_compare_DivMarg
    
    % edit explore_compare_DivMarg % this replaces several of the 5 files directly above.
    edit explore_compare_DistPairs
    edit scatter_plot_DistPairs
    edit plot64patches
    edit write_sbatch_loop_exp_comp_DistPairs
    edit write_sbatch_loop_exp_comp_DistPairsB
%     edit write_sbatch_loop_exp_comp_DivMarg
%     edit write_qsub_loop_exp_comp_DivMarg

    % edit find_range_of_variable_Loop % specifically for eigenvector, we are saying average range is 2/sqrt(N).  But I think this is an underestimate and we are giving eigenvectors an unfair advantage.
    edit compute_Dout_ideal
    edit write_sbatch_loop_compute_Dout_ideal
    edit combine_Dout_ideal
    edit calculate_Pmetric
    
end 


if(0)    
    edit compute_Eig_clusterPairROC_fixFiles
    edit write_sbatch_loop_Eig_clusterPairROC
    %edit write_sbatch_loop_Eig_cluster_meanNvarB
    edit compute_Kur_clusterPairROC_fixFiles
    edit write_sbatch_loop_Kur_clusterPairROC  
    %edit calc_ClusterMnNVars
    %edit calc_ClusterHitsNFas
    %edit calc_ClusterPairROC
end

    
    
    
if(1)    
    edit calc_ClusterPairROC_B
    edit calc_gTclusterSizeDist
    
    edit gen_ROC_AUC_scatter
    edit write_sbatch_loop_gen_ROC_AUC_scatter
    
    edit plot_Isotropic_Evecs
    edit plot_Isotropic_Kuramoto

    edit compare_methods_boundaryGradientMetric_Kur
    edit compare_methods_boundaryGradientMetric_Eig
    edit write_sbatch_loop_methods_boundaryGradientMetric_Kur
    edit write_sbatch_loop_methods_boundaryGradientMetric_Eig
    %
    edit compare_allParams_boundaryGradientMetric_Kur
    edit compare_allParams_boundaryGradientMetric_Eig
    edit write_sbatch_loop_allParams_boundaryGradientMetric_Kur
    edit write_sbatch_loop_allParams_boundaryGradientMetric_Eig
    %
    edit compareTotals_methods_boundaryGradientMetric_Kur
    edit compareTotals_methods_boundaryGradientMetric_Eig
    edit optimizeParameters_boundaryGradientMetric_Kur
    edit optimizeParameters_boundaryGradientMetric_Eig
    edit computeDiscriminability
    
    edit compute_Kur_bDgradMetric_fixFiles
    edit compute_Eig_bDgradMetric_fixFiles
    
    edit write_sbatch_loop_Kur_bDgradMetric_fixFiles
    edit write_sbatch_loop_Eig_bDgradMetric_fixFiles
    
    
    edit compute_Spatial_Gradient
    edit gradientB
    %edit compute_BoundaryGradientMetric
    %edit compute_BoundaryGradientPrecRec
    
    
    edit combine_boundaryGradientMetric_Kur
    edit write_sbatch_loop_combine_bGMK
    
    edit combine_boundaryGradientMetric_Eig
    edit write_sbatch_loop_combine_bGME
    
    edit addBoundaryToGT
    
    edit check_GT_disagreement

    %edit script2troubleshoot_AUC_ROC_results
    
    edit visualize_fullImage_Kuramoto_results
    %edit visKurPhase_inBone
    edit visKurPhase_inHSV
    edit plot_phaseAtClk_vs_Kscale
    
    
end
    
if(0) 
    edit visualize_HitsMetric
    
    
    edit compute_SensitivityIndex
    edit truncNormStats_MLEest
    % edit d_prime_analysis

    edit plot_RateDistCurves
    
    
    edit plot_2x2Dscatter_avgOvrImgPatches
    %edit visualize_pts_from_2x2D_scatter
    edit visualize_pts_from_2x2D_scatterB
    edit write_sbatch_loop_visualize_pts_from_2x2DscatterB
    
    edit compare_methods_mean_statistics
    edit plot_netMethods_vs_segMethods_best_perf
    
    edit histDistNorm
    
    edit makeMovie_phaseEvolution % self-contained function to make a video of phase relaxation dynamics of oscillator network.
                                  % will run this from visualize_pts_from_2x2D_scatter I think
    
    
    
    edit plot_DivMarg
    edit plot_DivMarg_Distrib
    
    edit compare_different_tscale
end



% Plot Statistics of Quality of Segmentation/Percept of Thresholded Cosine Distance of Oscillator System at different time points
if(0)
    edit analyze_cosDistSeg_stats % script (for now) to analyze statistics of thresholded pairwise cosine distance between oscillators
    edit plot_Ebars_KurHCN        % subfunction called in analyze_cosDistSeg_stats to make a subplot with errorbars and colors
    edit bar_KurHCN               % subfunction called in analyze_cosDistSeg_stats to make a subplot with bars and colors
end




if(0)
    
    edit ROCcompare_EigVKur             % Function to go through output data file saved from main_Kuramoto & main_EvecSeg
    edit write_sbatch_loop_ROCcompare_EigVKur  % Script to write sbatch files to call ROCcompare_EigVKur singularly for each directory (combination of N,C,Rmax,Pfar)
    edit ROC_KurEig_Stats_Stitch        % Script to take output matfiles from multiple runs using two functions above and stitch them together into one data structure.
    
    % edit calc_ROC_AUC         % Calculate Area Under ROC Curve for collection of parameter settings and return performance values
    edit plot_ROC_AUC_statsB  % Plot Segmentation Performance & Speed Measured from AUC
    % Note: plot_ROC_AUC_stats I left for posterity because I will change statsB a lot I think.
    edit plot_AUC_RmaxVsC_B     % Making plots to explore interaction between number/size of clusters (keeping N constant) and Rmax
    % Note: again I left plot_AUC_RmaxVsC for posterity.  I saved a copy as B and I am editing it pretty severely.

    edit plot_Evec_vs_params % script to take in data file saved from loop_main_Kuramoto that will produce eigenvector as function of parameter plots

end



%% Scripts and functions I am no longer using...
if(0)
    edit illustrate_ROC_calc % script to make plots that will illustrate ROC curve calculation and how it indicates quality of segmentation
                             % I am not using this anymore because it doesnt make sense to think about ground truth and ROC analysis of coupling matrix K.
    edit relativeMargin % a clustering-quality metric (ratio of distance between  oscillators within cluster vs. distance across clusters)
    edit Evaluate       % code downloaded from internet to do precision, recall, etc ...
end




edit makeMovie_KurPhaseRelax
edit makeMovie_KurPhaseRelax_compareKS
