function edit_all_ImageSeg()


%% This script to view functions for ImgSeg Project    
if(1)
    edit edit_all_ImageSeg
end

%% Functions to generate input image ensembles
if(0)
    edit GaussianBoxLoopGen         % Box inside field with pixels drawn from 2 different probability distributions.
        edit GaussianBox            %
    edit GradientBoxLoopGen         % Box inside field each containing gradient in pixel values going opposite directions.
        edit GradientBox            %
    edit HalfSplitLoopGen           % Simple image with half pixels light and half pixels dark.
        edit HalfSplit              %
    edit BSDSimsLoopGen_RandomPatch % Grabs patches randomly from Berkeley Segmentation Data Set images
    edit BSDSimsLoopGen_TilePatches % Grabs patches randomly from Berkeley Segmentation Data Set images
    
    edit stitchBSDSTilePatchesResults
        
end    


edit visualize_fullImage_Spectral_results



%% Functions to compute Image Segmentation
if(1)

    edit write_sbatch_loop_ImageSegB
    %edit write_sbatch_loop_ImageSeg
    %edit write_qsub_loop_ImageSeg
    edit Loop_ImgSegMethodsD
    edit SegmentMethod
    
    edit calc_weightsB
    edit compute_ModularityB
    edit compute_null_modelB
    edit compute_LaplacianB
    edit compute_AvgAssociationB

    %edit optimize_threshold2.m
end

%% Analysis of Gaussian Box segmentation results
if(0)
    edit SegMethodAnalyze_GaussianBox
    edit SegMethAnal_GaussianBox_Plots
    edit SegDataViz1
    edit SegDataViz2
    edit SegDataViz3
    edit SegDataViz3b
    edit SegDataViz4
    edit SubplotEvecNSeg
    edit onCluster % small function to check if on cluster and change save dir appropriately
end

%% Analysis of Gradient Box Segmentation results
if(0)
    % edit SegMethodAnalyze_GradientBox
    edit SegMethodAnalyze_GradientBox2
end

%% Analysis of BSDS Patches Segmentation results
if(0)
    edit SegMethodAnalyze_BSDS
end



%% Visualize Eigenvectors
if(0)
    edit Loop_plot_Evecs
    edit plot_ImgSeg_Evectors
    edit EvecVizF
    edit parse_Fname
end

%% Hillar's Maximum Entropy Code
if(0)
    edit compute_NM_MaxEnt
    edit nodePots2edgeWaits
    edit findMLE
    edit expectedDegree
end


%% Visualize Pixel Intensity Value vs Distance Plots for
if(0)
    edit PIV_vs_DiagDist
    edit PIV_vs_MaskDist
    %
    edit PIV_vs_DiagDistB
    edit PIV_vs_MaskDistB
end

%% Script to look at AA, GL and their normalized versions
if(0)
    edit troubleshoot_AAnGL.m
end

%% CODE TO DO SOMETHING...
if(0)
    edit extractPatch
    
    edit addBoundaryToGT
    
    edit tilePatches
    
    edit power_method
    edit Project_1s
end


%% to edit benchmark code included in Malik's BSDS package
if(1)
    
    [dirPre,sizeGoodIm] = onCluster;
    addpath([dirPre,'images/BSDS_images/BSR/bench/benchmarks/'])
    
    edit bench_bsds500
    edit write_sbatch_loop_bench_bsds500
    
    %edit bench_blur_bsds500
    %edit write_sbatch_loop_bench_blur_bsds500
    % Note: These are handled in bench_bsds500 with certain method if statments.
    
    edit allBench
    %edit test_benchs
    
    edit regionBench
    edit evaluation_reg_image
    edit collect_eval_reg
    
    edit boundaryBench
    %edit evaluation_bdry_image
    edit evaluation_bdry_imageB
    edit collect_eval_bdry
    
    edit match_segmentations
    edit match_segmentations2

    edit plot_bench_statsB
    edit calc_info

    edit plot_eval
    %edit plot_eval_compare_bsdsBench_Kur
    edit plot_eval_compare_bsdsBench_KurB
    %edit plot_eval_compare_bsdsBench_Eig
    edit plot_eval_compare_bsdsBench_EigB
    %edit plot_eval_compare_bsdsBench_blur
    edit plot_eval_compare_bsdsBench_blurB
    %
    edit plot_results_Kur_N_Eig_together
    
    edit visualize_single_imgPtch_compare_all_methods_optimized_Kur
    edit visualize_single_imgPtch_compare_all_methods_optimized_Eig
    edit visualize_single_imgPtch_compare_all_methods_optimized_KurNEig
    
    edit plot_single_imgPtch_compare_ImPix_Iso_SK_optimized
    edit plot_PRdelF_vector_field_Iso_SK_optimized
    
    edit construct_patch_groundTruth_4benchmark
    edit construct_pb_png_4benchmark
    
    edit write_sbatch_loop_pb_png_4benchmark
    
    edit correspondPixelsB
    
    
    
end



    edit make_gaussian_blur_pb_png_img_patches
    
    edit construct_gaussian_kernel
    edit construct_DoG_kernel




%% BSDS image & groundtruth access functions
if(0) 
    edit GtFilename   
    edit bsdsRoot
    edit imgFilename
    edit imgList
    edit imgRead
    edit readSeg
    edit readSegs
    edit seg2bmap
    edit segFilename
    edit writeSeg
end

