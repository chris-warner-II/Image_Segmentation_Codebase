Greetings. Here is an early sample of code from the first of my two Ph.D projects. Most of this project was done in Matlab with some bash scripting to interface with a computer cluster using sbatch. It is not the cleanest repo because this is a project that I worked on without collaboration and it was an exploratory project that developed organically. The second project I worked on in my Ph.D thesis in the Cell_Assembly_Codebase - done later, in Python, and with collaboration - is a better example of code. The purpose of this repo is to provide a sample of my coding capabilities - not a working project that can be picked up by others without significant effort. 

The paper associated with this work, entitled "A Model for Image Segmentation in Retina" can be found here: https://arxiv.org/abs/2005.02567

The goal of this project was to build a computational model as a proof of concept for visual grouping of information in retina by synchrony of firing in retinal ganglion cells (RGCs). The hypothesis being: A functional connectivity network is constructed on the fly, influenced by image contrast values within local receptive fields of RGCs, and electrical interactions within that network influence fine-timing of firing in RGCs without adjusting their spike rate. We model the fine-time change in cell firing relative to an underlying Gamma band oscillation by phase and mediate phase-phase interactions in the Kuramoto coupled oscillator model. With the phase relaxation, regions in the image with similar pixel contrast values attract each other in phase and fall into synchrony while regions across which contrast statistics differ, will repel in phase and desynchronize. The phase, when visualized, generally smooths out texture within objects and yet maintains crisp lines at boundaries between objects - resulting in a grouping in phase and what looks like a 'cartoonization' of the image. The phase representation can then be multiplexed into retinal spike train without effecting spike rates and a perceptual grouping of image locations or a coarse image segmentation results from this fast computation embedded in the circuitry and activity in retina. In this project, we compute probabilistic object boundaries based on spatial gradients in the phase distribution after Kuramoto relaxation and compare them to human segmentations provided in the Berkeley Segmentation Dataset (BSDS). We show quantitatively, using a precision - recall paradigm, that phase discontinuities across the 2D lattice of oscillators overlap with human drawn boundaries. We explore different ways of constructing phase interaction networks and show that Topographic Modularity performs best.  

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


Some useful functions to start in and some helpful hints:

	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

(1).	edit_all_imageSeg.m - a function that opens relevant m files in matlab editor.	

(2).	edit_all_Kuramoto.m - a function that opens relevant m files in matlab editor. 

(2). 	addProjPaths.m - adds the relevant paths to matlab path. Must run this first.

	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

(1). 	write_sbatch_loop_ImageSeg.m - Performs a parallelized grid search across model hyperparameters on 'Cortex' cluster running sbatch. Input is lists of hyperparameters. Writes two file types to scripts4cluster directory. First type is a bunch of sbatch script text files that can be called from the command line. Each calls the function Loop_ImgSegMethodsD.m with a different set of hyperparameter values. The second file type calls each of those files, launching a set of jobs to the cluster.

(2).	Loop_ImgSegMethodsD.m - Function that takes input that is hyperparameter set from the command line and translates them into inputs for and calls SegmentMethod.m function.

(3).	SegmentMethod.m - This is where the bulk of the calculation happens. Loops over image patches in the corpus and for each, calculates edge weights in network based on the input value for 'method' flags and some of the following functions (calc_weightsB.m, compute_LaplacianB.m, compute_ModularityB.m, compute_AvgAssociationB.m). Given the edge weights in matrix, it will compute eigenvectors and perform the kuramoto coupled oscillator phase relaxation simulation by running main_Kuramoto.m.

	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

(1).	write_sbatch_loop_pb_png_4benchmark.m

(2). 	construct_pb_png_4benchmark.m

	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

(1).	write_sbatch_loop_bench_bsds500.m - similar to other functions beginning write_sbatch_loop_<etc>, this function allows us to perform grid search across hyperparameters in parallel across multiple cluster machines by writing sbatch text files and a corresponding run file to call those files in order, submitting them to the queue. Here we compare boundaries found by humans to boundaries found by the algorithm by using bench_bsds500.m


(2).	bench_bsds500.m - for an input set of hyperparameters, compares probabilistic boundaries derived from spatial gradients in kuramoto phase distribution or reshaped eigenvector to boundaries drawn by humans segmenting BSDS images. Runs boundaryBench.m function in the code repository provided with BSDS, here: (https://www2.eecs.berkeley.edu/Research/Projects/CS/vision/bsds/)




	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

(1).	plot_eval_compare_bsdsBench_EigB.m \\
(2).	plot_eval_compare_bsdsBench_blurB.m \\
(3).	plot_eval_compare_bsdsBench_KurB.m - These functions run over Precision/Recall/F-measure metrics for all image patches for the various hyperparameter values, collect up various statistics into multidimensional arrays so they can be compared, analyzed, plotted and saved.


	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

(1).	visualize_single_imgPtch_compare_all_methods_optimized_Kur.m
(2).	visualize_single_imgPtch_compare_all_methods_optimized_Eig.m
(3).	visualize_single_imgPtch_compare_all_methods_optimized_KurNEig.m - These functions loop through various directories and grab the corresponding results files (some png, some txt) and will display the images and plot the data.  It will output ~500 jpg image  files that one can flip through to compare performance of different methods with optimized parameters and compare them. 


	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

evaluation_bdry_imageB.m


correspond_pixelsB.m


	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 



write_sbatch_loop_ bench_blur_bsds500.m

bench_blur_bsds500.m







	- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

NOT:
	main_ImageSeg.m, calc_weights.m, compute_Laplacian.m, compute_Modularity.m, compute_AvgAssociation.m
