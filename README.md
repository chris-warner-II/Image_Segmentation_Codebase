Greetings. Here is an early sample of code from the first of my two Ph.D projects. Most of this project was done in Matlab with some bash scripting to interface with a computer cluster using sbatch. It is not the cleanest repo because this is a project that I worked on without collaboration and it was an exploratory project that developed organically. The second project I worked on in my Ph.D thesis - done later, in Python, and with collaboration - is a better example of code. The purpose of this repo is to provide a sample of my coding capabilities - not a working project that can be picked up by others without significant effort. 

The paper associated with this work, entitled "A Model for Image Segmentation in Retina" can be found here: https://arxiv.org/abs/2005.02567

:elevator: :baseball: --> The goal of this project was to build a computational model as a proof of concept for visual grouping of information in retina by synchrony of firing in retinal ganglion cells (RGCs). The hypothesis being: A functional connectivity network is constructed on the fly, influenced by image contrast values within local receptive fields of RGCs, and electrical interactions within that network influence fine-timing of firing in RGCs without adjusting their spike rate. We model the fine-time change in cell firing relative to an underlying Gamma band oscillation by phase and mediate phase-phase interactions in the Kuramoto coupled oscillator model. With the phase relaxation, regions in the image with similar pixel contrast values attract each other in phase and fall into synchrony while regions across which contrast statistics differ, will repel in phase and desynchronize. The phase, when visualized, generally smooths out texture within objects and yet maintains crisp lines at boundaries between objects - resulting in a grouping in phase and what looks like a 'cartoonization' of the image. The phase representation can then be multiplexed into retinal spike train without effecting spike rates and a perceptual grouping of image locations or a coarse image segmentation results from this fast computation embedded in the circuitry and activity in retina. In this project, we compute probabilistic object boundaries based on spatial gradients in the phase distribution after Kuramoto relaxation and compare them to human segmentations provided in the Berkeley Segmentation Dataset. We show quantitatively, using a precision - recall paradigm, that phase discontinuities across the 2D lattice of oscillators overlap with human drawn boundaries. We explore different ways of constructing phase interaction networks and show that Topographic Modularity performs best.  






Some useful functions to start in and some helpful hints:

(1). 	write_sbatch_loop_ImageSeg.m - Performs a parallelized grid search across model hyperparameters on 'Cortex' cluster running sbatch. Input is lists of hyperparameters. Writes two file types to scripts4cluster directory. First type is a bunch of sbatch script text files that can be called from the command line. Each calls the function Loop_ImgSegMethodsD.m with a different set of hyperparameter values. The second file type calls each of those files, launching a set of jobs to the cluster.


(2).	Loop_ImgSegMethodsD.m - Function that takes input that is hyperparameter set from the command line and translates them into inputs for and calls SegmentMethod.m function.


(3). SegmentMethod.m - This is where the bulk of the calculation happens. Loops over image patches in the corpus and for each, calculates edge weights in network based on the input value for 'method' flags and some of the following functions (calc_weightsB.m, compute_LaplacianB.m, compute_ModularityB.m, compute_AvgAssociationB.m). Given the edge weights in matrix, it will compute eigenvectors and perform the kuramoto coupled oscillator phase relaxation simulation by running main_Kuramoto.m.


	NOT main_ImageSeg.m !!

(2). 

(3). addProjPaths.m is important to add the relevant paths to matlab. Run this first.






