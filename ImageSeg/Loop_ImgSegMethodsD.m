function [] = Loop_ImgSegMethodsD(method, radmax, sigpix, imdim,InputImages,imgNums2proc, SigW) %  Kscale,TiScale


% Note: imdim should be a cell where each entry is either [x,y] size of
% sub-image or is the word 'full'.  imdim = {[120,130],'full',...}
% Can loop thru different imdim values in write_sbatch_loop_ImageSeg too.
% So here, it is safe to assume that imdim is a 1D cell containing either
% 'full' or a 2D number vector [x,y].
%
% InputImages = 'BSDS_patch' or 'BSDS_tile' or 'BSDS_full'
%
% This script is supposed to run the image segmentation code for a list of 
% conditions Modularity Topo & Nontopo, MaxEnt & Heuristic, Graph
% Laplacian, Average Association, Normalized GL & AA, thresholding at 0 & 
% mean, allowing ME algorithm to converge & using scaling trick after 1 
% iteration.  And More.
%
% main_ImageSeg(im_st,im_fin,Inpt,ds_fctr,AvgAssoc,Modularity,maxEnt,topo,...
% maskFlg,GraphLap,NormCut,shiftEvals,neg,normlze,proj,meanTHevec,maxMEiter,...
% patches,sigma,sigdist,threshdist)


% EXAMPLE HOW TO RUN: Loop_ImgSegMethodsD(7,1,0.2,{[11,11]},'BSDS_patch',[1,5],0);

% [x,comp] = system('hostname')


%% General Flags ("Hyperparameters") that are the same for all the below function calls. They can be input by user.
%sigpix = [0.1, 0.2, 0.3, 0.4, 0.5];             % can be a vector
sigdist = inf; % [inf, radmax(1)/4];                    % can be a vector (sigdist=inf means 1 within Rmax and 0 beyond)
% Note: if radmax is a vector right now, this will act weird.  It expects a scalar rmax value.

shiftEvals = 1;             % flag to shift eigenvalues of Aux Matrix to be all positive (LEAVE THIS AS 1) ...
                            % Think I only need it for GL & NGL Eigenvectors.  But it does not change anything for anything else.
                            % Adds a positive multiple of identity matrix (diagonal matrix) to Q matrix.  Like adding self couplings that reinforce own phase.
                            
DoKurProcessing = 1; % flag to do kuramoto processing inside SegmentMethod function
DoEigProcessing = 1;
dataFileChk = 1;     % if output mat files exist (for Kuramoto or Eigen) do not rerun if flag=1.

% if ~exist('method','var')
%     method = 1:11;          % can be a vector
% end
% %
% if ~exist('radmax','var')
%     radmax = [1 2 3 4];     % can be a vector
% end
% %
% if ~exist('imdim','var')
%     imdim = [7,9,11,13,15,17,19,21,23,25,27,29,31]; % can be a vector % 
% end
% %
% if ~exist('inDir','var')
%     inDir = 'GaussianBox'; % 'GradientBox'
% end

% flags for saving plots and data
save_Evecs_Mat = 1;



for i=1:numel(imdim)
    
    
    if strcmp(imdim{i},'full')
        disp('Can not shift evals because full image matrices too large.')
        shiftEvals=0;
    end

    %% Just the structure of the calls that go into the main_ImageSeg function
    % SegmentMethod(InputImages,imgNums2proc,imdim,   
    %                   AvgAssoc,GraphLap,normlze,... 
    %                   	Modularity,maxEnt,topo,maskFlg,...
    %                       	NormCut,meanThresh,IsoDiff,...
    %                           	shiftEvals,maxMEiter,...
    %                               	sigpix,sigdist,radMax,...
    %                                       save_Evecs_Mat,dataFileChk, ...
    %                                           DoKurProcessing,DoEigProcessing, ...
    %                                               Kscale,SigW,TiScale)

    %% Threshold input image at mean pixel value.
        if ~isempty(find(method == 1,1))
            disp('Threshold at Mean Pixel Value')
            SegmentMethod(InputImages,imgNums2proc,imdim{i},  0,0,0,   0,0,0,0,  0,1,0,  shiftEvals,0,   sigpix,sigdist,radmax,   save_Evecs_Mat,dataFileChk, DoKurProcessing,DoEigProcessing, SigW); % Kscale,SigW,TiScale
        end
        
        %% Average Association
        if ~isempty(find(method == 2,1))
            disp('Average Association')
            SegmentMethod(InputImages,imgNums2proc,imdim{i},  1,0,0,   0,0,0,0,  0,0,0,  shiftEvals,0,   sigpix,sigdist,radmax,   save_Evecs_Mat,dataFileChk, DoKurProcessing,DoEigProcessing, SigW); % Kscale,SigW,TiScale
        end
        
        if ~isempty(find(method == 3,1))
            disp('Normalized AA')
            SegmentMethod(InputImages,imgNums2proc,imdim{i},  1,0,1,   0,0,0,0,  0,0,0,  shiftEvals,0,   sigpix,sigdist,radmax,   save_Evecs_Mat,dataFileChk, DoKurProcessing,DoEigProcessing, SigW); % Kscale,SigW,TiScale
            % This is not working.  Look at Weights distribution. I think Eig2 works.  Eig1 doesnt (kinda like GLnrm)
        end
        
        %% Graph Laplacian (Need 2nd Eigenvector of Negative GL)
        if ~isempty(find(method == 4,1))
            disp('Graph Laplacian')
            SegmentMethod(InputImages,imgNums2proc,imdim{i},  0,1,0,   0,0,0,1,  0,0,0,  shiftEvals,0,   sigpix,sigdist,radmax,   save_Evecs_Mat,dataFileChk, DoKurProcessing,DoEigProcessing, SigW); % Kscale,SigW,TiScale
        end
        
        if ~isempty(find(method == 5,1))
            disp('Normalized GL')
            SegmentMethod(InputImages,imgNums2proc,imdim{i},  0,1,1,   0,0,0,0,  0,0,0,  shiftEvals,0,   sigpix,sigdist,radmax,   save_Evecs_Mat,dataFileChk, DoKurProcessing,DoEigProcessing, SigW); % Kscale,SigW,TiScale
        end
        
        %% Modularity
        if ~isempty(find(method == 6,1))
            disp('Modularity N&G')
            SegmentMethod(InputImages,imgNums2proc,imdim{i},  0,0,0,   1,0,0,0,  0,0,0,  shiftEvals,0,   sigpix,sigdist,radmax,   save_Evecs_Mat,dataFileChk, DoKurProcessing,DoEigProcessing, SigW); % Kscale,SigW,TiScale
        end
        if ~isempty(find(method == 7,1))
            disp('Modularity SKH Adj')
            SegmentMethod(InputImages,imgNums2proc,imdim{i},  0,0,0,   1,0,1,0,  0,0,0,  shiftEvals,0,   sigpix,sigdist,radmax,   save_Evecs_Mat,dataFileChk, DoKurProcessing,DoEigProcessing, SigW); % Kscale,SigW,TiScale
        end
        if ~isempty(find(method == 8,1))
            disp('Modularity SKH Euc')
            SegmentMethod(InputImages,imgNums2proc,imdim{i},  0,0,0,   1,0,1,1,  0,0,0,  shiftEvals,0,   sigpix,sigdist,radmax,   save_Evecs_Mat,dataFileChk, DoKurProcessing,DoEigProcessing, SigW); % Kscale,SigW,TiScale
        end
        if ~isempty(find(method == 9,1))
            disp('Modularity SKH D&O')
            SegmentMethod(InputImages,imgNums2proc,imdim{i},  0,0,0,   1,0,1,2,  0,0,0,  shiftEvals,0,   sigpix,sigdist,radmax,   save_Evecs_Mat,dataFileChk, DoKurProcessing,DoEigProcessing, SigW); % Kscale,SigW,TiScale
        end
        if ~isempty(find(method == 10,1))
            disp('Modularity ME N&G')
            SegmentMethod(InputImages,imgNums2proc,imdim{i},  0,0,0,   1,1,0,0,  0,0,0,  shiftEvals,0,   sigpix,sigdist,radmax,   save_Evecs_Mat,dataFileChk, DoKurProcessing,DoEigProcessing, SigW); % Kscale,SigW,TiScale
        end
        if ~isempty(find(method == 11,1))
            disp('Modularity ME SKH')
            SegmentMethod(InputImages,imgNums2proc,imdim{i},  0,0,0,   1,1,1,0,  0,0,0,  shiftEvals,0,   sigpix,sigdist,radmax,   save_Evecs_Mat,dataFileChk, DoKurProcessing,DoEigProcessing, SigW); % Kscale,SigW,TiScale
        end
        
        %% Isotropic Diffusion: Just a uniform matrix constrained to have Rmax-sized Neighborhoods.
        if ~isempty(find(method == 12,1))
            disp('Isotropic Diffusion')
            SegmentMethod(InputImages,imgNums2proc,imdim{i},  0,0,0,   0,0,0,0,  0,0,1,  shiftEvals,0,   sigpix,sigdist,radmax,   save_Evecs_Mat,dataFileChk, DoKurProcessing,DoEigProcessing, SigW); % Kscale,SigW,TiScale
        end
    
end

%% Note:  After We Find Eigenvectors, We can project each pixel into high dimensional space by its value in each
% eigenvector (maybe we have to weight each one by its eigenvalue or something)