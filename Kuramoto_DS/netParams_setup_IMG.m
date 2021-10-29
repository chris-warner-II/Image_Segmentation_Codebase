%% Structure of parameters used to construct graph/network from image using Image Segmentation
%  Method (Graph Laplacian, Association, Modularity, etc.)

% Params extracted from network constructed from image using particular SegmentMethod

if ~exist('netParams','var')
    
    netParams.gT = gT;
    % netParams.bD = bD;
    
    netParams.im = im;
    
    try
        netParams.imBlurSig = sigB; % sigma std value on gaussian blurring kernel.
    catch
         netParams.imBlurSigC = sigC;
         netParams.imBlurSigS = sigS;
         netParams.imBlurKrat = Krat;
    end
    
    netParams.Qsparseness = Qsparseness;
    
    % Note: I do not want to save these matrices if they are not sparse
    % (i.e. in Mod_N&G, they are full matrices and the mat files are each
    % 750MB vs ~5MB with sparse matrices).  That fills up the cluster with
    % TB of data when I produce many files.
    if(Qsparseness < 0.5)
        netParams.Wdist = Wdist;
        netParams.Q = Q;
    end
                                      
    netParams.PconnFar = 0;
    netParams.Ndims = size(gT{1});
    netParams.N = numel(gT{1});
    
    for a = 1:numel(gT)
        netParams.C(a) = numel(unique(gT{a}));
    end
    
    netParams.sigP = sigP;
    netParams.sigD = sigD;
    netParams.Rmax = rmax;
    
else
    
    %disp('Not creating: ImgSeg - netParams already exists.')
    %netParams
    
end