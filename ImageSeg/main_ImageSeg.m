function main_ImageSeg(im_st,im_fin,Inpt,ds_fctr,patches,sigma,sigdist,radMax,distBinary,...
                        AvgAssoc,Modularity,maxEnt,maskFlg,GraphLap,NormCut,meanThresh,...
                        shiftEvals,neg,normlze,proj,topo,maxMEiter,evectors)

% main_ImageSeg(im_st,im_fin,Inpt,ds_fctr,patches,sigma,sigdist,radMax,distBinary,...
%                        AvgAssoc,Modularity,maxEnt,maskFlg,GraphLap,NormCut,...
%                        shiftEvals,neg,normlze,proj,topo,maxMEiter,evectors)
% This function is where i make a test image and then run it thru image segmentation algorithm.
%
% Currently, Graph Laplacian and Negative Graph Laplacian make sense.  The
% smallest and largest respectively eigenvectors are the all 1's vector and
% the second smallest/largest eigenvector provides a good segmentation.
%
% The normalized Graph Laplacian and negative normalized one do not make
% sense yet.  I expect the smallest and largest respectively eigenvectors
% to be the all 1's vector, but they are not.  I guess they are close, but
% there is some useful edge information in there.  And the 2nd smallest /
% largest eigenvectors are useless (all near zeros).



%% TODO:
%
% (1). Implement Projection Matrix in Graph Laplacian Calc                 [X]
% (2). Implement Normalized Avg Association                                [X]
% (3). Make sure Normalized GL & NGL are working properly                  [O]
% (3b.) Neg Norm GL Sparse Eigenvectors look funny...  
% (10).  Figure out why shifting is changing Eigenvectors (it shouldnt)    [ ]
% (4). Implement a way to visualize Power Method Dynamics (movie)          [X]
% (x). Look at Normalized Cut Algorithm to understand what it does.        [ ]
% (9). Compare smallest evec for Lap and largest for negLap.  One          [ ]
% subtracted from other should be zero.  Also, compare projected negLap
% largest evec to 2nd smallest of normal Lap.  Should be near zero diff.
% (2b). Implement Normalized Modularity if it exists (?)                   [ ]
% (5). Implement new test images (below) and test them                     [ ]
% (6). Compare different methods' performances on test images              [ ]
% (7). Compare different methods' performances on BSDS                     [ ]
% (8). Look at 8NN performance or 12 or ...                                [ ]
% (K). Plug in Hillar's MaxEnt code for Modularity Null Model              [ ] 

%% Random Number Generator for Repeatability. (Used in Input 1 - not sure if anywhere else)
rng(1234567); % Random Number Generator Seed.  Rand used in Ising Model.
% rng(316);

%% Flags of which method to run
if ~exist('AvgAssoc','var')
    AvgAssoc = 0;  % Just the W matrix.     % Largest Eigenvector
end


if ~exist('Modularity','var')
    Modularity = 1;      % Q = W - NM       % Largest Eigenvector
end 
if ~exist('maxEnt','var')
    maxEnt = 1;          % flag {0,1}  if 1, use MaxEnt NM. If 0, one of N&G approximations.
end
if ~exist('maxMEiter','var')
    maxMEiter = 1;       % [real integer] maximum number of iterations if MaxEnt algorithm doesnt converge.
end                      % (to look at RowSums scaling addition - may only need 1 iteration)
if ~exist('topo','var')
    topo = 1;            % flag {0,1} if 1, use topological or sparse null model.
end
if ~exist('maskFlg','var')
    maskFlg = 2;         % flag {0,1,2} 0 = Adj. Matrix Diagonal Distance PIV Means only  (BAD)
end                      % 1 = Img Euclidian Distance PIV dependent Means mask (MaskR)
                         % 2 = Img Distance & Orientation dependent PIV Means
                         %     mask (MaskDO) - Considers vertical & horizontal separately    
                         % Different Masks don't matter for MaxEnt Null Model calc.

if ~exist('GraphLap','var')
    GraphLap = 0;        % L = D - W          % 2nd smallest Eigenvector
end
if GraphLap
    neg = 1;        % flag {0,1} to use negative matrix caluclated by method (Q -> -Q)
    proj = 1;       % flag {0,1} to project out the Largest Eigenvector
else
    neg = 0;        
    proj = 0;
end
% if ~exist('neg','var')
%     neg = 1;        % flag {0,1} to use negative matrix caluclated by method (Q -> -Q)
% end
% if ~exist('proj','var')
%     proj = 1;       % flag {0,1} to project out the Largest Eigenvector 
% end                 % (for Neg Laplacian Dom Evec = 1's with Eval = 0 and rest <0)


if ~exist('NormCut','var')
    NormCut = 0;         % Run Shi-Malik Normalized Cut Algorithm (should verify this) 
end


if ~exist('normlze','var')
    normlze = 0;    % flag {0,1} to normalize matrix by degree of vertices
end
if ~exist('shiftEvals','var')
    shiftEvals = 0; % flag {0,1} to shift eigenvalues to be positive (does not matter if using MATLAB's eig)
end

if ~exist('meanThresh','var')
    meanThresh = 0; % flag {0,1} to Threshold input image at mean value (as a control).
end


%% Strings to tell user what specific method is was used in calculation
str = '';

if(shiftEvals) 
    str=[str,' Shft'];   % Descriptive strings are unnecesarily long.
end
if(normlze) 
    str=[str,' Norm'];
end
% if(proj)               % These are only and always with Graph Laplacian.
%     str=[str,' Proj'];
% end
% if(neg) 
%     str=[str,' Neg'];
% end
str = [str,' '];

 %% Images to run through the image segmentation code.
if ~exist('Inpt','var')
    Inpt = 2; % Can be {1:6} to choose differnt input images: see notes below
end
if ~exist('ds_fctr','var')
    ds_fctr= 1; % 8 or 10 works for BSDS.  5 for Kyoto. 10 for Waddle
end
if ~exist('patches','var') 
    patches=0; % a flag to extract patches from images instead of downsampling
end
if(patches && ~exist('numpatches','var') )
    numpatches = 3; % CHANGE THIS BACK!! TO WHAT!?!
end
if(~patches)
    numpatches = 1;
end


%% Set constants for defining weights matrix 
if ~exist('sigma','var')
    sigma = [0.4]; % spread on gaussian that determines weights from pixel intensity values (originally 0.1 - seems narrow!)
end
lambda = 1;  % weighting on null model for Modularity method. (Keep constant = 1.)
%
%
% These change size of radius of pixel interaction from Nearest Neighbor to
% all (using hard threshold or soft gaussian-like decay)
% Notes: [sig,th] = [1,0.55] for 4NN.
if ~exist('sigdist','var')
    sigdist = 1;  %[0.5, 1, 1.5]; % spread of gaussian distribution around difference of pixels intensity values (originally 1 - seems good)
end          % image. Determines smoothness of nonlinearity  Note:  Small sigma -> binary weights.
if ~exist('distBinary','var')
    distBinary = 0; % flag {0,1} - 1 to have distance dependence be binary, unweighted, 1 or 0
end                 %              0 to have distance dependence be real valued with a possible rmax cutoff
if ~exist('radMax','var')
    radMax = [1]; %[2,3,inf]; % euclidian maximum radius of pixels that can be connected. Vector that gets looped through below.
end
% passed through this thresholding step to produce binary weights (either 1 or 0) with
% weight = 1 for pixels close enough for tuned values of sig and thresh and 0 for 
% distances further than the value determined to be cutoff by tuning of the 2 values.

if ~exist('im_st','var')
    im_st = 1; 
end
if ~exist('im_fin','var')
    im_fin = 1;
end

if ~exist('evectors','var') % which eigenvectors to look at 1 for all but NGL
    evectors = [1];               
    %
    if (GraphLap && neg && proj) % 2 for Neg Graph Lap.
        evectors = [2];
    end
end % single value.  Change looping & saving code to handle multiple Evecs (easy addition).               

% flags for displaying and saving plots and data
save_Evecs = 1;
vizSlope = 1e-10; % slope of nonlinearity used to visualize Eigenvectors in EvecVizF (very small will binarize Evec) (used 1e-3 usually)
save_Mat = 1;
save_im = 1;     % plot and save only full image
save_im_sub = 0; % plot and save image and downsampled image or patch in image


%% Size of Rectangles for test images
if(Inpt>1 && Inpt <5)
    
%     % BIG IMAGE
%     ximg = 35; % size of image to be generated (Horizontal axis in imagesc)
%     yimg = 25; % (Vertical axis in imagesc)
% 
%     x2 = 10; % start and end of x&y of inner rectangle
%     x3 = 25;
%     y2 = 7;
%     y3 = 15;
    
    % SMALL IMAGE
    ximg = 20; % size of image to be generated (Horizontal axis in imagesc)
    yimg = 20; % (Vertical axis in imagesc)
    
    x2 = round(0.3*ximg); % start and end of x&y of inner rectangle
    x3 = round(0.7*ximg);
    y2 = round(0.3*yimg);
    y3 = round(0.7*yimg);
    
    inptParams.ximg = ximg; 
    inptParams.yimg = yimg; 
    inptParams.x2 = x2;
    inptParams.x3 = x3;
    inptParams.y2 = y2;
    inptParams.y3 = y3;
    
end

inptParams.Inpt = Inpt;

%% Test Image 1: Simplest Test I can think of (for sanity checks)
if(Inpt==1)
    disp(['Input Image is ',num2str(ximg),'x',num2str(yimg),' random box for sanity check.']);
    iids = struct('name','randBox    '); % note: spaces after are important
    im = rand(ximg,yimg);
    inptParams.groundTruth = 'There is none.';
    dirOut = ['rand_',num2str(ximg),'x',num2str(yimg),'_Box'];
    
end

%% Test Image 2: Binary rectangle inside an image with opposite background
if(Inpt==2)
    disp(['Input Image is a ',num2str(ximg),'x',num2str(yimg),' black(white) box in white(black) background.']);
    iids = struct('name','bwBox    '); % note: spaces after are important
    im = binary_box(ximg,yimg,x2,x3,y2,y3,0,0); % a box inside image
    inptParams.groundTruth = im; % where box actually is.  To be compared to segmentation results.
    dirOut = ['bw_',num2str(ximg),'x',num2str(yimg),'_Box'];
end

%% Test Image 3: Greyscale gradient rectangle inside an oppositely travelling greyscale gradient)
if(Inpt==3)
    disp(['Input Image is ',num2str(ximg),'x',num2str(yimg),' box with greyscale gradient inside and opposing gradient outside.']);
    iids = struct('name','gradBox    '); % note: spaces after are important
    outbeg = 0;
    outfin = 1;
    inbeg = 0;
    infin = 1;

    im = gradient_box(ximg,yimg,x2,x3,y2,y3,outbeg,outfin,inbeg,infin,0);
    dirOut = ['grad_',num2str(ximg),'x',num2str(yimg),'_Box'];
    
    inptParams.groundTruth = binary_box(ximg,yimg,x2,x3,y2,y3,0,0); % binary ground truth of where box is.
    inptParams.outbeg = outbeg;
    inptParams.outfin = outfin;
    inptParams.inbeg = inbeg;
    inptParams.infin = infin;
    
end


%% Test Image 4: Two separate Gaussian distributions of pixel intensity for inside and outside.
if(Inpt==4)
    disp(['Input Image is ',num2str(ximg),'x',num2str(yimg),' box with pixel intensity vals distribution different from surrounding.']);
    iids = struct('name','gaussBox    '); % note: spaces after are important
    mu1 = 0;       % Mean of Distribution of inside box
    sig1 = 0.01;      % Spread of Distribution of inside box
    mu2 = (1-mu1);   % Mean of Distribution of outside box
    sig2 = sig1;     % Spread of Distribution of outside box
    
    im = prob_dist_box(ximg,yimg,x2,x3,y2,y3,mu1,sig1,mu2,sig2,0);
    dirOut = ['gauss_',num2str(ximg),'x',num2str(yimg),'_Box_',num2str(mu1,'%0.5g'),'_',num2str(sig1,'%0.5g'),'_',num2str(mu2,'%0.5g'),'_',num2str(sig2,'%0.5g')];
    dirOut = dirOut(dirOut~='.');
    
    inptParams.groundTruth = binary_box(ximg,yimg,x2,x3,y2,y3,0,0);
    inptParams.mu1 = mu1;
    inptParams.sig1 = sig1;
    inptParams.mu2 = mu2;
    inptParams.sig2 = sig2;
    
    
    % INSIDE BOX HAS SAME MEAN BUT MUCH TIGHTER VARIANCE THAN OUTSIDE...
    % INSIDE BOX MOVES AROUND INSIDE OUTSIDE BOX TO INVESTIGATE EDGE EFFECTS...
    % MORE COMPLICATED EDGES AND CONTOURS THAN JUST BOX... (Curvature effects segmentation?)
    
end


%% Test Image 5++: Others & Future Ideas.
if(Inpt==5)
    disp('Not A valid input right now: Write Something...');
%
% (1). A gradient that runs diagonally (not vertical or horizontal)
% (2). Shapes other than rectangles.  Circles, Triangles, More complicated.
% (3). Some way to implement texture (Conditional probabilities - MRF's or CRF's)
end

%% Test Images 6: Use Berkeley Segmentation Data Set Images
if(Inpt==6)
    disp('Input Images are from BSDS500');
    iids = imgList('all'); % can be {'all','test','train','val'}  
    dirOut = 'BSDS';

    
%     
%     % code to find a given file name (where it is in the iids struct)
%     for i=1:numel(iids) 
%         if ~isempty(strfind(iids(i).name,'100075'))
%             i
%             iids(i).name
%         end
%     end
    
    iids_subset = [im_st:im_fin];
    iids = iids(iids_subset);
    inptParams.groundTruth = 'Given in BSDS data';
    
end

%% Test Image 7: Read in an image directly from its directory & filename
if(Inpt==7)
    disp('Loading Kyoto Image');
    iids = struct('name','Kyoto    '); % note: spaces after are important
	x = imread('./images/kyoto','jpg');
	im = double(x(:,:,1)); % make image grayscale.
    im = im./max(im(:)); % Normalize pixels in image to be between 0 & 1.
    dirOut = 'Kyoto';
    inptParams.groundTruth = 'Just Eye It up';
end


%% Test Image 8: Read in an image directly from its directory & filename
if(Inpt==8)
    disp('Loading Gray Waddle Image');
    iids = struct('name','waddle    '); % note: spaces after are important
    iids(2).name = 'waddle    ';
    iids(3).name = 'waddle    ';
	x = imread('./images/waddle','jpg');
	im = double(x(:,:,1)); % make image grayscale.
    im = im./max(im(:)); % Normalize pixels in image to be between 0 & 1.
    dirOut = 'Waddle';
    inptParams.groundTruth = 'Just Eye It up';
    
end

%% Create a structure called runFlags to save in mat file.
runFlags.im_st = im_st;
runFlags.im_fin = im_fin;
runFlags.ds_fctr = ds_fctr;
runFlags.patches = patches;
runFlags.distBinary = distBinary;
runFlags.AvgAssoc = AvgAssoc;
runFlags.Modularity = Modularity;
runFlags.maxEnt = maxEnt;
runFlags.maskFlg = maskFlg;
runFlags.GraphLap = GraphLap;
runFlags.NormCut = NormCut;
runFlags.shiftEvals = shiftEvals;
runFlags.neg = neg;
runFlags.normlze = normlze;
runFlags.proj = proj;
runFlags.topo = topo;
runFlags.maxMEiter = maxMEiter;
runFlags.evectors = evectors;
%
runFlags.save_Evecs = save_Evecs;
runFlags.vizSlope = vizSlope; 
runFlags.save_Mat = save_Mat;
runFlags.save_im = save_im;     
runFlags.save_im_sub = save_im_sub; 



%% Make a directories to put output images and data in
if ~exist(['./output/ImgSeg/simpleExamples/',dirOut],'dir')
    mkdir(['./output/ImgSeg/simpleExamples/',dirOut])
end
%
if ~exist(['./output/ImgSeg/simpleExamples/',dirOut,'/data from main'],'dir')
    mkdir(['./output/ImgSeg/simpleExamples/',dirOut,'/data from main']);
end
%
if ~exist(['./output/ImgSeg/simpleExamples/',dirOut,'/pics from main'],'dir')
    mkdir(['./output/ImgSeg/simpleExamples/',dirOut,'/pics from main']);
end
%
if ~exist(['./output/ImgSeg/simpleExamples/',dirOut,'/pics from main/evecs'],'dir')
    mkdir(['./output/ImgSeg/simpleExamples/',dirOut,'/pics from main/evecs']);
end

%% Plot some Gaussian Distributions for PIV of input and/or Weight dependencies on PIV & Distance
if(Inpt==4 && 1) % Plot 2 Gaussian Distributions of pixel values to view overlap. 
    x = linspace(0,1);
    G1 = (1/(sig1.*sqrt(2*pi))).*exp( -(x-mu1).^2 ./ (2.*sig1.^2) );
    G2 = (1/(sig2.*sqrt(2*pi))).*exp( -(x-mu2).^2 ./ (2.*sig2.^2) );
    G1 = G1./sum(G1); G2 = G2./sum(G2);
    hgd = figure; hold on, plot(x,G1,'r','LineWidth',2), plot(x,G2,'b','LineWidth',2)
    title(['Input PIV Distrib: \mu=[',num2str(mu1),',',num2str(mu2),'] \sigma=[',num2str(sig1),',',num2str(sig2),']'],'FontSize',20,'FontWeight','Bold');
    xlabel('Pixel Intensity Value','FontSize',16,'FontWeight','Bold')
    legend('Inside rect','Outside rect','Location','North')
    saveas(hgd,['./output/ImgSeg/simpleExamples/',dirOut,'/pics from main/evecs/',dirOut],'jpg');
    close
    
    %KL Divergence as measure of difference between distributions
    D_KL = sum(log(G2./G1).*G2); % = 0 if same; = large if different
    inptParams.D_KL = D_KL;

end


if(0) % Plot distance dependence gaussian and pixel itensity difference gaussians (for development & troubleshooting)
    h=figure;hold on,
    colour = ['rgkcmbrgkcmbrgkcmbrgkcmbrgkcmbrgkcmbrgkcmbrgkcmb'];
    
    % (1). plot Weight vs. Distance between pixels
    Xd = linspace(1,radMax);
    subplot(211), hold on
    for i=1:numel(sigdist)
        Wd = exp( -( Xd - 1 ) ./ (2*sigdist(i).^2) );
%         Wd(Xd>radMax)=0; % pixels beyond maximum radius have zero weight.
        plot(Xd,Wd,colour(i),'LineWidth',2)
        leg{i} = ['\sigma = ',num2str(sigdist(i))];
    end  
    legend(leg)
    xlabel('Pixel Separation, Euclidian Distance','FontSize',16,'FontWeight','Bold')
    ylabel('Weight Strength','FontSize',16,'FontWeight','Bold')
    title('Pixel Distance','FontSize',20,'FontWeight','Bold')
    text(4,0.8,['rmax = ',num2str(radMax)],'FontSize',16,'FontWeight','Bold')
    grid on

    % (2). plot Weight vs. Difference in pixel value
    Xd2 = linspace(0,1);
    subplot(212),hold on
    for i=1:numel(sigma)
        Wd2 = exp( -( Xd2 ) ./ (2*sigma(i).^2) );
        plot(Xd2,Wd2,colour(i),'LineWidth',2)
        leg2{i} = ['\sigma = ',num2str(sigma(i))];
    end
    legend(leg2)
    xlabel('Difference in Pixel Intensity Values','FontSize',16,'FontWeight','Bold')
    ylabel('Weight Strength','FontSize',16,'FontWeight','Bold')
    title('Pixel Similarity','FontSize',20,'FontWeight','Bold')
    grid on
    set(h,'Position',[1 1 1280 700]);
    saveas(h,['./output/ImgSeg/simpleExamples/',dirOut,'/Adjacency_Matrix_Gaussians'],'jpg');
    
    % (3). plot a color image in 2D of distance vs. pixel diff vs. weight
    figure, imagesc(Xd2,Xd,Wd'*Wd2),colorbar
    ylabel('Distance between Pixels','FontSize',16,'FontWeight','Bold')
    xlabel('Pixel Value Difference','FontSize',16,'FontWeight','Bold')
    title('Weight Strength','FontSize',20,'FontWeight','Bold')
    
    keyboard

end


%% Loop through input images
for iii = 1:numel(iids)
    
    iid = iids(iii).name;
    
    for j = 1:numpatches
    
        % read in image to be segmented if analyzing BSDS image.
        if(Inpt==6)
            try
                fname=iid;
                [im] = imgRead(fname,'gray'); % For running code on my computer
                load(GtFilename(fname(1:end-4)));
            catch
                fname=iid(3:end);
                [im] = imgRead(fname,'gray'); % Cluster does something funny -> ' ._name'
                load(GtFilename(fname(1:end-4)));
            end
            disp(['Image #',num2str(iii),' : ',fname(1:end-4)])
        else
            fname=iid;
    %         im = imm; % a hack because extracting a patch was making that patch replace entire image.
        end

        % extract a random patch from full resolution image (instead of downsampling it)
        if(patches)
            if ~exist('xpat','var')
                xpat=40;
            end
            if ~exist('ypat','var')
                ypat=30;
            end

    %             imm = im; % a hack

            [pach.patch,pach.xpbeg,pach.xpfin,pach.ypbeg,pach.ypfin,pach.pnum] = ...
                extractPatch(im,xpat,ypat,groundTruth,fname(1:end-4),dirOut,j); % extract random patch in image
            im = pach.patch;
            patstr=['patch',num2str(j)];
            pach.xpat = xpat;
            pach.ypat = ypat;
            inptParams.pach = pach;
        else
            patstr=''; % may not need this?
            inptParams.pach='Patch not used';
        end

        % downsample image
        imDS = imresize(im,1/ds_fctr);
        ximg = size(imDS,2);           % horizontal image dimension
        yimg = size(imDS,1);           % vertical image dimension
        disp(['Image size: ',num2str(ximg),'x',num2str(yimg)])

        inptParams.imDS = imDS;
        inptParams.im = im;
        inptParams.ximg = ximg;
        inptParams.yimg = yimg;
        inptParams.fname = fname;

        % display original image and downsampled image
        if(save_im_sub)
            hims = figure; subplot(121), imagesc(im), colormap('bone'), title([num2str(iii),' - ',fname])
            subplot(122), imagesc(imDS)
            set(gca,'xtick',[],'ytick',[]); % axis square
            saveas(hims,['./output/ImgSeg/simpleExamples/',dirOut,'/pics from main/evecs/'...
                            ,fname(1:end-4),' ',patstr],'jpg');
        end
        % display original image only (and save it if it is not already in directory)
        if(save_im)
            him = figure; imagesc(im), colormap('bone'), colorbar('SouthOutside')
            title([num2str(iii),' - ',fname],'FontSize',16,'FontWeight','Bold') 
            set(gca,'xtick',[],'ytick',[]); % axis square
            saveas(him,['./output/ImgSeg/simpleExamples/',dirOut,'/pics from main/evecs/'...
                            ,fname(1:end-4),' ',patstr],'jpg');
        end




        %% Do Network Computation to Segment Input Image
        for jjj = 1:numel(sigma) % loop through sigma dist values - fall off for pixel similarity dependence
        for kkk = 1:numel(sigdist) % loop through sigma dist values - fall off for distance dependence
        for lll = 1:numel(radMax) % loop through values of maximum radius for distance dependence

                sigP = sigma(jjj);
                sigD = sigdist(kkk);
                rmax = radMax(lll);
                
                % make appropriate directories for saving if they do not exist already 
                sD = num2str(sigD); sD(sD=='.')='p';
                sP = num2str(sigP); sP(sP=='.')='p'; % turning sigma and rmax parameters into strings for filenames
                rM = num2str(rmax); rM(rM=='.')='p';
                params = ['sP',sP,' sD',sD,' rM',rM];
                
%                 % get rid of one brightest pixel on bear's toenail.  I think its fucking things up.
%                 [xx,yy] = find(imDS == max(imDS(:)));
%                 imDS(xx,yy) = mean(imDS(:));          % come up with perm soln for this.
%                 % actually, things in Evec look mostly the same.

                % Calculate weights from Pixel Intensity Values and Closeness of Pixels
                disp(['Calculating Weights with sigma pix = ',num2str(sigP),'; sigma dist = ',num2str(sigD),'; rmax = ',num2str(rmax)])
                [W, Wconn, Wunc, Mask] = calc_weights(imDS,sigP,sigD,maskFlg,distBinary,rmax); % note: threshdist nolonger used.  using rmax


                % Method: Graph Laplacian (L = D - W).
                if(GraphLap)
                    method = ['GL'];
                    disp(['Method: ',method])
                    Q = compute_Laplacian(W,normlze,neg);
                end


                % Method: Modularity (Q = W - NM) where NM = null model or expected connectivity.
                if(Modularity)
                    method = ['Mod'];
                    %
                    disp(['Method: ',method])
                    [Q,NMiter] = compute_Modularity(W,Wconn,Wunc,Mask,ximg,yimg,ds_fctr,lambda,normlze,neg,maxEnt,maxMEiter,topo,maskFlg);
                    if(maxEnt)
                        if(topo) % if u want to use full (not sparse) MaxEnt null model
                            method = [method,' SKH ME']; % ,num2str(NMiter),'iter'
                        else
                            method = [method,' N&G ME']; % ,num2str(NMiter),'iter'
                        end
                    else
                        if(topo)
                            method = [method,' SKH'];
                            switch(maskFlg)
                            case(0)
                                method = [method,' Adj']; % diagonal in Adjacency
                            case(1)
                                method = [method,' Euc']; % mask Euclidian dist
                            case(2)
                                method = [method,' D&O']; % mask distance & orientation
                            end
                        else
                            method = [method,' N&G'];
                        end
                    end
                end



                % Method: Average Association (W) - Just the similarity, adjacency, or weights matrix
                if(AvgAssoc)
                    method = ['AA'];
                    disp(['Method: ', method]);
                    Q = compute_AvgAssociation(W,normlze);   
                end


                % Method: Normalized Cut (a la Shi & Malik)
                if(NormCut)
                    disp('That code from Shi''s website isnt working yet, fix it.')
                end
                
                % Method: Threshold at Average Pixel Value
                if(meanThresh)
                    method = ['THmean'];
                    disp(['Method: ',method]);
                    SegAtTh = imDS>mean(imDS(:));
                end
                
                method = [method,str];
                
                Meth.method = method;
                Meth.sD = sD;
                Meth.sP = sP;
                Meth.rM = rM;


                % Shift eigenvalues of Q matrix to make them all positive because power 
                % method pulls out eigenvector with largest absolute value eigenvalue.
                if(shiftEvals)
                    Q = Q + 100.*eye(size(Q)); % add C*(identity matrix) 
                end


                %% Calculate eigenvectors using MATLAB eig function
                if(1 && ~meanThresh)
                    tic
                    % Calculate Dominant Eigenvector using MATLAB's eig function for comparison
                    % This is just distance (Vertical & Horizontal mushed together)
                    % Diagonal Means only
                    [EigVec,EigVal] = eig(full(Q));
                    EigVal = diag(EigVal);
                    EigValMax = find(EigVal == max(EigVal));
                    
                    disp('First & Last Couple Eigenvalues are:')
                    [EigVal(1:5)',EigVal(end-5+1:end)']
                    
                    % plot eigenvectors in evector in one figure to see what info lower eigenvectors hold
                    if(save_Evecs)
                        xsub = round(sqrt(numel(evectors))); %+1 to plot image in too.
                        ysub = ceil(sqrt(numel(evectors)));  %+1
                        hev = figure;
%                         subplot(xsub,ysub,1)
%                         imagesc(imDS)
%                         xlabel(['Img In'],'FontSize',12,'FontWeight','Bold')
%                         title(method,'FontSize',16,'FontWeight','Bold')
%                         set(gca,'xtick',[],'ytick',[]); % axis square
%                         colorbar('SouthOutside')
                    end


                    for k = 1:numel(evectors) % Just showing largest and smallest but can show N largest and smallest

                        % k largest eigenvectors of Full matrix
                        % diagonal means only
                        EvecML = reshape(EigVec(:,size(EigVec,2)+1-evectors(k)),yimg,ximg);    % reshape to image size
                        EvecML = EvecML./max(abs(EvecML(:)));                                  % normalize
 
                        EvecMLs{k} = EvecML; % group eigenvectors together to save in case u are using more than one.


                        %% Find Best threshold by maximizing modularity scalar value over a line search of threshold values
                        % NOTE: THIS IS NOT WORKING 100% OF THE TIME YET. LOOK INTO IT!  NOT REALLY USING IT RIGHT NOW.
                        N=10; % number of thresholds to look through on lineSearch
                        [TH, THs, modu] = threshModuLineSearch(EvecML,W,imDS,N);
                        
                        THstr = num2str(TH,'%0.2g'); THstr(THstr=='.')='p';
%                         seg{k} = EvecML>TH;                % best segmentation (UNCOMMENT TO USE!!)

                        %% Plot Segmentation at Threshold that gives Max Modularity.
                        seg{k} = EvecML>mean(EvecML(:));                % segmentation at mean

                        if(0)
                            segA = imDS.*seg;
                            segB = imDS.*~seg;           % pixel overlays using segments as mask
                            %
                            h=figure; colormap('bone') % 0
                            subplot(223), imagesc(segA),title(['Best Segment A @ ',num2str(TH,'%0.2g')],'FontSize',20,'FontWeight','Bold')
                            xlabel(['(\sigma_p = ',num2str(sigP),') (\sigma_d = ',num2str(sigD),') ( rmax = ',num2str(rmax),')'],'FontSize',16,'FontWeight','Bold')
                            subplot(224), imagesc( segB ),title(['Segment B modu = ',num2str(max(modu),'%0.2g')],'FontSize',20,'FontWeight','Bold')
                            subplot(221), imagesc( imDS ),title(['Image: ',fname(1:end-4)],'FontSize',20,'FontWeight','Bold')
                            xlabel([method],'FontSize',16,'FontWeight','Bold')
                            subplot(222), imagesc( EvecVizF(EvecML,vizSlope) ),colorbar,title(['Log E-vec',num2str(evectors(k))],'FontSize',20,'FontWeight','Bold') 
                            % replace real(log(EvecML)) with steep sigmoid symmetric about zero that doesnt assymptotes but continues onto inf like tanh func or 1/(1 + e^kx) for big k
                            set(h,'Position',[1 1 1280 700]);

                            [THs;full(modu)]

                            % if I start using this, I should change the directory structure
                            if ~exist(['./output/ImgSeg/simpleExamples/',dirOut,'/pics from main/seg/',method,params],'dir')
                                mkdir(['./output/ImgSeg/simpleExamples/',dirOut,'/pics from main/seg/',method,params]);
                            end

                            % save plotted eigenvector and segmentation as a jpg image file
                            saveas(h,['./output/ImgSeg/simpleExamples/',dirOut,'/pics from main/seg/',method,params,...
                                '/Seg ' fname(1:end-4),' ',patstr,' ',' ev',num2str(evectors(k)),' at',THstr],'jpg'); %name is a string                            
                            close
                        end
                        
                        %% Plot the number of chosen eigenvector in 1 figure (in subplots)
                        if(save_Evecs)

                            % Calculate how modular or piecewise the eigenvector is by ratio (Larger # means more Piecewise or better)
                            % THIS ISNT THE BEST MEASURE BECAUSE IT DEPENDS ON THE SIZE OF SEGMENTS.  NOT USING NOW, JUST DISPLAYING
                            [H, hx] = hist(EvecML(:));
                            PWrat = max(H)./min(H(H~=0));
                            
                            % to plot H histogram and display max & min bins
                            if(0)
                                figure, plot(hx,H)
                                title(['Piecewise Ratio: max bin = ',num2str(max(H)),' | min bin = ',num2str(min(H(H~=0)))])
                            end

                            figure(hev), subplot(xsub,ysub,k)
                            imagesc(EvecVizF(EvecML,vizSlope))
                            title(method,'FontSize',16,'FontWeight','Bold')
                            xlabel(['Evec',num2str(evectors(k)),' (\lambda = ',num2str(EigVal(end-evectors(k)+1)),')'],'FontSize',16,'FontWeight','Bold') % ,' (PieceWise # = ',num2str(PWrat,'%0.1f'),')'
                            set(gca,'xtick',[],'ytick',[]); % axis square
%                             if ( k == evectors(round(xsub)/2) )
%                                 title([fname(1:end-4),' ',method,' ',params],'FontSize',12,'FontWeight','Bold')
%                             end
                            colorbar('SouthOutside')
                        end


                    end % Loop through to calculate different eigenvectors

                    disp('Time to find Eigenvectors:')
                    toc
                    
%                     keyboard
                    
                    if(save_Evecs)
                        % save plotted eigenvector and segmentation as a jpg image file
%                         set(hev,'Position',[1 1 1280 200]);
                        colormap('bone')
                        saveas(hev,['./output/ImgSeg/simpleExamples/',dirOut,'/pics from main/evecs/'...
                            ,fname(1:end-4),' ',patstr,' ev',num2str(evectors(1)),'_',num2str(evectors(end)),' ',method,'-',params],'jpg');    
                        close(hev)
                    end

                end % If statement to calculate eigenvectors using MATLAB eig function and not using mean threshold method.
                
                if(meanThresh)
                    SegAtTh = imDS>mean(imDS(:));
                    hst = figure; imagesc(SegAtTh), colormap('bone'), colorbar('SouthOutside')
                    title(['Segmenting Input Image at Mean Pixel Value'],'FontSize',16,'FontWeight','Bold')
                    set(gca,'xtick',[],'ytick',[]); % axis square
                    saveas(hst,['./output/ImgSeg/simpleExamples/',dirOut,'/pics from main/evecs/'...
                            ,fname(1:end-4),' ',patstr,method],'jpg');    
                    close(hst)
                end
                
                %%  Use Power Method to find 1st or 2nd Eigenvectors of L and Q and W
                if(0)
                    disp('Power Method Eigenvector Calc')
                    tic
                    
%                     if(proj)
%                         EvecTrue = EvecBig2;
%                         ev = '2';
%                     else
%                         EvecTrue = EvecBig1;
%                         ev = '1';
%                     end

                    % Works fine now... Even with Normalized Graph Laplacian & Average Association
                    threshPM = 1e-7; % if not small enough, Evector may not converge.
                    maxPMiter = 100000; % if too low, Evector doesnt converge and does not look right.
                    [EvecPM, ConvergencePM] = power_method(Q, ximg, yimg, EvecML, 'random', method, proj, threshPM, maxPMiter, 0,0,1,0);
                    toc

                    % Final Plot Dominant Eigenvector as found from Power Method
                    figure, imagesc(EvecPM); colorbar('FontSize',16,'FontWeight','Bold'); 
                    title(['Eigenvector #',ev,' PM ',method],'FontSize',20,'FontWeight','Bold')
                    
                else
                    
                    EvecPM = 'Power Method was not run.';

                end

                
            
            
            %% save an output mat file that I will then pipe through BSDS benchmark code.
            if(save_Mat & ~meanThresh)
                matOutDir = ['./output/ImgSeg/simpleExamples/',dirOut,'/data from main/'];
                matOutFname = [fname(1:end-4)];
                %
                save([matOutDir,matOutFname,' ',patstr,' ',method,'-',params,],'seg','EvecMLs','EvecPM','Meth','inptParams','runFlags')       
            end
%             saveas(h,['./output/ImgSeg/simpleExamples/',dirOut,...
%                                 '/Seg ' fname(1:end-4),' ',patstr,' ',method,params],'jpg'); %name is a string
            

        end % loop over max radius for distance dependence (weighting between pixels separated by more than rmax = 0)    
        end % loop over threshold and sigma pairs on distance dependence
        end % loop over different sigma values on pixel similarity weighting
        
        % close figure for image
        if(save_im)
           close(him) 
        end
        if(save_im_sub)
           close(hims) 
        end

    end % loop over patches chosen for BSDS images
    
    % display how much memory matlab is currently using
    s = whos;
    disp(['Matlab is using ~',num2str( (sum([s.bytes]) + 500e6) / 1e9 ),' GB of memory'])
    % s.bytes tells memory for each variable and matlab uses about 500MB on startup.
    
end % loop over images for BSDS images





%% Code from Saurabh meeting to do benchmark stuff...
% load('/Users/world7one/Desktop/Grad_School/Berkeley/Work/Fritz_Work/Projects/images/BSDS_images/BSR/bench/data/ucm2/2018.mat')
% reg = bwlabel(ucm2 < 0.5); % <--- BIGGER NUMBER MEANS BIGGER SEGMENTS!
% reg = reg(2:2:end, 2:2:end);
% figure, imagesc(reg)


