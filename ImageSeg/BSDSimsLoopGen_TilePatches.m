function [] = BSDSimsLoopGen_TilePatches(boxDims,ims2Proc,origin,ds_fctr)

% Create an ENSEMBLE OF IMAGES with an from the Berkeley Segmentation
% Dataset (BSDS) that consists of a number of patches and the groundtruth
% for each of the 500 images.  That data, stored in MAT files will be fed
% into the Segmentation Methods.

% 'origin' is where to start corner of first box (set to 0,0 for now).  May
% get different results depending in where edges of patches fall -
% especially if they do not communicate with eachother and allign phase
% offsets or whatever.

patches=1; % a flag to extract patches from images instead of downsampling

[dirPre,sizeGoodIm] = onCluster;

im_st = ims2Proc(1); 
im_fin = ims2Proc(2); %500; % can be 500

if ~exist('ds_fctr')
    ds_fctr = 1;
end

if ~isstr(boxDims)
    
    xpat = boxDims(1)-1; %40;
    ypat = boxDims(2)-1; %30;
    
    saveDir = [dirPre,'images/BSDS_tile/',num2str(xpat+1),'x',num2str(ypat+1),'_ds',num2str(ds_fctr),'/'];
    
else
    saveDir = [dirPre,'images/BSDS_full/ds',num2str(ds_fctr),'/'];    
end

% Make directories to house stacks of input images generated
    if ~exist(saveDir,'dir')
        mkdir(saveDir)
    end


disp('Input Images are from BSDS500');
iids = imgList('all') % can be {'all','test','train','val'}  
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



%% Loop through input images
for iii = 1:numel(iids)
    
    iid = iids(iii).name;
    disp(['Image #',num2str(iii),' of ',num2str(numel(iids)),' : ',iid])
    
    
    % read in image to be segmented if analyzing BSDS image.
%     try
        fname=iid;
        [im] = imgRead(fname,'gray'); % For running code on my computer
        load(GtFilename(fname(1:end-4)));
%     catch
%         disp('Cluster doing something funny.')
%         fname
%         fname=iid(3:end)
%         [im] = imgRead(fname,'gray'); % Cluster does something funny -> ' ._name'
%         load(GtFilename(fname(1:end-4)));
%     end

    % Save info/data on fullsize image
    imFull = im;                      % to save original size image
    gTfull = groundTruth;             % to save original size groundtruth data
    xFimg = size(imFull,2);           % horizontal image dimension
    yFimg = size(imFull,1);           % vertical image dimension
    
    
    if isstr(boxDims)
       xpat = xFimg;
       ypat = yFimg;
    end
    

    
    % Set up a grid to index into full image and extract patches
    xpbeg_grid = 1:(xpat+1):size(imFull,2);
    xpfin_grid = [(xpat+1):(xpat+1):size(imFull,2), size(imFull,2)];
    ypbeg_grid = 1:(ypat+1):size(imFull,1);
    ypfin_grid = [(ypat+1):(ypat+1):size(imFull,1), size(imFull,1)];
    xgrid = numel(xpbeg_grid);
    ygrid = numel(ypbeg_grid);
    
    % preallocate vectors to hold patch processing parameters.
    numpatches = xgrid*ygrid;
    pach.xpbeg = zeros(1,numpatches);
    pach.xpfin = zeros(1,numpatches);
    pach.ypbeg = zeros(1,numpatches);
    pach.ypfin = zeros(1,numpatches);
    pach.pnum = zeros(1,numpatches);
    
    
    % Loop through and extract patches at locations in image that tile it.
    k=1;
    for i = 1:xgrid
        for j = 1:ygrid

            pach.xpbeg(k) = xpbeg_grid(i);
            pach.xpfin(k) = xpfin_grid(i);
            pach.ypbeg(k) = ypbeg_grid(j);
            pach.ypfin(k) = ypfin_grid(j);
            pach.pnum(k) = k;
            %pach.xnum(k)
            k=k+1;

        end
    end
       
    
    
    
    
    
    
    % Chop up fill image into patches of size (xpat x ypat) as many as are
    % necessary and maybe dont worry about edge effects for now.
    for j = 1:numpatches

        % extract a random patch from full resolution image (instead of downsampling it)
        if(patches)

            [patch,gT] = tilePatches(imFull,groundTruth,pach,j); % extract patches that tile image
            
            im = patch;

            pach.xpat = xpat+1;
            pach.ypat = ypat+1;
            
        else
            pach='Patches not used';
        end

        im = imresize(im,1/ds_fctr); % downsample image (if ds_fctr = 1, no downsample)
        
        ximg = size(im,2);           % horizontal image dimension
        yimg = size(im,1);           % vertical image dimension
        

        % save each image patch separately stack to a mat file to be used as input later.
        fOut = [fname(1:end-4)];
        save([saveDir,fOut,'_ptch',num2str(j)],'im','gT','imFull','gTfull','pach','ds_fctr','ximg','yimg','xFimg','yFimg');
    
        
        

    end % loop over patches chosen for BSDS images
    
    
    
end % loop over images for BSDS images


