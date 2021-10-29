function [] = BSDSimsLoopGen_RandomPatch(boxDims,numpatches)

% Create an ENSEMBLE OF IMAGES with an from the Berkeley Segmentation
% Dataset (BSDS) that consists of a number of patches and the groundtruth
% for each of the 500 images.  That data, stored in MAT files will be fed
% into the Segmentation Methods.

[dirPre,sizeGoodIm] = onCluster;

rng(1234567); % Random Number Generator Seed.  For repeatability of "randomly chosen" patches.
% rng(316);

im_st = 1; 
im_fin = 500; % can be up to 500
ds_fctr = 1;

patches=1; % a flag to extract patches from images instead of downsampling
% numpatches = 10;
xpat = boxDims(1)-1; %40;
ypat = boxDims(2)-1; %30;


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

% Make directories to house stacks of input images generated
if ~exist([dirPre,'images/BSDS_patch/',num2str(xpat+1),'x',num2str(ypat+1),'_ds',num2str(ds_fctr),'/'],'dir')
    mkdir([dirPre,'images/BSDS_patch/',num2str(xpat+1),'x',num2str(ypat+1),'_ds',num2str(ds_fctr),'/'])
end



%% Loop through input images
for iii = 1:numel(iids)
    
    iid = iids(iii).name;
    disp(['Image #',num2str(iii),' of ',num2str(numel(iids)),' : ',iid])
    
    % preallocate vectors to hold patch processing parameters.
    pach.xpbeg = zeros(1,numpatches);
    pach.xpfin = zeros(1,numpatches);
    pach.ypbeg = zeros(1,numpatches);
    pach.ypfin = zeros(1,numpatches);
    pach.pnum = zeros(1,numpatches);
    
    
    % read in image to be segmented if analyzing BSDS image.
%     try
        fname=iid;
        [im] = imgRead(fname,'color'); % For running code on my computer
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
    
    for j = 1:numpatches

        % extract a random patch from full resolution image (instead of downsampling it)
        if(patches)

            [patch,gT,xpbeg,xpfin,ypbeg,ypfin,pnum] = ...
                extractPatch(imFull,xpat,ypat,groundTruth,fname(1:end-4),dirOut,j); % extract random patch in image
            
            if(numel(patch)==1)
                continue % move on to next patch.
            end
            
            im = patch;

%             pach.patch(:,:,j) = patch; % SAVED IN IMENS.

            pach.xpbeg = xpbeg;
            pach.xpfin = xpfin;
            pach.ypbeg = ypbeg;
            pach.ypfin = ypfin;
            pach.pnum = pnum;
            pach.xpat = xpat+1;
            pach.ypat = ypat+1;
            
        else
%             patstr=''; % may not need this?
            pach='Patches not used';
        end

        im = imresize(im,1/ds_fctr); % downsample image (if ds_fctr = 1, no downsample)
        
        ximg = size(im,2);           % horizontal image dimension
        yimg = size(im,1);           % vertical image dimension
        
%         disp(['Image size: ',num2str(ximg),'x',num2str(yimg)])
        
        %imEns(:,:,j) = im;
        
        
        % save each image patch separately stack to a mat file to be used as input later.
        fOut = [fname(1:end-4)];
        save([dirPre,'images/BSDS_patch/',num2str(xpat+1),'x',num2str(ypat+1),'_ds',num2str(ds_fctr),'/',fOut,'_ptch',num2str(j)],'im','gT','imFull','gTfull','pach','ds_fctr','ximg','yimg','xFimg','yFimg');
    
        
        

    end % loop over patches chosen for BSDS images
    
    
    
end % loop over images for BSDS images


