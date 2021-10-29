% This script is meant to create a connection between the results of our
% various Image Segmentation alrgorithms and the benchmark code provided
% with the Berkeley Segmentation Dataset. We have extracted patches from
% images and here we want to go into the full ground truth and extract
% those same patches from the Segmentations and Boundaries data structures.

[dirPre,sizeGoodIm] = onCluster;

% location of mat files with extracted image patches
imPtchDir = [dirPre,'images/BSDS_patch/101x101_ds1/'];

% location of full image provided with BSDS with expected groundTruth structure
gtFullDir = [dirPre,'images/BSDS_images/BSR/BSDS500/data/groundTruth/'];

% location to write out the patch groundTruth file with expected structure.
gtOutDir = [dirPre,'images/BSDS_patch/101x101_ds1/groundTruth/'];
if ~exist(gtOutDir,'dir')
    mkdir(gtOutDir)
end



files = dir([imPtchDir,'*.mat']);

for i = 1:numel(files)
    
    fname = files(i).name;
    
    % if groundTruth patch mat file already exists in gtOutDir, dont remake it.
    if exist([gtOutDir,fname],'file')
        disp(['Ground Truth Patch File: ',fname,' :already exists. Moving on.'])
        continue
    end
    
    disp([num2str(i),' / ',num2str(numel(files)),' : ', fname])
    
    % load in image patch file
    load([imPtchDir,fname]);
    
    % extract imgID from imgPtch file name
    brk = strfind(fname,'_ptch');
    imgID = fname(1:brk-1);
    
    % load in original full size groundTruth
    load([gtFullDir,imgID,'.mat']);
    
    % index into gtFull and pull out only parts pertaining to patch
    for j = 1:numel(groundTruth)
        groundTruthP{j}.Segmentation = groundTruth{j}.Segmentation(pach.ypbeg:pach.ypfin,pach.xpbeg:pach.xpfin);
        groundTruthP{j}.Boundaries = groundTruth{j}.Boundaries(pach.ypbeg:pach.ypfin,pach.xpbeg:pach.xpfin);
    end

    % replace the full groundTruth with the patch's groundTruth structure.
    groundTruth = groundTruthP;
    
    % save gtPtch into the gtOutDir
    save([gtOutDir,fname],'groundTruth','pach')
    
    clear groundTruthP groundTruth

end