function [] = combine_boundaryGradientMetric_Eig(method,RM)


% Script to combine all ~1500 mat files for individual image patches and
% combine their Dout results into a single mat file.  To make next step in
% processing quicker.  Do all the i/o now and only need to load 1 file
% later.


[dirPre, sizeGoodIm] = onCluster;

if(strmatch(method,'IsoDiff'))
    
    inDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/spectral/',method,'/'];
    files = dir([inDir,'Evecs_*rM',RM,'.mat']);

    
else % method = {AA,GL,NG,SK}
    
    inDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/spectral/',method,'/'];
    files = dir([inDir,'Evecs_*rM',RM,'_sDInf_sP0p2.mat']);
    
end

% Directory to Save Output File Into
saveDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/spectral/',method,'/'];
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end



% If output file already exists, break out and do not run this code ...
if exist([saveDir,'BoundaryGradientMetricEig_rM',RM,'_files',num2str(numel(files))],'file')
   
    disp(['The combined Boundary Gradient Metric File already exists at:'])
    disp([saveDir,'BoundaryGradientMetricEig_rM',RM,'_files',num2str(numel(files))])
    disp(['Not going to regenerate it.  Delete it and run again if u want to regenerate it.'])
    return
    
end


files

% Preallocate Arrays to Fill Up as We Loop Thru Files
meanGradientOnBoundary.im = zeros(7,numel(files));
stdGradientOnBoundary.im = zeros(7,numel(files));
meanGradientOffBoundary.im = zeros(7,numel(files));
stdGradientOffBoundary.im = zeros(7,numel(files));
boundaryDiscriminability.im = zeros(7,numel(files));
%
meanGradientOnBoundary.ev1 = zeros(7,numel(files));
stdGradientOnBoundary.ev1 = zeros(7,numel(files));
meanGradientOffBoundary.ev1 = zeros(7,numel(files));
stdGradientOffBoundary.ev1 = zeros(7,numel(files));
boundaryDiscriminability.ev1 = zeros(7,numel(files));
%
meanGradientOnBoundary.ev2o = zeros(7,numel(files));
stdGradientOnBoundary.ev2o = zeros(7,numel(files));
meanGradientOffBoundary.ev2o = zeros(7,numel(files));
stdGradientOffBoundary.ev2o = zeros(7,numel(files));
boundaryDiscriminability.ev2o = zeros(7,numel(files));
%
meanGradientOnBoundary.ev3o = zeros(7,numel(files));
stdGradientOnBoundary.ev3o = zeros(7,numel(files));
meanGradientOffBoundary.ev3o = zeros(7,numel(files));
stdGradientOffBoundary.ev3o = zeros(7,numel(files));
boundaryDiscriminability.ev3o = zeros(7,numel(files));
%
meanGradientOnBoundary.ev2 = zeros(7,numel(files));
stdGradientOnBoundary.ev2 = zeros(7,numel(files));
meanGradientOffBoundary.ev2 = zeros(7,numel(files));
stdGradientOffBoundary.ev2 = zeros(7,numel(files));
boundaryDiscriminability.ev2 = zeros(7,numel(files));
%
meanGradientOnBoundary.ev3 = zeros(7,numel(files));
stdGradientOnBoundary.ev3 = zeros(7,numel(files));
meanGradientOffBoundary.ev3 = zeros(7,numel(files));
stdGradientOffBoundary.ev3 = zeros(7,numel(files));
boundaryDiscriminability.ev3 = zeros(7,numel(files));
%
meanGradientOnBoundary.ev2w = zeros(7,numel(files));
stdGradientOnBoundary.ev2w = zeros(7,numel(files));
meanGradientOffBoundary.ev2w = zeros(7,numel(files));
stdGradientOffBoundary.ev2w = zeros(7,numel(files));
boundaryDiscriminability.ev2w = zeros(7,numel(files));
%
meanGradientOnBoundary.ev3w = zeros(7,numel(files));
stdGradientOnBoundary.ev3w = zeros(7,numel(files));
meanGradientOffBoundary.ev3w = zeros(7,numel(files));
stdGradientOffBoundary.ev3w = zeros(7,numel(files));
boundaryDiscriminability.ev3w = zeros(7,numel(files));


ImgPtchID = cell(1,numel(files));



% Loop through files
for i = 1:numel(files)
    
    disp(['File #',num2str(i),' / ',num2str(numel(files)),' : ',method,' : ',files(i).name ])
    
    
    ind = strfind(files(i).name,'_rM') - 1;
    ImgPtchID{i} = files(i).name(7:ind);
    
    try
        if(strmatch(method,'ImPix'))
            load([inDir,ImgPtchID{i},'.mat'])
        else
            load([inDir,files(i).name])
        end
    catch
        disp(['Something Wrong:  Can Not Load this File!  Maybe Corrupt! Look Into this.  Maybe Delete it.'])
    end
    
    meanGradientOnBoundary.im(:,i) = MC.im.M(:,1);
    stdGradientOnBoundary.im(:,i) = MC.im.M(:,2);
    meanGradientOffBoundary.im(:,i) = MC.im.S(:,1);
    stdGradientOffBoundary.im(:,i) = MC.im.S(:,2);
    boundaryDiscriminability.im(:,i) = MC.im.D;
    %
    meanGradientOnBoundary.ev1(:,i) = MC.ev1.M(:,1);
    stdGradientOnBoundary.ev1(:,i) = MC.ev1.M(:,2);
    meanGradientOffBoundary.ev1(:,i) = MC.ev1.S(:,1);
    stdGradientOffBoundary.ev1(:,i) = MC.ev1.S(:,2);
    boundaryDiscriminability.ev1(:,i) = MC.ev1.D;
    %
    meanGradientOnBoundary.ev2o(:,i) = MC.ev2o.M(:,1);
    stdGradientOnBoundary.ev2o(:,i) = MC.ev2o.M(:,2);
    meanGradientOffBoundary.ev2o(:,i) = MC.ev2o.S(:,1);
    stdGradientOffBoundary.ev2o(:,i) = MC.ev2o.S(:,2);
    boundaryDiscriminability.ev2o(:,i) = MC.ev2o.D;
    %
    meanGradientOnBoundary.ev3o(:,i) = MC.ev3o.M(:,1);
    stdGradientOnBoundary.ev3o(:,i) = MC.ev3o.M(:,2);
    meanGradientOffBoundary.ev3o(:,i) = MC.ev3o.S(:,1);
    stdGradientOffBoundary.ev3o(:,i) = MC.ev3o.S(:,2);
    boundaryDiscriminability.ev3o(:,i) = MC.ev3o.D;
    %
    meanGradientOnBoundary.ev2(:,i) = MC.ev2.M(:,1);
    stdGradientOnBoundary.ev2(:,i) = MC.ev2.M(:,2);
    meanGradientOffBoundary.ev2(:,i) = MC.ev2.S(:,1);
    stdGradientOffBoundary.ev2(:,i) = MC.ev2.S(:,2);
    boundaryDiscriminability.ev2(:,i) = MC.ev2.D;
    %
    meanGradientOnBoundary.ev3(:,i) = MC.ev3.M(:,1);
    stdGradientOnBoundary.ev3(:,i) = MC.ev3.M(:,2);
    meanGradientOffBoundary.ev3(:,i) = MC.ev3.S(:,1);
    stdGradientOffBoundary.ev3(:,i) = MC.ev3.S(:,2);
    boundaryDiscriminability.ev3(:,i) = MC.ev3.D;
    %
    meanGradientOnBoundary.ev2w(:,i) = MC.ev2w.M(:,1);
    stdGradientOnBoundary.ev2w(:,i) = MC.ev2w.M(:,2);
    meanGradientOffBoundary.ev2w(:,i) = MC.ev2w.S(:,1);
    stdGradientOffBoundary.ev2w(:,i) = MC.ev2w.S(:,2);
    boundaryDiscriminability.ev2w(:,i) = MC.ev2w.D;
    %
    meanGradientOnBoundary.ev3w(:,i) = MC.ev3w.M(:,1);
    stdGradientOnBoundary.ev3w(:,i) = MC.ev3w.M(:,2);
    meanGradientOffBoundary.ev3w(:,i) = MC.ev3w.S(:,1);
    stdGradientOffBoundary.ev3w(:,i) = MC.ev3w.S(:,2);
    boundaryDiscriminability.ev3w(:,i) = MC.ev3w.D;
    %

end

bDc_blurs_info = MC.bDc_blurs_info; % this assumes that it was the same for all individual files processed.



if(0)
    figure, errorbar(meanGradientOnBoundary,stdGradientOnBoundary)
    figure, errorbar(meanGradientOffBoundary,stdGradientOffBoundary)
end


% Save combined file
save([saveDir,'BoundaryGradientMetricEig_rM',RM,'_files',num2str(numel(files))],'ImgPtchID','meanGradientOnBoundary',...
        'stdGradientOnBoundary','meanGradientOffBoundary','stdGradientOffBoundary','boundaryDiscriminability','bDc_blurs_info',...
        'netParams','netflags')
    
    
    
if(0)
    figure, errorbar(mean(meanGradientOnBoundary,2),mean(stdGradientOnBoundary,2),'b')
    hold on, errorbar(mean(meanGradientOffBoundary,2),mean(stdGradientOffBoundary,2),'r')
    legend({'on boundary','off boundary'})
end