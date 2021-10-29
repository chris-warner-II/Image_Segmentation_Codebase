function [] = combine_boundaryGradientMetric_Kur(method,RM,KS)


% Script to combine all ~1500 mat files for individual image patches and
% combine their Dout results into a single mat file.  To make next step in
% processing quicker.  Do all the i/o now and only need to load 1 file
% later.


[dirPre, sizeGoodIm] = onCluster;

if(strmatch(method,'ImPix'))
    
    inDir = [dirPre,'images/BSDS_patch/101x101_ds1/'];
	%
    inDirB = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/IsoDiff/'];
    files = dir([inDirB,'KurMC_*rM3_NF_60_0_ksmid.mat']); % Why do I use this particular set of files?
    
    netParams = 'n/a';
    netflags = 'n/a';
    kurParams = 'n/a';
    kurflags = 'n/a';

    
elseif(strmatch(method,'IsoDiff'))
    
    inDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/',method,'/'];
    files = dir([inDir,'KurMC_*rM',RM,'_NF_60_0_ks',KS,'*.mat']);
    
else % method = {AA,GL,NG,SK}
    
    inDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/',method,'/'];
    files = dir([inDir,'KurMC_*rM',RM,'_sDInf_sP0p2_NF_60_0_ks',KS,'*.mat']);
    
end

% Directory to Save Output File Into
saveDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1/data/Kur_PIF_Fourier1/',method,'/'];
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end

saveDir
method

% If output file already exists, break out and do not run this code ...
if exist([saveDir,'BoundaryGradientMetricKur_rM',RM,'_KS',KS,'_files',num2str(numel(files)),'.mat'],'file')
   
    disp(['The combined Boundary Gradient Metric File already exists at:'])
    disp([saveDir,'BoundaryGradientMetricKur_rM',RM,'_KS',KS,'_files',num2str(numel(files))])
    disp(['Not going to regenerate it.  Delete it and run again if u want to regenerate it.'])
    return
    
end


% Preallocate Arrays to Fill Up as We Loop Thru Files
meanGradientOnBoundary = zeros(7,numel(files));
stdGradientOnBoundary = zeros(7,numel(files));
meanGradientOffBoundary = zeros(7,numel(files));
stdGradientOffBoundary = zeros(7,numel(files));
boundaryDiscriminability = zeros(7,numel(files));
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
    
    meanGradientOnBoundary(:,i) = MC.M(:,1);
    stdGradientOnBoundary(:,i) = MC.M(:,2);
    meanGradientOffBoundary(:,i) = MC.S(:,1);
    stdGradientOffBoundary(:,i) = MC.S(:,2);
    boundaryDiscriminability(:,i) = MC.D;
    
    

end

bDc_blurs_info = MC.bDc_blurs_info; % this assumes that it was the same for all individual files processed.





% Save combined file
save([saveDir,'BoundaryGradientMetricKur_rM',RM,'_KS',KS,'_files',num2str(numel(files))],'ImgPtchID','meanGradientOnBoundary',...
        'stdGradientOnBoundary','meanGradientOffBoundary','stdGradientOffBoundary','boundaryDiscriminability','bDc_blurs_info',...
        'netParams','netflags','kurParams','kurflags')
    
    
    
if(0)
    figure, errorbar(mean(meanGradientOnBoundary,2),mean(stdGradientOnBoundary,2),'b')
    hold on, errorbar(mean(meanGradientOffBoundary,2),mean(stdGradientOffBoundary,2),'r')
    legend({'on boundary','off boundary'})
end