function [] = exp_SepVPar_avgOvrImgPatches(fileType,fileSize,Tscale,sigP,Rmax,sigD,sigW,Kscale)

% This function will take in Kur_metaSummary mat files and Evec Rm sD sP 
% mat files that are both made in explore_Separation_vs_Parameters.  Each
% Kur mat file is for a single image patch and contains results for all
% parameter values inside.  Each Evec file contains results for single
% image patch and parameter combination.  This code will loop through these
% files and produce 2 x 2D scatter plots (Kur vs. Straw) & (Kur vs. Eig)
% Divisive Relative Margin plots.  Each pair of plots will be a single
% combination of parameters and each will contain one point for every image
% patch & ground truth pair (~25000 points after all data is processed).



% Directories to Input Images or Data :: Change these paths below to get a different images or patches
[dirPre,sizeGoodIm] = onCluster;
dirKur = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/Kur_PIF_Fourier1/Mod_SKHAdj/'];
dirEig = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/spectral/Mod_SKHAdj/'];

filesKur = dir([dirKur,'Kur_metaSummary*tscale',Tscale,'*']); % can loop through these later


% Directory to save out 2x2D scatter images into 
imgKur = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/imgs/Kur_PIF_Fourier1/Mod_SKHAdj/AvgAcrossImgs/'];
%
if ~exist(imgKur,'dir')
    mkdir(imgKur)
end

% turn numbers into strings for image filename.
sPs = num2str(sigP);
sPs(sPs=='.')='p';
rMs = num2str(Rmax);
rMs(rMs=='.')='p';
sDs = num2str(sigD);
sDs(sDs=='.')='p';
sWs = num2str(sigW);
sWs(sWs=='.')='p';
Kss = num2str(Kscale);
Kss(Kss=='.')='p';
% Note: Recall: sigD is either Inf or 1/4*Rmax




% Check if DivMarg_Scatter_2D_... file exists already.  If it does, bail out.
if exist([imgKur,'DivMarg_Scatter_2D_tscale',Tscale,'_sP',sPs,'_rM',rMs,'_sD',sDs,'_sW',sWs,'_Kscale',Kss,'.mat'],'file')

    disp(['Output mat file already exists: '])
    disp([imgKur,'DivMarg_Scatter_2D_tscale',Tscale,'_sP',sPs,'_rM',rMs,'_sD',sDs,'_sW',sWs,'_Kscale',Kss,'.mat'])
    disp(['Not regenerating it'])
    return
    
end


kkk = 0;         % counter when looping through filesKur & ground truths.
mis = 0;         % counter for misaligned Kur/Straw combinations.
misalignedID={}; % index into all vectors to be saved to find misaligned pairs
%
KurDivMarg = [];
EigDivMarg = [];
SMDivMarg  = [];

tic

for i = 1:numel(filesKur)

    KurFile = load([dirKur,filesKur(i).name]); % loading Kur_metaSummary files.

    disp(['Looking in file: ',filesKur(i).name])
    disp(['for tscale',Tscale,'_sP',num2str(sigP),'_rM',num2str(Rmax),'_sD',num2str(sigD),'_sW',num2str(sigW),'_Kscale',num2str(Kscale)])
    
    disp(['Image Patch # ',num2str(i),' / ',num2str(numel(filesKur))])
    indKur = find( KurFile.sigP == sigP & KurFile.sigD == sigD & KurFile.Rmax == Rmax & ...
                   KurFile.sigW == sigW & KurFile.Kscale == Kscale );
               
    if ~isempty(indKur)
        
        disp('found')

        if ~isfield(KurFile,'DivM_EV')
            KurFile.DivM_EV(indKur,:) = ones(size(KurFile.DivMarg_SM));
            disp('No Eigenvector Data')
            keyboard
        end
        
    
        % build up a data structure patch and ground truth info for each data point.
        KurDivMarg = [ KurDivMarg , KurFile.DivM_End(indKur,:) ];
        EigDivMarg = [ EigDivMarg , KurFile.DivM_EV(indKur,:) ];
        SMDivMarg  = [ SMDivMarg  , KurFile.DivMarg_SM ];

        gT = numel(KurFile.DivMarg_SM);

        beg = strfind(filesKur(i).name,fileType) + numel(fileType) + 1;
        fin = strfind(filesKur(i).name,'__') ;

        for j = 1:gT
            kkk = kkk + 1;
            ImgPtchID{kkk} = [filesKur(i).name(beg:fin)];
            ImgGtID(kkk) = j;
        end
    
    end

    % visualize the strawman & kuramoto to see if 2D plots make sense.
    % Much of this is copied from exp_SepVPar_avgOvrImgPatches.m file.
    [dirPre,sizeGoodIm] = onCluster;

    fileType = 'BSDS_patch';
    fileSize = '51x51_ds1';


    % directory to image patches extracted from BSDS.
    patchesDir = [dirPre,'images/',fileType,'/',fileSize,'/'];

    % directory to KurMC files and to Kur_metaSummary files (think I need the KurMC ones)
    dirKurMats = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/Kur_PIF_Fourier1/Mod_SKHAdj/'];

    % directory to the Evecs files
    dirEigMats = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/spectral/Mod_SKHAdj/'];



    % (1). Load in Image Patch to plot it and ground truth.
    I=numel(ImgPtchID);
    ImPtch = load([patchesDir,ImgPtchID{I}(1:end-1),'.mat']);




    % (2). Load in KurMC file and plot a.) ending solution and 
    %      b.) Kur DivMarg time evolution along with Strawman DivMarg.
    if(0)
        KurMetFileName = dir([dirKurMats,'Kur_metaSummary_',fileType,'_',ImgPtchID{I},'_tscale',Tscale,'_*files.mat']);
        KurMeta = KurMetFileName(end).name; % grab metaSummary file with most files contained in it (72 is max)
    end
    %
    % Note: Not really using the KurMeta file right now
    %
    KurMC = load([dirKurMats,'KurMC_',ImgPtchID{I},'rM',rMs,'_sD',sDs,'_sP',sPs,'_NF_60_',sWs,'_kscale',Kss,'_tscale',Tscale,'_runs1.mat']);



    % Check that Kur_MC file is looking at correct image patch...
    xxx = abs(ImPtch.im - KurMC.netParams.im);
    if mean(xxx(:)) % not = 0.  Meaning there is any difference,
        mis = mis+1;
        misalignedID{mis} = ImgPtchID(I);
    end

    % (3). Load in Evecs file to plot Eigenvector solution and Eig DivMarg.
    % EigDivMarg(I) % Note: Not dealing with Eig right now...
    
        

    
%     % probe to check memory used by Matlab (to determine how much to
%     % request in sbatch script when I submit job on Cluster)
%     check_memory_usage


end



% I could also maybe save vectors with 3 different RM results and which image, ground truth and patch the input belonged to in a mat file.
save([imgKur,'DivMarg_Scatter_2D_tscale',Tscale,'_sP',sPs,'_rM',rMs,'_sD',sDs,'_sW',sWs,'_Kscale',Kss],'KurDivMarg','EigDivMarg','SMDivMarg',...
    'ImgPtchID','ImgGtID','fileType','fileSize','Tscale','sigP','Rmax','sigD','sigW','Kscale','misalignedID')


toc
