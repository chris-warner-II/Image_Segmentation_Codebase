function [] = explore_compare_DivMarg(fileGeneral,fileSize,rM,sD,sP,sW,Ts,Ks)


% This function will take output mat file from main_Kuramoto function that is titled
% metaClusterKur_...  This will plot the phase of each oscillator at a 60Hz
% clock and will colorcode them by ground truth segment.
%
%
% Inputs to function must be of this form
%
% fileGeneral = 'GradientBox'; % 'BSDS_patch'; % 
% fileSize = '51x51_ds1';



plot2x2Dscatter = 1; % flag to make plot of 2x2D scatter Kur vs Straw & Eig vs Straw for each individual patch. 
                     % If this flag is zero, then this function is just creating the Kur_metaSummary files if they do not exist.

                     
                     
% Directories to Input Images or Data :: Change these paths below to get a different images or patches
[dirPre,sizeGoodIm] = onCluster;
dirKur = [dirPre,'output/Kuramoto/NetsFromImgs/',fileGeneral,'_',fileSize,'/data/Kur_PIF_Fourier1/Mod_SKHAdj/']; %Theta Scale ',tscale,'/'];
dirEig = [dirPre,'output/Kuramoto/NetsFromImgs/',fileGeneral,'_',fileSize,'/data/spectral/Mod_SKHAdj/'];
dirImg = [dirPre,'images/',fileGeneral,'/',fileSize,'/'];    % directory to image patches extracted from BSDS.



% Specify these with strings or with ''.
rM = num2str(rM); % '1';
sD = num2str(sD); % 'Inf';
sP = num2str(sP); % '0p2';
sP(sP=='.')='p';
Ts = num2str(Ts); % '1';
Ks = num2str(Ks); % '300';
sW = num2str(sW); % '0';

scatKur = 1; % This may only work when both scatKur & scatEig are 1.
scatEig = 1;

[dirKur,'KurMC','*rM',rM,'*sD',sD,'*sP',sP,'*NF_60_',sW,'*kscale',Ks,'*tscale',Ts,'*.mat']
[dirEig,'Evecs','*rM',rM,'_sD',sD,'_sP',sP,'*.mat']

filesKur = dir([dirKur,'KurMC','*rM',rM,'*sD',sD,'*sP',sP,'*NF_60_',sW,'*kscale',Ks,'*tscale',Ts,'*.mat']); % can loop through these later
filesEig = dir([dirEig,'Evecs','*rM',rM,'_sD',sD,'_sP',sP,'*.mat']); % will be in a different directory later.

% Directory to save images into (histograms, DivMarg Distributions for top & bottom parameter settings, etc...)
imgKur = [dirPre,'output/Kuramoto/NetsFromImgs/',fileGeneral,'_',fileSize,'/imgs/Kur_PIF_Fourier1/Mod_SKHAdj/AvgAcrossImgs/'];
if ~exist(imgKur,'dir')
    mkdir(imgKur)
end


% Directory to save 2D Scatter plots into
dirScat = [dirPre,'output/Kuramoto/NetsFromImgs/',fileGeneral,'_',fileSize,'/imgs/spectral/Mod_SKHAdj/'];
if ~exist(dirScat,'dir')
    mkdir(dirScat)
end

% Directory to save mat file containing location of each point along with identity of image patch where it came from.
dirDMmat = [dirPre,'output/Kuramoto/NetsFromImgs/',fileGeneral,'_',fileSize,'/data/Kur_PIF_Fourier1/Mod_SKHAdj/'];
if ~exist(dirDMmat,'dir')
    mkdir(dirDMmat)
end



% If Patch from BSDS images then this shows full image, patch and segmentation at pixel mean
if(0)
    
    fKur = filesKur(1).name;
    st = 7;
    nd = strfind(fKur,'_rM')-1;
    fImg = fKur(st:nd); 
    load([dirImg,fImg]);
    %
    try
        explore_BSDS_patches % figure handle in script is hIm
    catch
        hIm=figure; imagesc(im), colormap('bone'), title('Gradient Box Input Image')
    end
    
    % save this image
    saveGoodImg(hIm,[imgKur,'Input_Image'],sizeGoodIm)
    close(hIm)

end


    

    
disp('Time to loop through all metaCluster mat files:')
tic




% Preset data arrays to hold results across all BSDS image patches.
kkk = 0;         % counter when looping through filesKur & ground truths.
mis = 0;         % counter for misaligned Kur/Straw combinations.
misalignedID={}; % index into all vectors to be saved to find misaligned pairs
%
KurDivMarg = [];    % Divisive Margin from Kuramoto Coupled Oscillator final Phase Distribution
%
SMDivMarg_L  = [];    % Divisive Margin from StrawMan (Image Pixel Values)
SMDivMarg_C  = [];    % Divisive Margin from StrawMan (Image Pixel Values)
%
Eig1DivMarg_L = [];   % Divisive Margin from 1st Eigenvector (straight up)
Eig1vDivMarg_L = [];  % Divisive Margin from 1st Eigenvector with NonLinear Visualization applied.
%
Eig1DivMarg_C = [];   % Divisive Margin from 1st Eigenvector (straight up)
Eig1vDivMarg_C = [];  % Divisive Margin from 1st Eigenvector with NonLinear Visualization applied.
%
Eig2oDivMarg = [];  % Divisive Margin from 2nd Eigenvector only.
Eig2ovDivMarg = []; % Divisive Margin from 2nd Eigenvector only with NonLinear Visualization applied.
%
Eig3oDivMarg = [];  % Divisive Margin from 3rd Eigenvector only.
Eig3ovDivMarg = []; % Divisive Margin from 3rd Eigenvector only with NonLinear Visualization applied.
%  
Eig2DivMarg = [];   % Divisive Margin from Eigenvectors 1-2.
Eig2vDivMarg = [];  % Divisive Margin from Eigenvectors 1-2 with NonLinear Visualization applied.
Eig2wDivMarg = [];  % Divisive Margin from Eigenvectors 1-2 weighted by their Eigenvectors. 
Eig2wvDivMarg = []; % Divisive Margin from Eigenvectors 1-2 weighted by their Eigenvectors with NonLinear Visualization applied. 
%
Eig3DivMarg = [];   % Divisive Margin from Eigenvectors 1-3.
Eig3vDivMarg = [];  % Divisive Margin from Eigenvectors 1-3 with NonLinear Visualization applied.
Eig3wDivMarg = [];  % Divisive Margin from Eigenvectors 1-3 weighted by their Eigenvectors.
Eig3wvDivMarg = []; % Divisive Margin from Eigenvectors 1-3 weighted by their Eigenvectors with NonLinear Visualization applied. 

%




% Decide which files to loop through
if (scatEig & ~scatKur)
    
    numFiles2Loop = numel(filesEig);

else
    
    numFiles2Loop = numel(filesKur);

end


% numFiles2Loop = 500; % get rid of this.


% Preallocate arrays to hold avg difference between Linear and Circular DivMarg computation for Straw, Evec1 & Evec1viz.
Straw_LvC_diffMn = zeros(1,numFiles2Loop);
Evec1_LvC_diffMn = zeros(1,numFiles2Loop);
Evec1v_LvC_diffMn = zeros(1,numFiles2Loop);


for i = 1:numFiles2Loop % I dont know if this will work if I am not looping thru filesKur.

    disp([num2str(i),' / ',num2str(numFiles2Loop)])
    
    
    if (scatEig & scatKur)
        
        % Load file containing metaClustering Analysis results for Kuramoto & StrawMan for parameter combination
        try
            %[dirKur,filesKur(i).name]
            Kur = load([dirKur,filesKur(i).name]);
        catch
            disp('Kuramoto file is corrupted or not there: Skipping to next...')
            [dirKur,filesKur(i).name]
            %delete([dirKur,filesKur(i).name])
            continue
        end

        % Load in Eigenvector Segmentation results too so we can compare all 3 (Kuramoto, Eigenvector, Strawman)
        try
            %[dirEig,'Evecs_',Kur.kurflags.fname,'_rM',Kur.netflags.rM,'_sD',Kur.netflags.sD,'_sP',Kur.netflags.sP,'.mat']
            Eig = load([dirEig,'Evecs_',Kur.kurflags.fname,'_rM',Kur.netflags.rM,'_sD',Kur.netflags.sD,'_sP',Kur.netflags.sP,'.mat']);
        catch
            disp('Eigen file is corrupted or not there: Skipping to next...')
            [dirEig,'Evecs_',Kur.kurflags.fname,'_rM',Kur.netflags.rM,'_sD',Kur.netflags.sD,'_sP',Kur.netflags.sP,'.mat']
            %delete([dirEig,'Evecs_',kurflags.fname,'_rM',netflags.rM,' sD',netflags.sD,' sP',netflags.sP,'.mat'])
            continue
        end
    
    end

    


    % Loop through each human ground truth segmentation
    for g = 1:numel(Kur.netParams.gT)

        kkk = kkk + 1;
        ImgPtchID{kkk} = Kur.kurflags.fname; %[filesKur(i).name(beg:fin)];
        ImgGtID(kkk) = g;
        
    end
    %
    %
    %
    
    % Check that Kur_MC & Evecs files are looking at correct image patch...
    ImPtch = load([dirImg,ImgPtchID{end},'.mat']);
    %
    xxx = abs(ImPtch.im - Kur.netParams.im);
    yyy = abs(ImPtch.im - Eig.netParams.im);
    if ( mean(xxx(:)) | mean(yyy(:)) ) % not = 0.  Meaning there is any difference,
        mis = mis+1;
        misalignedID{mis} = ImgPtchID(end);
        
        figure, colormap('bone')
        subplot(131), imagesc(ImPtch.im), title('ImPtch')
        subplot(132), imagesc(Eig.netParams.im), title('Eig')
        subplot(133), imagesc(Kur.netParams.im), title('Kur')
        
        keyboard
        
        disp('Patches not aligned.  Throw away!. Break!')
        continue
        % Note: If I continue and dont compute all the stuff below, I dont
        % think I want to be saving this misaligned structure.
        % I DONT WANT ANY IMAGE PATCHES MISALIGNED!  - SEEMS FIXED NOW.
        close
    end
    
    
    %
    %
    %
    if iscell(Kur.MC)
        DivMarg_K = zeros( size(Kur.MC{1}.DistAvgPW,1), numel(Kur.netParams.gT) );
        DivMarg_SMa = zeros( 1, numel(Kur.netParams.gT));   
        for g = 1:numel(Kur.netParams.gT)
            DivMarg_K(:,g) =  Kur.MC{1}.DistAvgPW(:,1,g) ./ Kur.MC{1}.DistAvgPW(:,2,g);             % Divisive Margin for Coupled Osc. Model
            DivMarg_SMa(g) =  Kur.MC{1}.DistAvgPW(1,1,g) ./ Kur.MC{1}.DistAvgPW(1,2,g);             % Divisive Margin for Straw Man Model from Kur
        end
        KurDivMarg = [ KurDivMarg , DivMarg_K ];
        SMDivMarg_C  = [ SMDivMarg_C  , DivMarg_SMa];
    end
    %
    %
    %
    if isfield(Eig.MC,'im_l')
        DivMarg_SMb = zeros( 1, numel(Eig.netParams.gT)); 
        for g = 1:numel(Eig.netParams.gT)
            DivMarg_SMb(g) = Eig.MC.im_l{g}.DistAvgPW(1) ./ Eig.MC.im_l{g}.DistAvgPW(2);            % Divisive Margin for Straw Man Model from Eig using linear MCA
        end
        SMDivMarg_L  = [ SMDivMarg_L  , DivMarg_SMb];
    end
    if isfield(Eig.MC,'im_c_pi')
        DivMarg_SMc = zeros( 1, numel(Eig.netParams.gT)); 
        for g = 1:numel(Eig.netParams.gT)
            DivMarg_SMc(g) = Eig.MC.im_c_pi{g}.DistAvgPW(1) ./ Eig.MC.im_c_pi{g}.DistAvgPW(2);            % Divisive Margin for Straw Man Model from Eig using circular MCA
        end
    end
%     if isfield(Eig.MC,'im_c_asin')
%         DivMarg_SMd = zeros( 1, numel(Eig.netParams.gT)); 
%         for g = 1:numel(Eig.netParams.gT)
%             DivMarg_SMd(g) = Eig.MC.im_c_asin{g}.DistAvgPW(1) ./ Eig.MC.im_c_asin{g}.DistAvgPW(2);            % Divisive Margin for Straw Man Model from Eig using circular MCA
%         end
%     end
    %
    %
    %
    % Check that Strawman Divisive Margin is the same no matter which of 4 ways you compute it.
    disp('Strawman')
    Straw_LvC_diffMn(i) = mean(abs([DivMarg_SMa-DivMarg_SMb, DivMarg_SMa-DivMarg_SMc]));
    if ( mean(abs([DivMarg_SMa-DivMarg_SMb, DivMarg_SMa-DivMarg_SMc])) > 0.01 )
        disp('Straw man DivMarg doesnt match initial phase embedding DivMarg.')
        [DivMarg_SMa;DivMarg_SMb;DivMarg_SMc]
%         keyboard
    end
    %
    % Using 1st Eigenvector Only.
    %
    if isfield(Eig.MC,'ev1_l')
        DivMarg_EV1_l = zeros( 1, numel(Eig.netParams.gT)); 
        for g = 1:numel(Eig.netParams.gT)
            DivMarg_EV1_l(g) = Eig.MC.ev1_l{g}.DistAvgPW(1) ./ Eig.MC.ev1_l{g}.DistAvgPW(2);        % Divisive Margin for 1st Eigenvector (linear computation)
        end
        Eig1DivMarg_L = [ Eig1DivMarg_L , DivMarg_EV1_l ];
    end
    if isfield(Eig.MC,'ev1v_l')
        DivMarg_EV1v_l = zeros( 1, numel(Eig.netParams.gT)); 
        for g = 1:numel(Eig.netParams.gT)
            DivMarg_EV1v_l(g) = Eig.MC.ev1v_l{g}.DistAvgPW(1) ./ Eig.MC.ev1v_l{g}.DistAvgPW(2);     % Divisive Margin for 1st Eigenvector with nonlinear visualization (linear computation)
        end
        Eig1vDivMarg_L = [ Eig1vDivMarg_L , DivMarg_EV1v_l ];
    end
    if isfield(Eig.MC,'ev1_c_pi')
        DivMarg_EV1_c_pi = zeros( 1, numel(Eig.netParams.gT)); 
        for g = 1:numel(Eig.netParams.gT)
            DivMarg_EV1_c_pi(g) = Eig.MC.ev1_c_pi{g}.DistAvgPW(1) ./ Eig.MC.ev1_c_pi{g}.DistAvgPW(2);        % Divisive Margin for 1st Eigenvector (circular computation)
        end
        Eig1DivMarg_C = [ Eig1DivMarg_C , DivMarg_EV1_c_pi ];
    end
%     if isfield(Eig.MC,'ev1_c_asin')
%         DivMarg_EV1_c_asin = zeros( 1, numel(Eig.netParams.gT)); 
%         for g = 1:numel(Eig.netParams.gT)
%             DivMarg_EV1_c_asin(g) = Eig.MC.ev1_c_asin{g}.DistAvgPW(1) ./ Eig.MC.ev1_c_asin{g}.DistAvgPW(2);        % Divisive Margin for 1st Eigenvector (circular computation)
%         end
%     end
    if isfield(Eig.MC,'ev1v_c_pi')
        DivMarg_EV1v_c_pi = zeros( 1, numel(Eig.netParams.gT)); 
        for g = 1:numel(Eig.netParams.gT)
            DivMarg_EV1v_c_pi(g) = Eig.MC.ev1v_c_pi{g}.DistAvgPW(1) ./ Eig.MC.ev1v_c_pi{g}.DistAvgPW(2);     % Divisive Margin for 1st Eigenvector with nonlinear visualization (circular computation)
        end
        Eig1vDivMarg_C = [ Eig1vDivMarg_C , DivMarg_EV1_c_pi ];
    end
%     if isfield(Eig.MC,'ev1v_c_asin')
%         DivMarg_EV1v_c_asin = zeros( 1, numel(Eig.netParams.gT)); 
%         for g = 1:numel(Eig.netParams.gT)
%             DivMarg_EV1v_c_asin(g) = Eig.MC.ev1v_c_asin{g}.DistAvgPW(1) ./ Eig.MC.ev1v_c_asin{g}.DistAvgPW(2);     % Divisive Margin for 1st Eigenvector with nonlinear visualization (circular computation)
%         end
%     end
    %
    %
    %
    % Check that linear & circular eigenvector 1 divmarg computations match
    disp('Evec1 Viz. Pi Embedding.')
    Evec1v_LvC_diffMn(i) = mean(abs(DivMarg_EV1v_c_pi-DivMarg_EV1v_l));
    if ( mean(abs(DivMarg_EV1v_c_pi-DivMarg_EV1v_l)) > 0.01 )
        disp('DivMarg results on Evec1viz different using linear vs. circular computation.')
        [DivMarg_EV1v_c_pi;DivMarg_EV1v_l]
%         keyboard
    end
    %
%     disp('Evec1 Viz. ASin Embedding')
%     mean(abs(DivMarg_EV1v_c_asin-DivMarg_EV1v_l))
%     if ( mean(abs(DivMarg_EV1v_c_asin-DivMarg_EV1v_l)) > 0.01 )
%         disp('DivMarg results on Evec1viz different using linear vs. circular computation.')
%         %[DivMarg_EV1v_c_asin;DivMarg_EV1v_l]
%         %keyboard
%     end
    %
    
    
    
    
    disp('Evec1. Pi Embedding')
    Evec1_LvC_diffMn(i) = mean(abs(DivMarg_EV1_c_pi-DivMarg_EV1_l));
    if ( mean(abs(DivMarg_EV1_c_pi-DivMarg_EV1_l)) > 0.01 )
        disp('DivMarg results on Eigenvector1 different using linear vs. circular computation.')
       [DivMarg_EV1_c_pi;DivMarg_EV1_l]
%         figure, hist(Eig.EVecsML(:,1))
%         title('Histogram of values in Eigenvector 1.')         
%         keyboard
    end
    %
%     disp('Evec1. Asin Embedding')
%     mean(abs(DivMarg_EV1_c_asin-DivMarg_EV1_l))
%     if ( mean(abs(DivMarg_EV1_c_asin-DivMarg_EV1_l)) > 0.01 )
%         disp('DivMarg results on Eigenvector1 different using linear vs. circular computation.')
% %        [DivMarg_EV1_c_asin;DivMarg_EV1_l]
% %         figure, hist(Eig.EVecsML(:,1))
% %         title('Histogram of values in Eigenvector 1.')         
%         %keyboard
%     end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %
    % Using Only 2nd or Only 3rd Eigenvector.
    %
    if isfield(Eig.MC,'ev2o_l')
        DivMarg_EV2o = zeros( 1, numel(Eig.netParams.gT)); 
        for g = 1:numel(Eig.netParams.gT)
            DivMarg_EV2o(g) = Eig.MC.ev2o_l{g}.DistAvgPW(1) ./ Eig.MC.ev2o_l{g}.DistAvgPW(2);       % Divisive Margin for 2nd Eigenvector only
        end
        Eig2oDivMarg = [ Eig2oDivMarg , DivMarg_EV2o ];
    end
    if isfield(Eig.MC,'ev2ov_l')
        DivMarg_EV2ov = zeros( 1, numel(Eig.netParams.gT)); 
        for g = 1:numel(Eig.netParams.gT)
            DivMarg_EV2ov(g) = Eig.MC.ev2ov_l{g}.DistAvgPW(1) ./ Eig.MC.ev2ov_l{g}.DistAvgPW(2);    % Divisive Margin for 2nd Eigenvector only with nonlinear visualization
        end
        Eig2ovDivMarg = [ Eig2ovDivMarg , DivMarg_EV2ov ];
    end
    %
    if isfield(Eig.MC,'ev3o_l')
        DivMarg_EV3o = zeros( 1, numel(Eig.netParams.gT)); 
        for g = 1:numel(Eig.netParams.gT)
            DivMarg_EV3o(g) = Eig.MC.ev3o_l{g}.DistAvgPW(1) ./ Eig.MC.ev3o_l{g}.DistAvgPW(2);       % Divisive Margin for 3rd Eigenvector only
        end
        Eig3oDivMarg = [ Eig3oDivMarg , DivMarg_EV3o ];
    
    end
    if isfield(Eig.MC,'ev3ov_l')
        DivMarg_EV3ov = zeros( 1, numel(Eig.netParams.gT)); 
        for g = 1:numel(Eig.netParams.gT)
            DivMarg_EV3ov(g) = Eig.MC.ev3ov_l{g}.DistAvgPW(1) ./ Eig.MC.ev3ov_l{g}.DistAvgPW(2);    % Divisive Margin for 3rd Eigenvector only with nonlinear visualization
        end
        Eig3ovDivMarg = [ Eig3ovDivMarg , DivMarg_EV3ov ];
    end
    %
    % Using 1st & 2nd Eigenvector together.
    %
    if isfield(Eig.MC,'ev2_l')
        DivMarg_EV2 = zeros( 1, numel(Eig.netParams.gT)); 
        for g = 1:numel(Eig.netParams.gT)
            DivMarg_EV2(g) = Eig.MC.ev2_l{g}.DistAvgPW(1) ./ Eig.MC.ev2_l{g}.DistAvgPW(2);          % Divisive Margin for Eigenvectors 1&2
        end
        Eig2DivMarg = [ Eig2DivMarg , DivMarg_EV2 ];
    end
    if isfield(Eig.MC,'ev2v_l')
        DivMarg_EV2v = zeros( 1, numel(Eig.netParams.gT)); 
        for g = 1:numel(Eig.netParams.gT)
            DivMarg_EV2v(g) = Eig.MC.ev2v_l{g}.DistAvgPW(1) ./ Eig.MC.ev2v_l{g}.DistAvgPW(2);       % Divisive Margin for Eigenvectors 1&2 with nonlinear vizualization
        end
        Eig2vDivMarg = [ Eig2vDivMarg , DivMarg_EV2v ];
    end
    if isfield(Eig.MC,'ev2w_l')
        DivMarg_EV2w = zeros( 1, numel(Eig.netParams.gT)); 
        for g = 1:numel(Eig.netParams.gT)
            DivMarg_EV2w(g) = Eig.MC.ev2w_l{g}.DistAvgPW(1) ./ Eig.MC.ev2w_l{g}.DistAvgPW(2);       % Divisive Margin for Eigenvectors 1&2 weighted by eigenvectors
        end
        Eig2wDivMarg = [ Eig2wDivMarg , DivMarg_EV2w ];
    end
    if isfield(Eig.MC,'ev2wv_l')
        DivMarg_EV2wv = zeros( 1, numel(Eig.netParams.gT)); 
        for g = 1:numel(Eig.netParams.gT)
            DivMarg_EV2wv(g) = Eig.MC.ev2wv_l{g}.DistAvgPW(1) ./ Eig.MC.ev2wv_l{g}.DistAvgPW(2);    % Divisive Margin for Eigenvectors 1&2 with nonlinear vizualization weighted by eigenvectors
        end
        Eig2wvDivMarg = [ Eig2wvDivMarg , DivMarg_EV2wv ];
    end
    %
    % Using 1st thru 3rd Eigenvectors together.
    %
    if isfield(Eig.MC,'ev3_l')
        DivMarg_EV3 = zeros( 1, numel(Eig.netParams.gT)); 
        for g = 1:numel(Eig.netParams.gT)
            DivMarg_EV3(g) = Eig.MC.ev3_l{g}.DistAvgPW(1) ./ Eig.MC.ev3_l{g}.DistAvgPW(2);          % Divisive Margin for Eigenvectors 1-3
        end
        Eig3DivMarg = [ Eig3DivMarg , DivMarg_EV3 ];
    end
    if isfield(Eig.MC,'ev3v_l')
        DivMarg_EV3v = zeros( 1, numel(Eig.netParams.gT)); 
        for g = 1:numel(Eig.netParams.gT)
            DivMarg_EV3v(g) = Eig.MC.ev3v_l{g}.DistAvgPW(1) ./ Eig.MC.ev3v_l{g}.DistAvgPW(2);       % Divisive Margin for Eigenvectors 1-3 with nonlinear vizualization
        end
        Eig3vDivMarg = [ Eig3vDivMarg , DivMarg_EV3v ];
    end
    if isfield(Eig.MC,'ev3w_l')
        DivMarg_EV3w = zeros( 1, numel(Eig.netParams.gT)); 
        for g = 1:numel(Eig.netParams.gT)
            DivMarg_EV3w(g) = Eig.MC.ev3w_l{g}.DistAvgPW(1) ./ Eig.MC.ev3w_l{g}.DistAvgPW(2);       % Divisive Margin for Eigenvectors 1-3 weighted by eigenvectors
        end
        Eig3wDivMarg = [ Eig3wDivMarg , DivMarg_EV3w ];
    end
    if isfield(Eig.MC,'ev3wv_l')
        DivMarg_EV3wv = zeros( 1, numel(Eig.netParams.gT)); 
        for g = 1:numel(Eig.netParams.gT)
            DivMarg_EV3wv(g) = Eig.MC.ev3wv_l{g}.DistAvgPW(1) ./ Eig.MC.ev3wv_l{g}.DistAvgPW(2);    % Divisive Margin for Eigenvectors 1-3 with nonlinear vizualization weighted by eigenvectors
        end
        Eig3wvDivMarg = [ Eig3wvDivMarg , DivMarg_EV3wv ];
    end
    
    
    
    
    
    % Plot some stuff...
    if(0)
        
        figure
        
        subplot(221), imagesc(ImPtch.im), colormap('bone'),title('Image')
        set(gca,'XTick',[],'YTick',[])
        xlabel(['DM_L = ',num2str(mean(DivMarg_SMa),2),' // ','DM_C = ',num2str(mean(DivMarg_SMc),2)])
        %
        KurPhaseFinl = visKurPhase_inBone(Kur.netParams.im, reshape(Kur.metaCluster.phaseAtClk(:,end),Kur.netParams.Ndims(1),Kur.netParams.Ndims(2)));        
        subplot(222), imagesc(KurPhaseFinl), colormap('bone'),title('Kur')
        set(gca,'XTick',[],'YTick',[])
        xlabel(['DM = ',num2str(mean(DivMarg_K(end,:)),2)])
        %
        subplot(223), imagesc(reshape(Eig.EVecsML(:,1),Eig.netParams.Ndims(1),Eig.netParams.Ndims(2))), colormap('bone'),title('Eig')
        set(gca,'XTick',[],'YTick',[])
        xlabel(['DM = ',num2str(mean(DivMarg_EV1_l),2)])
        %
        subplot(224), imagesc( EvecVizF(  reshape(Eig.EVecsML(:,1),Eig.netParams.Ndims(1),Eig.netParams.Ndims(2)), 1e-16 ) ), colormap('bone'),title('EigViz')
        set(gca,'XTick',[],'YTick',[])
        xlabel(['DM = ',num2str(mean(DivMarg_EV1v_l),2)])
        
        
        keyboard
        
        
    end
    
    
    
    
    if(i==100)
        
        keyboard
        
    end
    
    
    
    
    
    

end % loop over files (either Kur ones or Eig ones depending on setting of scatKur & scatEig)









toc







% Visualize difference between Circular & Linear DivMarg Computations. Histograms
if(1)
    [ys,xs] = hist(Straw_LvC_diffMn,50);
    [ye,xe] = hist(Evec1_LvC_diffMn,50);
    [yv,xv] = hist(Evec1v_LvC_diffMn,50);
    %
    h=figure;
    subplot(131), bar(xs,ys), set(gca,'Yscale','log','FontSize',16,'FontWeight','Bold'), 
    title('Strawman')
    xlabel('Diff btwn Circ & Linr')
    ylabel('# Occurances / 5000')
    subplot(132), bar(xe,ye), set(gca,'Yscale','log','FontSize',16,'FontWeight','Bold'), 
    title('Eigenvector1')
    xlabel('Diff btwn Circ & Linr')
    %ylabel('# Occurances / 5000')
    subplot(133), bar(xv,yv), set(gca,'Yscale','log','FontSize',16,'FontWeight','Bold'), 
    title('Evec1 Viz')
    xlabel('Diff btwn Circ & Linr')
    %ylabel('# Occurances / 5000')
    
    annotation('textbox', [0 0.9 1 0.1],'String', ['# Pts= ',num2str(numFiles2Loop),', \sigmaP= ',sP,...
            ', Rmax: ',rM,', \sigmaD= ',sD,', \sigmaW= ',sW, ', Kscale: ',Ks,', Tscale: ',Ts,' '], ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')
        
   saveGoodImg(h,[dirScat,'Compare_Circ_vs_Linr_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts],sizeGoodIm)
   close(h)
    
end












KurDivMarg(isnan(KurDivMarg)) = 1;        % set to 1 if DivMarg = 0/0
SMDivMarg_L(isnan(SMDivMarg_L)) = 1;
SMDivMarg_C(isnan(SMDivMarg_C)) = 1;
%
Eig1DivMarg_L(isnan(Eig1DivMarg_L)) = 1;
Eig1DivMarg_C(isnan(Eig1DivMarg_C)) = 1;
Eig1vDivMarg_L(isnan(Eig1vDivMarg_L)) = 1;
Eig1vDivMarg_C(isnan(Eig1vDivMarg_C)) = 1;
%
Eig2oDivMarg(isnan(Eig2oDivMarg)) = 1;
Eig2ovDivMarg(isnan(Eig2ovDivMarg)) = 1;
%
Eig3oDivMarg(isnan(Eig3oDivMarg)) = 1;
Eig3ovDivMarg(isnan(Eig3ovDivMarg)) = 1;
%
Eig2DivMarg(isnan(Eig2DivMarg)) = 1;
Eig2vDivMarg(isnan(Eig2vDivMarg)) = 1;
Eig2wDivMarg(isnan(Eig2wDivMarg)) = 1;
Eig2wvDivMarg(isnan(Eig2wvDivMarg)) = 1;
%
Eig3DivMarg(isnan(Eig3DivMarg)) = 1;
Eig3vDivMarg(isnan(Eig3vDivMarg)) = 1;
Eig3wDivMarg(isnan(Eig3wDivMarg)) = 1;
Eig3wvDivMarg(isnan(Eig3wvDivMarg)) = 1;










% Go through misaligned patches and find their indecies in data vectors.
misaligned_indx = []; 
for i = 1:numel(misalignedID)
    misaligned_indx = [misaligned_indx, find(strcmp(misalignedID{i},ImgPtchID))];
end
%
good_indx = setxor(1:numel(ImgPtchID),misaligned_indx);








% Calculate avg. distance of point from diagonal or some other metric on this 2D scatter point cloud.
if ~isempty(KurDivMarg)
    metric_Kur_C = mean( SMDivMarg_C(good_indx) - KurDivMarg(end,good_indx) );
    metric_Kur_L = mean( SMDivMarg_L(good_indx) - KurDivMarg(end,good_indx) );
end
%
if ~isempty(Eig1DivMarg_L)
    metric_Eig1_L = mean( SMDivMarg_L(good_indx) - Eig1DivMarg_L(good_indx) );
end
if ~isempty(Eig1DivMarg_C)
    metric_Eig1_C = mean( SMDivMarg_C(good_indx) - Eig1DivMarg_C(good_indx) );
end
%
if ~isempty(Eig1vDivMarg_L)
    metric_Eig1v_L = mean( SMDivMarg_L(good_indx) - Eig1vDivMarg_L(good_indx) );
end
if ~isempty(Eig1vDivMarg_C)
    metric_Eig1v_C = mean( SMDivMarg_C(good_indx) - Eig1vDivMarg_C(good_indx) );
end


if(0)
    if ~isempty(Eig2oDivMarg)
        metric_Eig2o = mean( SMDivMarg_L(good_indx) - Eig2oDivMarg(good_indx) );
    end
    if ~isempty(Eig2ovDivMarg)
        metric_Eig2ov = mean( SMDivMarg_L(good_indx) - Eig2ovDivMarg(good_indx) );
    end
    %
    if ~isempty(Eig3oDivMarg)
        metric_Eig3o = mean( SMDivMarg_L(good_indx) - Eig3oDivMarg(good_indx) );
    end
    if ~isempty(Eig3ovDivMarg)
        metric_Eig3ov = mean( SMDivMarg_L(good_indx) - Eig3ovDivMarg(good_indx) );
    end

    if ~isempty(Eig2DivMarg)
        metric_Eig2 = mean( SMDivMarg_L(good_indx) - Eig2DivMarg(good_indx) );
    end
    if ~isempty(Eig2vDivMarg)
        metric_Eig2v = mean( SMDivMarg_L(good_indx) - Eig2vDivMarg(good_indx) );
    end
    if ~isempty(Eig2wDivMarg)
        metric_Eig2w = mean( SMDivMarg_L(good_indx) - Eig2wDivMarg(good_indx) );
    end
    if ~isempty(Eig2wvDivMarg)
        metric_Eig2wv = mean( SMDivMarg_L(good_indx) - Eig2wvDivMarg(good_indx) );
    end
    %
    if ~isempty(Eig3DivMarg)
        metric_Eig3 = mean( SMDivMarg_L(good_indx) - Eig3DivMarg(good_indx) );
    end
    if ~isempty(Eig3vDivMarg)
        metric_Eig3v = mean( SMDivMarg_L(good_indx) - Eig3vDivMarg(good_indx) );
    end
    if ~isempty(Eig3wDivMarg)
        metric_Eig3w = mean( SMDivMarg_L(good_indx) - Eig3wDivMarg(good_indx) );
    end
    if ~isempty(Eig3wvDivMarg)
        metric_Eig3wv = mean( SMDivMarg_L(good_indx) - Eig3wvDivMarg(good_indx) );
    end
end



%% Now Plot Scatter Plot of Straw vs Eigen & Straw vs. Kuramoto.
if(plot2x2Dscatter)

    % (1). Plot 2D scatter for 1st Eigenvector with and without Visualization Nonlinearity.
    if(1)
        hSc = figure;   
        %
        subplot(231)
        scatter( SMDivMarg_L(good_indx) , KurDivMarg(end,good_indx) ), hold on  
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel({'Linear','DivMarg Kur'},'FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        title('Kuramoto','FontSize',20,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Kur_L,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        subplot(232) 
        scatter( SMDivMarg_L(good_indx) , Eig1DivMarg_L(good_indx) ), hold on
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg Eig','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        title('Evecs 1','FontSize',20,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Eig1_L,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        subplot(233)
        scatter( SMDivMarg_L(good_indx) , Eig1vDivMarg_L(good_indx) ), hold on
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg Eig ','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        title(['w/ vizNonLin=',num2str(Eig.MC.vizNonLin,2)],'FontSize',20,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Eig1v_L,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        %
        subplot(234)
        scatter( SMDivMarg_C(good_indx) , KurDivMarg(end,good_indx) ), hold on  
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel({'Circular','DivMarg Kur'},'FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        title('Kuramoto','FontSize',20,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Kur_C,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        subplot(235) 
        scatter( SMDivMarg_C(good_indx) , Eig1DivMarg_C(good_indx) ), hold on
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg Eig','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        title('Evecs 1','FontSize',20,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Eig1_C,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        subplot(236)
        scatter( SMDivMarg_C(good_indx) , Eig1vDivMarg_C(good_indx) ), hold on
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg Eig ','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        title(['w/ vizNonLin=',num2str(Eig.MC.vizNonLin,2)],'FontSize',20,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Eig1v_C,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        annotation('textbox', [0 0.9 1 0.1],'String', ...
           ['On ',fileSize,' ',fileGeneral,' # Pts= ',num2str(numel(ImgPtchID)-numel(misaligned_indx)),' / ',num2str(numel(ImgPtchID)),', \sigmaP= ',num2str(Kur.netParams.sigP),', Rmax: ',num2str(Kur.netParams.Rmax),...
           ', \sigmaD= ',num2str(Kur.netParams.sigD),', \sigmaW= ',num2str(Kur.kurParams.sigW), ', Kscale: ',num2str(Kur.kurParams.Kscale),', Tscale: ',num2str(Kur.kurParams.TiScale),' '], ...
           'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')
        %
        saveGoodImg(hSc,[dirScat,'DivMarg_Scatter_2D_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'_Evecs1'],sizeGoodIm)
        close(hSc)
    end
    
    % (2). Plot 2D scatter for Eigenvector 2 with and without Visualization Nonlinearity.
    if(0)
        hSc = figure;   
        %
        subplot(131)
        scatter( SMDivMarg_C(good_indx) , KurDivMarg(end,good_indx) ), hold on  
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg Kur','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        title('Kuramoto','FontSize',20,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Kur,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        subplot(132)
        scatter( SMDivMarg_L(good_indx) , Eig2DivMarg(good_indx) ), hold on
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg Eig','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        title('Evecs 1-2','FontSize',20,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Eig2,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        subplot(133)
        scatter( SMDivMarg_L(good_indx) , Eig2vDivMarg(good_indx) ), hold on
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg Eig ','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        title(['w/ vizNonLin=',num2str(Eig.MC.vizNonLin,2)],'FontSize',20,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Eig2v,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        annotation('textbox', [0 0.9 1 0.1],'String', ...
           ['On ',fileSize,' ',fileGeneral,' # Pts= ',num2str(numel(ImgPtchID)-numel(misaligned_indx)),' / ',num2str(numel(ImgPtchID)),', \sigmaP= ',num2str(Kur.netParams.sigP),', Rmax: ',num2str(Kur.netParams.Rmax),...
           ', \sigmaD= ',num2str(Kur.netParams.sigD),', \sigmaW= ',num2str(Kur.kurParams.sigW), ', Kscale: ',num2str(Kur.kurParams.Kscale),', Tscale: ',num2str(Kur.kurParams.TiScale),' '], ...
           'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')
        %
        saveGoodImg(hSc,[dirScat,'DivMarg_Scatter_2D_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'_Evecs1-2'],sizeGoodIm)
        close(hSc)
    end
    
    % (3). Plot 2D scatter for Eigenvector 3 with and without Visualization Nonlinearity.
    if(0)
        hSc = figure;   
        %
        subplot(131)
        scatter( SMDivMarg_C(good_indx) , KurDivMarg(end,good_indx) ), hold on  
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg Kur','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        title('Kuramoto','FontSize',20,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Kur,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        subplot(132)
        scatter( SMDivMarg_L(good_indx) , Eig2DivMarg(good_indx) ), hold on
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg Eig','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        title('Evecs 1-2','FontSize',20,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Eig2,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        subplot(133)
        scatter( SMDivMarg_L(good_indx) , Eig2vDivMarg(good_indx) ), hold on
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg Eig ','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        title(['w/ vizNonLin=',num2str(Eig.MC.vizNonLin,2)],'FontSize',20,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Eig2v,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        annotation('textbox', [0 0.9 1 0.1],'String', ...
           ['On ',fileSize,' ',fileGeneral,' # Pts= ',num2str(numel(ImgPtchID)-numel(misaligned_indx)),' / ',num2str(numel(ImgPtchID)),', \sigmaP= ',num2str(Kur.netParams.sigP),', Rmax: ',num2str(Kur.netParams.Rmax),...
           ', \sigmaD= ',num2str(Kur.netParams.sigD),', \sigmaW= ',num2str(Kur.kurParams.sigW), ', Kscale: ',num2str(Kur.kurParams.Kscale),', Tscale: ',num2str(Kur.kurParams.TiScale),' '], ...
           'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')
        %
        saveGoodImg(hSc,[dirScat,'DivMarg_Scatter_2D_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'_Evecs1-2'],sizeGoodIm)
        close(hSc)
    end
    
    % (4). Plot 2D scatter for Eigenvectors 1-2 with and without Visualization Nonlinearity.
    if(0)
        hSc = figure;   
        %
        subplot(131)
        scatter( SMDivMarg_C(good_indx) , KurDivMarg(end,good_indx) ), hold on  
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg Kur','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        title('Kuramoto','FontSize',20,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Kur,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        subplot(132)
        scatter( SMDivMarg_L(good_indx) , Eig2DivMarg(good_indx) ), hold on
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg Eig','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        title('Evecs 1-2','FontSize',20,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Eig2,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        subplot(133)
        scatter( SMDivMarg_L(good_indx) , Eig2vDivMarg(good_indx) ), hold on
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg Eig ','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        title(['w/ vizNonLin=',num2str(Eig.MC.vizNonLin,2)],'FontSize',20,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Eig2v,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        annotation('textbox', [0 0.9 1 0.1],'String', ...
           ['On ',fileSize,' ',fileGeneral,' # Pts= ',num2str(numel(ImgPtchID)-numel(misaligned_indx)),' / ',num2str(numel(ImgPtchID)),', \sigmaP= ',num2str(Kur.netParams.sigP),', Rmax: ',num2str(Kur.netParams.Rmax),...
           ', \sigmaD= ',num2str(Kur.netParams.sigD),', \sigmaW= ',num2str(Kur.kurParams.sigW), ', Kscale: ',num2str(Kur.kurParams.Kscale),', Tscale: ',num2str(Kur.kurParams.TiScale),' '], ...
           'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')
        %
        saveGoodImg(hSc,[dirScat,'DivMarg_Scatter_2D_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'_Evecs1-2'],sizeGoodIm)
        close(hSc)
    end

    
    % (5). Plot 2D scatter for Eigenvectors 1-3 with and without Visualization Nonlinearity.
    if(0)
        hSc = figure;   
        %
        subplot(131)
        scatter( SMDivMarg(good_indx) , KurDivMarg(end,good_indx) ), hold on  
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg Kur','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        title('Kuramoto','FontSize',20,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Kur,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        subplot(132)
        scatter( SMDivMarg(good_indx) , Eig3DivMarg(good_indx) ), hold on
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg Eig','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        title('Evecs 1-3','FontSize',20,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Eig3,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        subplot(133) 
        scatter( SMDivMarg(good_indx) , Eig3vDivMarg(good_indx) ), hold on
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg Eig ','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        title(['w/ vizNonLin=',num2str(Eig.MC.vizNonLin,2)],'FontSize',20,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Eig3v,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        annotation('textbox', [0 0.9 1 0.1],'String', ...
           ['On ',fileSize,' ',fileGeneral,' # Pts= ',num2str(numel(ImgPtchID)-numel(misaligned_indx)),' / ',num2str(numel(ImgPtchID)),', \sigmaP= ',num2str(Kur.netParams.sigP),', Rmax: ',num2str(Kur.netParams.Rmax),...
           ', \sigmaD= ',num2str(Kur.netParams.sigD),', \sigmaW= ',num2str(Kur.kurParams.sigW), ', Kscale: ',num2str(Kur.kurParams.Kscale),', Tscale: ',num2str(Kur.kurParams.TiScale),' '], ...
           'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')
        %
        saveGoodImg(hSc,[dirScat,'DivMarg_Scatter_2D_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'_Evecs1-3'],sizeGoodIm)
        close(hSc)
    end
    
    
    % (6). Plot 2D scatter for Eigenvectors 1-2 weighted by their Eigenvalues with and without Visualization Nonlinearity.
    if(0)
        hSc = figure;   
        %
        subplot(131)
        scatter( SMDivMarg(good_indx) , KurDivMarg(end,good_indx) ), hold on  
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg Kur','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        title('Kuramoto','FontSize',20,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Kur,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        subplot(132)
        scatter( SMDivMarg(good_indx) , Eig2wDivMarg(good_indx) ), hold on
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg Eig','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        title('Weighted Evecs 1-2','FontSize',20,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Eig2w,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        subplot(133)
        scatter( SMDivMarg(good_indx) , Eig2wvDivMarg(good_indx) ), hold on
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg Eig ','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        title(['w/ vizNonLin=',num2str(Eig.MC.vizNonLin,2)],'FontSize',20,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Eig2wv,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        annotation('textbox', [0 0.9 1 0.1],'String', ...
           ['On ',fileSize,' ',fileGeneral,' # Pts= ',num2str(numel(ImgPtchID)-numel(misaligned_indx)),' / ',num2str(numel(ImgPtchID)),', \sigmaP= ',num2str(Kur.netParams.sigP),', Rmax: ',num2str(Kur.netParams.Rmax),...
           ', \sigmaD= ',num2str(Kur.netParams.sigD),', \sigmaW= ',num2str(Kur.kurParams.sigW), ', Kscale: ',num2str(Kur.kurParams.Kscale),', Tscale: ',num2str(Kur.kurParams.TiScale),' '], ...
           'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')
        %
        saveGoodImg(hSc,[dirScat,'DivMarg_Scatter_2D_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'_Evecs1-2w'],sizeGoodIm)
        close(hSc)
    end

    
    % (7). Plot 2D scatter for Eigenvectors 1-3 weighted by their Eigenvalues with and without Visualization Nonlinearity.
    if(0)
        hSc = figure;   
        %
        subplot(131)
        scatter( SMDivMarg(good_indx) , KurDivMarg(end,good_indx) ), hold on  
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg Kur','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        title('Kuramoto','FontSize',20,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Kur,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        subplot(132)
        scatter( SMDivMarg(good_indx) , Eig3wDivMarg(good_indx) ), hold on
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg Eig','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        title('Weighted Evecs 1-3','FontSize',20,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Eig3w,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        subplot(133) 
        scatter( SMDivMarg(good_indx) , Eig3wvDivMarg(good_indx) ), hold on
        plot([0 1.1], [0 1.1],'k--','LineWidth',1.5)
        axis([0 1.1 0 1.1])
        axis square
        set(gca,'FontSize',16,'FontWeight','Bold')
        ylabel('DivMarg Eig ','FontSize',18,'FontWeight','Bold')
        xlabel('DivMarg Img','FontSize',18,'FontWeight','Bold')
        title(['w/ vizNonLin=',num2str(Eig.MC.vizNonLin,2)],'FontSize',20,'FontWeight','Bold')
        text(0.9,0.1,'\color{green}{Better}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.1,0.9,'\color{red}{Worse}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.2,0.2,'\color{black}{Neutral}','FontSize',16,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        text(0.5,1.05,['\color{red}{<metric>=',num2str(metric_Eig3wv,2),'}'],'FontSize',18,'FontWeight','Bold','BackgroundColor','White','HorizontalAlignment','center')
        %
        annotation('textbox', [0 0.9 1 0.1],'String', ...
           ['On ',fileSize,' ',fileGeneral,' # Pts= ',num2str(numel(ImgPtchID)-numel(misaligned_indx)),' / ',num2str(numel(ImgPtchID)),', \sigmaP= ',num2str(Kur.netParams.sigP),', Rmax: ',num2str(Kur.netParams.Rmax),...
           ', \sigmaD= ',num2str(Kur.netParams.sigD),', \sigmaW= ',num2str(Kur.kurParams.sigW), ', Kscale: ',num2str(Kur.kurParams.Kscale),', Tscale: ',num2str(Kur.kurParams.TiScale),' '], ...
           'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')
        %
        saveGoodImg(hSc,[dirScat,'DivMarg_Scatter_2D_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'_Evecs1-3w'],sizeGoodIm)
        close(hSc)
    end

    
end
    


% Save a mat-file with important information - will be used in visualize_pts_from_2x2D_scatter function.
if(1)
    
    save([dirDMmat,'DivMarg_Results_allPatches_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts],...
        'KurDivMarg','SMDivMarg_L','SMDivMarg_C','Eig1DivMarg_L','Eig1DivMarg_C','Eig1vDivMarg_L',...
        'Eig1vDivMarg_C','Eig2oDivMarg','Eig2ovDivMarg','Eig3oDivMarg','Eig3ovDivMarg','Eig2DivMarg',...
        'Eig2vDivMarg','Eig2wDivMarg','Eig2wvDivMarg','Eig3DivMarg','Eig3vDivMarg','Eig3wDivMarg',...
        'Eig3wvDivMarg','misaligned_indx','misalignedID','good_indx','ImgPtchID','ImgGtID','sP','rM','sD','sW','Ks','Ts')

end
