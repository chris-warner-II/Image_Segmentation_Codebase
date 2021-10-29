function [] = explore_compare_DistPairs(fileGeneral,fileSize,methodType,rM,sD,sP,sW,Ts,Ks)


% This function will take output mat file from main_Kuramoto function that is titled
% metaClusterKur_...  This will plot the phase of each oscillator at a 60Hz
% clock and will colorcode them by ground truth segment.
%
%
% Inputs to function must be of this form
%
% fileGeneral = 'GradientBox'; % 'BSDS_patch'; % 
% fileSize = '51x51_ds1';
% methodType = 'GLnrm'; % 'Mod_SKHAdj'

plot2x2Dscatter = 0; % flag to make plot of 2x2D scatter Kur vs Straw & Eig vs Straw for each individual patch. 
                     % If this flag is zero, then this function is just creating the DistPW_Results_allPatches mat file.
                     % Note: I do plotting in the visualize_pts_from_2x2D_scatterB & scatter_plot_DistPairs functions

                     
                     
% Directories to Input Images or Data :: Change these paths below to get a different images or patches
[dirPre,sizeGoodIm] = onCluster;
dirKur = [dirPre,'output/Kuramoto/NetsFromImgs/',fileGeneral,'_',fileSize,'/data/Kur_PIF_Fourier1_fixed/',methodType,'/']; %Theta Scale ',tscale,'/'];
dirEig = [dirPre,'output/Kuramoto/NetsFromImgs/',fileGeneral,'_',fileSize,'/data/spectral_fixedT/',methodType,'/'];
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

% [dirKur,'KurMC','*rM',rM,'*sD',sD,'*sP',sP,'*NF_60_',sW,'*kscale',Ks,'*tscale',Ts,'*.mat']
% [dirEig,'Evecs','*rM',rM,'_sD',sD,'_sP',sP,'*.mat']

filesKur = dir([dirKur,'KurMC','*rM',rM,'_sD',sD,'_sP',sP,'_NF_60_',sW,'_kscale',Ks,'_tscale',Ts,'*.mat']); % can loop through these later
filesEig = dir([dirEig,'Evecs','*rM',rM,'_sD',sD,'_sP',sP,'*.mat']); % will be in a different directory later.



% Directory to save 2D Scatter plots into
dirScat = [dirPre,'output/Kuramoto/NetsFromImgs/',fileGeneral,'_',fileSize,'/imgs/DistPairs/',methodType,'/'];
if ~exist(dirScat,'dir')
    mkdir(dirScat)
end

% Directory to save mat file containing location of each point along with identity of image patch where it came from.
dirDMmat = [dirPre,'output/Kuramoto/NetsFromImgs/',fileGeneral,'_',fileSize,'/data/DistPairs/',methodType,'/'];
if ~exist(dirDMmat,'dir')
    mkdir(dirDMmat)
end



    

    
disp('Time to loop through all metaCluster mat files:')
tic




% Preset data arrays to hold results across all BSDS image patches.
kkk = 0;         % counter when looping through filesKur & ground truths.
mis = 0;         % counter for misaligned Kur/Straw combinations.
misalignedID={}; % index into all vectors to be saved to find misaligned pairs


somethingWeirdWithEig = 0;


%
% KurDistsPW = [];    % Divisive Margin from Kuramoto Coupled Oscillator final Phase Distribution
% %
% SMDistsPW_L  = [];    % Divisive Margin from StrawMan (Image Pixel Values)
% % SMDistsPW_C  = [];    % Divisive Margin from StrawMan (Image Pixel Values)
% %
% Eig1DistsPW_L = [];   % Divisive Margin from 1st Eigenvector (straight up)
% % Eig1DistsPW_C = [];   % Divisive Margin from 1st Eigenvector (straight up)
% %
% % Eig1vDistsPW_L = [];  % Divisive Margin from 1st Eigenvector with NonLinear Visualization applied.
% % Eig1vDistsPW_C = [];  % Divisive Margin from 1st Eigenvector with NonLinear Visualization applied.
% %
% Eig2oDistsPW = [];  % Divisive Margin from 2nd Eigenvector only.
% % Eig2ovDistsPW = []; % Divisive Margin from 2nd Eigenvector only with NonLinear Visualization applied.
% %
% Eig3oDistsPW = [];  % Divisive Margin from 3rd Eigenvector only.
% % Eig3ovDistsPW = []; % Divisive Margin from 3rd Eigenvector only with NonLinear Visualization applied.
% %  
% Eig2DistsPW = [];   % Divisive Margin from Eigenvectors 1-2.
% % Eig2vDistsPW = [];  % Divisive Margin from Eigenvectors 1-2 with NonLinear Visualization applied.
% Eig2wDistsPW = [];  % Divisive Margin from Eigenvectors 1-2 weighted by their Eigenvectors. 
% % Eig2wvDistsPW = []; % Divisive Margin from Eigenvectors 1-2 weighted by their Eigenvectors with NonLinear Visualization applied. 
% %
% Eig3DistsPW = [];   % Divisive Margin from Eigenvectors 1-3.
% % Eig3vDistsPW = [];  % Divisive Margin from Eigenvectors 1-3 with NonLinear Visualization applied.
% Eig3wDistsPW = [];  % Divisive Margin from Eigenvectors 1-3 weighted by their Eigenvectors.
% % Eig3wvDistsPW = []; % Divisive Margin from Eigenvectors 1-3 weighted by their Eigenvectors with NonLinear Visualization applied. 
% %

% mean location of each cluster in each ground truth
mnCluster_all.kur = []; 
mnCluster_all.im = [];
mnCluster_all.ev1 = [];
mnCluster_all.ev2o = [];
mnCluster_all.ev3o = [];
mnCluster_all.ev2 = [];
mnCluster_all.ev3 = [];
mnCluster_all.ev2w = [];
mnCluster_all.ev3w = [];

% Confidence Intervals for mean location of each cluster
mnClusterCI_all.kur = []; 
mnClusterCI_all.im = [];
mnClusterCI_all.ev1 = [];
mnClusterCI_all.ev2o = [];
mnClusterCI_all.ev3o = [];
mnClusterCI_all.ev2 = [];
mnClusterCI_all.ev3 = [];
mnClusterCI_all.ev2w = [];
mnClusterCI_all.ev3w = [];

% standard deviation of each cluster in each ground truth
stdCluster_all.kur = []; 
stdCluster_all.im = [];
stdCluster_all.ev1 = [];
stdCluster_all.ev2o = [];
stdCluster_all.ev3o = [];
stdCluster_all.ev2 = [];
stdCluster_all.ev3 = [];
stdCluster_all.ev2w = [];
stdCluster_all.ev3w = [];

% Confidence Intervals standard deviation of each cluster in each ground truth
stdClusterCI_all.kur = []; 
stdClusterCI_all.im = [];
stdClusterCI_all.ev1 = [];
stdClusterCI_all.ev2o = [];
stdClusterCI_all.ev3o = [];
stdClusterCI_all.ev2 = [];
stdClusterCI_all.ev3 = [];
stdClusterCI_all.ev2w = [];
stdClusterCI_all.ev3w = [];

% identity number and size of each cluster in each ground truth
idCluster_all = [];
sizeCluster_all = [];



% Decide which files to loop through
if (scatEig & ~scatKur)
    numFiles2Loop = numel(filesEig);
else
    numFiles2Loop = numel(filesKur);
end


% numFiles2Loop = 200; % get rid of this.




% If output mat file 'DistsPW_Results_allPatches' already exists in dirDMmat directory with same 
% number of files or more files as this run, break out.
f = dir([dirDMmat,'Clustering_Results_allPatches_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'*.mat']);
if ~isempty(f)
    for i = 1:numel(f)
        f(i).name
        beg = strfind(f(i).name,'files')+5;
        fin = strfind(f(i).name,'.mat')-1;
        numFilesPrev(i) = str2double(f(i).name(beg:fin));
    end
    %
    if( max(numFilesPrev) >= numFiles2Loop )
        disp('There already Exists an Output File. Not Rerunning...')
        return
    end
end



% % This will break out and not run if you are missing ANY KurMC or Evecs files.
% if( numel(filesEig)~=1500 || numel(filesKur)~=1500 ) 
%     disp('Not all 1500 KurMC & Evecs files produced yet.  Waiting on next step.')
%     return
% end



% % Preallocate arrays to hold avg difference between Linear and Circular DistsPW computation for Straw, Evec1 & Evec1viz.
% if(0)
%     Straw_LvC_diffMn = zeros(1,numFiles2Loop);
%     Evec1_LvC_diffMn = zeros(1,numFiles2Loop);
%     Evec1v_LvC_diffMn = zeros(1,numFiles2Loop);
% end




% Preallocate memory to compute range of each variable utelized.
%RangeUsedEvecs = zeros(numFiles2Loop,6);
%
% RangeHistogramPlotFlag = 0;

for i = 1:numFiles2Loop % I dont know if this will work if I am not looping thru filesKur.

    disp([num2str(i),' / ',num2str(numFiles2Loop)])
    
    
    if (scatEig & scatKur)
        
        % Load file containing metaClustering Analysis results for Kuramoto & StrawMan for parameter combination
        try
            %[dirKur,filesKur(i).name]
            Kur = load([dirKur,filesKur(i).name]);
            
            x=Kur.mnCluster;  % leave this here because i want it to fall into the catch clause if the mnCluster variable doesnt exist in the Kur structure.
            
%             %range! trying some stuff...
%             if(RangeHistogramPlotFlag)
%                 %
%                 numNodes2calcRange = 100;
%                 nodes2calcRange = unique(randi([1,Kur.netParams.N],1,numNodes2calcRange));
%                 Kur_dists = zeros(Kur.netParams.N,numel(nodes2calcRange));
%                 %
%                 for pp=1:numel(nodes2calcRange)
%                     Kur_dists(:,pp) = circ_dist(Kur.metaCluster.phaseAtClk(nodes2calcRange(pp),end),Kur.metaCluster.phaseAtClk(:,end));
%                 end
%                 %
%                 figure, hist(abs(Kur_dists(:)),100),
%                 title(['Hi']) % Kur.netParams.C
%                 
%                 %
%                 keyboard
%             end
            
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
            
            %RangeUsedEvecs(i,:) = range(Eig.EVecsML);
            
%             %range! trying some stuff
%             if(RangeHistogramPlotFlag)
%                 %
%                 numNodes2calcRange = 100;
%                 nodes2calcRange = unique(randi([1,Eig.netParams.N],1,numNodes2calcRange));
%                 E1_dists = zeros(Eig.netParams.N,numel(nodes2calcRange));
%                 E2_dists = zeros(Eig.netParams.N,numel(nodes2calcRange));
%                 E3_dists = zeros(Eig.netParams.N,numel(nodes2calcRange));
%                 %
%                 for pp=1:numel(nodes2calcRange)
%                     E1_dists(:,pp) = dist(   Eig.EVecsML(nodes2calcRange(pp),1),  Eig.EVecsML(:,1)  );
%                     E2_dists(:,pp) = dist(   Eig.EVecsML(nodes2calcRange(pp),2),  Eig.EVecsML(:,2)  );
%                     E3_dists(:,pp) = dist(   Eig.EVecsML(nodes2calcRange(pp),3),  Eig.EVecsML(:,3)  );
%                 end
%                 %
%                 figure, 
%                 subplot(311),hist(abs(E1_dists(:)),100)
%                 subplot(312),hist(abs(E2_dists(:)),100)
%                 subplot(313),hist(abs(E3_dists(:)),100)
%                 %
%                 keyboard
%             end
            
        catch
            disp('Eigen file is corrupted or not there: Skipping to next...')
            [dirEig,'Evecs_',Kur.kurflags.fname,'_rM',Kur.netflags.rM,'_sD',Kur.netflags.sD,'_sP',Kur.netflags.sP,'.mat']
            %delete([dirEig,'Evecs_',kurflags.fname,'_rM',netflags.rM,' sD',netflags.sD,' sP',netflags.sP,'.mat'])
            continue
        end
    
    end
    
    
    
    
%     % Compute Empirical Ranges for Eigenvectors and different combinations, for Kuramoto, for Image, etc...
%     EmpRange_Eig1 = range(Eig.EVecsML(:,1));
%     EmpRange_Eig2o = range(Eig.EVecsML(:,2));
%     EmpRange_Eig3o = range(Eig.EVecsML(:,3));
%     EmpRange_Eig2 = range(sum(Eig.EVecsML(:,1:2),2));
%     EmpRange_Eig3 = range(sum(Eig.EVecsML(:,1:3),2));
%     %
%     EmpRange_Eig1v = range(EvecVizF(Eig.EVecsML(:,1),Eig.MC.vizNonLin));
%     EmpRange_Eig2ov = range(EvecVizF(Eig.EVecsML(:,2),Eig.MC.vizNonLin));
%     EmpRange_Eig3ov = range(EvecVizF(Eig.EVecsML(:,3),Eig.MC.vizNonLin));
%     EmpRange_Eig2v = range(sum(EvecVizF(Eig.EVecsML(:,1:2),Eig.MC.vizNonLin),2));
%     EmpRange_Eig3v = range(sum(EvecVizF(Eig.EVecsML(:,1:3),Eig.MC.vizNonLin),2));
%     %
%     EmpRange_Kur = circ_EmpRange(Kur.metaCluster.phaseAtClk(:,end));
%     
% 
    % Loop through each human ground truth segmentation
    for g = 1:numel(Kur.netParams.gT)

        kkk = kkk + 1;
        ImgPtchID{kkk} = Kur.kurflags.fname; %[filesKur(i).name(beg:fin)];
        ImgGtID(kkk) = g;
        
        GtNumSegs(kkk) = numel(unique(Kur.netParams.gT{g}));
       
%         EmpRange.Eig1(kkk) = EmpRange_Eig1;
%         EmpRange.Eig2o(kkk) = EmpRange_Eig2o;
%         EmpRange.Eig3o(kkk) = EmpRange_Eig3o;
%         EmpRange.Eig2(kkk) = EmpRange_Eig2;
%         EmpRange.Eig3(kkk) = EmpRange_Eig3;
%         %
%         EmpRange.Eig1v(kkk) = EmpRange_Eig1v;
%         EmpRange.Eig2ov(kkk) = EmpRange_Eig2ov;
%         EmpRange.Eig3ov(kkk) = EmpRange_Eig3ov;
%         EmpRange.Eig2v(kkk) = EmpRange_Eig2v;
%         EmpRange.Eig3v(kkk) = EmpRange_Eig3v;
%         %
%         EmpRange.Kur(kkk) = EmpRange_Kur;
      
    end

    
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
   
   






    if ~isfield(Eig,'mnClusterT')
        disp('Think this was using the Empirical Mn & Std for Eig & Img...')
        Eig.mnClusterT = Eig.mnCluster;
        Eig.stdClusterT = Eig.stdCluster;
    end


    
    
    % NOTE:  Various checks to make sure that Kur and Eig files are
    % registered to be looking at same image patch and returning same
    % values when they should be.  A sanity check.  I can get rid of all or
    % parts of it because it may slow down computation considerably.
    for i = 1:numel(Kur.idCluster) % max([numel(Kur.idCluster),numel(Kur.idCluster)])

       try

           if any(Kur.idCluster{i}~=Eig.idCluster{i})
                disp('idCluster doesnt match up in Kur & Eig file')
                keyboard
           end
           %
           if any(Kur.sizeCluster{i}~=Eig.sizeCluster{i})
                disp('sizeCluster doesnt match up in Kur & Eig file')
                keyboard
           end


       catch

            % This try catch statement is because some of the Evec mat
            % files have this weird property of having the correct entries
            % and of having extra entries.  I dont know why.  I just want
            % to grab the early ones that add up to the number of
            % oscillators (N).
            disp('Try:Catch: Changing the length of Eig elements.')
%             i
%             Kur.idCluster{i}
%             Kur.sizeCluster{i}
%             Eig.idCluster{i}
%             Eig.sizeCluster{i}
%             % keyboard
%             somethingWeirdWithEig = somethingWeirdWithEig+1

%             iK = find( cumsum(Kur.sizeCluster{i}) == Kur.netParams.N );
            iE = find( cumsum(Eig.sizeCluster{i}) == Kur.netParams.N );

%             [iK, iE]

            Eig.idCluster{i} = Eig.idCluster{i}(1:iE);
            Eig.sizeCluster{i} = Eig.sizeCluster{i}(1:iE);

            Eig.mnClusterT.im{i} = Eig.mnClusterT.im{i}(1:iE);
            Eig.stdClusterT.im{i} = Eig.stdClusterT.im{i}(1:iE);
            
            Eig.mnClusterT.ev1{i} = Eig.mnClusterT.ev1{i}(1:iE);
            Eig.stdClusterT.ev1{i} = Eig.stdClusterT.ev1{i}(1:iE);
     
            Eig.mnClusterT.ev2o{i} = Eig.mnClusterT.ev2o{i}(1:iE);
            Eig.stdClusterT.ev2o{i} = Eig.stdClusterT.ev2o{i}(1:iE);
    
            Eig.mnClusterT.ev3o{i} = Eig.mnClusterT.ev3o{i}(1:iE);
            Eig.stdClusterT.ev3o{i} = Eig.stdClusterT.ev3o{i}(1:iE);

            Eig.mnClusterT.ev2{i} = Eig.mnClusterT.ev2{i}(1:iE);
            Eig.stdClusterT.ev2{i} = Eig.stdClusterT.ev2{i}(1:iE);
    
            Eig.mnClusterT.ev3{i} = Eig.mnClusterT.ev3{i}(1:iE);
            Eig.stdClusterT.ev3{i} = Eig.stdClusterT.ev3{i}(1:iE);

            Eig.mnClusterT.ev2w{i} = Eig.mnClusterT.ev2w{i}(1:iE);
            Eig.stdClusterT.ev2w{i} = Eig.stdClusterT.ev2w{i}(1:iE);
    
            Eig.mnClusterT.ev3w{i} = Eig.mnClusterT.ev3w{i}(1:iE);
            Eig.stdClusterT.ev3w{i} = Eig.stdClusterT.ev3w{i}(1:iE);



           if any(Kur.idCluster{i}~=Eig.idCluster{i})
                disp('idCluster STILL doesnt match up in Kur & Eig file')
                keyboard
           end
           %
           if any(Kur.sizeCluster{i}~=Eig.sizeCluster{i})
                disp('sizeCluster STILL doesnt match up in Kur & Eig file')
                keyboard
           end


       end


    end  
   
    
    

    
    
        
    idCluster_all = [idCluster_all, Kur.idCluster];
    sizeCluster_all = [sizeCluster_all, Kur.sizeCluster];
    
    mnCluster_all.kur = [ mnCluster_all.kur , Kur.mnCluster.kur];
    stdCluster_all.kur = [ stdCluster_all.kur , Kur.stdCluster.kur];
    
    mnCluster_all.im = [ mnCluster_all.im , Eig.mnClusterT.im];
    stdCluster_all.im = [ stdCluster_all.im , Eig.stdClusterT.im];

    
    mnCluster_all.ev1 = [mnCluster_all.ev1, Eig.mnClusterT.ev1];
    stdCluster_all.ev1 = [stdCluster_all.ev1, Eig.stdClusterT.ev1];
     
    
    mnCluster_all.ev2o = [mnCluster_all.ev2o, Eig.mnClusterT.ev2o];
    stdCluster_all.ev2o = [stdCluster_all.ev2o, Eig.stdClusterT.ev2o];
    
    
    mnCluster_all.ev3o = [mnCluster_all.ev3o, Eig.mnClusterT.ev3o];
    stdCluster_all.ev3o = [stdCluster_all.ev3o, Eig.stdClusterT.ev3o];

    
    mnCluster_all.ev2 = [mnCluster_all.ev2, Eig.mnClusterT.ev2];
    stdCluster_all.ev2 = [stdCluster_all.ev2, Eig.stdClusterT.ev2];
    
    
    mnCluster_all.ev3 = [mnCluster_all.ev3, Eig.mnClusterT.ev3];
    stdCluster_all.ev3 = [stdCluster_all.ev3, Eig.stdClusterT.ev3];

    
    mnCluster_all.ev2w = [mnCluster_all.ev2w, Eig.mnClusterT.ev2w];
    stdCluster_all.ev2w = [stdCluster_all.ev2w, Eig.stdClusterT.ev2w];
    
    
    mnCluster_all.ev3w = [mnCluster_all.ev3w, Eig.mnClusterT.ev3w];
    stdCluster_all.ev3w = [stdCluster_all.ev3w, Eig.stdClusterT.ev3w];
    
        

    
    
%     if isfield(Eig.mnCluster,'im')
% 
%         SMDistsPW_L  = [ SMDistsPW_L  , Eig.mnCluster.im];
%         
%     else
%         
%         SMDistsPW_L  = [ SMDistsPW_L  , nan( 2, numel(Eig.netParams.gT))];
%         
%     end
%     
%     
%     
% %     if isfield(Eig.MC,'im_c_pi')
% %         DistsPW_SMc = zeros( 2, numel(Eig.netParams.gT)); 
% %         for g = 1:numel(Eig.netParams.gT)
% %             DistsPW_SMc(:,g) = Eig.MC.im_c_pi{g}.DistAvgPW;         % Divisive Margin for Straw Man Model from Eig using circular MCA
% %         end
% %         
% %     else
% %         
% %         % Nothing
% %         
% %     end
%     %
%     %
%     %
%     % Check that Strawman Divisive Margin is the same no matter which of 4 ways you compute it.
% %     if(0)
% %         disp('Strawman')
% %         Straw_LvC_diffMn(i) = mean(abs([DistsPW_SMa-DistsPW_SMb, DistsPW_SMa-DistsPW_SMc]));
% %         if ( mean(abs([DistsPW_SMa-DistsPW_SMb, DistsPW_SMa-DistsPW_SMc])) > 0.01 )
% %             disp('Straw man DistsPW doesnt match initial phase embedding DistsPW.')
% %             [DistsPW_SMa;DistsPW_SMb;DistsPW_SMc]
% %     %         keyboard
% %         end
% %     end
%     %
%     % Using 1st Eigenvector Only.
%     %
%     if isfield(Eig.MC,'ev1_l')
%         DistsPW_EV1_l = zeros( 2, numel(Eig.netParams.gT)); 
%         for g = 1:numel(Eig.netParams.gT)
%             DistsPW_EV1_l(:,g) = Eig.MC.ev1_l{g}.DistAvgPW;        % Divisive Margin for 1st Eigenvector (linear computation)
%         end
%         Eig1DistsPW_L = [ Eig1DistsPW_L , DistsPW_EV1_l ];       
%     else
%         
%         Eig1DistsPW_L = [ Eig1DistsPW_L , nan( 2, numel(Eig.netParams.gT)) ];
%         
%     end
%     
% %     if isfield(Eig.MC,'ev1v_l')
% %         DistsPW_EV1v_l = zeros( 2, numel(Eig.netParams.gT)); 
% %         for g = 1:numel(Eig.netParams.gT)
% %             DistsPW_EV1v_l(:,g) = Eig.MC.ev1v_l{g}.DistAvgPW(1);     % Divisive Margin for 1st Eigenvector with nonlinear visualization (linear computation)
% %         end
% %         Eig1vDistsPW_L = [ Eig1vDistsPW_L , DistsPW_EV1v_l ];
% %         
% %     else
% %         
% %         Eig1vDistsPW_L = [ Eig1vDistsPW_L , nan( 2, numel(Eig.netParams.gT)) ];
% %         
% %     end
%     
%     %
%     %
%     %
%     %
%     %
%     % NOTE:  I am suspicious that the max value in Eig1DistsPW_C or
%     % Eig1vDistsPW_C is 1 and not pi or something at least >1.  I am not
%     % normalizing the quantities, but still they seem normalized. Is this
%     % also the case for Mod_SKHAdj?  
%     %
%     % In Mod_SKHAdj, they are >1.
%     %
%     %
%     % Maybe, I need to normalize 2D scatter space I am plotting in 
%     % visualize_pts_from_2x2D_scatterB by pi so that it goes from 0 to 1 so
%     % that I can compare it to Eigenvector scatter space and any other
%     % normalized scatter space.
%     %
%     % I am not using them in next step of processing right now, I will have
%     % to figure this out before I do.
%     % 
%     %
%     %
%     %
%     
% %     if isfield(Eig.MC,'ev1_c_pi')
% %         DistsPW_EV1_c_pi = zeros( 2, numel(Eig.netParams.gT)); 
% %         for g = 1:numel(Eig.netParams.gT)
% %             DistsPW_EV1_c_pi(:,g) = Eig.MC.ev1_c_pi{g}.DistAvgPW;        % Divisive Margin for 1st Eigenvector (circular computation)
% %         end
% %         Eig1DistsPW_C = [ Eig1DistsPW_C , DistsPW_EV1_c_pi ];
% %         
% %     else
% %         
% %         Eig1DistsPW_C = [ Eig1DistsPW_C , nan( 2, numel(Eig.netParams.gT)) ];
% %         
% %     end
% 
% %     if isfield(Eig.MC,'ev1v_c_pi')
% %         DistsPW_EV1v_c_pi = zeros( 2, numel(Eig.netParams.gT)); 
% %         for g = 1:numel(Eig.netParams.gT)
% %             DistsPW_EV1v_c_pi(:,g) = Eig.MC.ev1v_c_pi{g}.DistAvgPW;     % Divisive Margin for 1st Eigenvector with nonlinear visualization (circular computation)
% %         end
% %         Eig1vDistsPW_C = [ Eig1vDistsPW_C , DistsPW_EV1_c_pi ];
% %         
% %     else
% %         
% %         Eig1vDistsPW_C = [ Eig1vDistsPW_C , nan( 2, numel(Eig.netParams.gT)) ];
% %         
% %     end
%     %
%     %
%     %
%     % Check that linear & circular eigenvector 1 DistsPW computations match
%     if(0)
%         disp('Evec1 Viz. Pi Embedding.')
%         Evec1v_LvC_diffMn(i) = mean(abs(DistsPW_EV1v_c_pi-DistsPW_EV1v_l));
%         if ( mean(abs(DistsPW_EV1v_c_pi-DistsPW_EV1v_l)) > 0.01 )
%             disp('DistsPW results on Evec1viz different using linear vs. circular computation.')
%             [DistsPW_EV1v_c_pi;DistsPW_EV1v_l]
% %             keyboard
%         end
% 
%         disp('Evec1. Pi Embedding')
%         Evec1_LvC_diffMn(i) = mean(abs(DistsPW_EV1_c_pi-DistsPW_EV1_l));
%         if ( mean(abs(DistsPW_EV1_c_pi-DistsPW_EV1_l)) > 0.01 )
%             disp('DistsPW results on Eigenvector1 different using linear vs. circular computation.')
%            [DistsPW_EV1_c_pi;DistsPW_EV1_l]
% %            keyboard
%         end
%         
%     end
% 
% 
%     
%     
%     %
%     %% Using Only 2nd or Only 3rd Eigenvector.
%     %
%     if isfield(Eig.MC,'ev2o_l')
%         DistsPW_EV2o = zeros( 2, numel(Eig.netParams.gT)); 
%         for g = 1:numel(Eig.netParams.gT)
%             DistsPW_EV2o(:,g) = Eig.MC.ev2o_l{g}.DistAvgPW;       % Divisive Margin for 2nd Eigenvector only
%         end
%         Eig2oDistsPW = [ Eig2oDistsPW , DistsPW_EV2o ];
% 
%     else
%         
%         Eig2oDistsPW = [ Eig2oDistsPW , nan( 2, numel(Eig.netParams.gT)) ];
%         
%     end
%     
% %     if isfield(Eig.MC,'ev2ov_l')
% %         DistsPW_EV2ov = zeros( 2, numel(Eig.netParams.gT)); 
% %         for g = 1:numel(Eig.netParams.gT)
% %             DistsPW_EV2ov(:,g) = Eig.MC.ev2ov_l{g}.DistAvgPW;    % Divisive Margin for 2nd Eigenvector only with nonlinear visualization
% %         end
% %         Eig2ovDistsPW = [ Eig2ovDistsPW , DistsPW_EV2ov ];
% %         
% %     else
% %         
% %         Eig2ovDistsPW = [ Eig2ovDistsPW , nan( 2, numel(Eig.netParams.gT)) ];
% %         
% %     end
%     
%     
%     if isfield(Eig.MC,'ev3o_l')
%         DistsPW_EV3o = zeros( 2, numel(Eig.netParams.gT)); 
%         for g = 1:numel(Eig.netParams.gT)
%             DistsPW_EV3o(:,g) = Eig.MC.ev3o_l{g}.DistAvgPW;       % Divisive Margin for 3rd Eigenvector only
%         end
%         Eig3oDistsPW = [ Eig3oDistsPW , DistsPW_EV3o ];
%         
%     else
%         
%         Eig3oDistsPW = [ Eig3oDistsPW , nan( 2, numel(Eig.netParams.gT)) ];
%     
%     end
%     
% %     if isfield(Eig.MC,'ev3ov_l')
% %         DistsPW_EV3ov = zeros( 2, numel(Eig.netParams.gT)); 
% %         for g = 1:numel(Eig.netParams.gT)
% %             DistsPW_EV3ov(:,g) = Eig.MC.ev3ov_l{g}.DistAvgPW;    % Divisive Margin for 3rd Eigenvector only with nonlinear visualization
% %         end
% %         Eig3ovDistsPW = [ Eig3ovDistsPW , DistsPW_EV3ov ];
% %         
% %     else
% %         
% %         Eig3ovDistsPW = [ Eig3ovDistsPW , nan( 2, numel(Eig.netParams.gT)) ];
% %         
% %     end
%     %
%     %% Using 1st & 2nd Eigenvector together.
%     %
%     
%     if isfield(Eig.MC,'ev2_l')
%         DistsPW_EV2 = zeros( 2, numel(Eig.netParams.gT)); 
%         for g = 1:numel(Eig.netParams.gT)
%             DistsPW_EV2(:,g) = Eig.MC.ev2_l{g}.DistAvgPW;          % Divisive Margin for Eigenvectors 1&2
%         end
%         Eig2DistsPW = [ Eig2DistsPW , DistsPW_EV2 ];
%     else
%         
%         Eig2DistsPW = [ Eig2DistsPW , nan( 2, numel(Eig.netParams.gT)) ];
%         
%     end
%     
% %     if isfield(Eig.MC,'ev2v_l')
% %         DistsPW_EV2v = zeros( 2, numel(Eig.netParams.gT)); 
% %         for g = 1:numel(Eig.netParams.gT)
% %             DistsPW_EV2v(:,g) = Eig.MC.ev2v_l{g}.DistAvgPW;       % Divisive Margin for Eigenvectors 1&2 with nonlinear vizualization
% %         end
% %         Eig2vDistsPW = [ Eig2vDistsPW , DistsPW_EV2v ];
% %         
% %     else
% %         
% %         Eig2vDistsPW = [ Eig2vDistsPW , nan( 2, numel(Eig.netParams.gT)) ];
% %         
% %     end
% 
%     if isfield(Eig.MC,'ev2w_l')
%         DistsPW_EV2w = zeros( 2, numel(Eig.netParams.gT)); 
%         for g = 1:numel(Eig.netParams.gT)
%             DistsPW_EV2w(:,g) = Eig.MC.ev2w_l{g}.DistAvgPW;       % Divisive Margin for Eigenvectors 1&2 weighted by eigenvectors
%         end
%         Eig2wDistsPW = [ Eig2wDistsPW , DistsPW_EV2w ];
%         
%     else
%         
%         Eig2wDistsPW = [ Eig2wDistsPW , nan( 2, numel(Eig.netParams.gT)) ];
%         
%     end
%     
% %     if isfield(Eig.MC,'ev2wv_l')
% %         DistsPW_EV2wv = zeros( 2, numel(Eig.netParams.gT)); 
% %         for g = 1:numel(Eig.netParams.gT)
% %             DistsPW_EV2wv(:,g) = Eig.MC.ev2wv_l{g}.DistAvgPW;    % Divisive Margin for Eigenvectors 1&2 with nonlinear vizualization weighted by eigenvectors
% %         end
% %         normlz=1;
% %         Eig2wvDistsPW = [ Eig2wvDistsPW , DistsPW_EV2wv ];% ./normlz ];
% %         
% %     else
% %         
% %         Eig2wvDistsPW = [ Eig2wvDistsPW , nan( 2, numel(Eig.netParams.gT)) ];
% %         
% %     end
%     %
%     %% Using 1st thru 3rd Eigenvectors together.
%     %
%     
%     if isfield(Eig.MC,'ev3_l')
%         DistsPW_EV3 = zeros( 2, numel(Eig.netParams.gT)); 
%         for g = 1:numel(Eig.netParams.gT)
%             DistsPW_EV3(:,g) = Eig.MC.ev3_l{g}.DistAvgPW;          % Divisive Margin for Eigenvectors 1-3
%         end
%         Eig3DistsPW = [ Eig3DistsPW , DistsPW_EV3 ];
%         
%     else
%         
%         Eig3DistsPW = [ Eig3DistsPW , nan( 2, numel(Eig.netParams.gT)) ];
%         
%     end
%     
% %     if isfield(Eig.MC,'ev3v_l')
% %         DistsPW_EV3v = zeros( 2, numel(Eig.netParams.gT)); 
% %         for g = 1:numel(Eig.netParams.gT)
% %             DistsPW_EV3v(:,g) = Eig.MC.ev3v_l{g}.DistAvgPW;       % Divisive Margin for Eigenvectors 1-3 with nonlinear vizualization
% %         end
% %         Eig3vDistsPW = [ Eig3vDistsPW , DistsPW_EV3v ];
% %         
% %     else
% %         
% %         Eig3vDistsPW = [ Eig3vDistsPW , nan( 2, numel(Eig.netParams.gT)) ];
% %         
% %     end
% 
%     if isfield(Eig.MC,'ev3w_l')
%         DistsPW_EV3w = zeros( 2, numel(Eig.netParams.gT)); 
%         for g = 1:numel(Eig.netParams.gT)
%             DistsPW_EV3w(:,g) = Eig.MC.ev3w_l{g}.DistAvgPW;       % Divisive Margin for Eigenvectors 1-3 weighted by eigenvectors
%         end
%         Eig3wDistsPW = [ Eig3wDistsPW , DistsPW_EV3w ];
%         
%     else
%         
%         Eig3wDistsPW = [ Eig3wDistsPW , nan( 2, numel(Eig.netParams.gT)) ];
%         
%     end
%     
% %     if isfield(Eig.MC,'ev3wv_l')
% %         DistsPW_EV3wv = zeros( 2, numel(Eig.netParams.gT)); 
% %         for g = 1:numel(Eig.netParams.gT)
% %             DistsPW_EV3wv(:,g) = Eig.MC.ev3wv_l{g}.DistAvgPW;    % Divisive Margin for Eigenvectors 1-3 with nonlinear vizualization weighted by eigenvectors
% %         end
% %         Eig3wvDistsPW = [ Eig3wvDistsPW , DistsPW_EV3wv ];
% %         
% %     else
% %         
% %         Eig3wvDistsPW = [ Eig3wvDistsPW , nan( 2, numel(Eig.netParams.gT)) ];
% %         
% %     end
    
    
    
    
    
    % Plot some stuff...
    if(0)
        
        figure
        
        subplot(221), imagesc(ImPtch.im), colormap('bone'),title('Image')
        set(gca,'XTick',[],'YTick',[])
        xlabel(['DM_L = ',num2str(mean(DistsPW_SMa),2),' // ','DM_C = ',num2str(mean(DistsPW_SMc),2)])
        %
        KurPhaseFinl = visKurPhase_inBone(Kur.netParams.im, reshape(Kur.metaCluster.phaseAtClk(:,end),Kur.netParams.Ndims(1),Kur.netParams.Ndims(2)));        
        subplot(222), imagesc(KurPhaseFinl), colormap('bone'),title('Kur')
        set(gca,'XTick',[],'YTick',[])
        xlabel(['DM = ',num2str(mean(DistsPW_K(end,:)),2)])
        %
        subplot(223), imagesc(reshape(Eig.EVecsML(:,1),Eig.netParams.Ndims(1),Eig.netParams.Ndims(2))), colormap('bone'),title('Eig')
        set(gca,'XTick',[],'YTick',[])
        xlabel(['DM = ',num2str(mean(DistsPW_EV1_l),2)])
        %
        subplot(224), imagesc( EvecVizF(  reshape(Eig.EVecsML(:,1),Eig.netParams.Ndims(1),Eig.netParams.Ndims(2)), 1e-16 ) ), colormap('bone'),title('EigViz')
        set(gca,'XTick',[],'YTick',[])
        xlabel(['DM = ',num2str(mean(DistsPW_EV1v_l),2)])
        
        
        keyboard
        
        
    end
    
    
    


%     figure, 
%     subplot(131),imagesc(reshape(Eig.EVecsML(:,1),51,51)), colormap('bone'), colorbar
%     subplot(132),imagesc(reshape(Eig.EVecsML(:,2),51,51)), colormap('bone'), colorbar
%     subplot(133),imagesc(reshape(Eig.EVecsML(:,3),51,51)), colormap('bone'), colorbar
%     
%     keyboard

    

end % loop over files (either Kur ones or Eig ones depending on setting of scatKur & scatEig)




% keyboard


disp('How many files were there something weird with Eig numbers?')
somethingWeirdWithEig



% % Plot some Histograms of Pairwise Distances in & out of Clusters.
% % The Distance values are all normalized and so should fall between 0 & 1.
% if(0)
%     
%     % Kuramoto
%     H=histDistNorm(KurDistsPW,'Kuramoto');
% 
%     % Strawman L
%     H=histDistNorm(SMDistsPW_L,'Strawman L');
%     
%     % Strawman C
%     H=histDistNorm(SMDistsPW_C,'Strawman C');
%     
%     % Eigenvector 1 Linear
%     H=histDistNorm(Eig1DistsPW_L,'Eigenvector 1 Linear');
% 
%     % Eigenvector 1 Phase Embed
%     H=histDistNorm(Eig1DistsPW_C,'Eigenvector 1 Phase Embed');
%     
%     % Eigenvector 1 Viz Linear
%     H=histDistNorm(Eig1vDistsPW_L,'Eigenvector 1 Viz Linear');
% 
%     % Eigenvector 1 Phase Embed
%     H=histDistNorm(Eig1vDistsPW_C,'Eigenvector 1 Viz Phase Embed');
%     
%     % Eigenvector 2 Linear
%     H=histDistNorm(Eig2oDistsPW,'Eigenvector 2 Linear');
%     
%     % Eigenvector 2 Viz Linear
%     H=histDistNorm(Eig2ovDistsPW,'Eigenvector 2 Viz Linear');
%     
%     % Eigenvector 3 Linear
%     H=histDistNorm(Eig3oDistsPW,'Eigenvector 3 Linear');
%     
%     % Eigenvector 3 Viz Linear
%     H=histDistNorm(Eig3ovDistsPW,'Eigenvector 3 Viz Linear');
%     
%     % Eigenvectors 1-2 Linear
%     H=histDistNorm(Eig2DistsPW,'Eigenvectors 1-2 Linear');
%     
%     % Eigenvectors 1-2 Viz Linear
%     H=histDistNorm(Eig2vDistsPW,'Eigenvectors 1-2 Viz Linear');
%     
%     % Eigenvectors 1-3 Linear
%     H=histDistNorm(Eig3DistsPW,'Eigenvectors 1-3 Linear');
%     
%     % Eigenvectors 1-3 Viz Linear
%     H=histDistNorm(Eig3vDistsPW,'Eigenvectors 1-3 Viz Linear');
%     
%     keyboard
%     close all
% 
% end




% Eig2wDistsPW, Eig2wvDistsPW
% Eig3wDistsPW, Eig3wvDistsPW




toc



% keyboard





% % Visualize difference between Circular & Linear DistsPW Computations. Histograms
% if(0)
%     [ys,xs] = hist(Straw_LvC_diffMn,50);
%     [ye,xe] = hist(Evec1_LvC_diffMn,50);
%     [yv,xv] = hist(Evec1v_LvC_diffMn,50);
%     %
%     h=figure;
%     subplot(131), bar(xs,ys), set(gca,'Yscale','log','FontSize',16,'FontWeight','Bold'), 
%     title('Strawman')
%     xlabel('Diff btwn Circ & Linr')
%     ylabel('# Occurances / 5000')
%     subplot(132), bar(xe,ye), set(gca,'Yscale','log','FontSize',16,'FontWeight','Bold'), 
%     title('Eigenvector1')
%     xlabel('Diff btwn Circ & Linr')
%     %ylabel('# Occurances / 5000')
%     subplot(133), bar(xv,yv), set(gca,'Yscale','log','FontSize',16,'FontWeight','Bold'), 
%     title('Evec1 Viz')
%     xlabel('Diff btwn Circ & Linr')
%     %ylabel('# Occurances / 5000')
%     
%     annotation('textbox', [0 0.9 1 0.1],'String', ['# Pts= ',num2str(numFiles2Loop),', \sigmaP= ',sP,...
%             ', Rmax: ',rM,', \sigmaD= ',sD,', \sigmaW= ',sW, ', Kscale: ',Ks,', Tscale: ',Ts,' '], ...
%             'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')
%         
%    saveGoodImg(h,[dirScat,'Compare_Circ_vs_Linr_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts],sizeGoodIm)
%    close(h)
%     
% end














% % Go through misaligned patches and find their indecies in data vectors.
% misaligned_indx = []; 
% for i = 1:numel(misalignedID)
%     misaligned_indx = [misaligned_indx, find(strcmp(misalignedID{i},ImgPtchID))];
% end
% %
% good_indx = setxor(1:numel(ImgPtchID),misaligned_indx);





% MAYBE I WANT TO SAVE A SUMMARY FILE HERE AND COMPUTE SENSITIVIY IN A WHOLE OTHER STEP?




%% Now Plot Scatter Plot of Straw vs Eigen & Straw vs. Kuramoto.
disp('Computing Sensitivity Metrics')

% Note: For all RateDistSig ranges from [0 to 1] I think
RateDistSig = [0,0.01,0.1,0.3,0.5,0.7,0.9];

disp('d'' Weighted by Cluster Size')
disp('kur')
[d_prime_wtd.kur] = compute_SensitivityIndex(mnCluster_all.kur, stdCluster_all.kur, sizeCluster_all, RateDistSig, 1, 1);
disp('im')
[d_prime_wtd.im] = compute_SensitivityIndex(mnCluster_all.im, stdCluster_all.im, sizeCluster_all, RateDistSig, 0, 1);
disp('ev1')
[d_prime_wtd.ev1] = compute_SensitivityIndex(mnCluster_all.ev1, stdCluster_all.ev1, sizeCluster_all, RateDistSig, 0, 1);
disp('ev2o')
[d_prime_wtd.ev2o] = compute_SensitivityIndex(mnCluster_all.ev2o, stdCluster_all.ev2o, sizeCluster_all, RateDistSig, 0, 1);
disp('ev3o')
[d_prime_wtd.ev3o] = compute_SensitivityIndex(mnCluster_all.ev3o, stdCluster_all.ev3o, sizeCluster_all, RateDistSig, 0, 1);
disp('ev2')
[d_prime_wtd.ev2] = compute_SensitivityIndex(mnCluster_all.ev2, stdCluster_all.ev2, sizeCluster_all, RateDistSig, 0, 1);
disp('ev3')
[d_prime_wtd.ev3] = compute_SensitivityIndex(mnCluster_all.ev3, stdCluster_all.ev3, sizeCluster_all, RateDistSig, 0, 1);
disp('ev2w')
[d_prime_wtd.ev2w] = compute_SensitivityIndex(mnCluster_all.ev2w, stdCluster_all.ev2w, sizeCluster_all, RateDistSig, 0, 1);
disp('ev3w')
[d_prime_wtd.ev3w] = compute_SensitivityIndex(mnCluster_all.ev3w, stdCluster_all.ev3w, sizeCluster_all, RateDistSig, 0, 1);
%
disp('d'' Unweighted Clusters')
disp('kur')
[d_prime.kur] = compute_SensitivityIndex(mnCluster_all.kur, stdCluster_all.kur, sizeCluster_all, RateDistSig, 1, 0);
disp('im')
[d_prime.im] = compute_SensitivityIndex(mnCluster_all.im, stdCluster_all.im, sizeCluster_all, RateDistSig, 0, 0);
disp('ev1')
[d_prime.ev1] = compute_SensitivityIndex(mnCluster_all.ev1, stdCluster_all.ev1, sizeCluster_all, RateDistSig, 0, 0);
disp('ev2o')
[d_prime.ev2o] = compute_SensitivityIndex(mnCluster_all.ev2o, stdCluster_all.ev2o, sizeCluster_all, RateDistSig, 0, 0);
disp('ev3o')
[d_prime.ev3o] = compute_SensitivityIndex(mnCluster_all.ev3o, stdCluster_all.ev3o, sizeCluster_all, RateDistSig, 0, 0);
disp('ev2')
[d_prime.ev2] = compute_SensitivityIndex(mnCluster_all.ev2, stdCluster_all.ev2, sizeCluster_all, RateDistSig, 0, 0);
disp('ev3')
[d_prime.ev3] = compute_SensitivityIndex(mnCluster_all.ev3, stdCluster_all.ev3, sizeCluster_all, RateDistSig, 0, 0);
disp('ev2w')
[d_prime.ev2w] = compute_SensitivityIndex(mnCluster_all.ev2w, stdCluster_all.ev2w, sizeCluster_all, RateDistSig, 0, 0);
disp('ev3w')
[d_prime.ev3w] = compute_SensitivityIndex(mnCluster_all.ev3w, stdCluster_all.ev3w, sizeCluster_all, RateDistSig, 0, 0);






% % (1). Kuramoto vs Strawman. (note I am using Linear Strawman because DivMarg & Pmetric can handle it)
% goode1 = find(~isnan(KurDistsPW(1,good_indx))); goode2 = find(~isnan(SMDistsPW_L(1,good_indx))); goode = intersect(goode1,goode2);
% if ~isempty(goode)
%     [H, Pmtric.Kur, DivMarg.Kur] =  scatter_plot_DistPairs(KurDistsPW(:,goode),SMDistsPW_L(:,goode),pi,dirKur,dirImg,fileGeneral,fileSize,'Kuramoto',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID(goode),ImgGtID(goode),GtNumSegs(goode),plot2x2Dscatter);
%     if(H~=0)
%         saveGoodImg(H,[dirScat,'DistsPW_Scatter_2D_Kuramoto_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'_files',num2str(numFiles2Loop)],sizeGoodIm)
%         close(H)
%     end
% end
% 
% % (1b). Kuramoto using empirical range vs Strawman. (note I am using Linear Strawman because DivMarg & Pmetric can handle it)
% goode1 = find(~isnan(KurDistsPW(1,good_indx))); goode2 = find(~isnan(SMDistsPW_L(1,good_indx))); goode = intersect(goode1,goode2);
% if ~isempty(goode)
%     [H, Pmtric.KurEmp, DivMarg.KurEmp] =  scatter_plot_DistPairs(KurDistsPW(:,goode),SMDistsPW_L(:,goode),EmpRange.Kur(goode),dirKur,dirImg,fileGeneral,fileSize,'Kuramoto',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID(goode),ImgGtID(goode),GtNumSegs(goode),plot2x2Dscatter);
%     if(H~=0)
%         saveGoodImg(H,[dirScat,'DistsPW_Scatter_2D_Kuramoto_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'_files',num2str(numFiles2Loop)],sizeGoodIm)
%         close(H)
%     end
% end
% 
% % (2). Eigenvector 1 vs Strawman.
% goode1 = find(~isnan(Eig1DistsPW_L(1,good_indx))); goode2 = find(~isnan(SMDistsPW_L(1,good_indx))); goode = intersect(goode1,goode2);
% if ~isempty(goode)
%     [H, Pmtric.Eig1, DivMarg.Eig1] =  scatter_plot_DistPairs(Eig1DistsPW_L(:,goode),SMDistsPW_L(:,goode),EmpRange.Eig1(goode),dirEig,dirImg,fileGeneral,fileSize,'EigVec1',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID(goode),ImgGtID(goode),GtNumSegs(goode),plot2x2Dscatter);
%     if(H~=0)
%         saveGoodImg(H,[dirScat,'DistsPW_Scatter_2D_EigVec1_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'_files',num2str(numFiles2Loop)],sizeGoodIm)
%         close(H)
%     end
% end
% 
% % (3). Eigenvector 1 viz. vs Strawman.
% goode1 = find(~isnan(Eig1vDistsPW_L(1,good_indx))); goode2 = find(~isnan(SMDistsPW_L(1,good_indx))); goode = intersect(goode1,goode2);
% if ~isempty(goode)
%     [H, Pmtric.Eig1v, DivMarg.Eig1v] =  scatter_plot_DistPairs(Eig1vDistsPW_L(:,goode),SMDistsPW_L(:,goode),EmpRange.Eig1v(goode),dirEig,dirImg,fileGeneral,fileSize,'EigVec1 Viz',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID(goode),ImgGtID(goode),GtNumSegs(goode),plot2x2Dscatter);
%     if(H~=0)
%         saveGoodImg(H,[dirScat,'DistsPW_Scatter_2D_EigVec1v_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'_files',num2str(numFiles2Loop)],sizeGoodIm)
%         close(H)
%     end
% end
% 
% % (4). Eigenvector 2 vs Strawman.
% goode1 = find(~isnan(Eig2oDistsPW(1,good_indx))); goode2 = find(~isnan(SMDistsPW_L(1,good_indx))); goode = intersect(goode1,goode2);
% if ~isempty(goode)
%     [H, Pmtric.Eig2o, DivMarg.Eig2o] =  scatter_plot_DistPairs(Eig2oDistsPW(:,goode),SMDistsPW_L(:,goode),EmpRange.Eig2o(goode),dirEig,dirImg,fileGeneral,fileSize,'EigVec2',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID(goode),ImgGtID(goode),GtNumSegs(goode),plot2x2Dscatter);
%     if(H~=0)
%         saveGoodImg(H,[dirScat,'DistsPW_Scatter_2D_EigVec2_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'_files',num2str(numFiles2Loop)],sizeGoodIm)
%         close(H)
%     end
% end
% 
% % (5). Eigenvector 2 viz. vs Strawman.
% goode1 = find(~isnan(Eig2ovDistsPW(1,good_indx))); goode2 = find(~isnan(SMDistsPW_L(1,good_indx))); goode = intersect(goode1,goode2);
% if ~isempty(goode)
%     [H, Pmtric.Eig2ov, DivMarg.Eig2ov] =  scatter_plot_DistPairs(Eig2ovDistsPW(:,goode),SMDistsPW_L(:,goode),EmpRange.Eig2ov(goode),dirEig,dirImg,fileGeneral,fileSize,'EigVec2 Viz',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID(goode),ImgGtID(goode),GtNumSegs(goode),plot2x2Dscatter);
%     if(H~=0)
%         saveGoodImg(H,[dirScat,'DistsPW_Scatter_2D_EigVec2v_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'_files',num2str(numFiles2Loop)],sizeGoodIm)
%         close(H)
%     end
% end
% 
% % (6). Eigenvector 3 vs Strawman.
% goode1 = find(~isnan(Eig3oDistsPW(1,good_indx))); goode2 = find(~isnan(SMDistsPW_L(1,good_indx))); goode = intersect(goode1,goode2);
% if ~isempty(goode)
%     [H, Pmtric.Eig3o, DivMarg.Eig3o] =  scatter_plot_DistPairs(Eig3oDistsPW(:,goode),SMDistsPW_L(:,goode),EmpRange.Eig3o(goode),dirEig,dirImg,fileGeneral,fileSize,'EigVec3',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID(goode),ImgGtID(goode),GtNumSegs(goode),plot2x2Dscatter);
%     if(H~=0)
%         saveGoodImg(H,[dirScat,'DistsPW_Scatter_2D_EigVec3_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'_files',num2str(numFiles2Loop)],sizeGoodIm)
%         close(H)
%     end
% end
% 
% % (7). Eigenvector 3 viz. vs Strawman.
% goode1 = find(~isnan(Eig3ovDistsPW(1,good_indx))); goode2 = find(~isnan(SMDistsPW_L(1,good_indx))); goode = intersect(goode1,goode2);
% if ~isempty(goode)
%     [H, Pmtric.Eig3ov, DivMarg.Eig3ov] =  scatter_plot_DistPairs(Eig3ovDistsPW(:,goode),SMDistsPW_L(:,goode),EmpRange.Eig3ov(goode),dirEig,dirImg,fileGeneral,fileSize,'EigVec3 Viz',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID(goode),ImgGtID(goode),GtNumSegs(goode),plot2x2Dscatter);
%     if(H~=0)
%         saveGoodImg(H,[dirScat,'DistsPW_Scatter_2D_EigVec3v_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'_files',num2str(numFiles2Loop)],sizeGoodIm)
%         close(H)
%     end
% end
% 
% % (8). Eigenvectors 1-2 (unweighted, w/o nonlinear visualization)
% goode1 = find(~isnan(Eig2DistsPW(1,good_indx))); goode2 = find(~isnan(SMDistsPW_L(1,good_indx))); goode = intersect(goode1,goode2);
% if ~isempty(goode)
%     [H, Pmtric.Eig2, DivMarg.Eig2] =  scatter_plot_DistPairs(Eig2DistsPW(:,goode),SMDistsPW_L(:,goode),EmpRange.Eig2(goode),dirEig,dirImg,fileGeneral,fileSize,'EigVecs 1-2',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID(goode),ImgGtID(goode),GtNumSegs(goode),plot2x2Dscatter);
%     if(H~=0)
%         saveGoodImg(H,[dirScat,'DistsPW_Scatter_2D_EigVec1-2_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'_files',num2str(numFiles2Loop)],sizeGoodIm)
%         close(H)
%     end
% end
% 
% % (9). Eigenvectors 1-2 (unweighted, with nonlinear visualization)
% goode1 = find(~isnan(Eig2vDistsPW(1,good_indx))); goode2 = find(~isnan(SMDistsPW_L(1,good_indx))); goode = intersect(goode1,goode2);
% if ~isempty(goode)
%     [H, Pmtric.Eig2v, DivMarg.Eig2v] =  scatter_plot_DistPairs(Eig2vDistsPW(:,goode),SMDistsPW_L(:,goode),EmpRange.Eig2v(goode),dirEig,dirImg,fileGeneral,fileSize,'EigVecs 1-2 Viz',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID(goode),ImgGtID(goode),GtNumSegs(goode),plot2x2Dscatter);
%     if(H~=0)
%         saveGoodImg(H,[dirScat,'DistsPW_Scatter_2D_EigVec1-2v_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'_files',num2str(numFiles2Loop)],sizeGoodIm)
%         close(H)
%     end
% end
% 
% % % (10). Eigenvectors 1-2 (weighted by their eigenvalues, w/o nonlinear visualization)
% % goode1 = find(~isnan(Eig2wDistsPW(1,good_indx))); goode2 = find(~isnan(SMDistsPW_L(1,good_indx))); goode = intersect(goode1,goode2);
% % if ~isempty(goode)
% %     [H, Pmtric.Eig2w, DivMarg.Eig2w] =  scatter_plot_DistPairs(Eig2wDistsPW(:,goode),SMDistsPW_L(:,goode),dirEig,dirImg,fileGeneral,fileSize,'EigVecs 1-2 Weighted',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID(goode),ImgGtID(goode),GtNumSegs(goode),plot2x2Dscatter);
% %     if(H~=0)
% %         saveGoodImg(H,[dirScat,'DistsPW_Scatter_2D_EigVec1-2w_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'_files',num2str(numFiles2Loop)],sizeGoodIm)
% %         close(H)
% %     end
% % end
% 
% % % (11). Eigenvectors 1-2 (weighted by their eigenvalues, with nonlinear visualization)
% % goode1 = find(~isnan(Eig2wvDistsPW(1,good_indx))); goode2 = find(~isnan(SMDistsPW_L(1,good_indx))); goode = intersect(goode1,goode2);
% % if ~isempty(goode)
% %     [H, Pmtric.Eig2wv, DivMarg.Eig2wv] =  scatter_plot_DistPairs(Eig2wvDistsPW(:,goode),SMDistsPW_L(:,goode),dirEig,dirImg,fileGeneral,fileSize,'EigVecs 1-2 Weighted Viz',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID(goode),ImgGtID(goode),GtNumSegs(goode),plot2x2Dscatter);
% %     if(H~=0)
% %         saveGoodImg(H,[dirScat,'DistsPW_Scatter_2D_EigVec1-2wv_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'_files',num2str(numFiles2Loop)],sizeGoodIm)
% %         close(H)
% %     end
% % end
% 
% % (12). Eigenvectors 1-3 (unweighted, w/o nonlinear visualization)
% goode1 = find(~isnan(Eig3DistsPW(1,good_indx))); goode2 = find(~isnan(SMDistsPW_L(1,good_indx))); goode = intersect(goode1,goode2);
% if ~isempty(goode)
%     [H, Pmtric.Eig3, DivMarg.Eig3] =  scatter_plot_DistPairs(Eig3DistsPW(:,goode),SMDistsPW_L(:,goode),EmpRange.Eig3(goode),dirEig,dirImg,fileGeneral,fileSize,'EigVecs 1-3',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID(goode),ImgGtID(goode),GtNumSegs(goode),plot2x2Dscatter);
%     if(H~=0)
%         saveGoodImg(H,[dirScat,'DistsPW_Scatter_2D_EigVec1-3_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'_files',num2str(numFiles2Loop)],sizeGoodIm)
%         close(H)
%     end
% end
% 
% % (13). Eigenvectors 1-3 (unweighted, with nonlinear visualization)
% goode1 = find(~isnan(Eig3vDistsPW(1,good_indx))); goode2 = find(~isnan(SMDistsPW_L(1,good_indx))); goode = intersect(goode1,goode2);
% if ~isempty(goode)
%     [H, Pmtric.Eig3v, DivMarg.Eig3v] =  scatter_plot_DistPairs(Eig3vDistsPW(:,goode),SMDistsPW_L(:,goode),EmpRange.Eig3v(goode),dirEig,dirImg,fileGeneral,fileSize,'EigVecs 1-3 Viz',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID(goode),ImgGtID(goode),GtNumSegs(goode),plot2x2Dscatter);
%     if(H~=0)
%         saveGoodImg(H,[dirScat,'DistsPW_Scatter_2D_EigVec1-3v_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'_files',num2str(numFiles2Loop)],sizeGoodIm)
%         close(H)
%     end
% end
% 
% % % (14). Eigenvectors 1-3 (weighted by their eigenvalues, w/o nonlinear visualization)
% % goode1 = find(~isnan(Eig3wDistsPW(1,good_indx))); goode2 = find(~isnan(SMDistsPW_L(1,good_indx))); goode = intersect(goode1,goode2);
% % if ~isempty(goode)
% %     [H, Pmtric.Eig3w, DivMarg.Eig3w] =  scatter_plot_DistPairs(Eig3wDistsPW(:,goode),SMDistsPW_L(:,goode),dirEig,dirImg,fileGeneral,fileSize,'EigVecs 1-3 Weighted',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID(goode),ImgGtID(goode),GtNumSegs(goode),plot2x2Dscatter);
% %     if(H~=0)
% %         saveGoodImg(H,[dirScat,'DistsPW_Scatter_2D_EigVec1-3w_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'_files',num2str(numFiles2Loop)],sizeGoodIm)
% %         close(H)
% %     end
% % end
% 
% % % (15). Eigenvectors 1-3 (weighted by their eigenvalues, with nonlinear visualization)
% % goode1 = find(~isnan(Eig3wvDistsPW(1,good_indx))); goode2 = find(~isnan(SMDistsPW_L(1,good_indx))); goode = intersect(goode1,goode2);
% % if ~isempty(goode)
% %     [H, Pmtric.Eig3wv, DivMarg.Eig3wv] =  scatter_plot_DistPairs(Eig3wvDistsPW(:,goode),SMDistsPW_L(:,goode),dirEig,dirImg,fileGeneral,fileSize,'EigVecs 1-3 Weighted Viz',methodType,rM,sD,sP,sW,Ts,Ks,ImgPtchID(goode),ImgGtID(goode),GtNumSegs(goode),plot2x2Dscatter);
% %     if(H~=0)
% %         saveGoodImg(H,[dirScat,'DistsPW_Scatter_2D_EigVec1-3wv_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'_files',num2str(numFiles2Loop)],sizeGoodIm)
% %         close(H)
% %     end
% % end
 
    
    


% Save a mat-file with important information - will be used in visualize_pts_from_2x2D_scatter function.
if(1)
    
    disp('Saving Mat File.')
    
    save([dirDMmat,'Clustering_Results_allPatches_sP',sP,'_rM',rM,'_sD',sD,'_sW',sW,'_Ks',Ks,'_Ts',Ts,'_files',num2str(numFiles2Loop)],...
        'ImgPtchID','ImgGtID','GtNumSegs','sP','rM','sD','sW','Ks','Ts',...
        'd_prime','d_prime_wtd','RateDistSig','mnCluster_all','stdCluster_all','mnClusterCI_all','stdClusterCI_all','idCluster_all','sizeCluster_all')


        % 'misaligned_indx','misalignedID','good_indx',
        % 'KurDistsPW','SMDistsPW_L','SMDistsPW_C','Eig1DistsPW_L','Eig2oDistsPW','Eig3oDistsPW','Pmtric','DivMarg','EmpRange',...
        % 'Eig1vDistsPW_L','Eig2ovDistsPW','Eig3ovDistsPW','Eig2DistsPW','Eig3DistsPW','Eig2vDistsPW','Eig3vDistsPW',...
        
        % 
        % 'Eig1DistsPW_C', 'Eig1vDistsPW_C','Eig2wvDistsPW','Eig3wvDistsPW','Eig2wDistsPW','Eig3wDistsPW',
end

disp('This Function has successfully completed.')
clock
