function [rangeEvecs, H] = find_range_of_variable_Loop(fileGeneral,fileSize,methodType,rM,sD,sP,sW,Ts,Ks)


% This function will loop through all the 1500 or so Evecs mat file and
% compute the range of the 1st 6 eigenvectors.  It will then plot and/or
% return that information as a histogram or array.  This is an
% investigation I am doing when trying to understand our P-metric and how
% we should normalize it by the range when using Eigenvectors.  This
% basically shows the "empirical range".

plotRangeHistograms = 1; 

                     
                     
% Directories to Input Images or Data :: Change these paths below to get a different images or patches
[dirPre,sizeGoodIm] = onCluster;
dirEig = [dirPre,'output/Kuramoto/NetsFromImgs/',fileGeneral,'_',fileSize,'/data/spectral/',methodType,'/'];
% dirKur = [dirPre,'output/Kuramoto/NetsFromImgs/',fileGeneral,'_',fileSize,'/data/Kur_PIF_Fourier1/',methodType,'/'];


% Specify these with strings or with ''.
rM = num2str(rM); % '1';
sD = num2str(sD); % 'Inf';
sP = num2str(sP); % '0p2';
sP(sP=='.')='p';
Ts = num2str(Ts); % '1';
Ks = num2str(Ks); % '300';
sW = num2str(sW); % '0';

filesEig = dir([dirEig,'Evecs','*rM',rM,'_sD',sD,'_sP',sP,'*.mat']); % will be in a different directory later.

disp('Time to loop through all metaCluster mat files:')
tic

% Decide which files to loop through
numFiles2Loop = numel(filesEig);

% numFiles2Loop = 40; % get rid of this.

rangeEvecs = zeros(numFiles2Loop,6);


for i = 1:numFiles2Loop % I dont know if this will work if I am not looping thru filesKur.

    disp([num2str(i),' / ',num2str(numFiles2Loop)])
    
    try
        load([dirEig,filesEig(i).name]);
        rangeEvecs(i,:) = range(EVecsML);
    catch
        disp('File corrupt or missing.  Skipping:')
        [dirEig,filesEig(i).name]
    end
    
end
    

toc


numBars = 100;
if(plotRangeHistograms)

    H=figure;
    
    [y,x] = hist(rangeEvecs(:,1),numBars);
    subplot(311),bar(x,y),hold on,xlim([0 2]),xlabel('Eigenvector 1 Range','FontSize',18,'FontWeight','Bold')
    plot([2/sqrt(netParams.N) 2/sqrt(netParams.N)], [0 max(y)],'r--')
    text(2/sqrt(netParams.N), 1.1*max(y),'\color{red}{2/\sqrt(N)}')
    set(gca,'FontSize',16,'FontWeight','Bold')
    title([methodType,' : ',fileGeneral,' : ',fileSize,' : rM',rM,' sD',sD,' sP',sP],'FontSize',20,'FontWeight','Bold')
    %
    [y,x] = hist(rangeEvecs(:,2),numBars);
    subplot(312),bar(x,y),hold on,xlim([0 2]),xlabel('Eigenvector 2 Range','FontSize',18,'FontWeight','Bold')
    plot([2/sqrt(netParams.N) 2/sqrt(netParams.N)], [0 max(y)],'r--')
    text(2/sqrt(netParams.N), 1.1*max(y),'\color{red}{2/\sqrt(N)}')
    set(gca,'FontSize',16,'FontWeight','Bold')
    %
    [y,x] = hist(rangeEvecs(:,3),numBars);
    subplot(313),bar(x,y),hold on,xlim([0 2]),xlabel('Eigenvector 2 Range','FontSize',18,'FontWeight','Bold')
    plot([2/sqrt(netParams.N) 2/sqrt(netParams.N)], [0 max(y)],'r--')
    text(2/sqrt(netParams.N), 1.1*max(y),'\color{red}{2/\sqrt(N)}')
    set(gca,'FontSize',16,'FontWeight','Bold')
    
    
%     2/sqrt(netParams.N)
%     2

else
    
    H=0;
    
end

keyboard