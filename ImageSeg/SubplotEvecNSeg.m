% Script to make 3 figures:  Each one has subplots for each method and each
% values of sigma for pixel distribution.
%
%      (1).  Eigenvector used for segmentation
%      (2).  Segmentation when thresholding Evector at Mean
%      (3).  Segmentation when thresholding Evector at Optimal values

screen = get(0,'ScreenSize');
[dirPre] = onCluster;

InputImage = ['GaussianBox'];  

%% Load Input file (Can run for segmentation using threshold at mean pixel value or at optimal value).
switch SegAtIn
    case 'MeanPix'
        load([dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/data_from_segmentation/SegAnalInfo_R4_S10_D13_M11_segAtmean'])
    case 'Optimal'
        load([dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/data_from_segmentation/SegAnalInfo_R4_S10_D13_M11_segAtopt'])
end

%%
if ~exist([dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/pics_from_segmentation/SubplotEvecNSeg'],'dir')
    mkdir([dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/pics_from_segmentation/SubplotEvecNSeg']);
end

% %%
% line_color = {'y','y','g','g','c','r','b','m','c','r','k'};
% line_mrkr = {'o','x','+','*','s','d','^','<','>','v','h'};
% line_style = {'-','--','-','--','-','-','--','-.','-','--','-'};

%%
for f = 1:numel(imdimL)
        imdim = imdimL(f);
        fileIn = [num2str(imdim),'x',num2str(imdim),'_gauss_center'];
        ximg = imdim;
        yimg = imdim;

    for h = 1:numel(RmaxL) % Radius of extent of neighborhood for adjacency matrix calculation
        Rmax = RmaxL(h);

        hev = figure; % Figure to plot all the eigenvectors for different methods and input images.
        hse = figure; % Figure to plot segmentations for different methods and input images.
        
        disp(['ImDim = ',num2str(imdim),' & Rmax ',num2str(Rmax)])

        for i = 1:numel(Method)

            % Loop thru segmentation output MAT files in a directory...
            files = dir([dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/data_from_segmentation/',fileIn,'/*', ...
                num2str(ximg),'x',num2str(yimg),'*',Method{i},'-*','rM',num2str(Rmax),'*.mat']);

            % use 2nd Eigenvector if Graph Laplacian or Normalized Avg Association, 1st otherwise.
            if( ~isempty(strfind(Method{i},'GL')) || ~isempty(strfind(Method{i},'AAnrm')) )
                B = 2;
            elseif( ~isempty(strfind(Method{i},'N&G')) )
                B = 3; % Not sure why 3rd eigenvector Looks Best for nonTopo (N&G) Modu.
            else
                B = 1; % Dominant Eigenvector.
            end

            for j = 1:numel(files) % loop over files in directory
                fname = files(Id(j)).name;

                load([dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/data_from_segmentation/',fileIn,'/',fname]);

                switch segAt
                    case 'opt'
                        seg = segOpt; % to use segmentation that maximizes information (upper bound)
                    case 'mean'
                        seg = segMean; % to use segmentation thresholded at mean pixel value
                end

                % Plot Eigenvectors & Segmentation for all methods & input images
                if(i==1)
                    % Plot input Image
                    figure(hev), subplot(numel(Method)+2,numel(files),j)
                    imagesc(imEnsInfo.imEns(:,:,1)), title(['\sigma=',num2str(imEnsInfo.sig1)])
                    set(gca,'xtick',[],'ytick',[]);
                    if(j==1)
                        ylabel('img in')
                    end
                    %
                    figure(hse), subplot(numel(Method)+2,numel(files),j)
                    imagesc(imEnsInfo.imEns(:,:,1)), title(['\sigma=',num2str(imEnsInfo.sig1)])
                    set(gca,'xtick',[],'ytick',[]);
                    if(j==1)
                        ylabel('img in')
                    end
                    %
                    % Plot Segmentation at mean image pixel value
                    segAtMn = imEnsInfo.imEns(:,:,1)>mean(mean(imEnsInfo.imEns(:,:,1)));
                    figure(hev), subplot(numel(Method)+2,numel(files),numel(files)+j)
                    imagesc(segAtMn)
                    set(gca,'xtick',[],'ytick',[]);
                    if(j==1)
                        ylabel('Seg @ mn pix')
                    end
                    %
                    figure(hse), subplot(numel(Method)+2,numel(files),numel(files)+j)
                    imagesc(segAtMn)
                    set(gca,'xtick',[],'ytick',[]);
                    if(j==1)
                        ylabel('Seg @ mn pix')
                    end
                end
                figure(hev), subplot(numel(Method)+2,numel(files),numel(files)*(i+1)+j)
                imagesc(EvecsML{B}(:,:,1)), %xlabel(['Evec ',num2str(B),' (\lambda=',num2str(EvalsML(B)),')']) % colorbar,
                set(gca,'xtick',[],'ytick',[]);
                if(j==1)
                    ylabel([Method(i),' Ev',num2str(B)])
                end
                %
                figure(hse), subplot(numel(Method)+2,numel(files),numel(files)*(i+1)+j)
                imagesc(seg{B}(:,:,1)), %xlabel(['SegmEvec ',num2str(B),' (\lambda=',num2str(EvalsML(B)),')']) % colorbar,
                set(gca,'xtick',[],'ytick',[]);
                if(j==1)
                    ylabel([Method(i),' Ev',num2str(B)])
                end
                
            end % loop over files in the segmentation output directory
            
        end % loop over segmentation methods {AA, GL, Mod, ...}

        % SAVE HEV AND HSE IMAGES
        fileHev = [dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/pics_from_segmentation/SubplotEvecNSeg/',fileIn,'_Evec_rM',num2str(Rmax)];   
        set(hev,'Position',screen);
        hgexport(hev, [fileHev,'.jpg'], hgexport('factorystyle'), 'Format', 'jpeg');
        fileHse = [dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/pics_from_segmentation/SubplotEvecNSeg/',fileIn,'_Seg_',segAt,'_rM',num2str(Rmax)];   
        set(hse,'Position',screen);
        hgexport(hse, [fileHse,'.jpg'], hgexport('factorystyle'), 'Format', 'jpeg');
        
        close all

    end % Loop over Rmax = [1,2,3,4]

end % Loop over ImDim = [7,9,11,13,15,17,19,21, ... 31]