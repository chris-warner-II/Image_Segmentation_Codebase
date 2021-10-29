function SegMethodAnalyze_ProbDistBox(dirIn)

% This script (written specific for GaussianBox images for now) is going to
% loop through all MAT files containing output of LoopImgSegMethodsD & 
% SegmentMethod functions. It will take files using one method and
% calculate correctness of segmentation (still figuring out method) and
% plot it vs. KL-Divergence of the two Gaussian Distributions from which
% pixel values for inside and outside rectangles were drawn.

screen = get(0,'ScreenSize');

imdimL = [7 9 11 13 15 21]; %  31
RmaxL = [1 2 3 4];
segAtL = {'mean','opt'}; % ,'opt'
Method = {'AA','AAnrm','GL','GLnrm','N&G','SKH Adj','SKH Euc','SKH D&O','N&G ME','SKH ME'}; % 'THmean',   ,'max Info'
line_color = {'y','y','g','g','c','r','b','m','c','r','k'};
line_mrkr = {'o','x','+','*','s','d','^','<','>','v','h'};
line_style = {'-','--','-','--','-','-','--','-.','-','--','-'};

%% Loop through Methods and Different Images to Calculate Information Score
for f = 1:numel(imdimL)
    
    imdim = imdimL(f);
    fileIn = [num2str(imdim),'x',num2str(imdim),'_gauss_center'];
    dirIn = ['GaussianBox/',fileIn];              
    ximg = imdim;
    yimg = imdim;
    
    for g = 1:numel(segAtL) % Segmenting at mean or optimal threshold (or other). 
        segAt = segAtL{g};

        for h = 1:numel(RmaxL) % Radius of extent of neighborhood for adjacency matrix calculation
            Rmax = RmaxL(h);

            for i = 1:numel(Method)

                % Loop thru segmentation output MAT files in a directory...
                files = dir(['./output/ImgSeg/simpleExamples/',dirIn,'/data from segmentation/*',num2str(ximg),'x',num2str(yimg),'*',Method{i},'-*','rM',num2str(Rmax(h)),'*.mat']);
                disp(['Method: ',Method{i}])

                % use 2nd Eigenvector if Graph Laplacian or Normalized Avg Association, 1st otherwise.
                if( ~isempty(strfind(Method{i},'GL')) || ~isempty(strfind(Method{i},'AAnrm')) )
                    B = 2;
                elseif( ~isempty(strfind(Method{i},'N&G')) )
                    B = 3; % Not sure why 3rd eigenvector Looks Best.
                else
                    B = 1; % Dominant Eigenvector for Modularity Methods.
                end


                for j = 1:numel(files) % loop over files in directory
                    fname = files(j).name;
                    disp(['Image # ',num2str(j),'/',num2str(numel(files))])

                    load(['./output/ImgSeg/simpleExamples/',dirIn,'/data from segmentation/',fname]);

                    switch segAt
                        case 'opt'
                            seg = segOpt; % to use segmentation that maximizes information (upper bound)
                        case 'mean'
                            seg = segMean; % to use segmentation thresholded at mean pixel value
                    end

                    Io = zeros(1,size(seg{1},3));

                    for k = 1:size(seg{1},3)
                        % NOTE: Should I normalize by max info available in image or patch?
                        %       Maybe possible for binary ground truth (with just 2 segments)
                        Io(k) = calc_infoB(seg{B}(:,:,k), imEnsInfo.gndTruth);
                    end

                    % Pack into Matrices the important information for this figure
                    infoM(i,j) = mean(Io);
                    infoS(i,j) = std(Io);

                    if(j==1 && 0)
                        % Plot first 3 Eigenvectors for the near binary box in box case.
                        figure
                        subplot(numel(EvecsML)+1,1,1), imagesc(imEnsInfo.imEns(:,:,1)), colorbar, title(fname)
                        for n = 1:numel(EvecsML)
                            subplot(numel(EvecsML)+1,1,n+1), imagesc(EvecsML{n}(:,:,1)), colorbar, xlabel(['Evec ',num2str(n),' (\lambda=',num2str(EvalsML(n)),')'])
                        end

                        keyboard

                    end


                    if (i==numel(Method))
                        imgIns(:,:,j) = imEnsInfo.imEns(:,:,1);
                        D_KL(j) = imEnsInfo.D_KL;
                        mu1(j) = imEnsInfo.mu1;
                        mu2(j) = imEnsInfo.mu2;
                        sig1(j) = imEnsInfo.sig1;
                        sig2(j) = imEnsInfo.sig2;

                        % Segmentation by Thresholding at Mean Pixel Value
                        for m = 1:size(imEnsInfo.imEns,3)
                            segTH = (imEnsInfo.imEns(:,:,m) > mean(mean(imEnsInfo.imEns(:,:,m))));
                            IoTH(m) = calc_infoB(segTH, imEnsInfo.gndTruth);
                            infoM(i+1,j) = mean(IoTH);
                            infoS(i+1,j) = std(IoTH);
                        end 
                    end

                end % loop over files from given segmentation method in the segmentation output directory

%                 keyboard

            end % loop over methods


            %% Plotting Performance of Different Methods vs. increasing Distribution Overlap
            U = infoM + infoS;
            L = infoM - infoS;
            [Y,Id] = sort(D_KL,'descend');
            [Y,Ia] = sort(D_KL,'ascend');

            % Plot Performance of Different Methods.  Info vs. D_KL.
            if(0)
                figure, plot( sig1(Ia), D_KL(Ia), 'LineWidth',2 )
                title('\sigma vs. D_{KL}','FontWeight','Bold','FontSize',20)
                xlabel('\sigma of Gaussian prob Dists','FontWeight','Bold','FontSize',18)
                ylabel('KL-Divergence between 2 Dists','FontWeight','Bold','FontSize',18)
            end
            %
            % Plot information vs. Sigma (standard deviation) for each method.
            infoMax = infoM(end,Ia(end)); % This is segmenting at mean pixel value for well separated distributions
            his=figure; hold on
            axis([ 0 0.4 0 1.3*infoMax])
            for bb = 1:numel(Method)+1
                errorbar( sig1(Ia),infoM(bb,Ia), infoS(bb,Ia), 'LineWidth',2, 'LineStyle',line_style{bb}, 'Color',line_color{bb}, 'Marker',line_mrkr{bb}, 'MarkerSize',10 ) % 
            end
            plot(sig1(Ia), repmat(infoMax(1,1),1,numel(files)) ,'k--');
            legend([Method,'MeanTH'],'Location','SouthWest','FontSize',12,'FontWeight','Bold') %
            title(['Info vs. \sigma',' with  r_{max}=',num2str(Rmax(h)),' & TH ',segAt],'FontWeight','Bold','FontSize',20)
            xlabel('\sigma of Gaussian prob Dists','FontWeight','Bold','FontSize',18)
            ylabel('information','FontWeight','Bold','FontSize',18)

            % Plot Inset of Input Images
            axes('Position',[0.85  1.15*infoMax .04 .08]);
            box on
            imagesc(imgIns(:,:,Ia(1)))
            set(gca,'xtick',[],'ytick',[])
            %
            axes('Position',[0.68  1.15*infoMax .04 .08]);
            box on
            imagesc(imgIns(:,:,Ia(3)))
            set(gca,'xtick',[],'ytick',[])
            %
            axes('Position',[0.5 1.15*infoMax  .04 .08]);
            box on
            imagesc(imgIns(:,:,Ia(5)))
            set(gca,'xtick',[],'ytick',[])
            %
            axes('Position',[0.32  1.15*infoMax .04 .08]);
            box on
            imagesc(imgIns(:,:,Ia(7)))
            set(gca,'xtick',[],'ytick',[])
            %
            axes('Position',[0.15  1.15*infoMax .04 .08]);
            box on
            imagesc(imgIns(:,:,Ia(9)))
            set(gca,'xtick',[],'ytick',[])
            
            % SAVE THIS I VS SIG FIGURE (HIS)
            fileHis = ['./output/ImgSeg/simpleExamples/',dirIn,'/pics from segmentation/',fileIn,'_IvsSig_segAt_',segAt,'_rM',num2str(Rmax(h))];   
            set(his,'Position',screen);
            hgexport(his, [fileHis,'.jpg'], hgexport('factorystyle'), 'Format', 'jpeg');

            %% Plot Single Eigenvector used for each Input Box & each Method
            hev = figure; % Figure to plot all the eigenvectors for different methods and input images.
            hse = figure; % Figure to plot segmentations for different methods and input images.

            for i = 1:numel(Method)

                % Loop thru segmentation output MAT files in a directory...
                files = dir(['./output/ImgSeg/simpleExamples/',dirIn,'/data from segmentation/*',num2str(ximg),'x',num2str(yimg),'*',Method{i},'-*','rM',num2str(Rmax(h)),'*.mat']);
                disp(['Method: ',Method{i}])

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
                    disp(['Image # ',num2str(j),'/',num2str(numel(files))])

                    load(['./output/ImgSeg/simpleExamples/',dirIn,'/data from segmentation/',fname]);

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
            fileHev = ['./output/ImgSeg/simpleExamples/',dirIn,'/pics from segmentation/',fileIn,'_Evec_rM',num2str(Rmax(h))];   
            set(hev,'Position',screen);
            hgexport(hev, [fileHev,'.jpg'], hgexport('factorystyle'), 'Format', 'jpeg');
            fileHse = ['./output/ImgSeg/simpleExamples/',dirIn,'/pics from segmentation/',fileIn,'_Seg_',segAt,'_rM',num2str(Rmax(h))];   
            set(hse,'Position',screen);
            hgexport(hse, [fileHse,'.jpg'], hgexport('factorystyle'), 'Format', 'jpeg');

        end % loop over Rmax
    end % loop over segmentation thresholds {mean, optimal}
    
    keyboard
    
    % GONNA HAVE TO CLEAR SOME STUFF HERE
    clear imgIns
    close all
    
end % loop over directory (of files) for different image sizes 

