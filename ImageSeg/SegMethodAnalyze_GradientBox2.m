function SegMethodAnalyze_GradientBox2(imdimL,RmaxL)

% This script (written specific for GaussianBox images for now) is going to
% loop through all MAT files containing output of LoopImgSegMethodsD & 
% SegmentMethod functions. It will take files using one method and
% calculate correctness of segmentation (still figuring out method) and
% plot it vs. KL-Divergence of the two Gaussian Distributions from which
% pixel values for inside and outside rectangles were drawn.

screen = get(0,'ScreenSize');
[dirPre] = onCluster;

% imdimL = [7 9 11 13 15 21 23 25 27 29 31];
% RmaxL = [1 2 3 4];

Method = {'AA','AAnrm','GL','GLnrm','N&G','SKHAdj','SKHEuc','SKHD&O','N&G ME','SKH ME'}; % 'THmean',   ,'max Info'
segAtL = {'mean','opt'};
line_color = {'y','y','g','g','c','r','b','m','c','r','k'};
line_mrkr = {'o','x','+','*','s','d','^','<','>','v','h'};
line_style = {'-','--','-','--','-','-','--','-.','-','--','-'};

%% Preallocate (Note: 10 = number of files / pixel dist variances hardcoded for now)
infoMtot = zeros(numel(Method)+1, numel(imdimL), numel(RmaxL));
infoStot = zeros(numel(Method)+1, numel(imdimL),numel(RmaxL));
infoMaxTot = zeros(numel(imdimL),numel(RmaxL));

%% Loop through Methods and Different Images to Calculate Information Score
for g = 1:numel(segAtL) % Segmenting at mean or optimal threshold (or other). 
    segAt = segAtL{g};

    for f = 1:numel(imdimL)

        imdim = imdimL(f);            
        ximg = imdim;
        yimg = imdim;
        fileIn = [num2str(ximg),'x',num2str(yimg)];

        for h = 1:numel(RmaxL) % Radius of extent of neighborhood for adjacency matrix calculation 
            Rmax = RmaxL(h);

            for i = 1:numel(Method)
                Meth = Method{i};
                
                disp([segAt,' - ', num2str(imdim),' - ',num2str(Rmax),' - ',Meth])

                % Loop thru segmentation output MAT files in a directory...
                files = dir([dirPre,'output/ImgSeg/simpleExamples/GradientBox/data_from_segmentation/',fileIn,'/*',fileIn,'*',Meth,'-*','rM',num2str(Rmax),'*.mat']);

                % use 2nd Eigenvector if Graph Laplacian or Normalized Avg Association, 1st otherwise.
                if( ~isempty(strfind(Meth,'GL')) || ~isempty(strfind(Meth,'AAnrm')) )
                    B = 2;
                elseif( ~isempty(strfind(Meth,'N&G')) )
                    B = 3; % Not sure why 3rd eigenvector Looks Best for N&G.
                else
                    B = 1; % Dominant Eigenvector for Modularity Methods.
                end


                for j = 1:numel(files) % loop over files in directory
                    fname = files(j).name;

                    load([dirPre,'output/ImgSeg/simpleExamples/GradientBox/data_from_segmentation/',fileIn,'/',fname]);

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
                    
                    infoM(i,j) = mean(Io);
                    infoS(i,j) = std(Io);

                    if(j==1 && 0)
                        % Plot first 3 Eigenvectors for the near binary box in box case.
                        figure
                        subplot(numel(EvecsML)+1,1,1), imagesc(imEnsInfo.imEns(:,:,1)), colorbar, title(fname)
                        for n = 1:numel(EvecsML)
                            subplot(numel(EvecsML)+1,1,n+1), imagesc(EvecsML{n}(:,:,1)), colorbar, xlabel(['Evec ',num2str(n),' (\lambda=',num2str(EvalsML(n)),')'])
                        end
                        
                        figure,
                        subplot(131), imagesc(EvecVizF( EvecsML{B}, 0.01))
                        subplot(132), imagesc(segMean{B})
                        subplot(133), imagesc(segOpt{B})

                        keyboard

                    end


                    if (i==numel(Method))

                        % Segmentation by Thresholding at Mean Pixel Value
                        for m = 1:size(imEnsInfo.imEns,3)
                            segTH = (imEnsInfo.imEns(:,:,m) > mean(mean(imEnsInfo.imEns(:,:,m))));
                            IoTH(m) = calc_infoB(segTH, imEnsInfo.gndTruth);
                            infoM(i+1,j) = mean(IoTH);
                            infoS(i+1,j) = std(IoTH);
                        end 
                    end

                end % loop over files from given segmentation method in the segmentation output directory

            end % loop over methods


            %% Plot Single Eigenvector used for each Input Box & each Method
%             hev = figure; % Figure to plot all the eigenvectors for different methods and input images.
%             hse = figure; % Figure to plot segmentations for different methods and input images.

            %% Make a Bar Plot of Information Scores for each Method with Insets of Segmentations
%             [Y,I] = sort(infoM); % sort Methods by Information Score
%             figure, set(gcf,'Position',[0 0 1300 400]), bar(infoM(I),'w','LineWidth',2)
%             axis([0 numel(Method)+1 0 1])
%             set(gca,'XTickLabel',{Method{I}},'FontSize',14)
%             title(['Input Image: ',num2str(imEnsInfo.ximg),'x',num2str(imEnsInfo.yimg),' Gradient Box : SegAt ',segAt,' rM=',num2str(Rmax)],'FontWeight','Bold','Fontsize',20)
%             xlabel('Segmentation Method','FontWeight','Bold','Fontsize',18)
%             ylabel('Information (Segmented at Mean Pixel Value)','FontWeight','Bold','Fontsize',18)
%             % Gradient Box Input Image in inset intop corner
%             axes('Position',[.15 .60 .20 .30]);
%             box on
%             imagesc(imEnsInfo.imEns), colormap('bone')
%             set(gca,'xtick',[],'ytick',[])
%             % Plot Evector of Method in inset
%             for i = 1:numel(Method)%-1
% 
%                 % Loop thru segmentation output MAT files in a directory...
%                 files = dir([dirPre,'output/ImgSeg/simpleExamples/',dirIn,'/data from segmentation/*',num2str(imdim),'x',num2str(imdim),'*',Method{I(i)},'-*','rM',num2str(Rmax),'*.mat']);
%                 disp(['Method: ',Method{I(i)}])
% 
%                 % use 2nd Eigenvector if Graph Laplacian or Normalized Avg Association, 1st otherwise.
%                 if( ~isempty(strfind(Method{i},'GL')) || ~isempty(strfind(Method{i},'AAnrm')) )
%                     B = 2;
%                 elseif( ~isempty(strfind(Method{i},'N&G')) )
%                     B = 3; % Not sure why 3rd eigenvector Looks Best for N&G.
%                 else
%                     B = 1; % Dominant Eigenvector for Modularity Methods.
%                 end
% 
% 
%                 for j = 1:numel(files) % loop over files in directory
% 
%                     iid = files(j);
% 
%                     % Loop thru segmentation output MAT files in a directory...
%                     load([dirPre,'output/ImgSeg/simpleExamples/',dirIn,'/data from segmentation/',iid.name]);
% 
%                     % Gradient Box Input Image in inset intop corner (i+1)/(numel(Method)+2)
%                     if infoM(I(i)) > 0.5
%                         xx = infoM(I(i))-0.18;
%                         yy = infoM(I(i))-0.12;
%                     else
%                         xx = infoM(I(i))+0.12;
%                         yy = infoM(I(i))+0.18;
%                     end
%                     axes('Position',[(.1+i/(numel(infoM)+4)) xx .03 .05]);
%                     box on
%                     try
%                         imagesc(EvecVizF(EvecsML{B}(:,:,j),0.001)) %MethFlgs.vizSlope))
%                     catch
%                         imagesc(seg{B})
%                     end
%                     set(gca,'xtick',[],'ytick',[])
%                     % plot segmentation
%                     axes('Position',[(.1+i/(numel(infoM)+4)) yy .03 .05]);
%                     box on
%                     imagesc(seg{B})
%                     set(gca,'xtick',[],'ytick',[])
%                 end
% 
%             end
% 
% 
% %             keyboard
% 
% %             % SAVE HEV AND HSE IMAGES
% %             fileHev = ['./output/ImgSeg/simpleExamples/',dirIn,'/pics from segmentation/',fileIn,'_Evec_rM',num2str(Rmax)];   
% %             set(hev,'Position',screen);
% %             hgexport(hev, [fileHev,'.jpg'], hgexport('factorystyle'), 'Format', 'jpeg');
% %             fileHse = ['./output/ImgSeg/simpleExamples/',dirIn,'/pics from segmentation/',fileIn,'_Seg_',segAt,'_rM',num2str(Rmax)];   
% %             set(hse,'Position',screen);
% %             hgexport(hse, [fileHse,'.jpg'], hgexport('factorystyle'), 'Format', 'jpeg');


                    % Collect InfoM & InfoS information into 3-dim tensor to be plotted later.
                    infoMax = infoM(end); % This is segmenting at mean pixel value for well separated distributions
                    infoMtot(:,f,h) = infoM;
                    infoStot(:,f,h) = infoS;
                    infoMaxTot(f,h) = infoMax;
                    
                    
                    

        end % loop over Rmax

%         keyboard

        % GONNA HAVE TO CLEAR SOME STUFF HERE
        clear imgIns
        close all

    end % loop over directory (of files) for different image sizes 

    keyboard

        % Save info data Tensors
        save([dirPre,'output/ImgSeg/simpleExamples/GradientBox/data_from_segmentation/SegAnalInfo_R',num2str(numel(RmaxL)), ...
            '_D',num2str(numel(imdimL)),'_M',num2str(numel(Method)+1),'_segAt',segAt],...
            'infoMtot','infoStot','infoMaxTot','imdimL','RmaxL','Method','segAt')
    

end