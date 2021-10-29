function SegMethodAnalyze_GaussianBox(imdimL,RmaxL)

% This script (written specific for GaussianBox images for now) is going to
% loop through all MAT files containing output of LoopImgSegMethodsD & 
% SegmentMethod functions. It will take files using one method and
% calculate correctness of segmentation (still figuring out method) and
% plot it vs. KL-Divergence of the two Gaussian Distributions from which
% pixel values for inside and outside rectangles were drawn.

screen = get(0,'ScreenSize');
[dirPre] = onCluster;

%imdimL = [7 9 11 13 15 17 19 21 23 25 27 29 31];
%RmaxL = [1 2 3 4];

Method = {'AA','AAnrm','GL','GLnrm','N&G','SKHAdj','SKHEuc','SKHD&O','N&G ME','SKH ME'}; 
Methodb = [Method,'THmean'];
segAtL = {'mean','opt'};
line_color = {'y','y','g','g','c','r','b','m','c','r','k'};
line_mrkr = {'o','x','+','*','s','d','^','<','>','v','h'};
line_style = {'-','--','-','--','-','-','--','-.','-','--','-'};


%% Preallocate (Note: 10 = number of files / pixel dist variances hardcoded for now)
infoMnrmTot = zeros(numel(Method)+1, 10, numel(imdimL), numel(RmaxL));
infoSnrmTot = zeros(numel(Method)+1, 10, numel(imdimL),numel(RmaxL));
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

            for i = 1:numel(Method) % Add THmean to the methods ??
                Meth = Method{i};
                
                disp([segAt,' - ', num2str(imdim),' - ',num2str(Rmax),' - ',Meth])

                % Loop thru segmentation output MAT files in a directory...
                files = dir([dirPre,'output/ImgSeg/simpleExamples/GaussianBox/data_from_segmentation/',fileIn,'/','*',fileIn,'*',Meth,'-*','rM',num2str(Rmax),'*.mat']);

                % use 2nd Eigenvector if Graph Laplacian or Normalized Avg Association, 1st otherwise.
                if( ~isempty(strfind(Meth,'GL')) || ~isempty(strfind(Meth,'AAnrm')) )
                    B = 2;
                elseif( ~isempty(strfind(Meth,'N&G')) )
                    B = 3; % Not sure why 3rd eigenvector Looks Best.
                else
                    B = 1; % Dominant Eigenvector for Modularity Methods.
                end


                for j = 1:numel(files) % loop over files in directory
                    fname = files(j).name;

                    load([dirPre,'output/ImgSeg/simpleExamples/GaussianBox/data_from_segmentation/',fileIn,'/',fname]);

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

            end % loop over methods

            
            if(0)
                figure, plot( sig1(Ia), D_KL(Ia), 'LineWidth',2 )
                title('\sigma vs. D_{KL}','FontWeight','Bold','FontSize',20)
                xlabel('\sigma of Gaussian prob Dists','FontWeight','Bold','FontSize',18)
                ylabel('KL-Divergence between 2 Dists','FontWeight','Bold','FontSize',18)
            end
            
            
            %% Collect InfoM & InfoS information into 4-dim tensor to be plotted later.
            [Y,Id] = sort(D_KL,'descend');
            [Y,Ia] = sort(D_KL,'ascend');
            
            infoMax = infoM(end,Ia(end)); % This is segmenting at mean pixel value for well separated distributions
            infoMnrm = infoM./infoMax;
            infoSnrm = infoS./infoMax; % Scaling information mean and std by max possible info
            
            infoMnrmTot(:,:,f,h) = infoMnrm;
            infoSnrmTot(:,:,f,h) = infoSnrm; 
            infoMaxTot(f,h) = infoMax;
             

        end % loop over Rmax
        
        clear imgIns % because different size images lead to different dimensions
        
    end % loop over directory (of files) for different image sizes 
    
    
    % Save info data Tensors
    save([dirPre,'output/ImgSeg/simpleExamples/GaussianBox/data_from_segmentation/SegAnalInfo_R',num2str(numel(RmaxL)),'_S', ...
        num2str(numel(files)),'_D',num2str(numel(imdimL)),'_M',num2str(numel(Method)+1),'_segAt',segAt],...
        'infoMnrmTot','infoSnrmTot','infoMaxTot','imdimL','RmaxL','Method','mu1','mu2','sig1','sig2','D_KL','Ia','Id','segAt')
    
end % loop over segmentation thresholds {mean, optimal}   


