% SegDataViz1:
% Script to plot Information vs. Variance of Pixel Distributions for inside
% and outside box of "GaussianBox" images.  Plot lines with errorbars for
% all methods and make multiple plots for variations in 2 dimensions:
%     (1). Different Image dimensions [7x7, 9x9, ..., 31x31]
%     (2). Different Rmax values (rmax controls # nearest neighbors
%          considered in building graph from image.

screen = get(0,'ScreenSize');
[dirPre] = onCluster;

InputImage = 'GaussianBox';

%% Load Input file (Can run for segmentation using threshold at mean pixel value or at optimal value).
switch SegAtIn
    case 'MeanPix'
        load([dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/data_from_segmentation/SegAnalInfo_R4_S10_D13_M11_segAtmean'])
    case 'Optimal'
        load([dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/data_from_segmentation/SegAnalInfo_R4_S10_D13_M11_segAtopt'])
end

%%
if ~exist([dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/pics_from_segmentation/SegDataViz1'],'dir')
    mkdir([dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/pics_from_segmentation/SegDataViz1']);
end

%%
line_color = {'y','y','g','g','c','r','b','m','c','r','k'};
line_mrkr = {'o','x','+','*','s','d','^','<','>','v','h'};
line_style = {'-','--','-','--','-','-','--','-.','-','--','-'};

%%
for f = 1:numel(imdimL)
        imdim = imdimL(f);
        fileIn = [num2str(imdim),'x',num2str(imdim),'_gauss_center'];

    for h = 1:numel(RmaxL) % Radius of extent of neighborhood for adjacency matrix calculation
        Rmax = RmaxL(h);

        his=figure; hold on
        axis([ 0 0.4 0 1.3])
        for bb = 1:numel(Method)+1
            errorbar( sig1(Ia),infoMnrmTot(bb,Ia,f,h), infoSnrmTot(bb,Ia,f,h), 'LineWidth',2, 'LineStyle',line_style{bb}, ...
                'Color',line_color{bb}, 'Marker',line_mrkr{bb}, 'MarkerSize',10 ) % 
        end
        plot(sig1(Ia), ones(1,numel(sig1)) ,'k--');
        legend([Method,'MeanTH'],'Location','NorthEast','FontSize',12,'FontWeight','Bold') %
        title(['Info vs. \sigma',' with  r_{max}=',num2str(Rmax),' & TH ',segAt],'FontWeight','Bold','FontSize',20)
        xlabel('\sigma of Gaussian prob Dists','FontWeight','Bold','FontSize',18)
        ylabel(['Normalized Information (Max = ',num2str(infoMaxTot(f,h)),')'],'FontWeight','Bold','FontSize',18)

        % WAS DISPLAYING EXAMPLES OF INPUT IMAGES.  NEED TO FIX.
%         for j = 1:numel(files) % loop over files in directory
%             fname = files(j).name;
%             load([dirPre,'output/ImgSeg/simpleExamples/',dirIn,'/data_from_segmentation/',fname]);
%             imgIns(:,:,j) = imEnsInfo.imEns(:,:,1);
%         end
%         
% 
%         % Plot Inset of Input Images  % THIS ONE GOT IN THE WAY OF THE LEGEND
%         %             axes('Position',[0.85  0.8 .04 .08]);
%         %             box on
%         %             imagesc(imgIns(:,:,Ia(1)))
%         %             set(gca,'xtick',[],'ytick',[])
%         %
%         axes('Position',[0.68  0.8 .04 .08]);
%         box on
%         imagesc(imgIns(:,:,Ia(3)))
%         set(gca,'xtick',[],'ytick',[])
%         %
%         axes('Position',[0.5 0.8  .04 .08]);
%         box on
%         imagesc(imgIns(:,:,Ia(5)))
%         set(gca,'xtick',[],'ytick',[])
%         %
%         axes('Position',[0.32  0.8 .04 .08]);
%         box on
%         imagesc(imgIns(:,:,Ia(7)))
%         set(gca,'xtick',[],'ytick',[])
%         %
%         axes('Position',[0.15  0.8 .04 .08]);
%         box on
%         imagesc(imgIns(:,:,Ia(9)))
%         set(gca,'xtick',[],'ytick',[])

        % SAVE THIS I VS SIG FIGURE (HIS)
        fname = [fileIn,'_IvsSig_segAt',segAt,'_rM',num2str(Rmax)];
        fileHis = [dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/pics_from_segmentation/SegDataViz1/',fname];   
        set(his,'Position',screen);
        hgexport(his, [fileHis,'.jpg'], hgexport('factorystyle'), 'Format', 'jpeg');
        
        close(his)

    end
    
end