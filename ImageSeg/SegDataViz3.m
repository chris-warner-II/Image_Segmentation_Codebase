% SegDataViz3: (Building on SegDataViz2)  
%
% Put different Methods in a Subplot instead of plotting them individually.
% Script to plot a series of imagesc colormap plots.  Show image dimension
% vs. Rmax with Information on the color axis. Make multiple plots for 
% Different Variances on Gaussian Pixel Distribution.
%

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
if ~exist([dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/pics_from_segmentation/SegDataViz3'],'dir')
    mkdir([dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/pics_from_segmentation/SegDataViz3']);
end

Methodb = [Method,'Seg@Mn'];

%%
for j = 1:numel(sig1)
    
    his = figure;
    
    for i = 1:numel(Methodb)
        
        xx = reshape(infoMnrmTot(i,Ia(j),:,:),size(infoMnrmTot,3),size(infoMnrmTot,4));
        subplot(4,3,i); % hardcode 12 subplots since there are 11 Methods
        imagesc(xx',[0,1]), 
        set(gca,'XTick',[1:size(xx,1)],'XTickLabel',imdimL)
        %
        if (i==1)
            title([Methodb{i}],'FontWeight','Bold','FontSize',16) % ,' : \sigma = ',num2str(sig1(j))
        else
            title([Methodb{i}],'FontWeight','Bold','FontSize',18)
        end
        %
        if(i==10)
            ylabel('Rmax','FontWeight','Bold','FontSize',16)
            xlabel('Im Dim','FontWeight','Bold','FontSize',16)
        end

    end
    
    % Plot a colorbar in the 12th spot
    q = subplot(4,3,12,'Visible','off','CLim',[0 1]); 
    imagesc(zeros(size(xx')),'Visible','off')  
    c = colorbar('North');
    xlabel(c,'Normalized Information Colorbar','FontWeight','Bold','FontSize',16)
    axis('off')
    title(['\sigma = ',num2str(sig1(Ia(j)))],'FontWeight','Bold','FontSize',16)
    
    % SAVE THIS I VS SIG FIGURE (HIS)
    fname = ['imDim_v_Rmax_v_Info_segAt',segAt,'_',num2str(j),'sig',num2str(sig1(Ia(j)))];
    fname(fname=='.')='';
    fileHis = [dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/pics_from_segmentation/SegDataViz3/',fname];   
    set(his,'Position',screen);
    hgexport(his, [fileHis,'.jpg'], hgexport('factorystyle'), 'Format', 'jpeg');

    close(his)
        
end