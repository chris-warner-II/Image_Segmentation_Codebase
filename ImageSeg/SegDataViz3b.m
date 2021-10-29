% SegDataViz3b: (Building on SegDataViz3)  
%
% Same thing as SegDataViz3 with information value plotted Rmax vs. ImDim
% with subplots for different methods.  Difference here is that all
% information is plotted relative to the information gained by segmenting
% at the mean pixel value of the image.
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
if ~exist([dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/pics_from_segmentation/SegDataViz3b'],'dir')
    mkdir([dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/pics_from_segmentation/SegDataViz3b']);
end

Methodb = [Method,'Seg@Mn'];

%%
for j = 1:numel(sig1)
    
    his = figure;
    
    for i = 1:numel(Methodb)
      
        if(i==11) % plot actual results from segmenting at mean pixel value
            xx = reshape(infoMnrmTot(i,Ia(j),:,:),size(infoMnrmTot,3),size(infoMnrmTot,4));
        else
            xx = reshape(infoMnrmTot(i,Ia(j),:,:),size(infoMnrmTot,3),size(infoMnrmTot,4)) - ...
                reshape(infoMnrmTot(end,Ia(j),:,:),size(infoMnrmTot,3),size(infoMnrmTot,4));
        end
        subplot(4,3,i); % hardcode 12 subplots since there are 11 Methods
        imagesc(xx',[-1,1]), 
        set(gca,'XTick',[1:size(xx,1)],'XTickLabel',imdimL)
        %
        if(i==10)
            ylabel('Rmax','FontWeight','Bold','FontSize',16)
            xlabel('Im Dim','FontWeight','Bold','FontSize',16)
        end
        %
        if (i==1)
            title([Methodb{i}],'FontWeight','Bold','FontSize',16) % ,' : \sigma = ',num2str(sig1(j))
        else
            title([Methodb{i}],'FontWeight','Bold','FontSize',18)
        end
        
    end
    
    % Plot a colorbar in the 12th spot
    q = subplot(4,3,12,'Visible','off','CLim',[-1 1]); 
    imagesc(rand(size(xx')),'Visible','off')  
    c = colorbar('North');
    xlabel(c,{'Information beyond Segmenting', 'Image at Mean Pixel Value'},'FontWeight','Bold','FontSize',16)
    caxis([-1 1])
    axis('off')
    title(['\sigma = ',num2str(sig1(Ia(j)))],'FontWeight','Bold','FontSize',16)
    
    % SAVE THIS I VS SIG FIGURE (HIS)
    fname = ['imDim_v_Rmax_v_Info_segAt',segAt,'_',num2str(j),'sig',num2str(sig1(Ia(j)))];
    fname(fname=='.')='';
    fileHis = [dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/pics_from_segmentation/SegDataViz3b/',fname];   
    set(his,'Position',screen);
    hgexport(his, [fileHis,'.jpg'], hgexport('factorystyle'), 'Format', 'jpeg');

    close(his)
        
end