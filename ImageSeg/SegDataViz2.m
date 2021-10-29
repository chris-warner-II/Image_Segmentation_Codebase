% SegDataViz2:
% Script to plot a series of imagesc colormap plots.  Show image dimension
% vs. Rmax with Information on the color axis. Make multiple plots for 
% variations in 2 dimensions:
%     (1). Different Variances on Gaussian Pixel Distribution
%     (2). Different Methods {ThreshMean, AA, GL, Mod, ... }
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
if ~exist([dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/pics_from_segmentation/SegDataViz2'],'dir')
    mkdir([dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/pics_from_segmentation/SegDataViz2']);
end

Methodb = [Method,'Seg@Mn'];

%%
for i = 1:numel(Methodb)
    for j = 1:numel(sig1)
        xx = reshape(infoMnrmTot(i,Ia(j),:,:),size(infoMnrmTot,3),size(infoMnrmTot,4));
        his = figure; imagesc(xx,[0,1]), 
        title(['Normalized Information : ',Methodb{i},' - \sigma = ',num2str(sig1(Ia(j)))],'FontWeight','Bold','FontSize',20)
        set(gca,'YTick',[1:size(xx,1)],'YTickLabel',imdimL)
        xlabel('Max Radius of Influence','FontWeight','Bold','FontSize',18)
        ylabel('Image Dimensions','FontWeight','Bold','FontSize',18)
        colorbar
        
        % SAVE THIS I VS SIG FIGURE (HIS)
        fname = ['imDim_v_Rmax_v_Info_segAt',segAt,'_',Methodb{i},'_',num2str(j),'sig',num2str(sig1(Ia(j)))];
        fname(fname=='.')='';
        fileHis = [dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/pics_from_segmentation/SegDataViz2/',fname];   
%         set(his,'Position',screen);
        hgexport(his, [fileHis,'.jpg'], hgexport('factorystyle'), 'Format', 'jpeg');

        close(his)
        
    end
end