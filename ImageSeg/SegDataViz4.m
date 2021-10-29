% SegDataViz4: (Builds on all previous SegDataViz files)
%
% This calculates the area under the Info vs. Variance curve treating it
% kinda like a ROC curve.  It then reports one number for each Rmax & ImDim
% point.  I can distil all information down to one figure with a subplot
% for each Method and a grid of Rmax / ImDim points.

% SegAtIn = 'MeanPix'; % DELETE ME

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
if ~exist([dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/pics_from_segmentation/SegDataViz4'],'dir')
    mkdir([dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/pics_from_segmentation/SegDataViz4']);
end

Methodb = [Method,'Seg@Mn'];

%% %% %% %%

pp = sum(infoMnrmTot,2); % adding up area under info vs sigma curve.

his = figure;
    
for i = 1:numel(Methodb)

    xx = reshape(pp(i,:,:,:),size(infoMnrmTot,3),size(infoMnrmTot,4))./numel(sig1);
    subplot(4,3,i); % hardcode 12 subplots since there are 11 Methods
    imagesc(xx',[-1,1]), 
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
imagesc(rand(size(xx')),'Visible','off')  
c = colorbar('North');
xlabel(c,{'Area under I vs. \sigma curve', ['\sigma = ',num2str(sig1(Ia(1:2)))], num2str(sig1(Ia(3:6))), num2str(sig1(Ia(7:10)))},'FontWeight','Bold','FontSize',16)
caxis([0 1])
axis('off')

% SAVE THIS I VS SIG FIGURE (HIS)
fname = ['imDim_v_Rmax_v_AreaROCabs_segAt',segAt];
fname(fname=='.')='';
fileHis = [dirPre,'output/ImgSeg/simpleExamples/',InputImage,'/pics_from_segmentation/SegDataViz4/',fname];   
set(his,'Position',screen);
hgexport(his, [fileHis,'.jpg'], hgexport('factorystyle'), 'Format', 'jpeg');

close(his)