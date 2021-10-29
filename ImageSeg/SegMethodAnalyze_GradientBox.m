function SegMethodAnalyze_GradientBox(dirIn)
%
% USE SEGMETHODANALYZE_GRADEINTBOX2 !!!
%
% This script (written specific for GaussianBox images for now) is going to
% loop through all MAT files containing output of LoopImgSegMethodsD & 
% SegmentMethod functions. It will take files using one method and
% calculate correctness of segmentation (still figuring out method) and
% plot it vs. KL-Divergence of the two Gaussian Distributions from which
% pixel values for inside and outside rectangles were drawn.

dirIn = 'GradientBox';     
 
screen = get(0,'ScreenSize');

imdimL = [7 9 11 13 15 21]; %  31
RmaxL = [1 2 3 4];
segAtL = {'mean','opt'}; % ,'opt'
Method = {'AA','AAnrm','GL','GLnrm','N&G','SKH Adj','SKH Euc','SKH D&O','N&G ME','SKH ME'}; % 'THmean',   ,'max Info'
line_color = {'y','y','g','g','c','r','b','m','c','r','k'};
line_mrkr = {'o','x','+','*','s','d','^','<','>','v','h'};
line_style = {'-','--','-','--','-','-','--','-.','-','--','-'};

for h = 1:numel(Rmax)

    for i = 1:numel(Method)-1

        % Loop thru segmentation output MAT files in a directory...
        files = dir(['./output/ImgSeg/simpleExamples/',dirIn,'/data from segmentation/*',num2str(imSz),'x',num2str(imSz),'*',Method{i},'-*','rM',num2str(Rmax(h)),'*.mat']);
        disp(['Method: ',Method{i}])

        % use 2nd Eigenvector if Graph Laplacian, 1st otherwise.
        if( ~isempty(strfind(Method{i},'GL')) || ~isempty(strfind(Method{i},'AAnrm')) )
            B = 2;
        else
            B = 1;
        end


        for j = 1:numel(files) % loop over files in directory
            iid = files(j);
            fname = iid.name;
            disp(['Image # ',num2str(j),'/',num2str(numel(files))])

            load(['./output/ImgSeg/simpleExamples/',dirIn,'/data from segmentation/',iid.name]);
            Io = zeros(1,size(seg{1},3));


            for k = 1:size(seg{1},3)
                % NOTE: Should I normalize by max info available in image or patch?
                %       Maybe possible for binary ground truth (with just 2 segments)
                Io(k) = calc_infoB(seg{B}(:,:,k), imEnsInfo.gndTruth);
                
                
                % For Debugging Info Calculation, plot Segmentation, ground truth and Evector
                if(0)
                    figure, set(gcf,'Position',[0 0 1200 400]),
                    if (~ischar(EvecMLs))
                        x = EvecVizF(EvecMLs{B}(:,:,j),MethFlgs.vizSlope);
                        subplot(131),imagesc(x), colorbar  
                        title(num2str(mean(x(:))));
                    end
                    subplot(132),imagesc(seg{B}),
                    title(['Method: ',Method{i},'   Info: ',num2str(Io(k))])
                    subplot(133),imagesc(imEnsInfo.gndTruth{1}(:,:,j))

%                     keyboard
                end
                
            end

            % Pack into Matrices the important information for this figure
            infoM(i,j) = mean(Io);
            infoS(i,j) = std(Io);
            
%             figure, imagesc(seg{B}),title([Method{i},' ',num2str(infoM(i,j))])



        end % loop over files from given segmentation method in the segmentation output directory

    end
    
    
    
    %% Make a Bar Plot of Information Scores for each Method with Insets of Segmentations
    [Y,I] = sort(infoM); % sort Methods by Information Score
    figure, set(gcf,'Position',[0 0 1300 400]), bar(infoM(I),'w','LineWidth',2)
    axis([0 numel(Method) 0 1])
    set(gca,'XTickLabel',{Method{I}},'FontSize',14)
    title(['Input Image: ',num2str(imEnsInfo.ximg),'x',num2str(imEnsInfo.yimg),' Gradient Box'],'FontWeight','Bold','Fontsize',20)
    xlabel('Segmentation Method','FontWeight','Bold','Fontsize',18)
    ylabel('Information (Segmented at Mean Pixel Value)','FontWeight','Bold','Fontsize',18)
    % Gradient Box Input Image in inset intop corner
    axes('Position',[.15 .60 .20 .30]);
    box on
    imagesc(imEnsInfo.imEns), colormap('bone')
    set(gca,'xtick',[],'ytick',[])
    % Plot Evector of Method in inset
    for i = 1:numel(Method)-1
        
        % Loop thru segmentation output MAT files in a directory...
        files = dir(['./output/ImgSeg/simpleExamples/',dirIn,'/data from segmentation/*',num2str(imSz),'x',num2str(imSz),'*',Method{I(i)},'-*','rM',num2str(Rmax(h)),'*.mat']);
        disp(['Method: ',Method{I(i)}])

        % use 2nd Eigenvector if Graph Laplacian, 1st otherwise.
        if( ~isempty(strfind(Method{I(i)},'GL')) || ~isempty(strfind(Method{I(i)},'AAnrm')) )
            B = 2;
        else
            B = 1;
        end


        for j = 1:numel(files) % loop over files in directory
            
            iid = files(j);

            % Loop thru segmentation output MAT files in a directory...
            load(['./output/ImgSeg/simpleExamples/',dirIn,'/data from segmentation/',iid.name]);


            % Gradient Box Input Image in inset intop corner (i+1)/(numel(Method)+2)
            if infoM(I(i)) > 0.5
                xx = infoM(I(i))-0.12;
                yy = infoM(I(i))-0.18;
            else
                xx = infoM(I(i))+0.12;
                yy = infoM(I(i))+0.18;
            end
            axes('Position',[(.1+i/(numel(infoM)+4)) xx .03 .05]);
            box on
            try
                imagesc(EvecVizF(EvecMLs{B}(:,:,j),MethFlgs.vizSlope))
            catch
                imagesc(seg{:})
            end
            set(gca,'xtick',[],'ytick',[])
            % plot segmentation
            axes('Position',[(.1+i/(numel(infoM)+4)) yy .03 .05]);
            box on
            imagesc(seg{B})
            set(gca,'xtick',[],'ytick',[])
        end
        
    end

    

end

keyboard




% I EXISTED BEFORE YOU ENTERED THE SCENE, AND I WILL PERSIST AFTER !!