function SegMethodAnalyze_BSDS(dirIn)

% This script (written specific for GaussianBox images for now) is going to
% loop through all MAT files containing output of LoopImgSegMethodsD & 
% SegmentMethod functions. It will take files using one method and
% calculate correctness of segmentation (still figuring out method) and
% plot it vs. KL-Divergence of the two Gaussian Distributions from which
% pixel values for inside and outside rectangles were drawn.

dirIn = 'BSDS_proc';              
 
Rmax = [1];

Method = {'THmean','AA','AAnrm','GL','GLnrm','Mod N&G','Mod SKH Adj','Mod SKH Euc','Mod SKH D&O','Mod N&G ME','Mod SKH ME','max Info'};
line_color = {'k','y','y','g','g','c','r','b','m','c','r'};
line_mrkr = {'o','x','+','*','s','d','^','<','>','v','h'};
line_style = {'-','-','--','-','--','-','-','--','-.','-','--'};

for h = 1:numel(Rmax)

    for i = 1:numel(Method)-1

        % Loop thru segmentation output MAT files in a directory...
        files = dir(['./output/ImgSeg/simpleExamples/',dirIn,'/data from segmentation/*',Method{i},'-*','rM',num2str(Rmax(h)),'*.mat']);
        disp(['Method: ',Method{i}])

        % use 2nd Eigenvector if Graph Laplacian, 1st otherwise.
        if( ~isempty(strfind(Method{i},'GL')) || ~isempty(strfind(Method{i},'AAnrm')) )
            B = 2;
        else
            B = 1;
        end
    
        
        

        for j = 1:10 %numel(files) % loop over files in directory
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
                if(1 & mod(k,3)==0 & mod(j,4)==0)
                    figure, set(gcf,'Position',[0 0 1200 400]),
                    subplot(141),imagesc(imEnsInfo.imEns(:,:,k)), colorbar  
                    if (~ischar(EvecMLs))
                        x = EvecVizF(EvecMLs{B}(:,:,k),MethFlgs.vizSlope);
                        subplot(142),imagesc(x), colorbar  
                        title(['Evec - mean ',num2str(mean(x(:)))]);
                    end
                    subplot(143),imagesc(seg{B}(:,:,k)),
                    title(['Method: ',Method{i},'   Info: ',num2str(Io(k))])
                    subplot(144),imagesc(imEnsInfo.gndTruth{k,1}) % Only showing single groundtruth for now.
                    title(['Ground Truth'])

                    pause
                    
                    close
                end
                
                
            end

            % Pack into Matrices the important information for this figure
            infoM(i,j) = mean(Io);
            infoS(i,j) = std(Io);
            
            if (i==1)
                imgIns(:,:,j) = imEnsInfo.imEns(:,:,1);
            end
            
        end % loop over files from given segmentation method in the segmentation output directory

    end

    % Upper and Lower bounds for Errorbars
    U = infoM + infoS;
    L = infoM - infoS;

    [Y,I] = sort(mean(D_KL));

    % Plot Performance of Different Methods.  Info vs. D_KL.
    figure, plot( mean(sig1(:,I)),mean(D_KL(:,I)),'LineWidth',2 )
    title('\sigma vs. D_{KL}','FontWeight','Bold','FontSize',20)
    xlabel('\sigma of Gaussian prob Dists','FontWeight','Bold','FontSize',18)
    ylabel('KL-Divergence between 2 Dists','FontWeight','Bold','FontSize',18)
    %
    figure, hold on
    axis([ 0 0.4 0 1])
    for bb = 1:numel(Method)-1
        errorbar( sig1(bb,I),infoM(bb,I), infoS(bb,I), 'LineWidth',2, 'LineStyle',line_style{bb}, 'Color',line_color{bb}, 'Marker',line_mrkr{bb}, 'MarkerSize',5 ) % 
    end
    infoMax = infoM(1,I(end)); % This is segmenting at mean pixel value for well separated distributions
    plot(sig1(bb,I), repmat(infoMax(1,1),1,numel(files)) ,'k--');
    legend(Method,'Location','SouthWest','FontSize',12,'FontWeight','Bold') %
    title(['Info vs. \sigma',' with  r_{max}=',num2str(Rmax(h))],'FontWeight','Bold','FontSize',20)
    xlabel('\sigma of Gaussian prob Dists','FontWeight','Bold','FontSize',18)
    ylabel('information','FontWeight','Bold','FontSize',18)
    
    % Plot Inset of Input Images
    axes('Position',[0.85  0.82 .04 .08]);
    box on
    imagesc(imgIns(:,:,I(1)))
    set(gca,'xtick',[],'ytick',[])
    %
    axes('Position',[0.5  0.82 .04 .08]);
    box on
    imagesc(imgIns(:,:,I(5)))
    set(gca,'xtick',[],'ytick',[])
    %
    axes('Position',[0.15  0.82 .04 .08]);
    box on
    imagesc(imgIns(:,:,I(9)))
    set(gca,'xtick',[],'ytick',[])

end

keyboard




% I EXISTED BEFORE YOU ENTERED THE SCENE, AND I WILL PERSIST AFTER !!