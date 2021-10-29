% This script/function written by CW 11/15 to plot multiple precision
% recall curves on the same isoF plot.  We can use this to compare
% different methods or different parameter settings for the same method.


[dirPre,sizeGoodIm] = onCluster;
addpath([dirPre,'images/BSDS_images/BSR/bench/benchmarks/'])


method = {'AAnrm','GLnrm','Mod_SKHAdj','Mod_N&G'}; % 'IsoDiff',


rM = {'1','3','5','10'};
shapes = {'d-','o-','s-','^-'};
         % {'r','m','b','c'};

ev = {'ev1','ev2o','ev3o','ev2','ev3','ev2w','ev3w'}; 
colors = {'r','m','b','c'}; % {'r','g','b','c','m','y','k'};

plot_max = 1; % set to 1 to plot average of Precision & Recall leading to largest F, regardless of threshold
              % set to 0 to plot average of Precision & Recall at each threshold value for each image

              
          
              
% Flag to use a gaussian kernel (sig=1) to preblur image before running network Kuramoto computation.
blur_flg=0;
if(blur_flg)
    blur_tag_M = '_blur_sig1';
    blur_tit = 'w/ Pre-Blurring (\sigma=1)';
else
    blur_tag_M = ''; % if we are not blurring.
    blur_tit = 'w/ No Pre-Blurring.';
end
blur_tag_I = 'blur_sz13_sig1/';
    
F_all = zeros(numel(rM),numel(ev),numel(method));
F_max = zeros(numel(method)+2,1);

mean_maxF = zeros(numel(rM),numel(ev),numel(method));
std_maxF = zeros(numel(rM),numel(ev),numel(method));

numFiles = zeros(numel(rM),numel(ev),numel(method));


% matrix of parameter values for nice labeling later.
param_matrix = cell(numel(rM),numel(ev));
for i = 1:numel(rM)
    for j = 1:numel(ev)
        param_matrix{i,j} = ['rM',rM{i},':',ev{j}];
    end   
end             
              
              
        

plot_IsoF_flag = 0; % no longer plotting this for now.

              
              
              
for A = 1:numel(method)              
              



    % Set up one figure with multiple subplot (one for each eigenvector combination)
    if(plot_IsoF_flag)
        h1 = openfig('isoF.fig','reuse'); % open figure
        ax1 = gca; % get handle to axes of figure
        fig1 = get(ax1,'children'); %get handle to all the children in the figure
        H = figure; %create new figure
        s1 = subplot(2,4,1); axis square %create and get handle to the subplot axes
        s2 = subplot(2,4,2); axis square
        s3 = subplot(2,4,3); axis square
        s8 = subplot(2,4,4); axis off
        s4 = subplot(2,4,5); axis square
        s5 = subplot(2,4,6); axis square
        s6 = subplot(2,4,7); axis square
        s7 = subplot(2,4,8); axis square
        %
        copyobj(fig1,s1); colormap(bone); %copy children to new parent axes i.e. the subplot axes
        copyobj(fig1,s2); colormap(bone); 
        copyobj(fig1,s3); colormap(bone);  
        copyobj(fig1,s4); colormap(bone);   
        copyobj(fig1,s5); colormap(bone);   
        copyobj(fig1,s6); colormap(bone);   
        copyobj(fig1,s7); colormap(bone);  

        close(h1);

    end


    for j = 1:numel(ev)
        for i = 1:numel(rM)

            if strcmp(method{A},'IsoDiff')
                evalDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/spectral/',method{A},'/benchmark_results/rM',rM{i},'/',ev{j},'/'];
            else
                evalDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/spectral/',method{A},'/benchmark_results/rM',rM{i},'/sDInf/sP0p2/',ev{j},'/'];
            end

            
            if(plot_IsoF_flag)
                [a, evalRes] = plot_eval(evalDir,[colors{i},shapes{i}],eval(['s',num2str(j)]),plot_max); 
            end
           
            
            
            
            % Loop through each image patch and grab maxF. So later I can compute their mean and std.
            files = dir([evalDir,'*_ev1.txt']);
            
            maxF = zeros(1,numel(files));
            
            numFiles(i,j,A) = numel(files);
            
            for k = 1:numel(files)
                
                filename = fullfile(evalDir,files(k).name);
                AA  = dlmread(filename);
                %cntR = AA(:, 2);
                sumR = AA(:, 3);
                %cntP = AA(:, 4);
                sumP = AA(:, 5);
                %
                cntR0 = AA(:, 6);                      % CW: added these new ways to compute pixel correspondence in evaluation_bdry_image.
                cntP0 = AA(:, 7);

                R = cntR0 ./ (sumR + (sumR==0)); % Note: this (sumR==0) in the denominator is just to ensure we do not
                P = cntP0 ./ (sumP + (sumP==0)); %       divide by zero. If sumR/P == 0, then we divide by 1. (CW).
                
                F = 2*P.*R./(P+R+((P+R)==0));    % F-measure.

                maxF(k) = max(F);
                
                k
                
            end % loop over image patches
            
            % found maxF at each image patch. Now get mean and std across image patches.
            mean_maxF(i,j,A) = mean(maxF);
            std_maxF(i,j,A) = std(maxF);
            
            
            
            
            
            
            

        end % loop over rM parameter
        
        
        
        
        if(plot_IsoF_flag)
            % Plot all against the Precision Recall curve of boundaries in the Raw Image Pixel Spatial Gradients
            evalDir = [dirPre,'images/BSDS_patch/101x101_ds1/benchmark_results/'];
            [a, evalRes] = plot_eval(evalDir,'ko--',eval(['s',num2str(j)]),plot_max);

            % Plot all against the Precision Recall curve of boundaries in the Blurred Image Pixel Spatial Gradients
            evalDir = [dirPre,'images/BSDS_patch/101x101_ds1/',blur_tag_I,'benchmark_results/'];
            [a, evalResB] =  plot_eval(evalDir,'go--',eval(['s',num2str(j)]),plot_max);


            title(ev{j})
            set(gca,'XTick',[0:0.2:1],'YTick',[0:0.2:1])
            if j==4
                xlabel('Recall')
                ylabel('Precision')
            end
            grid on
        end

        
        
        
        

    end % loop over ev parameter
    
    
    
    
    
    % Sort average F-measure results (for different parameter settings)
%     x = F_all(:,:,A);
%     [p1,q1] = sort(x(:),'descend');
    
    y = mean_maxF(:,:,A);
    z = std_maxF(:,:,A);
    [p2,q2] = sort(y(:),'descend');
    
    
    
    
    
    % Now, Get mean & std of F-measure for same image patches using just raw image pixels and using optimal gaussian blurring.
    imPixDir = [dirPre,'images/BSDS_patch/101x101_ds1/benchmark_results/'];
    imBlurDir = [dirPre,'images/BSDS_patch/101x101_ds1/',blur_tag_I,'benchmark_results/'];
    %
    for k = 1:numel(files)
                
        filename = fullfile(imPixDir,files(k).name);
        AA  = dlmread(filename);
        %cntR = AA(:, 2);
        sumR = AA(:, 3);
        %cntP = AA(:, 4);
        sumP = AA(:, 5);
        %
        cntR0 = AA(:, 6);                      % CW: added these new ways to compute pixel correspondence in evaluation_bdry_image.
        cntP0 = AA(:, 7);

        R = cntR0 ./ (sumR + (sumR==0)); % Note: this (sumR==0) in the denominator is just to ensure we do not
        P = cntP0 ./ (sumP + (sumP==0)); %       divide by zero. If sumR/P == 0, then we divide by 1. (CW).

        F = 2*P.*R./(P+R+((P+R)==0));    % F-measure.

        maxF_imPix(k) = max(F);

        k

    end
    %
    for k = 1:numel(files)
                
        filename = fullfile(imBlurDir,files(k).name);
        AA  = dlmread(filename);
        %cntR = AA(:, 2);
        sumR = AA(:, 3);
        %cntP = AA(:, 4);
        sumP = AA(:, 5);
        %
        cntR0 = AA(:, 6);                      % CW: added these new ways to compute pixel correspondence in evaluation_bdry_image.
        cntP0 = AA(:, 7);

        R = cntR0 ./ (sumR + (sumR==0)); % Note: this (sumR==0) in the denominator is just to ensure we do not
        P = cntP0 ./ (sumP + (sumP==0)); %       divide by zero. If sumR/P == 0, then we divide by 1. (CW).

        F = 2*P.*R./(P+R+((P+R)==0));    % F-measure.

        maxF_imBlur(k) = max(F);

        k

    end
    
    
    
    
    
    F_max2(A) = max(max(mean_maxF(:,:,A)));                                % using our average F-measure (mean & std) over image patches
    [rM_max2(A),ks_max2(A)] = find( mean_maxF(:,:,A) == F_max2(A) );
    std_max2(A) = std_maxF(rM_max2(A),ks_max2(A),A);
    
    



    if(plot_IsoF_flag)
        % Make a legend
        axes(s8)

        hold on
        %
        text(0.5,1,{'\color{red}{rM1}','\color{blue}{rM3}','\color{cyan}{rM5}','\color{magenta}{rM10}'},'VerticalAlignment','top','HorizontalAlignment','center','FontSize',16,'FontWeight','Bold')
        %

        if(plot_max)
            % open vs. filled
            scatter(0.5,0.7,200,'ko','LineWidth',2)
            text(0.55,0.7,'= const THs','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
            %
            scatter(0.5,0.6,200,'ko','filled','LineWidth',2)
            text(0.55,0.6,'= vary THs','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
        else
            scatter(0.5,0.6,200,'ko','LineWidth',2)
            text(0.55,0.6,'= max F','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
        end
        %
        plot([0.5,0.4],[0.5,0.5],'g-','LineWidth',1.5)
        text(0.5,0.5,' \color{gray}{= Iso F}','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
        %
        plot([0.5,0.4],[0.4,0.4],'ko--','LineWidth',1.5)
        text(0.5,0.4,'= ImPix','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
        %
        plot([0.5,0.3],[0.3,0.3],'go--','LineWidth',1.5)
        text(0.5,0.3,'= ImBlur','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
        % 
        axis([0 1 0.2 1])


        % Make a title
        method_tag = method{A};
        method_tag(method_tag=='_')=' ';

        annotation('textbox', [0 0.9 1 0.1],'String',[method_tag,' : Spectral Results'], ...
                    'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold')


        if(plot_max)
            max_tag = 'best&max';
        else
            max_tag = 'curve';
        end

        saveGoodImg(H,[dirPre,'../Documentation/Cosyne_2016/Overall_PR_curve_results_thinpbOFF/',method{A},'_',max_tag,'_Eig_allParams.jpg'],sizeGoodIm)
    end
    
    
    
    
    % Plot mean F-measure (mean & std) for each method for different parameters. Also plot imPix & imBlur F-measure values.
    H1 = figure; hold on
    errorbar(mean_maxF(:,:,A)',std_maxF(:,:,A)','LineStyle','none','Marker','o','LineWidth',2);
    %
    errorbar(8, mean(maxF_imPix), std(maxF_imPix),'LineStyle','none','Color','green','Marker','o','LineWidth',2)
    errorbar(9, mean(maxF_imBlur), std(maxF_imBlur),'LineStyle','none','Color','red','Marker','o','LineWidth',2)
    %
    plot([1 8],[mean(maxF_imPix) mean(maxF_imPix)],'g--','LineWidth',2)
    plot([1 9],[mean(maxF_imBlur) mean(maxF_imBlur)],'r--','LineWidth',2)
    %
    plot([1 8],[mean(maxF_imPix)+std(maxF_imPix) mean(maxF_imPix)+std(maxF_imPix)],'g--')
    plot([1 9],[mean(maxF_imBlur)+std(maxF_imBlur) mean(maxF_imBlur)+std(maxF_imBlur)],'r--')
    %
    bar(mean_maxF(:,:,A)')
    bar(8, mean(maxF_imPix),'g')
    bar(9, mean(maxF_imBlur),'r')
    %
    set(gca,'XTick',1:9,'XTickLabel',[ev(:);'ImPix';'ImBlur'],'FontSize',12,'FontWeight','Bold')
    title([method{A},' ',blur_tit,' :  Spectral Methods '],'FontSize',18,'FontWeight','Bold')
    xlabel('Parameters','FontSize',18,'FontWeight','Bold')
    ylabel(['Average F-measure (across ~ ',num2str(round(mean(mean(numFiles(:,:,A))))),' image patches)'],'FontSize',18,'FontWeight','Bold')
    grid on
    saveGoodImg(H1,[dirPre,'../Documentation/Cosyne_2016/Overall_PR_curve_results_thinpbOFF/Fmeasure_ErrBar_',method{A},blur_tag_M,'_Eig_allParams.jpg'],sizeGoodIm)
    close(H1)
    
    
    
    
end % loop over method



