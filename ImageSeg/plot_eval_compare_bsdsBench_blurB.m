% This script/function written by CW 11/15 to plot multiple precision
% recall curves on the same isoF plot.  We can use this to compare
% different methods or different parameter settings for the same method.


[dirPre,sizeGoodIm] = onCluster;
addpath([dirPre,'images/BSDS_images/BSR/bench/benchmarks/'])




sig = {'0p5','1','1p5','2','2p5','5','10'};
colors =  {'r','b','y','c','k','m','g'};

sz =  {'7','13','17','21','27','49','95'}; 
shapes = {'v-','d-','v-','^-','v-','s-','o-'};


%H=open('isoF.fig');

plot_max = 1; % set to 1 to plot average of Precision & Recall leading to largest F, regardless of threshold
              % set to 0 to plot average of Precision & Recall at each threshold value for each image

num_cPdist = 4; % number of distances allowed in correspondPixelsB function
cPd_colors = {'r','g','b','k'};


%% Preallocate memory for arrays to hold statistics (across image patches) for each different blur.

% justMethod - Absolute value of performance. Not relative to any baseline or null model.

% (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
justMethod.maxF_newMean_meanAccImgs = zeros(num_cPdist,numel(sig)+1);
justMethod.maxF_newMean_stdAccImgs = zeros(num_cPdist,numel(sig)+1);
justMethod.maxF_newStd_meanAccImgs = zeros(num_cPdist,numel(sig)+1);
%
justMethod.maxR_newMean_meanAccImgs = zeros(num_cPdist,numel(sig)+1);
justMethod.maxR_newMean_stdAccImgs = zeros(num_cPdist,numel(sig)+1);
justMethod.maxR_newStd_meanAccImgs = zeros(num_cPdist,numel(sig)+1);
%
justMethod.maxP_newMean_meanAccImgs = zeros(num_cPdist,numel(sig)+1);
justMethod.maxP_newMean_stdAccImgs = zeros(num_cPdist,numel(sig)+1);
justMethod.maxP_newStd_meanAccImgs = zeros(num_cPdist,numel(sig)+1);


% (4). New way to compute P-R-F - (Max value using best single ground truther.). 
justMethod.maxF_newMax_meanAccImgs = zeros(num_cPdist,numel(sig)+1);
justMethod.maxF_newMax_stdAccImgs = zeros(num_cPdist,numel(sig)+1);
%
justMethod.maxR_newMax_meanAccImgs = zeros(num_cPdist,numel(sig)+1);
justMethod.maxR_newMax_stdAccImgs = zeros(num_cPdist,numel(sig)+1);
%
justMethod.maxP_newMax_meanAccImgs = zeros(num_cPdist,numel(sig)+1);
justMethod.maxP_newMax_stdAccImgs = zeros(num_cPdist,numel(sig)+1);



% relImPix - Performance measures relative to Raw image Pixels.

% (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
relImPix.maxF_newMean_meanAccImgs = zeros(num_cPdist,numel(sig)+1);
relImPix.maxF_newMean_stdAccImgs = zeros(num_cPdist,numel(sig)+1);
relImPix.maxF_newStd_meanAccImgs = zeros(num_cPdist,numel(sig)+1);
%
relImPix.maxR_newMean_meanAccImgs = zeros(num_cPdist,numel(sig)+1);
relImPix.maxR_newMean_stdAccImgs = zeros(num_cPdist,numel(sig)+1);
relImPix.maxR_newStd_meanAccImgs = zeros(num_cPdist,numel(sig)+1);
%
relImPix.maxP_newMean_meanAccImgs = zeros(num_cPdist,numel(sig)+1);
relImPix.maxP_newMean_stdAccImgs = zeros(num_cPdist,numel(sig)+1);
relImPix.maxP_newStd_meanAccImgs = zeros(num_cPdist,numel(sig)+1);


% (4). New way to compute P-R-F - (Max value using best single ground truther.). 
relImPix.maxF_newMax_meanAccImgs = zeros(num_cPdist,numel(sig)+1);
relImPix.maxF_newMax_stdAccImgs = zeros(num_cPdist,numel(sig)+1);
%
relImPix.maxR_newMax_meanAccImgs = zeros(num_cPdist,numel(sig)+1);
relImPix.maxR_newMax_stdAccImgs = zeros(num_cPdist,numel(sig)+1);
%
relImPix.maxP_newMax_meanAccImgs = zeros(num_cPdist,numel(sig)+1);
relImPix.maxP_newMax_stdAccImgs = zeros(num_cPdist,numel(sig)+1);




numFiles = zeros(1,numel(sig)+1);



%% Do analysis on raw image pixels first (before doing blurring by different sigma)
%  NOTE: To be able to find measures relative to ImPix, you have to do this
%  inside the blur loop, I think.  No other way... Is there?
%
% Plot all against the Precision Recall curve of boundaries in the Unblurred Raw Image Pixel Spatial Gradients
evalDirIm = [dirPre,'images/BSDS_patch/101x101_ds1/benchmark_results/'];

% i=1;
% [H, evalRes] = plot_eval(evalDirIm,[colors{i},shapes{i}],H,plot_max); % This uses ev1.txt files

for d = 1:num_cPdist
    
    disp('Distance tolerance')
    disp(d)
    i=1;


    % Loop through each image patch and grab maxF. So later I can compute their mean and std.
    files = dir([evalDirIm,'*_d',num2str(d),'_ev2.txt']);


    maxF_new_max_imPix = zeros(1,numel(files));    % For these next 3, we compute P/R/F for each groundtruth associated with an
    maxR_new_max_imPix = zeros(1,numel(files));    % image patch.  Here we just grab the max value.  What is the algorithm's
    maxP_new_max_imPix = zeros(1,numel(files));    % best match to any single human?

    maxF_new_mean_imPix = zeros(1,numel(files));   % For a single image patch, the average P/R/F across the different groundtruths
    maxR_new_mean_imPix = zeros(1,numel(files));
    maxP_new_mean_imPix = zeros(1,numel(files));

    maxF_new_std_imPix = zeros(1,numel(files));    % Standard Deviation of P/R/F across different groundtruths for single image ptch
    maxR_new_std_imPix = zeros(1,numel(files));
    maxP_new_std_imPix = zeros(1,numel(files));

    numFiles(i) = numel(files);              % number of image patches processed for a given method or blur

    
    disp(['Looping thru ',num2str(numel(files)),' files.'])
    for k = 1:numel(files)

        filename = fullfile(evalDirIm,files(k).name);
        AA  = dlmread(filename);
        thr = AA(:, 1); 
        cP_d = AA(:, 2);
        %
        Rmean = AA(:, 3);
        Rstd = AA(:, 4);
        Rmax = AA(:, 5);
        %
        Pmean = AA(:, 6);
        Pstd = AA(:, 7);
        Pmax = AA(:, 8);
        %
        Fmean = AA(:, 9);
        Fstd = AA(:, 10);
        Fmax = AA(:, 11);
        %
        num_gTs = AA(:, 12);
        Fmax_whichGTs = AA(:, 13);
        
        
        maxF_new_max_imPix(k) = max(Fmax);
        ind = find(Fmax==max(Fmax));
        ind = ind(1);
        maxR_new_max_imPix(k) = Rmax(ind);
        maxP_new_max_imPix(k) = Pmax(ind);
        maxF_bestGT_imPix(k) = Fmax_whichGTs(ind);
        thr_new_max_imPix(k) = thr(ind);

        maxF_new_mean_imPix(k) = max(Fmean);
        ind = find(Fmean==max(Fmean));
        ind = ind(1);
        maxR_new_mean_imPix(k) = Rmean(ind);
        maxP_new_mean_imPix(k) = Pmean(ind);
        maxF_new_std_imPix(k) = Fstd(ind);   % std across groundtruths as threshold that gave max avgF value.
        thr_new_mean_imPix(k) = thr(ind);

        
        %k

    end % looping through image patches


    % justMethod - Take Mean & Std across image patches of Precision,Recall,F-measure computed in different ways.


    % (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
    justMethod.maxF_newMean_meanAccImgs(d,i) = mean(maxF_new_mean_imPix);
    justMethod.maxF_newMean_stdAccImgs(d,i) = std(maxF_new_mean_imPix);
    justMethod.maxF_newStd_meanAccImgs(d,i) = mean(maxF_new_std_imPix);
    %
    justMethod.maxR_newMean_meanAccImgs(d,i) = mean(maxR_new_mean_imPix);
    justMethod.maxR_newMean_stdAccImgs(d,i) = std(maxR_new_mean_imPix);
    justMethod.maxR_newStd_meanAccImgs(d,i) = mean(maxR_new_std_imPix);
    %
    justMethod.maxP_newMean_meanAccImgs(d,i) = mean(maxP_new_mean_imPix);
    justMethod.maxP_newMean_stdAccImgs(d,i) = std(maxP_new_mean_imPix);
    justMethod.maxP_newStd_meanAccImgs(d,i) = mean(maxP_new_std_imPix);


    % (4). New way to compute P-R-F - (Max value using best single ground truther.). 
    justMethod.maxF_newMax_meanAccImgs(d,i) = mean(maxF_new_max_imPix);
    justMethod.maxF_newMax_stdAccImgs(d,i) = std(maxF_new_max_imPix);
    %
    justMethod.maxR_newMax_meanAccImgs(d,i) = mean(maxR_new_max_imPix);
    justMethod.maxR_newMax_stdAccImgs(d,i) = std(maxR_new_max_imPix);
    %
    justMethod.maxP_newMax_meanAccImgs(d,i) = mean(maxP_new_max_imPix);
    justMethod.maxP_newMax_stdAccImgs(d,i) = std(maxP_new_max_imPix);


    % relImPix - calculate P/R/F different ways of method relative to raw image pixels.

    % (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
    relImPix.maxF_newMean_meanAccImgs(d,i) = 0;
    relImPix.maxF_newMean_stdAccImgs(d,i) = 0;
    relImPix.maxF_newStd_meanAccImgs(d,i) = 0;
    %
    relImPix.maxR_newMean_meanAccImgs(d,i) = 0;
    relImPix.maxR_newMean_stdAccImgs(d,i) = 0;
    relImPix.maxR_newStd_meanAccImgs(d,i) = 0;  % these are zero because it is imPix relative to imPix.
    %
    relImPix.maxP_newMean_meanAccImgs(d,i) = 0;
    relImPix.maxP_newMean_stdAccImgs(d,i) = 0;
    relImPix.maxP_newStd_meanAccImgs(d,i) = 0;


    % (4). New way to compute P-R-F - (Max value using best single ground truther.). 
    relImPix.maxF_newMax_meanAccImgs(d,i) = 0;
    relImPix.maxF_newMax_stdAccImgs(d,i) = 0;
    %
    relImPix.maxR_newMax_meanAccImgs(d,i) = 0;
    relImPix.maxR_newMax_stdAccImgs(d,i) = 0;
    %
    relImPix.maxP_newMax_meanAccImgs(d,i) = 0;
    relImPix.maxP_newMax_stdAccImgs(d,i) = 0;






    %% Repeat analysis in For Loop for each different Blurring amount

    for i = 2:numel(sig)+1
        
        disp(['Different blurring values ', num2str(i),' / ', num2str(numel(sig)+1)])
        disp(['blur_sz',num2str(sz{i-1}),'_sig',num2str(sig{i-1})])

        evalDir = [dirPre,'images/BSDS_patch/101x101_ds1/blur_sz',sz{i-1},'_sig',sig{i-1},'/benchmark_results/'];

        % [H, evalRes] = plot_eval(evalDir,[colors{i-1},shapes{i-1}],H,plot_max); % 




        % Loop through each image patch and grab maxF. So later I can compute their mean and std.
        files = dir([evalDir,'*_d',num2str(d),'_ev2.txt']);

        maxF_new_max = zeros(1,numel(files));    % For these next 3, we compute P/R/F for each groundtruth associated with an
        maxR_new_max = zeros(1,numel(files));    % image patch.  Here we just grab the max value.  What is the algorithm's
        maxP_new_max = zeros(1,numel(files));    % best match to any single human?

        maxF_new_mean = zeros(1,numel(files));   % For a single image patch, the average P/R/F across the different groundtruths
        maxR_new_mean = zeros(1,numel(files));
        maxP_new_mean = zeros(1,numel(files));

        maxF_new_std = zeros(1,numel(files));    % Standard Deviation of P/R/F across different groundtruths for single image ptch
        maxR_new_std = zeros(1,numel(files));
        maxP_new_std = zeros(1,numel(files));

        
        numFiles(i) = numel(files);              % number of image patches processed for a given method or blur
        disp(['Looping thru ',num2str(numel(files)),' files.'])
        for k = 1:numel(files)

            filename = fullfile(evalDir,files(k).name);
            AA  = dlmread(filename);
            thr = AA(:, 1); 
            cP_d = AA(:, 2);
            %
            Rmean = AA(:, 3);
            Rstd = AA(:, 4);
            Rmax = AA(:, 5);
            %
            Pmean = AA(:, 6);
            Pstd = AA(:, 7);
            Pmax = AA(:, 8);
            %
            Fmean = AA(:, 9);
            Fstd = AA(:, 10);
            Fmax = AA(:, 11);
            %
            num_gTs = AA(:, 12);
            Fmax_whichGTs = AA(:, 13);

            maxF_new_max(k) = max(Fmax);
            ind = find(Fmax==max(Fmax));
            ind = ind(1);
            maxR_new_max(k) = Rmax(ind);
            maxP_new_max(k) = Pmax(ind);
            maxF_bestGT(k) = Fmax_whichGTs(ind);
            thr_new_max(k) = thr(ind);

            maxF_new_mean(k) = max(Fmean);
            ind = find(Fmean==max(Fmean));
            ind = ind(1);
            maxR_new_mean(k) = Rmean(ind);
            maxP_new_mean(k) = Pmean(ind);
            maxF_new_std(k) = Fstd(ind);   % std across groundtruths as threshold that gave max avgF value.
            thr_new_mean(k) = thr(ind);

            %k

        end % looping through image patches



        % justMethod - Take Mean & Std across image patches of Precision,Recall,F-measure computed in different ways.

        % (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
        justMethod.maxF_newMean_meanAccImgs(d,i) = mean(maxF_new_mean);
        justMethod.maxF_newMean_stdAccImgs(d,i) = std(maxF_new_mean);
        justMethod.maxF_newStd_meanAccImgs(d,i) = mean(maxF_new_std);
        %
        justMethod.maxR_newMean_meanAccImgs(d,i) = mean(maxR_new_mean);
        justMethod.maxR_newMean_stdAccImgs(d,i) = std(maxR_new_mean);
        justMethod.maxR_newStd_meanAccImgs(d,i) = mean(maxR_new_std);
        %
        justMethod.maxP_newMean_meanAccImgs(d,i) = mean(maxP_new_mean);
        justMethod.maxP_newMean_stdAccImgs(d,i) = std(maxP_new_mean);
        justMethod.maxP_newStd_meanAccImgs(d,i) = mean(maxP_new_std);


        % (4). New way to compute P-R-F - (Max value using best single ground truther.). 
        justMethod.maxF_newMax_meanAccImgs(d,i) = mean(maxF_new_max);
        justMethod.maxF_newMax_stdAccImgs(d,i) = std(maxF_new_max);
        %
        justMethod.maxR_newMax_meanAccImgs(d,i) = mean(maxR_new_max);
        justMethod.maxR_newMax_stdAccImgs(d,i) = std(maxR_new_max);
        %
        justMethod.maxP_newMax_meanAccImgs(d,i) = mean(maxP_new_max);
        justMethod.maxP_newMax_stdAccImgs(d,i) = std(maxP_new_max);


        % relImPix - Different ways of computing P/R/F all relative to raw image pixels performance.


        % (3). New way to compute P-R-F - (Mean & STD across different ground truthers). 
        relImPix.maxF_newMean_meanAccImgs(d,i) = mean(maxF_new_mean-maxF_new_mean_imPix);
        relImPix.maxF_newMean_stdAccImgs(d,i) = std(maxF_new_mean-maxF_new_mean_imPix);
        relImPix.maxF_newStd_meanAccImgs(d,i) = mean(maxF_new_std-maxF_new_std_imPix);
        %
        relImPix.maxR_newMean_meanAccImgs(d,i) = mean(maxR_new_mean-maxR_new_mean_imPix);
        relImPix.maxR_newMean_stdAccImgs(d,i) = std(maxR_new_mean-maxR_new_mean_imPix);
        relImPix.maxR_newStd_meanAccImgs(d,i) = mean(maxR_new_std-maxR_new_std_imPix);
        %
        relImPix.maxP_newMean_meanAccImgs(d,i) = mean(maxP_new_mean-maxP_new_mean_imPix);
        relImPix.maxP_newMean_stdAccImgs(d,i) = std(maxP_new_mean-maxP_new_mean_imPix);
        relImPix.maxP_newStd_meanAccImgs(d,i) = mean(maxP_new_std-maxP_new_std_imPix);


        % (4). New way to compute P-R-F - (Max value using best single ground truther.). 
        relImPix.maxF_newMax_meanAccImgs(d,i) = mean(maxF_new_max-maxF_new_max_imPix);
        relImPix.maxF_newMax_stdAccImgs(d,i) = std(maxF_new_max-maxF_new_max_imPix);
        %
        relImPix.maxR_newMax_meanAccImgs(d,i) = mean(maxR_new_max-maxR_new_max_imPix);
        relImPix.maxR_newMax_stdAccImgs(d,i) = std(maxR_new_max-maxR_new_max_imPix);
        %
        relImPix.maxP_newMax_meanAccImgs(d,i) = mean(maxP_new_max-maxP_new_max_imPix);
        relImPix.maxP_newMax_stdAccImgs(d,i) = std(maxP_new_max-maxP_new_max_imPix);




        % Plotting histograms of all different P/R/F calculations for each blur.
        % SOME STUFF NOT DEFINED IN HERE.
        if(0)
            disp('Plotting histograms of all different P/R/F calculations for each blur.')
            F_hist_bins = linspace(0,1,20);

            N_maxF = hist(maxF,F_hist_bins);
            N_maxF_old = hist(maxF_old,F_hist_bins);

            N_maxF_new_max =  hist(maxF_new_max,F_hist_bins);
            N_maxF_new_mean =  hist(maxF_new_mean,F_hist_bins);

            N_maxF_unionGT =  hist(maxF_unionGT,F_hist_bins);



            figure, hold on,
            plot(F_hist_bins, N_maxF,'k')
            plot(F_hist_bins, N_maxF_old,'r--')

            plot(F_hist_bins, N_maxF_new_max,'b')
            plot(F_hist_bins, N_maxF_new_mean,'g')

            plot(F_hist_bins, N_maxF_unionGT,'c')

            xlabel('F-measure')
            ylabel('counts')
            title(['Blurring \sigma = ',sig{i-1}])

            legend('bench','bench2','best gT','mean across gT','union gT')
        end




%         % Note: evalRes = [bestT,bestR,bestP,bestF,R_max,P_max,F_max,Area_PR]
%         disp(['sig:',sig{i-1},' sz:',sz{i-1}])
%         disp(['      Same TH -- R:',num2str(evalRes(2),2),' P:',num2str(evalRes(3),2),' F:',num2str(evalRes(4),2)])
%         disp(['       Any TH -- R:',num2str(evalRes(5),2),' P:',num2str(evalRes(6),2),' F:',num2str(evalRes(7),2)])


    end % Looping over different Gaussian Blur sigma values (i)


end % Looping over cP_distances (d)





%% Plot mean & std of different P/R/F measures for justMethod.
H1=figure; hold on
for d = 1:num_cPdist
    errorbar([1:numel(sig)+1]+0.1*(d-1), justMethod.maxF_newMax_meanAccImgs(d,:), justMethod.maxF_newMax_stdAccImgs(d,:),[cPd_colors{d},'s--'],'LineWidth',2)
end
legend('d_t=0', 'd_t=1','d_t=1.4', 'd_t=2')% 'maxGT (d=0)','meanGT (d=0)','maxGT (d=1)','meanGT (d=1)','maxGT (d=1.4)','meanGT (d=1.4)','maxGT (d=2)','meanGT (d=2)'
for d = 1:num_cPdist
    errorbar([1:numel(sig)+1]+0.1*(d-1), justMethod.maxF_newMean_meanAccImgs(d,:), justMethod.maxF_newMean_stdAccImgs(d,:),[cPd_colors{d},'s-'],'LineWidth',2)
end
set(gca,'XTick',1:numel(sig)+1,'XTickLabel',{'0',sig{:}},'FontSize',16,'FontWeight','Bold')
grid on
xlabel('\sigma Gaussian RF Kernel','FontSize',18,'FontWeight','Bold')
ylabel('Absolute F-measure','FontSize',18,'FontWeight','Bold')
title('Fmeasure of Gaussian RFs','FontSize',20,'FontWeight','Bold')


saveGoodImg(H1,[dirPre,'../Documentation/Cosyne_2016/GaussianRF_1.jpg'],sizeGoodIm)
plot2svg([dirPre,'../Documentation/Cosyne_2016/GaussianRF_1.svg'],H1)
close(H1)


%% Plot mean & std of different P/R/F measures for relImPix.
H2=figure; hold on

for d = 1:num_cPdist
    errorbar([1:numel(sig)+1]+0.1*(d-1), relImPix.maxF_newMax_meanAccImgs(d,:), justMethod.maxF_newMax_stdAccImgs(d,:),[cPd_colors{d},'s--'],'LineWidth',2)
end
legend('d_t=0', 'd_t=1','d_t=1.4', 'd_t=2')% 'maxGT (d=0)','meanGT (d=0)','maxGT (d=1)','meanGT (d=1)','maxGT (d=1.4)','meanGT (d=1.4)','maxGT (d=2)','meanGT (d=2)'
for d = 1:num_cPdist
    errorbar([1:numel(sig)+1]+0.1*(d-1), relImPix.maxF_newMean_meanAccImgs(d,:), justMethod.maxF_newMean_stdAccImgs(d,:),[cPd_colors{d},'s-'],'LineWidth',2)
end
plot([0 numel(sig)+2], [0 0],'k--', 'LineWidth',1.5)
set(gca,'XTick',1:numel(sig)+1,'XTickLabel',{'0',sig{:}},'FontSize',16,'FontWeight','Bold')
grid on
xlabel('\sigma Gaussian RF Kernel','FontSize',18,'FontWeight','Bold')
ylabel('\Delta F_i','FontSize',18,'FontWeight','Bold')
title('Fmeasure of Gaussian RFs','FontSize',20,'FontWeight','Bold')


saveGoodImg(H2,[dirPre,'../Documentation/Cosyne_2016/GaussianRF_2.jpg'],sizeGoodIm)
plot2svg([dirPre,'../Documentation/Cosyne_2016/GaussianRF_2.svg'],H2)
close(H2)







% Make up my own legend here...
if(0)
    figure(H);
    hold on
    %
    text(0.87,1,{'\color{red}{\sigma=0.5}','\color{blue}{\sigma=1}','\color{yellow}{\sigma=1.5}','\color{cyan}{\sigma=2}','\color{black}{\sigma=2.5}','\color{magenta}{\sigma=5}','\color{green}{\sigma=10}'},'VerticalAlignment','top','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')

    if(plot_max)
        % open vs. filled
        scatter(0.65,0.76,200,'ko','LineWidth',2)
        text(0.67,0.76,'= const THs','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
        %
        scatter(0.65,0.74,200,'ko','filled','LineWidth',2)
        text(0.67,0.74,'= vary THs','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
    else
        scatter(0.65,0.70,200,'ko','LineWidth',2)
        text(0.67,0.70,'= max F','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
    end

    %
    plot([0.63,0.67],[0.65,0.65],'g-','LineWidth',1.5)
    text(0.67,0.65,' \color{green}{= Iso F}','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')
    %
    plot([0.63,0.87],[0.60,0.60],'ko--','LineWidth',1.5)
    text(0.67,0.60,'= ImPix','VerticalAlignment','middle','HorizontalAlignment','left','FontSize',16,'FontWeight','Bold')

    if(plot_max)
        max_tag = 'best&max';
    else
        max_tag = 'curve';
    end


    title(['Gaussian Blurred Image. ',max_tag],'FontSize',18,'FontWeight','Bold')

    hold off


    saveGoodImg(H,[dirPre,'../Documentation/Cosyne_2016/Overall_PR_curve_results_thinpbOFF/GaussianBlur_',max_tag,'_Kur_allParams.jpg'],sizeGoodIm)
end