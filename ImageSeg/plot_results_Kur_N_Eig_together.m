% This script takes in the output mat files from
% plot_eval_compare_bsdsBench_EigB and plot_eval_compare_bsdsBench_KurB
% that are run separately to make plots about F-measure segmentation
% performance of Spectral and Kuramoto versions of image segmentation. This
% takes those two results and plots them side by side.  There is likely a
% lot of extraneous stuff in this code.

[dirPre,sizeGoodIm] = onCluster;

do_maxGT = 1;
do_meanGT = 1;

method_colorE = {'g','y','r','b','k'};
method_colorK = {'k','y','g','r','b','c'};
E_to_K_translate = [3, 2, 5, 6, 1]

which_errbars = 'sem';
blur_tag_M = '_blur_sig1';

% E = load(['Eig_plot_results_sem']); % HOW IT WAS. REPLACE IF RESULTS DISAPPEAR.
dirr = './output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1_blur_sig1/data/EigVKur/';
E = load([dirr,'Eig_plot_results', blur_tag_M, '_', which_errbars]);
K = load([dirr,'Kur_plot_results', blur_tag_M, '_', which_errbars]);

outImgDir = dirPre;










%% Plot F-measure of bestEig and bestKur side by-side in subplots with
%           (meanGT & maxGT) for different values of d_t.
Hoo = figure; % figure to hold best Eig and best Kur performances.
subplot(121), hold on


switch E.which_errbars
    case 'sem'
        maxY = 0.45;
        minY = 0;
    case 'std'
        maxY = 0.55;
        minY = 0;
end

which_eig_combo = zeros(1,numel(E.method_short));
for m = 1:numel(E.method_short)
    which_eig_combo(m) =  find( mean(squeeze(E.justMethod.maxF_newMax_meanAccImgs(4,:,:,m)),2) == max(mean(squeeze(E.justMethod.maxF_newMax_meanAccImgs(4,:,:,m)),2)) );
    best_legend{m} = [E.method_short{m},' (',E.ev{which_eig_combo(m)},')'];
end

% subplot(ha(4)), hold on
if(do_maxGT)
    % performance using just imPix
    errorbar([1:4]+0.1*(0), E.justImPix.maxF_newMax_meanAccImgs, E.justImPix.maxF_newMax_semAccImgs, 'LineWidth',3,'LineStyle','-','Color','magenta')
    %
    % performance using just imBlur
    errorbar([1:4]+0.1*(0), E.justImBlur.maxF_newMax_meanAccImgs, E.justImBlur.maxF_newMax_semAccImgs, 'LineWidth',3,'LineStyle','-','Color','cyan')
    %
    % Show maxGT for methods, imPix and imBlur
    for m = 1:numel(E.method_short) % {GLnrm AAnrm Mod_SKHAdj Mod_N&G IsoDiff}
        errorbar([1:4]+0.1*(1-1), squeeze(E.justMethod.maxF_newMax_meanAccImgs(4,which_eig_combo(m),:,m)), squeeze(E.justMethod.maxF_newMax_semAccImgs(4,which_eig_combo(m),:,m)), 'LineWidth',3,'LineStyle','-','Color',method_colorE{m})
    end
end


% if show meanGT in addition to maxGT results.
if(do_meanGT)
    % performance using just imPix
    errorbar([1:4]+0.1*(0), E.justImPix.maxF_newMean_meanAccImgs, E.justImPix.maxF_newMean_semAccImgs, 'LineWidth',3,'LineStyle','--','Color','magenta')
    %
    % performance using just imBlur
    errorbar([1:4]+0.1*(0), E.justImBlur.maxF_newMean_meanAccImgs, E.justImBlur.maxF_newMean_semAccImgs, 'LineWidth',3,'LineStyle','--','Color','cyan')
    %
    for m = 1:numel(E.method_short) % {GLnrm AAnrm Mod_SKHAdj Mod_N&G IsoDiff}
        errorbar([1:4]+0.1*(1-1), squeeze(E.justMethod.maxF_newMean_meanAccImgs(4,which_eig_combo(m),:,m)), squeeze(E.justMethod.maxF_newMean_semAccImgs(4,which_eig_combo(m),:,m)), 'LineWidth',3,'LineStyle','--','Color',method_colorE{m})
    end
    %
end

legend('Pix','Blur',best_legend{:},'Location','SouthEast')

xlabel('d_t','FontSize',18,'FontWeight','Bold')
ylabel('F','FontSize',18,'FontWeight','Bold')

title(['Best Eig'],'FontSize',18,'FontWeight','Bold')
ylim([minY maxY])
grid on
set(gca,'XTick',1:4,'XTickLabel',{'0','1','1.4','2'},'YTick',0:0.05:0.45,'FontSize',16,'FontWeight','Bold')

% % % % % % %

subplot(122), hold on

    
which_F_computation = {'meanGT','meanGT','meanGT','maxGT','maxGT','maxGT'};
relative_to_what = {'justMethod','relImPix','relImBlur','justMethod','relImPix','relImBlur'};


for i = 1:1 % half the size. trying to plot maxGT & meanGT on same plot.
         % numel(which_F_computation)
              % relative_to_what


    % NOTE: TO DO: IF IT IS JUST METHOD, PLOT IMPIX (MAGENTA) & IMBLUR (CYAN) TOO FOR COMPARISON.
    switch relative_to_what{i}
        case 'justMethod'
            errorbar( [1:K.num_cPdist] +0.01*(-1),   K.justImPix.maxF_newMax_meanAccImgs, K.justImPix.maxF_newMax_semAccImgs,['m-'], 'LineWidth', 4)
            errorbar( [1:K.num_cPdist] +0.01*(-1),   K.justImPix.maxF_newMean_meanAccImgs, K.justImPix.maxF_newMean_semAccImgs,['m--'], 'LineWidth', 4)
            %
            errorbar( [1:K.num_cPdist] +0.01*(-2),   K.justImBlur.maxF_newMax_meanAccImgs, K.justImBlur.maxF_newMax_semAccImgs,['c-'], 'LineWidth', 4)
            errorbar( [1:K.num_cPdist] +0.01*(-2),   K.justImBlur.maxF_newMean_meanAccImgs, K.justImBlur.maxF_newMean_semAccImgs,['c--'], 'LineWidth', 4)
        case 'relImPix'
            errorbar( [1:K.num_cPdist] +0.01*(-1), [0 0 0 0], [0 0 0 0], ['m-'], 'LineWidth', 4)
            errorbar( [1:K.num_cPdist] +0.01*(-1), [0 0 0 0], [0 0 0 0], ['m--'], 'LineWidth', 4)
            %
            errorbar( [1:K.num_cPdist] +0.01*(-2), K.imPixVsImBlur.maxF_newMax_meanAccImgs, K.imPixVsImBlur.maxF_newMax_semAccImgs, ['c-'], 'LineWidth', 4)
            errorbar( [1:K.num_cPdist] +0.01*(-2), K.imPixVsImBlur.maxF_newMean_meanAccImgs, K.imPixVsImBlur.maxF_newMean_semAccImgs, ['c--'], 'LineWidth', 4)

        case 'relImBlur'
            errorbar( [1:K.num_cPdist] +0.01*(-1), -K.imPixVsImBlur.maxF_newMax_meanAccImgs, K.imPixVsImBlur.maxF_newMax_semAccImgs, ['m-'], 'LineWidth', 4)
            errorbar( [1:K.num_cPdist] +0.01*(-1), -K.imPixVsImBlur.maxF_newMean_meanAccImgs, K.imPixVsImBlur.maxF_newMean_semAccImgs, ['m--'], 'LineWidth', 4)
            %
            errorbar( [1:K.num_cPdist] +0.01*(-2), [0 0 0 0], [0 0 0 0], ['c-'], 'LineWidth', 4)
            errorbar( [1:K.num_cPdist] +0.01*(-2), [0 0 0 0], [0 0 0 0], ['c--'], 'LineWidth', 4)
    end




    for A = 1:numel(K.method)
        if(do_maxGT)
            errorbar( [1:K.num_cPdist] +0.01*A,   K.maxF_meanAccImgs(:,3+i,A), K.maxF_semAccImgs(:,3+i,A), ... % 'maxGT'
                [method_colorK{A},'-'], 'LineWidth', 3)
        end

        if(do_meanGT)
            errorbar( [1:K.num_cPdist] +0.01*A,   K.maxF_meanAccImgs(:,i,A), K.maxF_semAccImgs(:,i,A), ...     % 'meanGT'
                [method_colorK{A},'--'], 'LineWidth', 3)
        end
    end

    xlabel(['d_t'],'FontSize',18,'FontWeight','Bold')
    ylabel(['\Delta F  '],'FontSize',18,'FontWeight','Bold')
    set(gca,'XTick',[1:K.num_cPdist],'XTickLabel',[0 1 1.4 2],'FontSize',16,'FontWeight','Bold')

    % plot( [0.5 K.num_cPdist+0.5], [0 0],'k.-.','LineWidth',2)
    grid on

    % Hardcode y-axis limits.
    switch i
        case 1
            if ~isempty(strmatch(K.which_errbars,'sem'))
                ylim([0 0.45]); % justMethod w/ sem.
            else
                ylim([0 0.60]); % justMethod w/ std.
            end
        case 2
            if ~isempty(strmatch(K.which_errbars,'sem'))
                ylim([-0.01 0.08]); % relImPix w/ sem. 
            else
                ylim([-0.10 0.18]); % relImPix w/ std.
            end

        case 3
            if ~isempty(strmatch(K.which_errbars,'sem'))
                ylim([-0.05 0.05]); % relImBlur w/ sem. 
            else
                ylim([-0.11 0.11]); % relImBlur w/ std.
            end
    end


end

title(['Best Kur'],'FontSize',18,'FontWeight','Bold')
    

% 
disp([outImgDir,'compare_BestEig_N_BestKur_errbars_',which_errbars,blur_tag_M,'.svg'])
% SOMETHING OUTDATED ABOUT PLOT2SVG. MAYBE FIX IT UP?
saveas(Hoo,[outImgDir,'compare_BestEig_N_BestKur_errbars_',which_errbars,blur_tag_M,'.svg'],'svg')
%
saveGoodImg(Hoo,[outImgDir,'compare_BestEig_N_BestKur_errbars_',which_errbars,blur_tag_M,'.jpg'],sizeGoodIm)
close(Hoo)
% 
 













%% Plot F-measure of bestEig vs. bestKur (meanGT & maxGT) for different values of d_t.
Hii = figure; % figure to hold best Eig and best Kur performances.
hold on

i=1 %trying to plot maxGT & meanGT on same plot.
switch E.which_errbars
    case 'sem'
        maxY = 0.45;
        minY = 0;
    case 'std'
        maxY = 0.55;
        minY = 0;
end
%
% %
%
if(do_maxGT)
    % performance using just imPix
    herrorbar(E.justImPix.maxF_newMax_meanAccImgs, K.justImPix.maxF_newMax_meanAccImgs, ...
        E.justImPix.maxF_newMax_semAccImgs, 'm.')%,'LineWidth',3,'LineStyle','none')
    errorbar(E.justImPix.maxF_newMax_meanAccImgs, K.justImPix.maxF_newMax_meanAccImgs, ...
        K.justImPix.maxF_newMax_semAccImgs, 'm.')%, 'LineWidth',3,'LineStyle','none','Color','magenta')
    scatter(E.justImPix.maxF_newMax_meanAccImgs, K.justImPix.maxF_newMax_meanAccImgs, ...
        70, 'mo','LineWidth',3)
    %
    % performance using just imBlur
    herrorbar(E.justImBlur.maxF_newMax_meanAccImgs, K.justImBlur.maxF_newMax_meanAccImgs, ...
        E.justImBlur.maxF_newMax_semAccImgs, 'c.')%, 'LineWidth',3,'LineStyle','none')
    errorbar(E.justImBlur.maxF_newMax_meanAccImgs, K.justImBlur.maxF_newMax_meanAccImgs, ...
        K.justImBlur.maxF_newMax_semAccImgs,'c.') %, 'LineWidth',3,'LineStyle','none','Color','cyan')
    scatter(E.justImBlur.maxF_newMax_meanAccImgs, K.justImBlur.maxF_newMax_meanAccImgs, ...
        70, 'mo','LineWidth',3)
    %
    % Show maxGT for methods, imPix and imBlur
    for m = 1:numel(E.method_short) % {GLnrm AAnrm Mod_SKHAdj Mod_N&G IsoDiff}
        herrorbar(squeeze(E.justMethod.maxF_newMax_meanAccImgs(4,E.which_eig_combo(m),:,m)), K.maxF_meanAccImgs(:,3+i,E_to_K_translate(m)), ...
            squeeze(E.justMethod.maxF_newMax_semAccImgs(4,E.which_eig_combo(m),:,m)), [method_colorE{m},'.'])% ,'LineWidth',3,'LineStyle','none')
        errorbar(squeeze(E.justMethod.maxF_newMax_meanAccImgs(4,E.which_eig_combo(m),:,m)), K.maxF_meanAccImgs(:,3+i,E_to_K_translate(m)), ...
            K.maxF_semAccImgs(:,3+i,E_to_K_translate(m)), [method_colorE{m},'.']) % 'LineWidth',3,'LineStyle','none','Color',method_colorE{m})
        scatter(squeeze(E.justMethod.maxF_newMax_meanAccImgs(4,E.which_eig_combo(m),:,m)), K.maxF_meanAccImgs(:,3+i,E_to_K_translate(m)), ...
            70, [method_colorE{m},'o'],'LineWidth',3)
    end
end
%
% %
%
if(do_meanGT)
    % performance using just imPix
    herrorbar(E.justImPix.maxF_newMean_meanAccImgs, K.justImPix.maxF_newMean_meanAccImgs, ...
        E.justImPix.maxF_newMean_semAccImgs, 'm.')%,'LineWidth',3,'LineStyle','none')
    errorbar(E.justImPix.maxF_newMean_meanAccImgs, K.justImPix.maxF_newMean_meanAccImgs, ...
        K.justImPix.maxF_newMean_semAccImgs, 'm.')%, 'LineWidth',3,'LineStyle','none','Color','magenta')
    scatter(E.justImPix.maxF_newMean_meanAccImgs, K.justImPix.maxF_newMean_meanAccImgs, ...
        70, 'mx','LineWidth',3)
    %
    % performance using just imBlur
    herrorbar(E.justImBlur.maxF_newMean_meanAccImgs, K.justImBlur.maxF_newMean_meanAccImgs, ...
        E.justImBlur.maxF_newMax_semAccImgs, 'c.')%, 'LineWidth',3,'LineStyle','none')
    errorbar(E.justImBlur.maxF_newMean_meanAccImgs, K.justImBlur.maxF_newMean_meanAccImgs, ...
        K.justImBlur.maxF_newMean_semAccImgs,'c.') %, 'LineWidth',3,'LineStyle','none','Color','cyan')
    scatter(E.justImBlur.maxF_newMean_meanAccImgs, K.justImBlur.maxF_newMean_meanAccImgs, ...
        70, 'cx','LineWidth',3)
    %
    % Show maxGT for methods, imPix and imBlur
    for m = 1:numel(E.method_short) % {GLnrm AAnrm Mod_SKHAdj Mod_N&G IsoDiff}
        herrorbar(squeeze(E.justMethod.maxF_newMean_meanAccImgs(4,E.which_eig_combo(m),:,m)), K.maxF_meanAccImgs(:,i,E_to_K_translate(m)), ...
            squeeze(E.justMethod.maxF_newMax_semAccImgs(4,E.which_eig_combo(m),:,m)), [method_colorE{m},'x'])% ,'LineWidth',3,'LineStyle','none')
        errorbar(squeeze(E.justMethod.maxF_newMean_meanAccImgs(4,E.which_eig_combo(m),:,m)), K.maxF_meanAccImgs(:,i,E_to_K_translate(m)), ...
            K.maxF_semAccImgs(:,i,E_to_K_translate(m)), [method_colorE{m},'x']) % 'LineWidth',3,'LineStyle','none','Color',method_colorE{m})
        scatter(squeeze(E.justMethod.maxF_newMean_meanAccImgs(4,E.which_eig_combo(m),:,m)), K.maxF_meanAccImgs(:,i,E_to_K_translate(m)), ...
            70, [method_colorE{m},'x'],'LineWidth',3)
    end
end
%
% %
%
plot([0 0.5],[0 0.5],'k--')
xticks([0, .1, .2, .3, .4, .5])
yticks([0, .1, .2, .3, .4, .5])
ylabel('Kuramoto Net','FontSize',18,'FontWeight','Bold')
xlabel('Eigen-method','FontSize',18,'FontWeight','Bold')

title(['Segmentation F-measure'],'FontSize',20,'FontWeight','Bold')
grid on
set(gca,'FontSize',16,'FontWeight','Bold')
% 
disp([outImgDir,'compare_BestEigVKur_',which_errbars,blur_tag_M,'.svg'])
saveas(Hii,[outImgDir,'compare_BestEigVKur_',which_errbars,blur_tag_M,'.svg'],'svg')
%plot2svg([outImgDir,'compare_BestEig_N_BestKur_errbars_',which_errbars,blur_tag_M,'.svg'],Hii)
%
saveGoodImg(Hii,[outImgDir,'compare_BestEigVKur_',which_errbars,blur_tag_M,'.jpg'],[0 0 .7 1])
close(Hii)
