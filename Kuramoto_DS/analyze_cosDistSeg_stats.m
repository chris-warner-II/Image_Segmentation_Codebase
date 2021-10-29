% script (for now) to analyze statistics of thresholded pairwise cosine distance between oscillators

Cl = '4';
NF = '0p6';

% Directory to save output plots into...
imgsDir = ['/Users/world7one/Desktop/Grad_School/Berkeley/Work/Fritz_Work/',...
    'Projects/output/Kuramoto/HandCookedNetwork/imgs/SegStatsParamSearch/'];
if ~exist(imgsDir,'dir') 
    mkdir(imgsDir)
end


% Directory to look for data output from Kuramoto main to analyze here
dataDir = ['/Users/world7one/Desktop/Grad_School/Berkeley/Work/Fritz_Work/',...
    'Projects/output/Kuramoto/HandCookedNetwork/data/PIF_Fourier1_N1x12_C',num2str(Cl),'_rmaxInf/'];


files = dir([dataDir,'*',NF,'*']);
k=0; % a counter to keep track of number of files.

%% Below: Loop through output files saved from main_Kuramoto code and plot statistics
for i = 1:numel(files)


    % Load data file saved from main_Kuramoto_HandCookedNetwork
    if strfind(files(i).name,'.mat')
        x = load([dataDir,files(i).name]);
        k=k+1;
    else
        continue
    end
    

    % Extract some variables from the x data structure
    muW = x.runParams.muW;
    rmax = x.runParams.rmax;
    sigSg = x.runParams.sigSg;
    sigWk = x.runParams.sigWk;
    negWts = x.runParams.negWts;
    runs = x.runParams.runs;
    Tsec = x.runParams.Tsec;
    spp = x.runParams.spp;
    tau = x.runParams.tau;
    T = x.runParams.T;
    PIFlg = x.runParams.PIFlg;
    PIFparams = x.runParams.PIFparams;
    numInt = x.runParams.numInt;
    gndTruth = x.runParams.gndTruth;
    C = x.runParams.C;
    Ndims = x.runParams.Ndims;
    N = x.runParams.N;
    sigW = x.runParams.sigW;
    Weak = x.runParams.Weak;
    Strng = x.runParams.Strng;

    paramsCell{k} = ['Win=',num2str(Strng,2),' : Wout=',num2str(Weak,2)];
    sigNF_k(k) = sigW; 
    Win_k(k) = Strng;
    Wout_k(k) = Weak;


    % Statistics from occurances where Quickest Onset, Longest Duration and Lowest Error coincide.
    all3_coincide = ( x.QuickestOn.onset_t == x.LowestErr.onset_t) & ( x.QuickestOn.onset_t == x.LongestDur.onset_t);
    all3_num(:,k) = sum( all3_coincide, 2 );

    % What are mean and std of onset_t, duration, false_pos & missed conns when
    % all3 (LongestDur, QuickestOn, Lowest Error) coincide?
    for j = 1:numel(x.thresh_cosDist_seg)
        all3_onset_t_mean(j,k) = mean(x.LongestDur.onset_t(j,all3_coincide(j,:)));
        all3_dur_mean(j,k) = mean(x.LongestDur.dur(j,all3_coincide(j,:)));
        all3_false_pos_mean(j,k) = mean(x.LongestDur.false_pos(j,all3_coincide(j,:)));
        all3_missed_conns_mean(j,k) = mean(x.LongestDur.missed_conns(j,all3_coincide(j,:)));
        %
        all3_onset_t_std(j,k) = std(x.LongestDur.onset_t(j,all3_coincide(j,:)));
        all3_dur_std(j,k) = std(x.LongestDur.dur(j,all3_coincide(j,:)));
        all3_false_pos_std(j,k) = std(x.LongestDur.false_pos(j,all3_coincide(j,:)));
        all3_missed_conns_std(j,k) = std(x.LongestDur.missed_conns(j,all3_coincide(j,:)));
    end



    % Find statistics from 3 different cases individually.
    %
    % for Longest Duration of Stable Percept (Segmentation)
    LD_nan = isnan(x.LongestDur.onset_t);
    LongestDur_numNans(:,k) = sum(LD_nan');
    for j = 1:numel(x.thresh_cosDist_seg)
        indx = find(LD_nan(j,:)==0);
        LongestDur_onset_t_mean(j,k) = mean(x.LongestDur.onset_t(j,indx));
        LongestDur_dur_mean(j,k) = mean(x.LongestDur.dur(j,indx));
        LongestDur_false_pos_mean(j,k) = mean(x.LongestDur.false_pos(j,indx));
        LongestDur_missed_conns_mean(j,k) = mean(x.LongestDur.missed_conns(j,indx));
        %
        LongestDur_onset_t_std(j,k) = std(x.LongestDur.onset_t(j,indx));
        LongestDur_dur_std(j,k) = std(x.LongestDur.dur(j,indx));
        LongestDur_false_pos_std(j,k) = std(x.LongestDur.false_pos(j,indx));
        LongestDur_missed_conns_std(j,k) = std(x.LongestDur.missed_conns(j,indx));
    end
    %
    % for Quickest Onset of Stable Percept (Segmentation)
    QO_nan = isnan(x.QuickestOn.onset_t);
    QuickestOn_numNans(:,k) = sum(QO_nan');
    for j = 1:numel(x.thresh_cosDist_seg)
        indx = find(QO_nan(j,:)==0);
        QuickestOn_onset_t_mean(j,k) = mean(x.QuickestOn.onset_t(j,indx));
        QuickestOn_dur_mean(j,k) = mean(x.QuickestOn.dur(j,indx));
        QuickestOn_false_pos_mean(j,k) = mean(x.QuickestOn.false_pos(j,indx));
        QuickestOn_missed_conns_mean(j,k) = mean(x.QuickestOn.missed_conns(j,indx));
        %
        QuickestOn_onset_t_std(j,k) = std(x.QuickestOn.onset_t(j,indx));
        QuickestOn_dur_std(j,k) = std(x.QuickestOn.dur(j,indx));
        QuickestOn_false_pos_std(j,k) = std(x.QuickestOn.false_pos(j,indx));
        QuickestOn_missed_conns_std(j,k) = std(x.QuickestOn.missed_conns(j,indx));
    end
    %
    % for Lowest Error of Stable Percept (Segmentation)
    LE_nan = isnan(x.LowestErr.onset_t);
    LowestErr_numNans(:,k) = sum(LE_nan');
    for j = 1:numel(x.thresh_cosDist_seg)
        indx = find(QO_nan(j,:)==0);
        LowestErr_onset_t_mean(j,k) = mean(x.LowestErr.onset_t(j,indx));
        LowestErr_dur_mean(j,k) = mean(x.LowestErr.dur(j,indx));
        LowestErr_false_pos_mean(j,k) = mean(x.LowestErr.false_pos(j,indx));
        LowestErr_missed_conns_mean(j,k) = mean(x.LowestErr.missed_conns(j,indx));
        %
        LowestErr_onset_t_std(j,k) = std(x.LowestErr.onset_t(j,indx));
        LowestErr_dur_std(j,k) = std(x.LowestErr.dur(j,indx));
        LowestErr_false_pos_std(j,k) = std(x.LowestErr.false_pos(j,indx));
        LowestErr_missed_conns_std(j,k) = std(x.LowestErr.missed_conns(j,indx));
    end
    %


end





%% Plot


% Order Parameters So Plots will look more coherent [Hierarchical Sorting]
% 1st, NF.  2nd, Win. 3rd, Wout.
indx = []; % Note:  This indx and all ones below are different from ones above.
[V_NF,I_NF] = sort(sigNF_k);
NFs = unique(V_NF);
for a = 1:numel(NFs)
    
    ind_NF = find(V_NF == NFs(a) );
    
    [V_Win,I_Win] = sort(Win_k(ind_NF));
    WinS = unique(V_Win);
    
    for b = 1:numel(WinS)
        
        ind_Win = find(V_Win == WinS(b) );
        
        [V_Wout,I_Wout] = sort(Wout_k(ind_NF(ind_Win)));
        WoutS = unique(V_Wout);
        
        indx = [indx, (ind_NF(ind_Win(I_Wout)))];  % Figure out ordering.  Why is Win = 1 before Win = 5 ?
%         paramsCell{indx}
%         keyboard
        
        
    end
    
end

% paramsCell{indx}
% keyboard

for j = 1:numel(indx)
    paramsCell2{j} = paramsCell{indx(j)};
end



colors = colormap('jet');
num_colors = size(colors,1);
num_Wout = numel(I_Wout);
cstep = floor(num_colors/num_Wout);


% Plot Mean with Std ErrorBars these statistics
%
% 
for j = 1:numel(x.thresh_cosDist_seg)
    
    
    % When all 3 measurements (Quickest Onset, Longest Duration, Lowest Error) coincide at a time point.
    ha3 = figure; 
    subplot(611), hold on % t_onset statistics
    title(['Statistics when Stable Percept is Earliest Onset, Longest Duration & Minimum Error with  \sigma_{NF} = ',num2str(sigW,2),...
        ' # Clusters = ',num2str(C),' & Cos Dist Threshold = ',num2str(x.thresh_cosDist_seg(j),2)],'FontSize',20,'FontWeight','Bold')
    plot_Ebars_KurHCN(indx, all3_onset_t_mean(j,:), all3_onset_t_std(j,:), I_Wout, 'Onset Time');
    %
    subplot(612), hold on % duration of stable percept
    plot_Ebars_KurHCN(indx, all3_dur_mean(j,:), all3_dur_std(j,:), I_Wout, 'Duration');
    %
    subplot(613), hold on % false positives
    plot_Ebars_KurHCN(indx, all3_false_pos_mean(j,:), all3_false_pos_std(j,:), I_Wout, 'False Pos');
    %
    subplot(614), hold on % missed connections
    plot_Ebars_KurHCN(indx, all3_missed_conns_mean(j,:), all3_missed_conns_std(j,:), I_Wout, 'Missed Conns');
    %
    subplot(615), hold on % number of failures.
    bar_KurHCN(indx, all3_num(j,:), I_Wout, '# Coincidences', paramsCell2)
    %
    saveGoodImg(ha3,[imgsDir,'SegStats_all3_NF',NF,'_C',num2str(C),'_thr',num2str(x.thresh_cosDist_seg(j))],[0 0 1 1])
    close(ha3) 


    %
    %
    % Longest Duration of Stable Percept
    hld = figure; 
    subplot(611), hold on % t_onset statistics
    title(['Statistics when Stable Percept is Longest Duration with \sigma_{NF} = ',num2str(sigW,2),' # Clusters = ',num2str(C),...
        ' & Cos Dist Threshold = ',num2str(x.thresh_cosDist_seg(j),2)],'FontSize',20,'FontWeight','Bold')
    plot_Ebars_KurHCN(indx, LongestDur_onset_t_mean(j,:), LongestDur_onset_t_std(j,:), I_Wout, 'Onset Time');
    %
    subplot(612), hold on % duration of stable percept
    plot_Ebars_KurHCN(indx, LongestDur_dur_mean(j,:), LongestDur_dur_std(j,:), I_Wout, 'Duration');
    %
    subplot(613), hold on % false positives
    plot_Ebars_KurHCN(indx, LongestDur_false_pos_mean(j,:), LongestDur_false_pos_std(j,:), I_Wout, 'False Pos');
    %
    subplot(614), hold on % missed connections
    plot_Ebars_KurHCN(indx, LongestDur_missed_conns_mean(j,:), LongestDur_missed_conns_std(j,:), I_Wout, 'Missed Conns');
    %
    subplot(615), hold on % number of failures
    bar_KurHCN(indx, LongestDur_numNans(j,:), I_Wout, '# Failures', paramsCell2)
    %
    saveGoodImg(hld,[imgsDir,'SegStats_LongestDur_NF',NF,'_thr',num2str(x.thresh_cosDist_seg(j))],[0 0 1 1])
    close(hld) 
    
    
    
    %
    %
    % Quickest Onset of a Stable Percept
    hqo = figure;
    subplot(611), hold on % t_onset statistics
    title(['Statistics when Stable Percept is Earliest Onset with \sigma_{NF} = ',num2str(sigW,2),' # Clusters = ',num2str(C),...
        ' & Cos Dist Threshold = ',num2str(x.thresh_cosDist_seg(j),2)],'FontSize',20,'FontWeight','Bold')
    plot_Ebars_KurHCN(indx, QuickestOn_onset_t_mean(j,:), QuickestOn_onset_t_std(j,:), I_Wout, 'Onset Time');
    %
    subplot(612), hold on % duration of stable percept
    plot_Ebars_KurHCN(indx, QuickestOn_dur_mean(j,:), QuickestOn_dur_std(j,:), I_Wout, 'Duration');
    %
    subplot(613), hold on % false positives
    plot_Ebars_KurHCN(indx, QuickestOn_false_pos_mean(j,:), QuickestOn_false_pos_std(j,:), I_Wout, 'False Pos');
    %
    subplot(614), hold on % missed connections
    plot_Ebars_KurHCN(indx, QuickestOn_missed_conns_mean(j,:), QuickestOn_missed_conns_std(j,:), I_Wout, 'Onset Time');
    %
    subplot(615), hold on % number of failures.
    bar_KurHCN(indx, QuickestOn_numNans(j,:), I_Wout, '# Failures', paramsCell2)
    %
    saveGoodImg(hqo,[imgsDir,'SegStats_QuickestOn_NF',NF,'_thr',num2str(x.thresh_cosDist_seg(j))],[0 0 1 1])
    close(hqo) 
    
    
    
    
    
    
    %
    %
    % Lowest Error (discrepancy from ground truth) of a Stable Percept
    hle = figure;
    subplot(611), hold on % t_onset statistics
    title(['Statistics when Stable Percept is Minimum Error with \sigma_{NF} = ',num2str(sigW,2),' # Clusters = ',num2str(C),...
        ' & Cos Dist Threshold = ',num2str(x.thresh_cosDist_seg(j),2)],'FontSize',20,'FontWeight','Bold')
    plot_Ebars_KurHCN(indx, LowestErr_onset_t_mean(j,:), LowestErr_onset_t_std(j,:), I_Wout, 'Onset Time');
    %
    subplot(612), hold on % duration of stable percept
    plot_Ebars_KurHCN(indx, LowestErr_dur_mean(j,:), LowestErr_dur_std(j,:), I_Wout, 'Duration');
    %
    subplot(613), hold on % false positives
    plot_Ebars_KurHCN(indx, LowestErr_false_pos_mean(j,:), LowestErr_false_pos_std(j,:), I_Wout, 'False Pos');
    %
    subplot(614), hold on % missed connections
    plot_Ebars_KurHCN(indx, LowestErr_missed_conns_mean(j,:), LowestErr_missed_conns_std(j,:), I_Wout, 'Onset Time');
    %
    subplot(615), hold on % number of failures.
    bar_KurHCN(indx, LowestErr_numNans(j,:), I_Wout, '# Failures', paramsCell2)
    %
    saveGoodImg(hle,[imgsDir,'SegStats_LowestErr_NF',NF,'_thr',num2str(x.thresh_cosDist_seg(j))],[0 0 1 1])
    close(hle) 
    
    
    
    
    
end

% TODO:
%   Sort different files in order by Wout  !!!
%   Make 3 series of plots for 3 different thresholds
%   Make different series of plots for different number of clusters
