dirPre = onCluster;
dataDistDir = [dirPre,'output/Kuramoto/HandCookedNetwork/data/ROC_KurEig_Distil/'];

files2stitch = dir([dataDistDir,'*.mat']);


%% Preallocate Data Structure Space.
% parameters
DTot.params.Win = [];
DTot.params.Wout = [];
DTot.params.Rmax = [];
DTot.params.Pfar = [];
DTot.params.N = [];
DTot.params.C = [];
DTot.params.MuW = [];
DTot.params.SigW = [];
DTot.params.runParams = [];
DTot.params.runflags = [];
%
% Eigenvector Segmentation Stats
DTot.Eigen.AUCe = [];
DTot.Eigen.robustoE = [];
%
% Kuramoto Run-Avg'd Results
DTot.acrossRuns.AUCk_mn =[];
DTot.acrossRuns.AUCk_std =[];
DTot.acrossRuns.robustoK_mn =[];
DTot.acrossRuns.robustoK_std =[];
%
% Kuramoto Single Run Statistics
DTot.singleRuns.maxAUC_mn = [];
DTot.singleRuns.pctTimeOn_mn = [];
DTot.singleRuns.FirstTimeOn_mn = [];
DTot.singleRuns.meanOnDuration = [];
DTot.singleRuns.meanOffDuration = [];




%% Loop through mat file and save their contents into Dtot.
for i = 1:numel(files2stitch)
    
    
    load([dataDistDir,files2stitch(i).name])  

    
    DTot.params.Win = [DTot.params.Win, Distil.params.Win];
    DTot.params.Wout = [DTot.params.Wout, Distil.params.Wout];
    DTot.params.Rmax = [DTot.params.Rmax, Distil.params.Rmax];
    DTot.params.Pfar = [DTot.params.Pfar, Distil.params.Pfar];
    DTot.params.N = [DTot.params.N, Distil.params.N];
    DTot.params.C = [DTot.params.C, Distil.params.C];
    DTot.params.MuW = [DTot.params.MuW, Distil.params.MuW];
    DTot.params.SigW = [DTot.params.SigW, Distil.params.SigW];
    DTot.params.runParams = [DTot.params.runParams, Distil.params.runParams];
    DTot.params.runflags = [DTot.params.runflags, Distil.params.runflags];
    %
    % Eigenvector Segmentation Stats
    DTot.Eigen.AUCe = [DTot.Eigen.AUCe, Distil.Eigen.AUCe];
    DTot.Eigen.robustoE = [DTot.Eigen.robustoE, Distil.Eigen.robustoE];
    %
    % Kuramoto Run-Avg'd Results
    DTot.acrossRuns.AUCk_mn =[DTot.acrossRuns.AUCk_mn, Distil.acrossRuns.AUCk_mn];
    DTot.acrossRuns.AUCk_std =[DTot.acrossRuns.AUCk_std, Distil.acrossRuns.AUCk_std];
    DTot.acrossRuns.robustoK_mn =[DTot.acrossRuns.robustoK_mn, Distil.acrossRuns.robustoK_mn];
    DTot.acrossRuns.robustoK_std =[DTot.acrossRuns.robustoK_std, Distil.acrossRuns.robustoK_std];
    %
    % Kuramoto Single Run Statistics
    DTot.singleRuns.maxAUC_mn = [DTot.singleRuns.maxAUC_mn, Distil.singleRuns.maxAUC_mn];
    DTot.singleRuns.pctTimeOn_mn = [DTot.singleRuns.pctTimeOn_mn, Distil.singleRuns.pctTimeOn_mn];
    DTot.singleRuns.FirstTimeOn_mn = [DTot.singleRuns.FirstTimeOn_mn, Distil.singleRuns.FirstTimeOn_mn];
    DTot.singleRuns.meanOnDuration = [DTot.singleRuns.meanOnDuration, Distil.singleRuns.meanOnDuration];
    DTot.singleRuns.meanOffDuration = [DTot.singleRuns.meanOffDuration, Distil.singleRuns.meanOffDuration];
    
    
end




%% Save Total Distil Data Structure
save([dataDistDir,'Stitch_KurEig_ROC_Distil_1_',num2str(numel(files2stitch))],'DTot','files2stitch'); 

