% This script will take in the different metaSummary_Kur mat files saved
% from explore_Separation_vs_Parameters.  We ran that code for different
% parameter settings (rM, sD, sP, sW, Kscale) for different values of
% Tscale. Note: Tscale is the initial spread of theta values - phase of
% oscillators at initialization.  We saved the results into 3 different
% directories sorted by Tscale value.  This script will input the
% metaCluster_Kur mat file and directly compare the segmentation results
% for a set of parameters at the 3 different Tscale values.

dirPre = '/Users/world7one/Desktop/Grad_School/Berkeley/Work/Fritz_Work/Projects/output/Kuramoto/NetsFromImgs/GradientBox_21x21_ds1/data/Kur_PIF_Fourier1/Mod SKHAdj/';

Ts1 = load([dirPre,'Kur_metaSummary_GradientBox__tscale1_72files.mat']);
Ts0p5 = load([dirPre,'Kur_metaSummary_GradientBox__tscale0p5_72files.mat']);
Ts0p1 = load([dirPre,'Kur_metaSummary_GradientBox__tscale0p1_72files.mat']);


rM = [1,2,3,4,inf];
sD = 1; % Do below inside for loops.
sW = [0, 0.6];
sP = [0.1, 0.2, 0.3, 0.4];
kS = [300]; % not looking over Kscale right now.

for i = 1:numel(rM)
    
    sD = [rM(i)./4, inf];
    
    for j = 1:numel(sD)
        
        for k = 1:numel(sW)
        
            for L = 1:numel(sP)
            
                for M = 1:numel(kS)
            
                    % Find file for each tscale setting that matches other parameter settings
                    Ts1_ind = find( Ts1.Rmax==rM(i) & Ts1.sigD==sD(j) & Ts1.sigW==sW(k) & Ts1.sigP==sP(L) & Ts1.Kscale==kS(M) );
                    Ts0p1_ind = find( Ts0p1.Rmax==rM(i) & Ts0p1.sigD==sD(j) & Ts0p1.sigW==sW(k) & Ts0p1.sigP==sP(L) & Ts0p1.Kscale==kS(M) );
                    Ts0p5_ind = find( Ts0p5.Rmax==rM(i) & Ts0p5.sigD==sD(j) & Ts0p5.sigW==sW(k) & Ts0p5.sigP==sP(L) & Ts0p5.Kscale==kS(M) );
                    
                    % Look more closely at individual files (not the metaCluster Summary file) to compare trajectory of
                    % oscillators (if they are same for different tscale values)
                    rM_str = num2str(rM(i));
                    sD_str = num2str(sD(j));
                    sD_str(sD_str=='.')='p';
                    sW_str = num2str(sW(k));
                    sW_str(sW_str=='.')='p';
                    sP_str = num2str(sP(L));
                    sP_str(sP_str=='.')='p';
                    kS_str = num2str(kS(M));
                    
                    figure,
                    
                    try
                        %Ts1b = 
                        load([dirPre,'KurMC_GradientBox_0101_rM',rM_str,'_sD',sD_str,'_sP',sP_str,'_NF_60_',sW_str,'_kscale',kS_str,'_tscale1_runs1.mat']);
                        subplot(2,3,1), plot_phaseAtClk(metaCluster,MC,netParams,kurParams,1,1)
                        title('Tscale = 1')
                        subplot(2,3,4), plot_DivMarg(MC,kurParams,kurflags,netflags)
                        title(['AUC = ',num2str( Ts1.AUC(Ts1_ind), 3 )])
                        %plot_DivMarg_Distrib
                    catch
                        
                    end
                    try
                        %Ts0p5b = 
                        load([dirPre,'KurMC_GradientBox_0101_rM',rM_str,'_sD',sD_str,'_sP',sP_str,'_NF_60_',sW_str,'_kscale',kS_str,'_tscale0p5_runs1.mat']);
                        subplot(2,3,2), plot_phaseAtClk(metaCluster,MC,netParams,kurParams,1,1)
                        title('Tscale = 0.5')
                        subplot(2,3,5), plot_DivMarg(MC,kurParams,kurflags,netflags)
                        title(['AUC = ',num2str( Ts0p5.AUC(Ts0p5_ind), 3 )])
                        %plot_DivMarg_Distrib
                    catch
                        
                    end
                    try
                        %Ts0p1b = 
                        load([dirPre,'KurMC_GradientBox_0101_rM',rM_str,'_sD',sD_str,'_sP',sP_str,'_NF_60_',sW_str,'_kscale',kS_str,'_tscale0p1_runs1.mat']);
                        subplot(2,3,3), plot_phaseAtClk(metaCluster,MC,netParams,kurParams,1,1)
                        title('Tscale = 0.1')
                        subplot(2,3,6), plot_DivMarg(MC,kurParams,kurflags,netflags)
                        title(['AUC = ',num2str( Ts0p1.AUC(Ts0p1_ind), 3 )])
                        %plot_DivMarg_Distrib
                    catch
                        
                    end
                    
                    keyboard
                    
                end
            end
        end
    end
end
        