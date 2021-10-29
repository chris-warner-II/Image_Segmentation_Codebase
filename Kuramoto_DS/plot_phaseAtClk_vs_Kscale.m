function plot_phaseAtClk_vs_Kscale(netMethod, img_ptch, fileSpec)

[dirPre,sizeGoodIm] = onCluster;

if ~exist([dirPre,'Kuramoto_Phase_Vs_Kscale_Single_Image_Experiment_Plots/KurPhase_vs_SimTime_vary_Kscale_',img_ptch,'_',netMethod,'_',fileSpec,'.jpg'],'file')

    % Script to look at Kuramoto Phases when using Isotropic Diffusion - uninformed by image. 

    fileType = 'BSDS_patch';
    fileSize = '51x51_ds1';

    %netMethod = 'IsoDiff'; % ,'AAnrm','Mod_N&G','Mod_SKHAdj'}; {,'Mod_SKHAdj'

    matKurFilesDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/Kur_PIF_Fourier1/',netMethod,'/'];
    matEigFilesDir = [dirPre,'output/Kuramoto/NetsFromImgs/',fileType,'_',fileSize,'/data/spectral/',netMethod,'/'];


    phase_time_vector = [1:19];

    %netMethod = 'IsoDiff';
    %img_ptch = '10081_ptch1';
    %fileSpec = 'rM1_';

    files = dir([matKurFilesDir,'*',img_ptch,'*',fileSpec,'*.mat']); % different Kscale values

    LogPhaseChange = zeros(numel(files),18);



    h=figure;





    for k = 1:numel(files)


        nf = strfind(files(k).name,'_NF');
        params_deets = files(k).name(7:nf-1);


        kur = load([matKurFilesDir,files(k).name]);
        eig = load([matEigFilesDir,'Evecs_',params_deets,'.mat']);


        % How much is the phase of oscillators changing with simulation time?
        x = diff(kur.metaCluster.phaseAtClk,1,2);
        LogPhaseChange(k,:) = log(mean(abs(x)));




       % Compute mean pairwise Area under ROC Curve to put on figure.
       for j = 1:numel(kur.netParams.gT)
            kur1D(k,j) = mean( kur.AUC_ROC_1D.kur{j}(kur.AUC_ROC_1D.kur{j}>0) );
            kurUB(k,j) = mean( kur.AUC_ROC_GS.kur{j}(kur.AUC_ROC_GS.kur{j}>0) );
       end



       for j = 1:numel(phase_time_vector)
           figure(h), subplot(numel(files)+2,numel(phase_time_vector),numel(phase_time_vector)*(k-1)+j), 
           imagesc( visKurPhase_inHSV( kur.netParams.im, reshape(kur.metaCluster.phaseAtClk(:,j),kur.netParams.Ndims) ) ), 
           set(gca,'XTick',[],'YTick',[]), axis square, colormap('hsv'), freezeColors
           if(k==1); title(['(t=',num2str(phase_time_vector(j)),')']); end
           if(j==1); ylabel(['(Ks=',num2str(kur.kurParams.Kscale,2),')']); end
       end



        for j = 1:numel(eig.netParams.gT)
            im(j) = mean( eig.AUC_ROC_1D.im{j}(eig.AUC_ROC_1D.im{j}>0) );
        end

    end    




    % Plot Phase Change vs. Kuramoto Simulation Time-Step for Different Kscales.
    subplot(numel(files)+2,1,numel(files)+1),
    plot(LogPhaseChange')
    xlabel('Simulation Time Step')
    ylabel('Log \Delta in Phase')
    legend({'1','2','3'})



    subplot(numel(files)+2,6,6*(numel(files)+1)+1), 
    imagesc(eig.netParams.im), set(gca,'XTick',[],'YTick',[]), 
    axis square, colormap('bone'),freezeColors
    ylabel(['Img'])
    xlabel({ 'im pix', 'ks=1', 'ks=2', 'ks=3'})
    %
    for j = 1:numel(eig.netParams.gT)
        subplot(numel(files)+2,(numel(eig.netParams.gT)+1),(numel(eig.netParams.gT)+1)*(numel(files)+1)+j+1), 
        imagesc(eig.netParams.gT{j}), set(gca,'XTick',[],'YTick',[]), 
        axis square, colormap('jet'), freezeColors
        ylabel(['gT',num2str(j)])
        xlabel({num2str(im(j),2),num2str(kur1D(:,j),2)})
    end


    annotation('textbox', [0 0.9 1 0.1],'String',[netMethod,' : ',params_deets,' : Kuramoto Phase vs. Simulation Time (varying Kscale)'], ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',16,'FontWeight','Bold')


    saveGoodImg(h,[dirPre,'Kuramoto_Phase_Vs_Kscale_Single_Image_Experiment_Plots/KurPhase_vs_SimTime_vary_Kscale_',img_ptch,'_',netMethod,'_',fileSpec],sizeGoodIm)
    close(h)

end
