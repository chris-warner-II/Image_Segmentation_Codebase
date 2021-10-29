% function makeMovie_KurPhaseRelax()

%
% Function to make a movie of phase evolution of oscillators
%
% Movie making works fine on my machine.  On the cluster, however, it just
% produces a black movie - no content.  Not sure why.  Look into it.
% 



% THIS IS RUNNING LIKE A SCRIPT. LOAD THE APPROPRIATE MAT FILE. AND RUN
% THIS.


[dirPre,sizeGoodIm] = onCluster;

% A Picture of a colorwheel I will use as a colorbar for oscillator phase
colorwheel = imread([dirPre,'images/HSV_colorwheel.jpeg']);


method = 'Mod_SKHAdj'; 
rM = '5';
ks = {'sml','mid','lrg'};

gTdir = [dirPre,'images/BSDS_patch/101x101_ds1/groundTruth/'];

matDirIn = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1_blur_sig1/data/Kur_PIF_Fourier1/',method,'/'];

imgDirOut = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1_blur_sig1/imgs/Kur_PIF_Fourier1/',method,'/'];
%
if ~exist(imgDirOut,'dir')
   mkdir(imgDirOut) 
end


files = dir([matDirIn,'KurMC*_rM',num2str(rM),'_*_ks',ks{1},'.mat']);




for F = 407:numel(files)

    % load in mat file from SegmentMethod and Loop code.
    load([matDirIn,files(F).name]);
    kssml = load([matDirIn,files(F).name]);
    ksmid = load([matDirIn,files(F).name(1:end-7),ks{2},'.mat']);
    kslrg = load([matDirIn,files(F).name(1:end-7),ks{3},'.mat']);
    
    gTfile = [gTdir,kurflags.fname,'.mat'];


    % Set up and save a movie of time evolution of coupled oscillator system
    movName=[imgDirOut,'/PhaseRelax_compareKS_',netflags.imageIn,'_',netflags.fname,kurflags.KurParamsTag];
    writerObj = VideoWriter([movName,'.avi']);
    writerObj.FrameRate = 1;
    open(writerObj);

    disp(['File # ',num2str(F)])
    movName
    disp('Making a movie of oscilllator phases at 2\pi clock intervals.')




    % (1). Show image pixels first and PB

    pb_im = compute_Spatial_Gradient(netParams.im, 0);
    pb_im = pb_im./max(max(pb_im));

    pp = figure;
    set(gcf,'units','normalized','position',[0 0 1 1]) % [from left, from bottom, width, height]
    % 
    subplot(121), imagesc(netParams.im), colormap(bone), freezeColors, 
    axis square
    set(gca,'XTick',[],'YTick',[])
    title('Input Image Pixels \in (0,1)','FontSize',18,'FontWeight','Bold')
    subplot(122), imagesc(pb_im), colormap(bone), freezeColors, 
    axis square
    set(gca,'XTick',[],'YTick',[])
    title('Probabalistic Boundaries \in (0,1)','FontSize',18,'FontWeight','Bold')

    % Show Phase map color wheel
    hx2=axes('parent',pp,'position',[0.42 0.18 0.19 0.19]); % normalizedunits are used for position
    axes(hx2),imagesc(colorwheel) % this is the magnified subplot
    axis off square
    title('phase','FontSize',16,'FontWeight','Bold')
    
    % Show Ground Truth Boundariesclr
    
    hx3=axes('parent',pp,'position',[0.42 0.65 0.19 0.19]); % normalizedunits are used for position
    axes(hx3),imagesc(netParams.bD), colormap(bone), freezeColors % this is the magnified subplot
    axis off square
    title('GT','FontSize',16,'FontWeight','Bold')

    annotation('textbox', [0 0.05 1 0.1],'String',['Simulation Initialization'], ...
               'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')

    % Save Frame as part of the movie file
    mov = getframe(pp);

    writeVideo(writerObj,mov);
    close(pp)


    % (2). Loop thru time steps in simulation and show oscillator phases and PB
    phasemap_sml = zeros([netParams.Ndims,size(metaCluster.phaseAtClk,2)]);
    pb_sml = zeros([netParams.Ndims,size(metaCluster.phaseAtClk,2)]);
    %
    phasemap_mid = zeros([netParams.Ndims,size(metaCluster.phaseAtClk,2)]);
    pb_mid = zeros([netParams.Ndims,size(metaCluster.phaseAtClk,2)]);
    %
    phasemap_lrg = zeros([netParams.Ndims,size(metaCluster.phaseAtClk,2)]);
    pb_lrg = zeros([netParams.Ndims,size(metaCluster.phaseAtClk,2)]);
    %
    for i = 1:size(metaCluster.phaseAtClk,2) 

        disp([num2str(i),'/',num2str(size(metaCluster.phaseAtClk,2))])

        pp = figure;
        set(gcf,'units','normalized','position',[0 0 1 1]) % [from left, from bottom, width, height]
        %
        
        
        
        
        
        phasemap_sml(:,:,i) = reshape(kssml.metaCluster.phaseAtClk(:,i),kssml.netParams.Ndims);
        pb_sml(:,:,i) = compute_Spatial_Gradient(phasemap_sml(:,:,i), 1);
        pb_sml(:,:,i) = pb_sml(:,:,i)./max(max(pb_sml(:,:,i)));
        % 
        
        phasemap_mid(:,:,i) = reshape(ksmid.metaCluster.phaseAtClk(:,i),ksmid.netParams.Ndims);
        pb_mid(:,:,i) = compute_Spatial_Gradient(phasemap_mid(:,:,i), 1);
        pb_mid(:,:,i) = pb_mid(:,:,i)./max(max(pb_mid(:,:,i)));
        % 
        
        phasemap_lrg(:,:,i) = reshape(kslrg.metaCluster.phaseAtClk(:,i),kslrg.netParams.Ndims);
        pb_lrg(:,:,i) = compute_Spatial_Gradient(phasemap_lrg(:,:,i), 1);
        pb_lrg(:,:,i) = pb_lrg(:,:,i)./max(max(pb_lrg(:,:,i)));
        % 
        
        
        
        
        
        
        % write these to pb_png images to be used by evaluation_bdry_imageB function to compute F-measure.
        imwrite(pb_sml(:,:,i),'./temp_pb_sml.png','PNG','BitDepth',8);
        imwrite(pb_mid(:,:,i),'./temp_pb_mid.png','PNG','BitDepth',8);
        imwrite(pb_lrg(:,:,i),'./temp_pb_lrg.png','PNG','BitDepth',8);
        
        
        

        % compute R,P,F and save them in temporary text files
        nthr = 10;
        dt = 2; % in pixels
        xxx = evaluation_bdry_imageB('./temp_pb_sml.png',gTfile,'temporary_sml_prd_file',nthr,dt);
        xxx = evaluation_bdry_imageB('./temp_pb_mid.png',gTfile,'temporary_mid_prd_file',nthr,dt);
        xxx = evaluation_bdry_imageB('./temp_pb_lrg.png',gTfile,'temporary_lrg_prd_file',nthr,dt);
        
        
        
        % read in from those temporary text files saved just above (the contents is:) 
        %          1        2     3      4        5      6      7        8      9      10      11      12       13 
        vars = ['thresh', 'dt', 'mnR', 'stdR', 'Rmax', 'mnP', 'stdP', 'Pmax', 'mnF', 'stdF', 'Fmax', '#GT', 'BestGT'];
        for j = 1:4
            eval(['Asml(:,:,j) = dlmread(''./temporary_sml_pd',num2str(j),'_rd_file'');'])
            eval(['Amid(:,:,j) = dlmread(''./temporary_mid_pd',num2str(j),'_rd_file'');'])
            eval(['Alrg(:,:,j) = dlmread(''./temporary_lrg_pd',num2str(j),'_rd_file'');'])
        end
        
        
        
        % Find best Fmax & Fmean for each {ks,dt} combination
        for j = 1:4
            best_Fmax_kssml(j) = max(Asml(:,11,j));
            best_Fmean_kssml(j) = max(Asml(:,9,j));
            %
            best_Fmax_ksmid(j) = max(Amid(:,11,j));
            best_Fmean_ksmid(j) = max(Amid(:,9,j));
            %
            best_Fmax_kslrg(j) = max(Alrg(:,11,j));
            best_Fmean_kslrg(j) = max(Alrg(:,9,j));
        end
        

        
        
        subplot(321), imagesc(phasemap_sml(:,:,i)), caxis([0 2*pi]), colormap(hsv), freezeColors, 
        axis square
        set(gca,'XTick',[],'YTick',[])
        title('Model Phase Map','FontSize',18,'FontWeight','Bold')
        ylabel('ks small','FontSize',18,'FontWeight','Bold')
        %
        subplot(323), imagesc(phasemap_mid(:,:,i)), caxis([0 2*pi]), colormap(hsv), freezeColors, 
        axis square
        set(gca,'XTick',[],'YTick',[])
        ylabel('ks mid','FontSize',18,'FontWeight','Bold')
        %title('Model Phase Map','FontSize',18,'FontWeight','Bold')
        %
        subplot(325), imagesc(phasemap_lrg(:,:,i)), caxis([0 2*pi]), colormap(hsv), freezeColors, 
        axis square
        set(gca,'XTick',[],'YTick',[])
        ylabel('ks large','FontSize',18,'FontWeight','Bold')
        %title('Model Phase Map','FontSize',18,'FontWeight','Bold')
        
        
        subplot(322), imagesc(pb_sml(:,:,i)), colormap(bone), freezeColors, 
        axis square
        set(gca,'XTick',[],'YTick',[])
        title('Probabalistic Boundaries \in (0,1)','FontSize',18,'FontWeight','Bold')
        xlabel(['F(d_t) = (',num2str(best_Fmean_kssml,2),')'],'FontSize',16,'FontWeight','Bold')
        %
        subplot(324), imagesc(pb_mid(:,:,i)), colormap(bone), freezeColors, 
        axis square
        set(gca,'XTick',[],'YTick',[])
        xlabel(['F(d_t) = (',num2str(best_Fmean_ksmid,2),')'],'FontSize',16,'FontWeight','Bold')
        %title('Probabalistic Boundaries \in (0,1)','FontSize',18,'FontWeight','Bold')
        %
        subplot(326), imagesc(pb_lrg(:,:,i)), colormap(bone), freezeColors, 
        axis square
        set(gca,'XTick',[],'YTick',[])
        xlabel(['F(d_t) = (',num2str(best_Fmean_kslrg,2),')'],'FontSize',16,'FontWeight','Bold')
        %title('Probabalistic Boundaries \in (0,1)','FontSize',18,'FontWeight','Bold')
        
        
        

        
        
        
        
        
        % Show Phase map color wheel
        hx2=axes('parent',pp,'position',[0.42 0.18 0.19 0.19]); % normalizedunits are used for position
        axes(hx2),imagesc(colorwheel) % this is the magnified subplot
        axis off square
        title('phase','FontSize',14,'FontWeight','Bold')
        
        % Show Ground Truth Boundaries
        hx3=axes('parent',pp,'position',[0.42 0.65 0.19 0.19]); % normalizedunits are used for position
        axes(hx3),imagesc(netParams.bD), colormap(bone), freezeColors % this is the magnified subplot
        axis off square
        title('GT','FontSize',14,'FontWeight','Bold')

        annotation('textbox', [0 0.05 1 0.1],'String',['Simulation Time = ',num2str((kurParams.spp*(i-1))/kurParams.T*kurParams.Tsec,2)], ...
                   'EdgeColor', 'none', 'HorizontalAlignment', 'center','FontSize',20,'FontWeight','Bold')

        % Save Frame as part of the movie file
        mov = getframe(pp);

        writeVideo(writerObj,mov);
        close(pp)
    end

    close(writerObj);


end

