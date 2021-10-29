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
ks = 'lrg';

% NOTE: Here, add the ability to look at all 3 ks values side by side in a single movie.

matDirIn = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1_blur_sig1/data/Kur_PIF_Fourier1/',method,'/'];

imgDirOut = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1_blur_sig1/imgs/Kur_PIF_Fourier1/',method,'/'];
%
if ~exist(imgDirOut,'dir')
   mkdir(imgDirOut) 
end


files = dir([matDirIn,'KurMC*_rM',num2str(rM),'_*_ks',ks,'.mat']);




for F = 200:numel(files)

    % load in mat file from SegmentMethod and Loop code.
    load([matDirIn,files(F).name]);


    % Set up and save a movie of time evolution of coupled oscillator system
    movName=[imgDirOut,'/PhaseRelax_',netflags.imageIn,'_',netflags.fname,kurflags.KurParamsTag];
    writerObj = VideoWriter([movName,'.avi']);
    writerObj.FrameRate = 1;
    open(writerObj);

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
    phasemap = zeros([netParams.Ndims,size(metaCluster.phaseAtClk,2)]);
    pb = zeros([netParams.Ndims,size(metaCluster.phaseAtClk,2)]);
    %
    for i = 1:size(metaCluster.phaseAtClk,2) 

        disp([num2str(i),'/',num2str(size(metaCluster.phaseAtClk,2))])

        pp = figure;
        set(gcf,'units','normalized','position',[0 0 1 1]) % [from left, from bottom, width, height]
        %
        phasemap(:,:,i) = reshape(metaCluster.phaseAtClk(:,i),netParams.Ndims);
        pb(:,:,i) = compute_Spatial_Gradient(phasemap(:,:,i), 1);
        pb(:,:,i) = pb(:,:,i)./max(max(pb(:,:,i)));
        % 
        subplot(121), imagesc(phasemap(:,:,i)), caxis([0 2*pi]), colormap(hsv), freezeColors, 
        axis square
        set(gca,'XTick',[],'YTick',[])
        title('Model Phase Map','FontSize',18,'FontWeight','Bold')
        % colorbar, cbfreeze
        subplot(122), imagesc(pb(:,:,i)), colormap(bone), freezeColors, 
        axis square
        set(gca,'XTick',[],'YTick',[])
        title('Probabalistic Boundaries \in (0,1)','FontSize',18,'FontWeight','Bold')
        % colorbar, cbfreeze

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

