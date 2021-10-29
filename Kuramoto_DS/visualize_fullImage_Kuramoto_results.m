[dirPre,sizeGoodIm] = onCluster;

dirVertIn =  [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_full_321x481_ds1/data/Kur_PIF_Fourier1/Mod_SKHAdj/'];
dirVertOut = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_full_321x481_ds1/imgs/Kur_PIF_Fourier1/Mod_SKHAdj/'];
%
dirHorzIn =  [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_full_481x321_ds1/data/Kur_PIF_Fourier1/Mod_SKHAdj/'];
dirHorzOut = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_full_481x321_ds1/imgs/Kur_PIF_Fourier1/Mod_SKHAdj/'];

dirPickIn = dirHorzIn;   % dirVertIn or dirHorzIn.
dirPickOut = dirHorzOut; % should be same as above (with Out)

filesInPick = dir(dirPickIn);

if ~exist(dirPickOut,'dir')
    mkdir(dirPickOut)
end

for i = 1:numel(filesInPick)
    
    if exist([dirPickOut,filesInPick(i).name(1:end-4),'.jpg'],'file')
       disp(['Output Image file already exists:'])
       [dirPickOut,filesInPick(i).name(1:end-4),'.jpg']
       disp(['Moving on to next...'])
       continue
        
    end


    try
        load([dirPickIn,filesInPick(i).name])
    catch
        disp(['Iteration # ',num2str(i),': Well that didnt work, maybe it was a directory.'])
        dirPickIn
        filesInPick(i).name
        continue
    end


    disp(['Iteration #',num2str(i)])


    % Compute Divisive Margin
    if iscell(MC)
        DivMarg = squeeze(MC{1}.DistAvgPW(:,1,:))./squeeze(MC{1}.DistAvgPW(:,2,:));
    end
    
    
    

    % Convert final phase solution back into linear variable to compare to pixels (just for visualization).
    imageOrig = netParams.im;
    [one, two] = find(imageOrig == min(imageOrig(:))); % find brightest pixel.
    pixMax(1) = one(1);
    pixMax(2) = two(1);

    phaseFinal = reshape(metaCluster.phaseAtClk(:,end),netParams.Ndims);
    phaseOffset = pi - phaseFinal(pixMax(1),pixMax(2));           % global phase to add to all phases to make phase at brightest pixel = pi
    phaseFinal = mod( phaseFinal + phaseOffset , 2*pi);           % rotate all phases by phase offset
    %phaseFinal = abs(phaseFinal - pi) ./ pi;                      % fold circular variabls down to [0 1] number line - ie. pixel space

    phaseInit = reshape(metaCluster.phaseAtClk(:,1),netParams.Ndims);
    phaseOffset = pi - phaseInit(pixMax(1),pixMax(2));            % global phase to add to all phases to make phase at brightest pixel = pi
    phaseInit = mod( phaseInit + phaseOffset , 2*pi);             % rotate all phases by phase offset
    %phaseInit = abs(phaseInit - pi) ./ pi;                       % fold circular variabls down to [0 1] number line - ie. pixel space





    % Plot these things: 
    h=figure; 
    subplot(4,2,[1,3]), imagesc(imageOrig), set(gca,'Xtick',[],'Ytick',[]), title('Original Image & Initial Phase Embedding','FontSize',20,'FontWeight','Bold'),
    colormap('bone'), caxis([0 1]), colorbar, cbfreeze, freezeColors
    subplot(4,2,[2,4]), imagesc(phaseFinal), set(gca,'Xtick',[],'Ytick',[]), title(['Final Phase after Simulation (t=0.3s) '],'FontSize',20,'FontWeight','Bold'),
    colormap('hsv'), caxis([0 2*pi]), colorbar, cbfreeze, freezeColors
    %
    for g = 1:numel(netParams.gT)
        subplot(4,numel(netParams.gT),2*numel(netParams.gT)+g), imagesc(netParams.gT{g}), colormap('jet'), set(gca,'Xtick',[],'Ytick',[]), freezeColors
    end
    subplot(4,numel(netParams.gT),2*numel(netParams.gT)+1), ylabel('gnd truth','FontSize',18,'FontWeight','Bold')
    %
    if iscell(MC)
        subplot(414),plot(kurParams.tau*kurParams.spp*(0:size(DivMarg,1)-1) , DivMarg)
        xlabel('Time (sec)','FontSize',18,'FontWeight','Bold')
        ylabel('Divisive Margin','FontSize',18,'FontWeight','Bold')
        axis([0 kurParams.Tsec 0 1])
    end


    
    
    % Save image comparing image pixels to Final Oscillator Phase Distribution
    saveGoodImg(h,[dirPickOut,filesInPick(i).name(1:end-4),'.jpg'],sizeGoodIm)
    close(h);



end







% Can most likely borrow a lot from this code.
%edit visualize_pts_from_2x2D_scatter
