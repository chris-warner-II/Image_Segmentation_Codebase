[dirPre,sizeGoodIm] = onCluster;


% dirVertIn =  [dirPre,'output/Kuramoto/NetsFromImgs/GradientBox_51x51_ds1/data/spectral/Mod_SKHAdj/'];
% dirVertOut = [dirPre,'output/Kuramoto/NetsFromImgs/GradientBox_51x51_ds1/imgs/spectral/Mod_SKHAdj/'];

dirVertIn =  [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_full_321x481_ds1/data/spectral/Mod_SKHAdj/'];
dirVertOut = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_full_321x481_ds1/imgs/spectral/Mod_SKHAdj/'];
%
dirHorzIn =  [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_full_481x321_ds1/data/spectral/Mod_SKHAdj/'];
dirHorzOut = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_full_481x321_ds1/imgs/spectral/Mod_SKHAdj/'];

dirPickIn = dirVertIn;   % dirVertIn or dirHorzIn.
dirPickOut = dirVertOut; % should be same as above (with Out)

filesInPick = dir(dirPickIn);

if ~exist(dirPickOut,'dir')
    mkdir(dirPickOut)
end

% Nonlinearity applied to Eigenvectors to visualize them.
vizNonlinScale=1e-15;

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
    
    fileDeets = filesInPick(i).name(1:end-4);
    fileDeets(fileDeets=='_') = ' ';
    fileDeets = [fileDeets,' vizNonlin ',num2str(vizNonlinScale)];

    % Plot these things: 
    h=figure; colormap('bone')
    subplot(1,2,1), imagesc(netParams.im), set(gca,'Xtick',[],'Ytick',[]), axis square, title({'Original Image',fileDeets},'FontSize',20,'FontWeight','Bold')
    %
    subplot(2,6,4),imagesc(reshape( EvecVizF(EigVec(:,1),vizNonlinScale), netParams.Ndims(1),netParams.Ndims(2))), set(gca,'Xtick',[],'Ytick',[]), axis square, title('Evec1')
    subplot(2,6,5),imagesc(reshape( EvecVizF(EigVec(:,2),vizNonlinScale), netParams.Ndims(1),netParams.Ndims(2))), set(gca,'Xtick',[],'Ytick',[]), axis square, title('Evec2')
    subplot(2,6,6),imagesc(reshape( EvecVizF(EigVec(:,3),vizNonlinScale), netParams.Ndims(1),netParams.Ndims(2))), set(gca,'Xtick',[],'Ytick',[]), axis square, title('Evec3')
    %
    subplot(2,6,10),imagesc(reshape( EvecVizF(EigVec(:,4),vizNonlinScale), netParams.Ndims(1),netParams.Ndims(2))), set(gca,'Xtick',[],'Ytick',[]), axis square, title('Evec4')
    subplot(2,6,11),imagesc(reshape( EvecVizF(EigVec(:,5),vizNonlinScale), netParams.Ndims(1),netParams.Ndims(2))), set(gca,'Xtick',[],'Ytick',[]), axis square, title('Evec5')
    subplot(2,6,12),imagesc(reshape( EvecVizF(EigVec(:,6),vizNonlinScale), netParams.Ndims(1),netParams.Ndims(2))), set(gca,'Xtick',[],'Ytick',[]), axis square, title('Evec6')
            
    %keyboard

    % Save image comparing image pixels to Final Oscillator Phase Distribution
    fileDeets(fileDeets==' ') = '_';
    saveGoodImg(h,[dirPickOut,fileDeets,'.jpg'],sizeGoodIm)
    close(h);



end







% Can most likely borrow a lot from this code.
%edit visualize_pts_from_2x2D_scatter
