% Script to combine all ~1500 mat files for individual image patches and
% combine their Dout results into a single mat file.  To make next step in
% processing quicker.  Do all the i/o now and only need to load 1 file
% later.


[dirPre, sizeGoodIm] = onCluster;
patchesDir = [dirPre,'images/BSDS_patch/51x51_ds1_wDopt/'];

% Directory to save patches and result from this optimization into
saveDir = [dirPre,'images/BSDS_patch/51x51_ds1_wDopt/'];
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end


files = dir([patchesDir,'*.mat']);

% Preallocate Arrays to Fill Up as We Loop Thru Files
Dout_BestC_tot = [];
Dout_BestL_tot = [];
Dout_BestE_tot = [];
Dout_DiffC_tot = [];
Dout_DiffL_tot = [];
Dout_DiffE_tot = [];
ImgGtID = [];
ImgPtchID = {};
GtSegSize = {};
k = 0; % counter for cells


% Loop through files

for i = 1:numel(files)
    
    disp(['File #',num2str(i),' / ',num2str(numel(files))])
    load([patchesDir,files(i).name])
    
    Dout_BestC_tot = [Dout_BestC_tot, Dout_BestC];
    Dout_BestL_tot = [Dout_BestL_tot, Dout_BestL];
    Dout_BestE_tot = [Dout_BestE_tot, Dout_BestE];
    Dout_DiffC_tot = [Dout_DiffC_tot, Dout_DiffC];
    Dout_DiffL_tot = [Dout_DiffL_tot, Dout_DiffL];
    Dout_DiffE_tot = [Dout_DiffE_tot, Dout_DiffE];
    ImgGtID = [ImgGtID, 1:numel(gT)];
    
    for j = 1:numel(gT)
        k=k+1;
        ImgPtchID{k} = files(i).name(1:end-4);
        GtSegSize{k} = clusterSize{j};
    end
    
end


% Save combined file
save([saveDir,'DoutIdeal_allPatches'],'Dout_BestC_tot','Dout_BestL_tot','Dout_BestE_tot',...
    'Dout_DiffC_tot','Dout_DiffL_tot','Dout_DiffE_tot','ImgGtID','ImgPtchID','GtSegSize')