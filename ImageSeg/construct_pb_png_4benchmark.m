function construct_pb_png_4benchmark(method,rM,blur_flg)

% This script is meant to create a connection between the results of our
% various Image Segmentation alrgorithms and the benchmark code provided
% with the Berkeley Segmentation Dataset. We have saved output mat files
% (KurMC & Evecs) which contain spatial gradient maps. This will take those
% spatial gradient maps, normalize them to be between 0 & 1 and save them 
% as png images.  These probabalistic boundary (pb) png images will be used
% as input to the benchmark code which will use them to compute precision &
% recall & f-measure.


% method = 'IsoDiff';


% if ~exist('rM','var')
%     
% end


[dirPre,sizeGoodIm] = onCluster;

%;
if(blur_flg==1)
    blur_tag_M = '_blur_sig1';
    blur_tag_I = 'blur_sz13_sig1/';
elseif(blur_flg==2)
    blur_tag_M = '_blur_sigC1_S8_Kr0p01';
    blur_tag_I = 'blur_sigC1_S8_Kr0p01/';
else
    blur_tag_M = ''; % if we are not blurring.
    blur_tag_I = '';
end



%% (1). This is for Image file (just spatial gradients in image).

if ~blur_flg % can do this here but it is already being handled in make_gaussian_blur_pb_png_img_patches.m function.
    
    disp(['Creating PB_png image files for spatial gradients in actual image.'])

    % location of mat files with extracted image patches
    matInDir = [dirPre,'images/BSDS_patch/101x101_ds1/'];
    matOutDir = [dirPre,'images/BSDS_patch/101x101_ds1/',blur_tag_I];

    % location to write out the probabalistic boundaries png image file (input to benchmark code).
    pbOutDir = [matInDir,'pb_png/'];
    if ~exist(pbOutDir,'dir')
        mkdir(pbOutDir)
    end



    disp(['Creating Probabalistic Boundary PB image files for Image Spatial Gradients'])
    
    
    files = dir([matInDir,'*.mat']);
    disp(['Number of files: ',numel(files)])
    

    % Loop through all the files to convert F matrix (of spatial gradients) to png pb image file.
    for i = 1:numel(files)

        fname = files(i).name;
        disp([num2str(i),' : ',fname])

        % do not rerun this if the pb_png files already exist.
        if exist([pbOutDir,fname(1:end-4),'.png'],'file')
            disp([num2str(i),' Probabalistic Boundary PNG File: ',fname(1:end-4),'.png',' :already exists. Moving on.'])
            continue
        end

        disp([num2str(i),' / ',num2str(numel(files)),' : ', fname(1:end-4),'.png'])


        % load KurMC mat file
        load([matInDir,fname])


        if(blur_flg)

            % TODO: Maybe...
            %
            % can do this here but it is already being handled in
            % make_gaussian_blur_pb_png_img_patches.m function &
            % make_DoG_filter_pb_png_img_patches.m function.

        end




        pb = MC.F ./ max(MC.F(:)); % turn spatial gradients into "probabalistic boundaries"

        imwrite(pb,[pbOutDir,fname(1:end-4),'.png'],'PNG','BitDepth',8)


    end

end



%% (2). This is for KurMC file.

disp(['Creating PB_png image files for spatial gradients Phase of Kuramoto Simulation.'])

% location of mat files with extracted image patches
matInDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/Kur_PIF_Fourier1/',method,'/'];

% 

files = dir([matInDir,'KurMC_*rM',rM,'_*.mat']);

disp(['Number of files: ',numel(files)])

% Loop through all the files to convert F matrix (of spatial gradients) to png pb image file.
for i = 1:numel(files)
    
    fname = files(i).name;
    disp([num2str(i),' : ',fname])
    
    % parse fname into it's constituent parts
    rM_loc = strfind(fname,'_rM');
    sD_loc = strfind(fname,'_sD');
    sP_loc = strfind(fname,'_sP');
    NF_loc = strfind(fname,'_NF');
    ks_loc = strfind(fname,'_ks');
    
    if strcmp(method,'IsoDiff') 
        fn_tag = fname(7:rM_loc-1);
        rM_tag = fname(rM_loc+1:NF_loc-1);
        NF_tag = fname(NF_loc+1:ks_loc-1);
        ks_tag = fname(ks_loc+1:end-4);
        % location to write out the probabalistic boundaries png image file (input to benchmark code).
        pbOutDir = [matInDir,'pb_png/',rM_tag,'/',NF_tag,'/',ks_tag,'/'];
    else
        fn_tag = fname(7:rM_loc-1);
        rM_tag = fname(rM_loc+1:sD_loc-1);
        sD_tag = fname(sD_loc+1:sP_loc-1);
        sP_tag = fname(sP_loc+1:NF_loc-1);
        NF_tag = fname(NF_loc+1:ks_loc-1);
        ks_tag = fname(ks_loc+1:end-4);
        % location to write out the probabalistic boundaries png image file (input to benchmark code).
        pbOutDir = [matInDir,'pb_png/',rM_tag,'/',sD_tag,'/',sP_tag,'/',NF_tag,'/',ks_tag,'/'];
        
    end
    
    
    
    % Make output directory if it does not already exist.
    if ~exist(pbOutDir,'dir')
        mkdir(pbOutDir)
    end
    

    
    % do not rerun this if the pb_png files already exist.
    if exist([pbOutDir,fn_tag,'.png'],'file')
        disp([num2str(i),' Probabalistic Boundary PNG File: ',fname(7:end-4),' :already exists. Moving on.'])
        continue
    end
    
    disp([num2str(i),' / ',num2str(numel(files)),' : ', fname(7:end-4)])
    
    % load KurMC mat file
    try
        load([matInDir,fname])
    catch
        disp('Something went wrong.  File must be corrupt. Skipping:')
        [matInDir,fname]
        continue
    end
    

    
    pb = MC.F ./ max(MC.F(:)); % turn spatial gradients into "probabalistic boundaries"
    
    imwrite(pb,[pbOutDir,fn_tag,'.png'],'PNG','BitDepth',8)
    
end



%% (3). This is for Evecs files.


disp(['Creating PB_png image files for spatial gradients in Eigenvectors.'])


% location of mat files with extracted image patches
matInDir = [dirPre,'output/Kuramoto/NetsFromImgs/BSDS_patch_101x101_ds1',blur_tag_M,'/data/spectral/',method,'/'];

% location to write out the probabalistic boundaries png image file (input to benchmark code).
pbOutDir = [matInDir,'pb_png/'];
if ~exist(pbOutDir,'dir')
    mkdir(pbOutDir)
end


% 
if strcmp(method,'IsoDiff')
    files = dir([matInDir,'Evecs_*_rM',rM,'.mat']);
else
    files = dir([matInDir,'Evecs_*_rM',rM,'_*.mat']);
end

disp(['Number of files: ',numel(files)])


% Loop through all the files to convert F matrix (of spatial gradients) to png pb image file.
for i = 1:numel(files)
    
    fname = files(i).name;
    disp([num2str(i),' : ',fname])
    
    % parse fname into it's constituent parts
    rM_loc = strfind(fname,'_rM');
    sD_loc = strfind(fname,'_sD');
    sP_loc = strfind(fname,'_sP');
    %
    fn_tag = fname(7:rM_loc-1);
    rM_tag = fname(rM_loc+1:sD_loc-1);
    sD_tag = fname(sD_loc+1:sP_loc-1);
    sP_tag = fname(sP_loc+1:end-4);
    
    
    if strcmp(method,'IsoDiff') 
        fn_tag = fname(7:rM_loc-1);
        rM_tag = fname(rM_loc+1:end-4);
        % location to write out the probabalistic boundaries png image file (input to benchmark code).
        pbOutDir = [matInDir,'pb_png/',rM_tag,'/'];
    else
        fn_tag = fname(7:rM_loc-1);
        rM_tag = fname(rM_loc+1:sD_loc-1);
        sD_tag = fname(sD_loc+1:sP_loc-1);
        sP_tag = fname(sP_loc+1:end-4);
        % location to write out the probabalistic boundaries png image file (input to benchmark code).
        pbOutDir = [matInDir,'pb_png/',rM_tag,'/',sD_tag,'/',sP_tag,'/'];
        
    end
    
    
    
    % Make directories to hold output pb.png files if they dont exist already.
    if ~exist([pbOutDir,'ev1/'],'dir')
        mkdir([pbOutDir,'ev1/'])
    end
    %
    if ~exist([pbOutDir,'ev2o/'],'dir')
        mkdir([pbOutDir,'ev2o/'])
    end
    %
    if ~exist([pbOutDir,'ev3o/'],'dir')
        mkdir([pbOutDir,'ev3o/'])
    end
    %
    if ~exist([pbOutDir,'ev2/'],'dir')
        mkdir([pbOutDir,'ev2/'])
    end
    %
    if ~exist([pbOutDir,'ev3/'],'dir')
        mkdir([pbOutDir,'ev3/'])
    end
    %
    if ~exist([pbOutDir,'ev2w/'],'dir')
        mkdir([pbOutDir,'ev2w/'])
    end
    %
    if ~exist([pbOutDir,'ev3w/'],'dir')
        mkdir([pbOutDir,'ev3w/'])
    end
    

    
    
    % do not rerun this if the pb_png files already exist (check just ev1 assuming others will be there if it is).
    if exist([pbOutDir,'ev1/',fn_tag,'.png'],'file')
        disp([num2str(i),' Probabalistic Boundary PNG File: ',fname(7:end-4),' :already exists. Moving on.'])
        continue
    end
    
    
    
    
    
    disp([num2str(i),' / ',num2str(numel(files)),' : ', fname(7:end-4),'.png'])
    
    % load Evecs mat file
    try
        load([matInDir,fname])
    catch
        disp('Something went wrong.  File must be corrupt. Skipping:')
        [matInDir,fname]
        continue
    end
    
    pb = MC.ev1.F ./ max(MC.ev1.F(:)); % turn spatial gradients into "probabalistic boundaries"
    imwrite(pb,[pbOutDir,'ev1/',fn_tag,'.png'],'PNG','BitDepth',8)
    
    pb = MC.ev2o.F ./ max(MC.ev2o.F(:)); % turn spatial gradients into "probabalistic boundaries"
    imwrite(pb,[pbOutDir,'ev2o/',fn_tag,'.png'],'PNG','BitDepth',8)
    
    pb = MC.ev3o.F ./ max(MC.ev3o.F(:)); % turn spatial gradients into "probabalistic boundaries"
    imwrite(pb,[pbOutDir,'ev3o/',fn_tag,'.png'],'PNG','BitDepth',8)
    
    pb = MC.ev2.F ./ max(MC.ev2.F(:)); % turn spatial gradients into "probabalistic boundaries"
    imwrite(pb,[pbOutDir,'ev2/',fn_tag,'.png'],'PNG','BitDepth',8)
    
    pb = MC.ev3.F ./ max(MC.ev3.F(:)); % turn spatial gradients into "probabalistic boundaries"
    imwrite(pb,[pbOutDir,'ev3/',fn_tag,'.png'],'PNG','BitDepth',8)
    
    pb = MC.ev2w.F ./ max(MC.ev2w.F(:)); % turn spatial gradients into "probabalistic boundaries"
    imwrite(pb,[pbOutDir,'ev2w/',fn_tag,'.png'],'PNG','BitDepth',8)
    
    pb = MC.ev3w.F ./ max(MC.ev3w.F(:)); % turn spatial gradients into "probabalistic boundaries"
    imwrite(pb,[pbOutDir,'ev3w/',fn_tag,'.png'],'PNG','BitDepth',8)


end

