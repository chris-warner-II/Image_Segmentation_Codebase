method = [7,8]; % [3,5,6,7,12]; %  % 3,5,6,7,% 5=GLnrm, 6=Mod_N&G or 7=Mod_SKHAdj 8=Mod_SKHEuc or 3=AAnrm or 12=IsoDiff
radmax = [1,3,5,10]; % [1:10,inf]; ,3,5,10
sigpix = [0.2]; % [0.1,0.2,0.3,0.4];
imdim =  {[101,101]}; % {'full'}; %
InputImages = 'BSDS_patch'; % 'BSDS_full'; %'BSDS_full'; % 'GradientBox'; %  'BSDS_tile'; % ;


L=1;
if ischar(imdim{L})
    imdim_cell = ['{''',imdim{L},'''}'];
else
    imdim_cell = ['{[',num2str(imdim{L}(1)),',',num2str(imdim{L}(2)),']}'];
end


% Put this here so I can now parallelize by patches extracted from single image if I like.
PatchOffset = 0;
numPatches = 500; % 1500; % 1000 possible for the 101x101 kind.
numP4loop = 5; % 100; % number of jobs to run within the same matlab call.

imgNums2procLoop(:,1) = PatchOffset + [1:numP4loop:numPatches]; %1:70; % start image patch number
imgNums2procLoop(:,2) = PatchOffset + [numP4loop:numP4loop:numPatches]; %1:70; % end image patch number



jobs_per_cpu = 1; % 44 combinations of rM & sP
numCPUs = 4; % is 4 on cluster I think.
% 1500 script_ImgSeg & taskfile_ImgSeg files


SigW = 0; %[0.15*60]; % 0.01*60,0.03*60,0.05*60,0.1*60,



numel(method)
numel(radmax)
numel(sigpix)
numel(SigW)

rM = zeros( 1,size(imgNums2procLoop,1)*numel(radmax)*numel(sigpix) );
sP = zeros( 1,size(imgNums2procLoop,1)*numel(radmax)*numel(sigpix) );
% sW = zeros( 1,size(imgNums2procLoop,1)*numel(radmax)*numel(sigpix) );
ims = zeros( size(imgNums2procLoop,1)*numel(radmax)*numel(sigpix), 2 );

st=1;
st2=1;
for i = 1:numel(radmax)
    nd = st - 1 + size(imgNums2procLoop,1)*numel(sigpix);
    rM(st:nd) = radmax(i);
    st = nd + 1;
    
    
    for j = 1:numel(sigpix)
        nd2 = st2 - 1 + size(imgNums2procLoop,1);
        sP(st2:nd2) = sigpix(j);
        ims(st2:nd2,:) = imgNums2procLoop;
        st2 = nd2 + 1;
    end
    
end


% duplicate everything for multiple methods.
rM = repmat(rM,1,numel(method));
sP = repmat(sP,1,numel(method));
ims = repmat(ims,numel(method),1);

mth = [];
for i = 1:numel(method)
    mth = [mth,repmat(method(i),1,size(imgNums2procLoop,1)*numel(radmax)*numel(sigpix))];
end


% just a plot for visualization to see if things are lined up right.
if(1)
    figure, hold on,
    plot(1:numel(rM),rM,'b')
    plot(1:numel(sP),sP./max(sP),'r')
    plot(1:size(ims,1),ims./max(ims(:)),'g')
    plot(1:numel(mth),mth,'c')
end






if ~exist('scripts4cluster','dir')
   mkdir('scripts4cluster') 
end



% Delete Old Versionss of scripts, tasks & task files and repopulate them.
!rm scripts4cluster/script_ImgSeg*
!rm scripts4cluster/taskfile_ImgSeg*
!rm scripts4cluster/task_ImgSeg*.sh


for k = 1:(numel(rM)./(numCPUs*jobs_per_cpu) )


    %% (#1). Build up script_ImgSeg file that calls helper_ht & 
    fid = fopen(['scripts4cluster/script_ImgSeg',num2str(k)], 'w+');
    fprintf(fid,['#!/bin/bash\n']);
    fprintf(fid,['# Job name:\n']);
    fprintf(fid,['#SBATCH --job-name=ImSeg',num2str(k),'\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['# Partition:\n']);
    fprintf(fid,['#SBATCH --partition=cortex\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['# Processors:\n']);
    fprintf(fid,['#SBATCH --ntasks=4\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['# Constrain Nodes:\n']);
    fprintf(fid,['#SBATCH --constraint=cortex_nogpu\n']);
%     fprintf(fid,['#\n']);
%     fprintf(fid,['# Exclude the Following Nodes:\n']);
%     fprintf(fid,['#SBATCH -x n0000.cortex0,n0001.cortex0\n']); % ,n0007.cortex0,n0008.cortex0,n0009.cortex0,n0010.cortex0,n0011.cortex0,,n0012.cortex0,n0013.cortex0
    fprintf(fid,['#\n']);
    fprintf(fid,['# Memory:\n']);
    fprintf(fid,['#SBATCH --mem-per-cpu=3500M\n']); % note: 3.5GB good for 51x51 patches.  7.5GB good for full images?
    fprintf(fid,['#\n']);
    fprintf(fid,['# Wall clock limit:\n']);
    fprintf(fid,['#SBATCH --time=0-12:0:00\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['ht_helper.sh -t taskfile_ImgSeg',num2str(k),' -n1 -s60 -dvk']);
    fclose(fid);
    

    
    %% (#2). Make the taskfile_ImgSeg# file that pertains to each script_ImgSeg# file
    fid = fopen(['scripts4cluster/taskfile_ImgSeg',num2str(k)], 'w+');
    for i = 1:(numCPUs*jobs_per_cpu)
        fprintf(fid,['task_ImgSeg',num2str( (k-1)*(numCPUs*jobs_per_cpu) + i  ),'.sh\n']);
    end
    fclose(fid);
                            
end                          
                            
                            
                            
for i = 1: numel(rM)                            
                            
    %% (#2.) Build up each individual task_ImgSeg#.sh file
    fid = fopen(['scripts4cluster/task_ImgSeg',num2str(i),'.sh'], 'w+');
    fprintf(fid,['#!/bin/bash\n']);
    fprintf(fid,['module load matlab/R2013a\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['echo Kuramoto Image Segmentation',num2str(i),'\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['matlab -nosplash -nodesktop -noFigureWindows << EOF\n']);   % was matlab -nodisplay but that stopped working on cluster.
    fprintf(fid,['\n']);
    fprintf(fid,['cd /global/home/users/cwarner/Projects\n']);
    fprintf(fid,['addpath(genpath(pwd))\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['Loop_ImgSegMethodsD(',num2str(mth(i)),',',num2str(rM(i)),',',num2str(sP(i)),',',imdim_cell,...
                ',''',InputImages,''',[',num2str(ims(i,1)),',',num2str(ims(i,2)),'],',num2str(SigW),');\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['exit\n']);
    fprintf(fid,['EOF\n']);
    fclose(fid);
                        


end





% Make a run file that will call sbatch script_ImgSeg##
% Note: can not just run file because I do not have permission on cluster.
%       But I can just copy and paste its contents into command line. Works.
fid = fopen(['scripts4cluster/Run_scripts_ImgSeg'], 'w+');
for i=1:( numel(rM)./(numCPUs*jobs_per_cpu))
    fprintf(fid,['sbatch script_ImgSeg',num2str(i),'\n']);
    fprintf(fid,['sleep 1 \n']);
end

fprintf(fid,['squeue\n']);
fclose(fid);
% 
% 
% 
% % Make a run file that will display last few lines of all *.out files to look thru quickly
% fid = fopen(['scripts4cluster/View_outs_ImSeg'], 'w+');
% for i=1:k
%     fprintf(fid,['echo --------------------------------------------------------------------------\n']);
%     fprintf(fid,['echo ImSeg',num2str(i),'.out\n']);
%     fprintf(fid,['tail -n 30 ImSeg',num2str(i),'.out\n']);
%     fprintf(fid,['\n']);
% end
% 
% fclose(fid);
% 
% 
% % Make a run file that will display last few lines of all *.err files to look thru quickly
% fid = fopen(['scripts4cluster/View_errs_ImSeg'], 'w+');
% for i=1:k
%     fprintf(fid,['echo --------------------------------------------------------------------------\n']);
%     fprintf(fid,['echo ImSeg',num2str(i),'.err\n']);
%     fprintf(fid,['tail -n 20 ImSeg',num2str(i),'.err\n']);
%     fprintf(fid,['\n']);
% end
% 
% fclose(fid);



