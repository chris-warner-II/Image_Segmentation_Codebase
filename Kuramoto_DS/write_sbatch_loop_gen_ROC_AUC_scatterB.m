

% Input: Loop through these
netMethod = {'IsoDiff','AAnrm','GLnrm','Mod_SKHAdj','Mod_N&G'}; % 
radMax = [1, 3, 5, 10]; % [1:10,inf];
sigPix = [0.2]; %[0.1, 0.2, 0.3, 0.4];

sigDist = inf;
sigW = 0;
% Kscale = 1; %300;  
% Tscale = 1;

imdim =  '51x51_ds1';
InputImages = 'BSDS_patch'; % 'GradientBox'; %  'BSDS_tile'; % 'BSDS_full'; %

KSstr = {'sml','mid','lrg'};


jobs_per_cpu = 1; % 44 combinations of rM & sP
numCPUs = 1; % can be up to 4 if you want to cluster jobs on same machine.




% Set up vector with RadiusMax values to loop over
rM = zeros( 1,numel(radMax)*numel(sigPix) );
st=1;
for i = 1:numel(radMax)
    nd = st - 1 + numel(sigPix);
    rM(st:nd) = radMax(i);
    st = nd + 1;    
end
rM = repmat(rM,1,numel(netMethod)); % duplicate everything for multiple methods.





% Set up vector with SigmaPixels values to loop over
sP = repmat(sigPix,1,numel(radMax));
sP = repmat(sP,1,numel(netMethod)); % duplicate everything for multiple methods.






% Set up cell array of NetMethod values to loop over.
Meth = cell( 1,numel(radMax)*numel(sigPix) );
st=1;
for i = 1:numel(netMethod)
    nd = st - 1 + numel(sigPix)*numel(radMax);
    for j = st:nd
        Meth{j} = netMethod{i};
    end
    st = nd + 1;    
end








if ~exist('scripts4cluster','dir')
   mkdir('scripts4cluster') 
end



% Delete Old Versionss of scripts, tasks & task files and repopulate them.
!rm scripts4cluster/script_genROCAUC*
!rm scripts4cluster/taskfile_genROCAUC*
!rm scripts4cluster/task_genROCAUC*.sh



for k = 1:( numel(rM)./(numCPUs*jobs_per_cpu) )


    %% (#1). Build up script_ImgSeg file that calls helper_ht & 
    fid = fopen(['scripts4cluster/script_genROCAUC',num2str(k)], 'w+');
    fprintf(fid,['#!/bin/bash\n']);
    fprintf(fid,['# Job name:\n']);
    fprintf(fid,['#SBATCH --job-name=AUC',num2str(k),'\n']);
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
    fprintf(fid,['#SBATCH --mem-per-cpu=3500M\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['# Wall clock limit:\n']);
    fprintf(fid,['#SBATCH --time=5:0:00\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['ht_helper.sh -t taskfile_genROCAUC',num2str(k),' -n1 -s60 -dvk']);
    fclose(fid);
    

    
    %% (#2). Make the taskfile_ImgSeg# file that pertains to each script_ImgSeg# file
    fid = fopen(['scripts4cluster/taskfile_genROCAUC',num2str(k)], 'w+');
    for i = 1:(numCPUs*jobs_per_cpu)
        fprintf(fid,['task_genROCAUC',num2str( (k-1)*(numCPUs*jobs_per_cpu) + i  ),'.sh\n']);
    end
    fclose(fid);
                            
end                          
                            
                            
                            
for i = 1: numel(rM)                            
                            
    %% (#2.) Build up each individual task_ImgSeg#.sh file
    fid = fopen(['scripts4cluster/task_genROCAUC',num2str(i),'.sh'], 'w+');
    fprintf(fid,['#!/bin/bash\n']);
    fprintf(fid,['module load matlab/R2013a\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['echo generating ROC_AUC scatter plots and data files',num2str(i),'\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['matlab -nosplash -nodesktop -noFigureWindows << EOF\n']);   % was matlab -nodisplay but that stopped working on cluster.
    fprintf(fid,['\n']);
    fprintf(fid,['cd /global/home/users/cwarner/Projects\n']);
    fprintf(fid,['addpath(genpath(pwd))\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['gen_ROC_AUC_scatter(''',InputImages,''',''',imdim,''',{''',Meth{i},'''},',num2str(rM(i)),',',...
        num2str(sigDist),',',num2str(sP(i)),',',num2str(sigW),');\n']); % ',',num2str(Tscale),',',num2str(Kscale),
    fprintf(fid,['\n']);
    fprintf(fid,['exit\n']);
    fprintf(fid,['EOF\n']);
    fclose(fid);
                        


end



% Make a run file that will call sbatch script_ImgSeg##
% Note: can not just run file because I do not have permission on cluster.
%       But I can just copy and paste its contents into command line. Works.
% or run chmod 777 Run_scripts_... before running it by ./Run_scripts_...
fid = fopen(['scripts4cluster/Run_scripts_genROCAUC'], 'w+');
for i=1:( numel(rM)./(numCPUs*jobs_per_cpu) )
    fprintf(fid,['sbatch script_genROCAUC',num2str(i),'\n']);
    fprintf(fid,['sleep 1 \n']);
end

fprintf(fid,['squeue\n']);
fclose(fid);



