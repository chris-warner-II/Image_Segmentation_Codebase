
% imdim =  {[101,101]}; % {'full'}; %
InputImages = 'BSDS_patch';
MethodType = {'Mod_SKHAdj'}; % ,'AAnrm','GLnrm','Mod_N&G'

PTCH = {'1'}; % ,'2','3'
RM = {'1'}; % ,'2','3','4','5','6','7','8','9','10','Inf'
SP = {'0p1'}; % ,'0p3','0p2','0p4'


if ~exist('scripts4cluster','dir')
   mkdir('scripts4cluster') 
end


% Delete Old Versions of scripts, tasks & task files and repopulate them.
!rm scripts4cluster/script_ClusterEigPairROC*
!rm scripts4cluster/Run_scripts_ClusterEigPairROC
!rm scripts4cluster/View_outs_ClusterEigPairROC
!rm scripts4cluster/View_errs_ClusterEigPairROC

k=0;

for g = 1:numel(MethodType)
    for i = 1:numel(PTCH)
        for j = 1:numel(RM)
            for h = 1:numel(SP)

                k=k+1;

                fid = fopen(['scripts4cluster/script_ClusterEigPairROC',num2str(k)], 'w+');
                fprintf(fid,['#!/bin/bash\n']);
                fprintf(fid,['# Job name:\n']);
                fprintf(fid,['#SBATCH --job-name=EROC',num2str(k),'\n']);
                fprintf(fid,['#\n']);
                fprintf(fid,['# Partition:\n']);
                fprintf(fid,['#SBATCH --partition=cortex\n']);
                fprintf(fid,['#\n']);
                fprintf(fid,['# Processors:\n']);
                fprintf(fid,['#SBATCH --ntasks=4\n']);
                fprintf(fid,['#\n']);
                fprintf(fid,['# Constrain Nodes:\n']);
                fprintf(fid,['#SBATCH --constraint=cortex_nogpu\n']);
                %fprintf(fid,['#\n']);
                %fprintf(fid,['# Exclude the Following Nodes:\n']);
                %fprintf(fid,['#SBATCH -x n0000.cortex0,n0001.cortex0,n0012.cortex0,n0013.cortex0\n']); % n0007.cortex0,n0009.cortex0,n0010.cortex0,n0011.cortex0,n0001.cortex0,n0012.cortex0,n0013.cortex0,
                fprintf(fid,['#\n']);
                fprintf(fid,['# Memory:\n']);
                fprintf(fid,['#SBATCH --mem-per-cpu=2G\n']);
                fprintf(fid,['#\n']);
                fprintf(fid,['# Wall clock limit:\n']);
                fprintf(fid,['#SBATCH --time=4:0:00\n']);
                fprintf(fid,['#\n']);
                fprintf(fid,['#\n']);
                fprintf(fid,['#SBATCH -o EROC',num2str(k),'.out\n']);
                fprintf(fid,['#\n']);
                fprintf(fid,['#SBATCH -e EROC',num2str(k),'.err\n']);
                fprintf(fid,['\n']);
                fprintf(fid,['module load matlab/R2013a\n']);
                fprintf(fid,['\n']);
                fprintf(fid,['echo Computing Eigen Cluster Pairwise ROC Area under Curve for all files',num2str(k),'\n']);
                fprintf(fid,['\n']);
                fprintf(fid,['matlab -nosplash -nodesktop -noFigureWindows << EOF\n']);   % was matlab -nodisplay but that stopped working on cluster.
                fprintf(fid,['\n']);
                fprintf(fid,['cd /global/home/users/cwarner/Projects\n']);
                fprintf(fid,['addpath(genpath(pwd))\n']);
                fprintf(fid,['\n']);
                fprintf(fid,['compute_Eig_clusterPairROC_fixFiles(''',InputImages,''',''',MethodType{g},''',''',PTCH{i},''',''',RM{j},''',''',SP{h},''');\n']);
                fprintf(fid,['\n']);
                fprintf(fid,['exit\n']);
                fprintf(fid,['EOF\n']);
                fclose(fid);

            end
        end
    end
end






% Make a run file that will call sbatch script_ImgSeg##
% Note: can not just run file because I do not have permission on cluster.
%       But I can just copy and paste its contents into command line. Works.
fid = fopen(['scripts4cluster/Run_scripts_ClusterEigPairROC'], 'w+');
for i=1:k
    fprintf(fid,['sbatch script_ClusterEigPairROC',num2str(i),'\n']);
    fprintf(fid,['sleep 1 \n']);
end

fprintf(fid,['squeue\n']);
fclose(fid);



% Make a run file that will display last few lines of all *.out files to look thru quickly
fid = fopen(['scripts4cluster/View_outs_ClusterEigPairROC'], 'w+');
for i=1:k
    fprintf(fid,['echo --------------------------------------------------------------------------\n']);
    fprintf(fid,['echo EROC',num2str(i),'.out\n']);
    fprintf(fid,['tail -n 30 EROC',num2str(i),'.out\n']);
    fprintf(fid,['\n']);
end

fclose(fid);


% Make a run file that will display last few lines of all *.err files to look thru quickly
fid = fopen(['scripts4cluster/View_errs_ClusterEigPairROC'], 'w+');
for i=1:k
    fprintf(fid,['echo --------------------------------------------------------------------------\n']);
    fprintf(fid,['echo EROC',num2str(i),'.err\n']);
    fprintf(fid,['tail -n 20 EROC',num2str(i),'.err\n']);
    fprintf(fid,['\n']);
end

fclose(fid);



