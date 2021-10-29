% Make sbatch file to call loop_main_Kuramoto_HandCookedNetwork with some
% parameters to loop through


if ~exist('scripts4cluster','dir')
   mkdir('scripts4cluster') 
end



% Do this in bash window to see how many directories there are...
%  ls -l ./SegResEig | wc -l
% 
%  Or run ./count_SegResFiles from Projects folder.
%
%  Then loop from one to the value it returns.


for i = 1:113

    fid = fopen(['scripts4cluster/script_ROCcompare_EigVKur',num2str(i)], 'w+');

    fprintf(fid,['#!/bin/bash\n']);
    fprintf(fid,['# Job name:\n']);
    fprintf(fid,['#SBATCH --job-name=KROC',num2str(i),'\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['# Partition:\n']);
    fprintf(fid,['#SBATCH --partition=cortex\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['# Processors:\n']);
    fprintf(fid,['#SBATCH --ntasks=1\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['# Memory:\n']);
    fprintf(fid,['#SBATCH --mem-per-cpu=8G\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['# Wall clock limit:\n']);
    fprintf(fid,['#SBATCH --time=2-00:0:00\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['#SBATCH -o KROC',num2str(i),'.out\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['#SBATCH -e KROC',num2str(i),'.err\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['module load matlab/R2013a\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['echo Kuramoto / Eigenvector ROC Compare',num2str(i),'\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['matlab -nosplash -nodesktop -noFigureWindows << EOF\n']);   % was matlab -nodisplay but that stopped working on cluster.
    fprintf(fid,['\n']);
    fprintf(fid,['cd /global/home/users/cwarner/Projects\n']);
    fprintf(fid,['addpath(genpath(pwd))\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['ROCcompare_EigVKur(',num2str(i),');\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['exit\n']);
    fprintf(fid,['EOF\n']);
%     fprintf(fid,['fi\n']);
    fclose(fid);
    
end







% make a run file that will call sbatch script_HAWSamp##
% Note: can not just run file because I do not have permission on cluster.
%       But I can just copy and paste its contents into command line. Works.
fid = fopen(['scripts4cluster/Run_ROCcompare'], 'w+');
for i=1:113
    fprintf(fid,['sbatch script_ROCcompare_EigVKur',num2str(i),'\n']);
    fprintf(fid,['sleep 3 \n']);
end

fprintf(fid,['squeue\n']);
fclose(fid);