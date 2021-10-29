%start_loop = 96:5:490;
%finish_loop = 100:5:490;

start_loop =  [99,  110, 114, 118, 123, 127, 132, 137, 144, 149, 154, 159, 164, 168, 173, 178, 183, 187, 192, 197, 202];
finish_loop = [100, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205];


if ~exist('scripts4cluster','dir')
   mkdir('scripts4cluster') 
end


for i = 1:numel(start_loop) 


    fid = fopen(['scripts4cluster/script_cleanup_Kur_metaSummary',num2str(i)], 'w+');

    fprintf(fid,['#!/bin/bash\n']);
    fprintf(fid,['# Job name:\n']);
    fprintf(fid,['#SBATCH --job-name=KurCln',num2str(i),'\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['# Partition:\n']);
    fprintf(fid,['#SBATCH --partition=cortex\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['# Processors:\n']);
    fprintf(fid,['#SBATCH --nodes=1\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['# Memory:\n']);
    fprintf(fid,['#SBATCH --mem-per-cpu=1G\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['# Wall clock limit:\n']);
    fprintf(fid,['#SBATCH --time=2-00:0:00\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['#SBATCH -o KurCln',num2str(i),'.out\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['#SBATCH -e KurCln',num2str(i),'.err\n']);
    fprintf(fid,['\n']);
    %
    fprintf(fid,['\n']);
    fprintf(fid,['module load matlab/R2013a\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['echo Cleaning out duplicate Kur_metaSummary files ',num2str(i),'\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['matlab -nosplash -nodesktop -noFigureWindows << EOF\n']);   % was matlab -nodisplay but that stopped working on cluster.
    fprintf(fid,['\n']);
    fprintf(fid,['cd /global/home/users/cwarner/Projects\n']);
    fprintf(fid,['addpath(genpath(pwd))\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['cleanup_Kur_metaSummary_files(',num2str(start_loop(i)),',',num2str(finish_loop(i)),');\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['exit\n']);
    fprintf(fid,['EOF\n']);
    fclose(fid);


end







% make a run file that will call sbatch script_HAWSamp##
% Note: can not just run file because I do not have permission on cluster.
%       But I can just copy and paste its contents into command line. Works.
fid = fopen(['scripts4cluster/Run_scripts_cleanup_Kur_metaSummary'], 'w+');
for j=1:i
    fprintf(fid,['sbatch script_cleanup_Kur_metaSummary',num2str(j),'\n']);
    fprintf(fid,['sleep 1 \n']);
end

fprintf(fid,['squeue\n']);
fclose(fid);






% Make a run file that will display last few lines of all *.out files to look thru quickly
fid = fopen(['scripts4cluster/View_outs_KurCln'], 'w+');
for j=1:i
    fprintf(fid,['echo --------------------------------------------------------------------------\n']);
    fprintf(fid,['echo KurCln',num2str(j),'.out\n']);
    fprintf(fid,['tail -n 30 KurCln',num2str(j),'.out\n']);
    fprintf(fid,['\n']);
end

fclose(fid);


% Make a run file that will display last few lines of all *.err files to look thru quickly
fid = fopen(['scripts4cluster/View_errs_KurCln'], 'w+');
for j=1:i
    fprintf(fid,['echo --------------------------------------------------------------------------\n']);
    fprintf(fid,['echo KurCln',num2str(j),'.err\n']);
    fprintf(fid,['tail -n 20 KurCln',num2str(j),'.err\n']);
    fprintf(fid,['\n']);
end

fclose(fid);









