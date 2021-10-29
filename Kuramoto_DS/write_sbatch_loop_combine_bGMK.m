
%% Look in Ground Truths Directory to find what image patches
[dirPre,sizeGoodIm] = onCluster;


method =  {'Mod_N&G'}; % % 'IsoDiff','AAnrm','Mod_SKHAdj','GLnrm', 'ImPix'
RM = {'1','3','5','10'}; 
KS = {'sml','mid','lrg'}; 

C = 0; % counter for nested for loops


if ~exist('scripts4cluster','dir')
   mkdir('scripts4cluster') 
end


for i = 1:numel(method)
    for j = 1:numel(RM)
        for k = 1:numel(KS)
        
            C = C+1;

            fid = fopen(['scripts4cluster/script_combine_bGMK',num2str(C)], 'w+');

            fprintf(fid,['#!/bin/bash\n']);
            fprintf(fid,['# Job name:\n']);
            fprintf(fid,['#SBATCH --job-name=cbGMK',num2str(C),'\n']);
            fprintf(fid,['#\n']);
            fprintf(fid,['# Partition:\n']);
            fprintf(fid,['#SBATCH --partition=cortex\n']);
            fprintf(fid,['#\n']);
            fprintf(fid,['# Constrain Nodes:\n']);
            fprintf(fid,['#SBATCH --constraint=cortex_nogpu\n']);
            fprintf(fid,['#\n']);
            fprintf(fid,['# Processors:\n']);
            fprintf(fid,['#SBATCH --ntasks=4\n']);
            %fprintf(fid,['#\n']);
            %fprintf(fid,['# Exclude the Following Nodes:\n']);
            %fprintf(fid,['#SBATCH -x n0000.cortex0,n0001.cortex0,n0012.cortex0,n0013.cortex0,n0007.cortex0,n0008.cortex0,n0009.cortex0,n0010.cortex0,n0011.cortex0\n']); % ,n0001.cortex0,n0012.cortex0,n0013.cortex0,
            fprintf(fid,['#\n']);
            fprintf(fid,['# Memory:\n']);
            fprintf(fid,['#SBATCH --mem-per-cpu=7500M\n']);
            fprintf(fid,['#\n']);
            fprintf(fid,['# Wall clock limit:\n']);
            fprintf(fid,['#SBATCH --time=20:0:00\n']);
            fprintf(fid,['#\n']);
            fprintf(fid,['#\n']);
            fprintf(fid,['#SBATCH -o cbGMK',num2str(C),'.out\n']);
            fprintf(fid,['#\n']);
            fprintf(fid,['#SBATCH -e cbGMK',num2str(C),'.err\n']);
            fprintf(fid,['\n']);
            %
            fprintf(fid,['\n']);
            fprintf(fid,['module load matlab/R2013a\n']);
            fprintf(fid,['\n']);
            fprintf(fid,['echo "Combining Boundary Gradient Metric Data Kuramoto ',num2str(C),'"\n']);
            fprintf(fid,['date\n']);
            fprintf(fid,['\n']);
            fprintf(fid,['matlab -nosplash -nodesktop -noFigureWindows << EOF\n']);   % was matlab -nodisplay but that stopped working on cluster.
            fprintf(fid,['\n']);
            %fprintf(fid,['echo "Made it past"\n']);
            fprintf(fid,['datestr(clock)\n']);
            fprintf(fid,['\n']);
            fprintf(fid,['cd /global/home/users/cwarner/Projects\n']);
            fprintf(fid,['\n']);
            fprintf(fid,['addpath(genpath(pwd))\n']);
            fprintf(fid,['\n']);
            fprintf(fid,['combine_boundaryGradientMetric_Kur(''',method{i},''',''',RM{j},''',''',KS{k},''');\n']);
            fprintf(fid,['\n']);
            fprintf(fid,['exit\n']);
            fprintf(fid,['EOF\n']);
            fclose(fid);
            

        end
    end
end



% Make a run file that will call sbatch script_ImgSeg##
% Note: can not just run file because I do not have permission on cluster.
%       But I can just copy and paste its contents into command line. Works.
fid = fopen(['scripts4cluster/Run_scripts_combine_bGMK'], 'w+');
for j=1:C
    fprintf(fid,['sbatch script_combine_bGMK',num2str(j),'\n']);
    fprintf(fid,['sleep 1 \n']);
end

fprintf(fid,['squeue\n']);
fclose(fid);



% Make a run file that will display last few lines of all *.out files to look thru quickly
fid = fopen(['scripts4cluster/View_outs_combine_bGMK'], 'w+');
for j=1:C
    fprintf(fid,['echo --------------------------------------------------------------------------\n']);
    fprintf(fid,['echo cbGMK',num2str(j),'.out\n']);
    fprintf(fid,['tail -n 30 cbGMK',num2str(j),'.out\n']);
    fprintf(fid,['\n']);
end

fclose(fid);


% Make a run file that will display last few lines of all *.err files to look thru quickly
fid = fopen(['scripts4cluster/View_errs_combine_bGMK'], 'w+');
for j=1:C
    fprintf(fid,['echo --------------------------------------------------------------------------\n']);
    fprintf(fid,['echo cbGMK',num2str(j),'.err\n']);
    fprintf(fid,['tail -n 20 cbGMK',num2str(j),'.err\n']);
    fprintf(fid,['\n']);
end

fclose(fid);



