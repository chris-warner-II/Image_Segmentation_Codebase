
% imdim =  {[101,101]}; % {'full'}; %
% InputImages = 'BSDS_patch'; % 'GradientBox'; %  'BSDS_tile'; % 'BSDS_full'; %

% Put this here so I can now parallelize by patches extracted from single image if I like.
imgNums2procLoop(:,1) = [1:20:1500]; %1:70; % start image patch number
imgNums2procLoop(:,2) = [20:20:1500]; %1:70; % end image patch number

numIter = 100;
plots = 0; % dont make this 1 unless I write up how to save them inside the function!
outFileCheck = 1; % flag to check if output file exists before regenerating it.

if ~exist('scripts4cluster','dir')
   mkdir('scripts4cluster') 
end


% Delete Old Versions of scripts, tasks & task files and repopulate them.
!rm scripts4cluster/script_DoutIdeal*
!rm scripts4cluster/Run_scripts_DoutIdeal
!rm scripts4cluster/View_outs_DoutIdeal
!rm scripts4cluster/View_errs_DoutIdeal

for k = 1:size(imgNums2procLoop,1)

    imgNums2proc = imgNums2procLoop(k,:);

    fid = fopen(['scripts4cluster/script_DoutIdeal',num2str(k)], 'w+');
    fprintf(fid,['#!/bin/bash\n']);
    fprintf(fid,['# Job name:\n']);
    fprintf(fid,['#SBATCH --job-name=Dout',num2str(k),'\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['# Partition:\n']);
    fprintf(fid,['#SBATCH --partition=cortex\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['# Processors:\n']);
    fprintf(fid,['#SBATCH --nodes=1\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['# Exclude the Following Nodes:\n']);
    fprintf(fid,['#SBATCH -x n0000.cortex0,n0001.cortex0,n0012.cortex0,n0013.cortex0\n']); % n0007.cortex0,n0009.cortex0,n0010.cortex0,n0011.cortex0,n0001.cortex0,n0012.cortex0,n0013.cortex0,
    fprintf(fid,['#\n']);
    fprintf(fid,['# Memory:\n']);
    fprintf(fid,['#SBATCH --mem-per-cpu=3500M\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['# Wall clock limit:\n']);
    fprintf(fid,['#SBATCH --time=3:0:00\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['#SBATCH -o Dout',num2str(k),'.out\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['#SBATCH -e Dout',num2str(k),'.err\n']);
    fprintf(fid,['\n']);
    %
    fprintf(fid,['\n']);
    fprintf(fid,['module load matlab/R2013a\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['echo Computing Dout Ideal with Optimization',num2str(k),'\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['matlab -nosplash -nodesktop -noFigureWindows << EOF\n']);   % was matlab -nodisplay but that stopped working on cluster.
    fprintf(fid,['\n']);
    fprintf(fid,['cd /global/home/users/cwarner/Projects\n']);
    fprintf(fid,['addpath(genpath(pwd))\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['compute_Dout_ideal(',num2str(imgNums2proc(1)),',',num2str(imgNums2proc(2)),',',num2str(numIter),',',num2str(plots),',',num2str(outFileCheck),');\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['exit\n']);
    fprintf(fid,['EOF\n']);
    fclose(fid);

end







% Make a run file that will call sbatch script_ImgSeg##
% Note: can not just run file because I do not have permission on cluster.
%       But I can just copy and paste its contents into command line. Works.
fid = fopen(['scripts4cluster/Run_scripts_DoutIdeal'], 'w+');
for i=1:k
    fprintf(fid,['sbatch script_DoutIdeal',num2str(i),'\n']);
    fprintf(fid,['sleep 1 \n']);
end

fprintf(fid,['squeue\n']);
fclose(fid);



% Make a run file that will display last few lines of all *.out files to look thru quickly
fid = fopen(['scripts4cluster/View_outs_DoutIdeal'], 'w+');
for i=1:k
    fprintf(fid,['echo --------------------------------------------------------------------------\n']);
    fprintf(fid,['echo ImSeg',num2str(i),'.out\n']);
    fprintf(fid,['tail -n 30 Dout',num2str(i),'.out\n']);
    fprintf(fid,['\n']);
end

fclose(fid);


% Make a run file that will display last few lines of all *.err files to look thru quickly
fid = fopen(['scripts4cluster/View_errs_DoutIdeal'], 'w+');
for i=1:k
    fprintf(fid,['echo --------------------------------------------------------------------------\n']);
    fprintf(fid,['echo ImSeg',num2str(i),'.err\n']);
    fprintf(fid,['tail -n 20 Dout',num2str(i),'.err\n']);
    fprintf(fid,['\n']);
end

fclose(fid);



