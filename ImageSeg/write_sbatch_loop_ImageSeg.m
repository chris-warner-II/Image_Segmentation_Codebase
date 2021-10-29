
method = 7; % Mod SKHAdj
radmax = [1:6,inf];
sigpix = [0.1, 0.2, 0.3, 0.4, 0.5];
imdim = 21;
InputImages = 'GradientBox';


if ~exist('scripts4cluster','dir')
   mkdir('scripts4cluster') 
end


k=0;
for i = 1:numel(radmax)
    for j = 1:numel(sigpix)
        for L = 1:numel(imdim)
        
            k=k+1;

            fid = fopen(['scripts4cluster/script_ImgSeg',num2str(k)], 'w+');

            fprintf(fid,['#!/bin/bash\n']);
            fprintf(fid,['# Job name:\n']);
            fprintf(fid,['#SBATCH --job-name=ImSeg',num2str(k),'\n']);
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
            fprintf(fid,['#SBATCH -o ImSeg',num2str(k),'.out\n']);
            fprintf(fid,['#\n']);
            fprintf(fid,['#SBATCH -e ImSeg',num2str(k),'.err\n']);
            fprintf(fid,['\n']);
            %
            fprintf(fid,['\n']);
            fprintf(fid,['module load matlab/R2013a\n']);
            fprintf(fid,['\n']);
            fprintf(fid,['echo Kuramoto Image Segmentation',num2str(k),'\n']);
            fprintf(fid,['\n']);
            fprintf(fid,['matlab -nosplash -nodesktop -noFigureWindows << EOF\n']);   % was matlab -nodisplay but that stopped working on cluster.
            fprintf(fid,['\n']);
            fprintf(fid,['cd /global/home/users/cwarner/Projects\n']);
            fprintf(fid,['addpath(genpath(pwd))\n']);
            fprintf(fid,['\n']);
            fprintf(fid,['Loop_ImgSegMethodsD(',num2str(method),',',num2str(radmax(i)),',',num2str(sigpix(j)),',',num2str(imdim(L)),',',InputImages,');\n']);
            fprintf(fid,['\n']);
            fprintf(fid,['exit\n']);
            fprintf(fid,['EOF\n']);
            fclose(fid);
        
        end
    end
end







% make a run file that will call sbatch script_HAWSamp##
% Note: can not just run file because I do not have permission on cluster.
%       But I can just copy and paste its contents into command line. Works.
fid = fopen(['scripts4cluster/Run_scripts_ImgSeg'], 'w+');
for i=1:k
    fprintf(fid,['sbatch script_ImgSeg',num2str(i),'\n']);
    fprintf(fid,['sleep 3 \n']);
end

fprintf(fid,['squeue\n']);
fclose(fid);








