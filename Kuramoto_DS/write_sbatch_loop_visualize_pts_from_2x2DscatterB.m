
imdim =  {[51,51]}; % {'full'}; %
InputImages = 'BSDS_patch'; % 'GradientBox'; %  'BSDS_tile'; % 'BSDS_full'; %
MethodType = {'Mod_SKHAdj'}; % ,'Mod_N&G', 'GLnrm'
ImPtchSize = [num2str(imdim{1}(1)),'x',num2str(imdim{1}(2)),'_ds1'];


Rmax = [9]; % [1:10,inf];
SigP = [0.1]; % [0.1, 0.2, 0.3, 0.4];

Ks = 300;  
SigW = 0; % [0, 0.01*60]; 

SigD = Inf;

Ts = 1; % [1, 0.5, 0.1];


if ~exist('scripts4cluster','dir')
   mkdir('scripts4cluster') 
end


% Delete Old Versions of scripts, tasks & task files and repopulate them.
!rm scripts4cluster/script_vis_2Dscatter*
!rm scripts4cluster/Run_scripts_visScat
!rm scripts4cluster/View_outs_visScat
!rm scripts4cluster/View_errs_visScat


k=0; 
for h = 1:numel(MethodType)
    for i = 1:numel(Rmax)
        for j = 1:numel(SigD)
            for L = 1:numel(SigP) 
                for m = 1:numel(SigW)
                    for n = 1:numel(Ts)
                        for o = 1:numel(Ks)

                            k=k+1;

                            fid = fopen(['scripts4cluster/script_vis_2Dscatter',num2str(k)], 'w+');

                            fprintf(fid,['#!/bin/bash\n']);
                            fprintf(fid,['# Job name:\n']);
                            fprintf(fid,['#SBATCH --job-name=vSct',num2str(k),'\n']);
                            fprintf(fid,['#\n']);
                            fprintf(fid,['# Partition:\n']);
                            fprintf(fid,['#SBATCH --partition=cortex\n']);
                            fprintf(fid,['#\n']);
                            fprintf(fid,['# Processors:\n']);
                            fprintf(fid,['#SBATCH --nodes=1\n']);
                            fprintf(fid,['#\n']);
                            fprintf(fid,['# Exclude the Following Nodes:\n']);
                            fprintf(fid,['#SBATCH -x n0000.cortex0,n0001.cortex0,n0012.cortex0,n0013.cortex0\n']); % n0007.cortex0,n0008.cortex0,n0009.cortex0,n0010.cortex0,n0011.cortex0,n0001.cortex0,n0012.cortex0,n0013.cortex0,
                            fprintf(fid,['#\n']);
                            fprintf(fid,['# Memory:\n']);
                            fprintf(fid,['#SBATCH --mem-per-cpu=1G\n']);
                            fprintf(fid,['#\n']);
                            fprintf(fid,['# Wall clock limit:\n']);
                            fprintf(fid,['#SBATCH --time=5:0:00\n']);
                            fprintf(fid,['#\n']);
                            fprintf(fid,['#\n']);
                            fprintf(fid,['#SBATCH -o visScat',num2str(k),'.out\n']);
                            fprintf(fid,['#\n']);
                            fprintf(fid,['#SBATCH -e visScat',num2str(k),'.err\n']);
                            fprintf(fid,['\n']);
                            %
                            fprintf(fid,['\n']);
                            fprintf(fid,['module load matlab/R2013a\n']);
                            fprintf(fid,['\n']);
                            fprintf(fid,['echo Visualize Points from 2x2D Scatter Plots ',num2str(k),'\n']);
                            fprintf(fid,['\n']);
                            fprintf(fid,['matlab -nosplash -nodesktop -noFigureWindows << EOF\n']);   % was matlab -nodisplay but that stopped working on cluster.
                            fprintf(fid,['\n']);
                            fprintf(fid,['cd /global/home/users/cwarner/Projects\n']);
                            fprintf(fid,['addpath(genpath(pwd))\n']);
                            fprintf(fid,['\n']);
                            fprintf(fid,['visualize_pts_from_2x2D_scatterB(''',InputImages,''',''',ImPtchSize,''',''',MethodType{h},''',',num2str(Rmax(i)),',',num2str(SigD(j)),...
                                            ',',num2str(SigP(L)),',',num2str(SigW(m)),',',num2str(Ts(n)),',',num2str(Ks(o)),');\n']);
                            fprintf(fid,['\n']);
                            fprintf(fid,['exit\n']);
                            fprintf(fid,['EOF\n']);
                            fclose(fid);

                        end
                    end
                end
            end
        end
    end
end



% Make a run file that will call sbatch script_ImgSeg##
% Note: can not just run file because I do not have permission on cluster.
%       But I can just copy and paste its contents into command line. Works.
fid = fopen(['scripts4cluster/Run_scripts_visScat'], 'w+');
for i=1:k
    fprintf(fid,['sbatch script_vis_2Dscatter',num2str(i),'\n']);
    fprintf(fid,['sleep 1 \n']);
end

fprintf(fid,['squeue\n']);
fclose(fid);



% Make a run file that will display last few lines of all *.out files to look thru quickly
fid = fopen(['scripts4cluster/View_outs_visScat'], 'w+');
for i=1:k
    fprintf(fid,['echo --------------------------------------------------------------------------\n']);
    fprintf(fid,['echo visScat',num2str(i),'.out\n']);
    fprintf(fid,['tail -n 30 visScat',num2str(i),'.out\n']);
    fprintf(fid,['\n']);
end

fclose(fid);


% Make a run file that will display last few lines of all *.err files to look thru quickly
fid = fopen(['scripts4cluster/View_errs_visScat'], 'w+');
for i=1:k
    fprintf(fid,['echo --------------------------------------------------------------------------\n']);
    fprintf(fid,['echo visScat',num2str(i),'.err\n']);
    fprintf(fid,['tail -n 20 visScat',num2str(i),'.err\n']);
    fprintf(fid,['\n']);
end

fclose(fid);



