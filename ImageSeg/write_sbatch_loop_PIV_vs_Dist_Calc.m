%% Look in Ground Truths Directory to find what image patches
[dirPre,sizeGoodIm] = onCluster;

rmin     = [31];
rmax     = [40];
r_range  = [rmin; rmax];

im_st    = 1:50:1500;
im_fin   = 50:50:1500;
im_range = [im_st; im_fin];

ds  = num2str(1);
sz  = num2str(51);
sP  = num2str(0.2);
sD  = num2str(inf);
plt = num2str(0);


if ~exist('scripts4cluster','dir')
   mkdir('scripts4cluster') 
end

C = 0; % loop variable

%% Loop thru r_range and im_range and make scripts to run on cluster nodes. 
for j = 1:size(r_range,2)
    for i = 1:size(im_range,2)
        
        C = C+1;

        fid = fopen(['scripts4cluster/script_PIV_vs_Dist',num2str(C)], 'w+');

        fprintf(fid,['#!/bin/bash\n']);
        fprintf(fid,['# Job name:\n']);
        fprintf(fid,['#SBATCH --job-name=PvD',num2str(C),'\n']);
        fprintf(fid,['#\n']);
        fprintf(fid,['# Partition:\n']);
        fprintf(fid,['#SBATCH --partition=cortex\n']);
        fprintf(fid,['#\n']);
%         fprintf(fid,['# Constrain Nodes:\n']);
%         fprintf(fid,['#SBATCH --constraint=cortex_nogpu\n']);
        fprintf(fid,['#\n']);
        fprintf(fid,['# Processors:\n']);
        fprintf(fid,['#SBATCH --ntasks=1\n']);
        %fprintf(fid,['#\n']);
        %fprintf(fid,['# Exclude the Following Nodes:\n']);
        %fprintf(fid,['#SBATCH -x n0000.cortex0,n0001.cortex0,n0012.cortex0,n0013.cortex0,n0007.cortex0,n0008.cortex0,n0009.cortex0,n0010.cortex0,n0011.cortex0\n']); % ,n0001.cortex0,n0012.cortex0,n0013.cortex0,
        fprintf(fid,['#\n']);
        fprintf(fid,['# Memory:\n']);
        fprintf(fid,['#SBATCH --mem-per-cpu=15G\n']);
        fprintf(fid,['#\n']);
        fprintf(fid,['# Wall clock limit:\n']);
        fprintf(fid,['#SBATCH --time=76:0:00\n']);
        fprintf(fid,['#\n']);
        fprintf(fid,['#\n']);
        fprintf(fid,['#SBATCH -o PvD',num2str(C),'.out\n']);
        fprintf(fid,['#\n']);
        fprintf(fid,['#SBATCH -e PvD',num2str(C),'.err\n']);
        fprintf(fid,['\n']);
        %
        fprintf(fid,['\n']);
        fprintf(fid,['module load matlab/R2013a\n']);
        fprintf(fid,['\n']);
        fprintf(fid,['echo "Calculating Pixel (PIV) Similarity vs Separation Distance fer image patches ',...
                        num2str(im_range(1,i)),'-',num2str(im_range(2,i)),' and r = [',num2str(r_range(1,j)),'-',num2str(r_range(2,j)),']"\n']);
        fprintf(fid,['date\n']);
        fprintf(fid,['\n']);
        fprintf(fid,['matlab -nosplash -nodesktop -noFigureWindows << EOF\n']);   % was matlab -nodisplay but that stopped working on cluster.
        fprintf(fid,['\n']);
        fprintf(fid,['disp(''Made it into Matlab'')\n']);
        fprintf(fid,['datestr(clock)\n']);
        fprintf(fid,['\n']);
        fprintf(fid,['cd /global/home/users/cwarner/Projects\n']);
        fprintf(fid,['\n']);
        fprintf(fid,['addpath(genpath(pwd))\n']);
        fprintf(fid,['\n']);
        fprintf(fid,['BSDS_PIV_vs_Dist_Stats_Calc(',num2str(im_range(1,i)),',',num2str(im_range(2,i)),',',...
            ds,',',sz,',',sP,',',sD,',','[',num2str(r_range(1,j)),',',num2str(r_range(2,j)),'],',plt,');\n']);
        fprintf(fid,['\n']);
        fprintf(fid,['exit\n']);
        fprintf(fid,['EOF\n']);
        fclose(fid);

    end
end



% Make a run file that will call sbatch script_ImgSeg##
% Note: can not just run file because I do not have permission on cluster.
%       But I can just copy and paste its contents into command line. Works.
fid = fopen(['scripts4cluster/Run_scripts_PIV_vs_Dist'], 'w+');
for j=1:C
    fprintf(fid,['sbatch script_PIV_vs_Dist',num2str(j),'\n']);
    fprintf(fid,['sleep 1 \n']);
end

fprintf(fid,['squeue\n']);
fclose(fid);



% Make a run file that will display last few lines of all *.out files to look thru quickly
fid = fopen(['scripts4cluster/View_outs_PIV_vs_Dist'], 'w+');
for j=1:C
    fprintf(fid,['echo --------------------------------------------------------------------------\n']);
    fprintf(fid,['echo PvD',num2str(j),'.out\n']);
    fprintf(fid,['tail -n 30 PvD',num2str(j),'.out\n']);
    fprintf(fid,['\n']);
end

fclose(fid);


% Make a run file that will display last few lines of all *.err files to look thru quickly
fid = fopen(['scripts4cluster/View_errs_PIV_vs_Dist'], 'w+');
for j=1:C
    fprintf(fid,['echo --------------------------------------------------------------------------\n']);
    fprintf(fid,['echo PvD',num2str(j),'.err\n']);
    fprintf(fid,['tail -n 20 PvD',num2str(j),'.err\n']);
    fprintf(fid,['\n']);
end

fclose(fid);



