
% Input: Loop through these
netMethod = {'Mod_N&G'}; % 'IsoDiff','AAnrm','GLnrm','Mod_SKHAdj',
KSstr = {'sml','mid','lrg'};

radmax = [5]; %[1,3,5,10];
sigpix = [0.2]; %[0.1, 0.2, 0.3, 0.4];
imdim =  {[101,101]}; %[51,51], {'full'}; % 
InputImages = 'BSDS_patch'; %'BSDS_full'; %  'GradientBox'; %  'BSDS_tile'; % 

SigW = 0; %[0.15*60]; % 0.01*60,0.03*60,0.05*60,0.1*60,
sigDist = inf;


if ~exist('scripts4cluster','dir')
   mkdir('scripts4cluster') 
end


k=0;

for H = 1:numel(netMethod)
    for i = 1:numel(radmax)
        for j = 1:numel(sigpix)
            for L = 1:numel(imdim)

                imdim_cell = [num2str(imdim{L}(1)),'x',num2str(imdim{L}(2)),'_ds1'];

                for m = 1:numel(KSstr)
                    for n = 1:numel(SigW)

                        k=k+1;

                        fid = fopen(['scripts4cluster/script_genROCAUC',num2str(k)], 'w+');

                        fprintf(fid,['#!/bin/bash\n']);
                        fprintf(fid,['# Job name:\n']);
                        fprintf(fid,['#SBATCH --job-name=AUC',num2str(k),'\n']);
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
                        fprintf(fid,['#SBATCH --time=10:0:00\n']);
                        fprintf(fid,['#\n']);
                        fprintf(fid,['#\n']);
                        fprintf(fid,['#SBATCH -o AUC',num2str(k),'.out\n']);
                        fprintf(fid,['#\n']);
                        fprintf(fid,['#SBATCH -e AUC',num2str(k),'.err\n']);
                        fprintf(fid,['\n']);
                        %
                        fprintf(fid,['\n']);
                        fprintf(fid,['module load matlab/R2013a\n']);
                        fprintf(fid,['\n']);
                        fprintf(fid,['echo ROC AUC accretion',num2str(k),'\n']);
                        fprintf(fid,['date\n']);
                        fprintf(fid,['\n']);
                        fprintf(fid,['matlab -nosplash -nodesktop -noFigureWindows << EOF\n']);   % was matlab -nodisplay but that stopped working on cluster.
                        fprintf(fid,['\n']);
                        fprintf(fid,['echo Made it past\n']);
                        fprintf(fid,['datestr(clock)\n']);
                        fprintf(fid,['\n']);
                        fprintf(fid,['cd /global/home/users/cwarner/Projects\n']);
                        fprintf(fid,['\n']);
                        fprintf(fid,['addpath(genpath(pwd))\n']);
                        fprintf(fid,['\n']);
                        fprintf(fid,['gen_ROC_AUC_scatter(''',InputImages,''',''',imdim_cell,''',{''',netMethod{H},'''},',num2str(radmax(i)),',',...
                                    num2str(sigDist),',',num2str(sigpix(j)),',',num2str(SigW(n)),',''',KSstr{m},''');\n']);
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



% Make a run file that will call sbatch script_ImgSeg##
% Note: can not just run file because I do not have permission on cluster.
%       But I can just copy and paste its contents into command line. Works.
fid = fopen(['scripts4cluster/Run_scripts_genROCAUC'], 'w+');
for i=1:k
    fprintf(fid,['sbatch script_genROCAUC',num2str(i),'\n']);
    fprintf(fid,['sleep 1 \n']);
end

fprintf(fid,['squeue\n']);
fclose(fid);



% Make a run file that will display last few lines of all *.out files to look thru quickly
fid = fopen(['scripts4cluster/View_outs_genROCAUC'], 'w+');
for i=1:k
    fprintf(fid,['echo --------------------------------------------------------------------------\n']);
    fprintf(fid,['echo AUC',num2str(i),'.out\n']);
    fprintf(fid,['tail -n 30 AUC',num2str(i),'.out\n']);
    fprintf(fid,['\n']);
end

fclose(fid);


% Make a run file that will display last few lines of all *.err files to look thru quickly
fid = fopen(['scripts4cluster/View_errs_genROCAUC'], 'w+');
for i=1:k
    fprintf(fid,['echo --------------------------------------------------------------------------\n']);
    fprintf(fid,['echo AUC',num2str(i),'.err\n']);
    fprintf(fid,['tail -n 20 AUC',num2str(i),'.err\n']);
    fprintf(fid,['\n']);
end

fclose(fid);



