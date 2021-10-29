method = [6,7,8]; %[3,5,6,7,12]; % 3=AAnrm, 5=GLnrm, 6=Mod_N&G, 7=Mod_SKHAdj, 12=IsoDiff, 8=Mod_SKHEuc
radmax = [1,3,5,10]; %[1:7,inf]; % [1:4,inf];
sigpix = [0.2]; %[0.1, 0.2, 0.3, 0.4];
imdim = {[101,101]}; % ,[101,101]   % {'full'}; %
InputImages = 'BSDS_patch'; %  'GradientBox'; %  'BSDS_tile'; % 'BSDS_full'; %

% Put this here so I can now parallelize by patches extracted from single image if I like.
imgNums2procLoop(:,1) = [1:5:500]; % start image patch number
imgNums2procLoop(:,2) = [1:5:500]+4; % end image patch number


%Kscale = 300;  
SigW = 0; % [0.01 0.03 0.05 0.1].*60; 

% TiScale = 1; % [1, 0.5, 0.1];


if ~exist('scripts4cluster','dir')
   mkdir('scripts4cluster') 
end


k=0;

for H = 1:numel(method)
    for i = 1:numel(radmax)
        for j = 1:numel(sigpix)
            for L = 1:numel(imdim)

                if ischar(imdim{L})
                    imdim_cell = ['{''',imdim{L},'''}'];
                else
                    imdim_cell = ['{[',num2str(imdim{L}(1)),',',num2str(imdim{L}(2)),']}'];
                end

                %for m = 1:numel(Kscale)
                    for n = 1:numel(SigW)
                        %for o = 1:numel(TiScale)
                            for p = 1:size(imgNums2procLoop,1)

                                imgNums2proc = imgNums2procLoop(p,:);

                                k=k+1;

                                fid = fopen(['scripts4cluster/script_ImgSeg',num2str(k)], 'w+');

                                fprintf(fid,['#!/bin/bash\n']);
                                fprintf(fid,['# Job name:\n']);
                                fprintf(fid,['#SBATCH --job-name=ImSeg',num2str(k),'\n']);
                                fprintf(fid,['#\n']);
                                fprintf(fid,['# Partition:\n']);
                                fprintf(fid,['#SBATCH --partition=cortex\n']);
                                fprintf(fid,['#\n']);
                                fprintf(fid,['# Constrain Nodes:\n']);
                                fprintf(fid,['#SBATCH --constraint=cortex_nogpu\n']);
                                fprintf(fid,['#\n']);
                                fprintf(fid,['# Processors:\n']);
                                fprintf(fid,['#SBATCH --ntasks=4\n']);
                                fprintf(fid,['#\n']);
                                %fprintf(fid,['# Exclude the Following Nodes:\n']); % n0000.cortex0, ,n0012.cortex0,n0013.cortex0,n0007.cortex0,n0008.cortex0,n0009.cortex0,n0010.cortex0,n0011.cortex0\n
                                %fprintf(fid,['#SBATCH -x n0002.cortex0']); % ,n0001.cortex0,n0012.cortex0,n0013.cortex0,
                                fprintf(fid,['#\n']);
                                fprintf(fid,['# Memory:\n']);
                                fprintf(fid,['#SBATCH --mem-per-cpu=3500M\n']);
                                fprintf(fid,['#\n']);
                                fprintf(fid,['# Wall clock limit:\n']);
                                fprintf(fid,['#SBATCH --time=3-12:0:00\n']);
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
                                fprintf(fid,['Loop_ImgSegMethodsD(',num2str(method(H)),',',num2str(radmax(i)),',',num2str(sigpix(j)),',',imdim_cell,...
                                            ',''',InputImages,''',[',num2str(imgNums2proc(1)),',',num2str(imgNums2proc(2)),'],',num2str(SigW(n)),');\n']);
                                fprintf(fid,['\n']);
                                fprintf(fid,['exit\n']);
                                fprintf(fid,['EOF\n']);
                                fclose(fid);

                            end
                        %end
                    end
                %end
            end
        end
    end
end







% Make a run file that will call sbatch script_ImgSeg##
% Note: can not just run file because I do not have permission on cluster.
%       But I can just copy and paste its contents into command line. Works.
fid = fopen(['scripts4cluster/Run_scripts_ImgSeg'], 'w+');
for i=1:k
    fprintf(fid,['sbatch script_ImgSeg',num2str(i),'\n']);
    %fprintf(fid,['sleep 3 \n']);
end

fprintf(fid,['squeue\n']);
fclose(fid);



% Make a run file that will display last few lines of all *.out files to look thru quickly
fid = fopen(['scripts4cluster/View_outs_ImSeg'], 'w+');
for i=1:k
    fprintf(fid,['echo --------------------------------------------------------------------------\n']);
    fprintf(fid,['echo ImSeg',num2str(i),'.out\n']);
    fprintf(fid,['tail -n 30 ImSeg',num2str(i),'.out\n']);
    fprintf(fid,['\n']);
end

fclose(fid);


% Make a run file that will display last few lines of all *.err files to look thru quickly
fid = fopen(['scripts4cluster/View_errs_ImSeg'], 'w+');
for i=1:k
    fprintf(fid,['echo --------------------------------------------------------------------------\n']);
    fprintf(fid,['echo ImSeg',num2str(i),'.err\n']);
    fprintf(fid,['tail -n 20 ImSeg',num2str(i),'.err\n']);
    fprintf(fid,['\n']);
end

fclose(fid);



