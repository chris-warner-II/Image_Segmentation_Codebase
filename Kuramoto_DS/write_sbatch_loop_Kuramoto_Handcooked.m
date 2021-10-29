% Make sbatch file to call loop_main_Kuramoto_HandCookedNetwork with some
% parameters to loop through

runParams.muW = 60;         % mean of oscillator resonant phase distribution (Hz) 
% %
% loop_sigW = 0.02*runParams.muW; % [0.001*runParams.muW 0.01*runParams.muW 0.1*runParams.muW];
% loop_strng = [10 5 1];
loop_weak = [1:1:10]; %[-100 -50 -20 -10 -1 -0.1 0.1 1]; % [-10 -1 -0.1 -0.01 -0.001 0 0.001 0.01 0.1 1];

                  
% N = 48
loop_Gnd_Truth = {'{[ones(1,24), 2*ones(1,24)]}', ... % C = 2
                  '{[ones(1,16), 2*ones(1,16), 3*ones(1,16)]}', ... % C = 3
                  '{[ones(1,12), 2*ones(1,12), 3*ones(1,12), 4*ones(1,12)]}', ... % C = 4
                  '{[ones(1,8), 2*ones(1,8), 3*ones(1,8), 4*ones(1,8), 5*ones(1,8), 6*ones(1,8)]}', ... % C = 6
                  '{[ones(1,6), 2*ones(1,6), 3*ones(1,6), 4*ones(1,6), 5*ones(1,6), 6*ones(1,6), 7*ones(1,6), 8*ones(1,6)]}'}; %, ... % C = 8
                     
% NOTE:  I WILL NEED TO FIX THIS TO BE A CELL OF CELLS. I THINK IT WILL ERROR WHEN I TRY TO RUN IT ON THE CLUSTER.  JUST BE AWARE.
                      
%                       '{[ones(1,4), 2*ones(1,4), 3*ones(1,4), 4*ones(1,4), 5*ones(1,4), 6*ones(1,4), 7*ones(1,4), 8*ones(1,4), 9*ones(1,4), 10*ones(1,4), 11*ones(1,4), 12*ones(1,4)]}', ... % C = 12
%                       '{[ones(1,3), 2*ones(1,3), 3*ones(1,3), 4*ones(1,3), 5*ones(1,3), 6*ones(1,3), 7*ones(1,3), 8*ones(1,3), 9*ones(1,3), 10*ones(1,3), 11*ones(1,3), 12*ones(1,3), 13*ones(1,3), 14*ones(1,3), 15*ones(1,3), 16*ones(1,3)]}', ... % C = 16
%                       '{[ones(1,2), 2*ones(1,2), 3*ones(1,2), 4*ones(1,2), 5*ones(1,2), 6*ones(1,2), 7*ones(1,2), 8*ones(1,2), 9*ones(1,2), 10*ones(1,2), 11*ones(1,2), 12*ones(1,2), 13*ones(1,2), 14*ones(1,2), 15*ones(1,2), 16*ones(1,2), 17*ones(1,2), 18*ones(1,2), 19*ones(1,2), 20*ones(1,2), 21*ones(1,2), 22*ones(1,2), 23*ones(1,2), 24*ones(1,2)]}', ... % C = 24
%                        };    
                  



if ~exist('scripts4cluster','dir')
   mkdir('scripts4cluster') 
end


k=0;
for i = 1:numel(loop_Gnd_Truth)
    for j = 1:numel(loop_weak)
        
        k=k+1;

        fid = fopen(['scripts4cluster/script_Kuramoto_HandCookedNet',num2str(k)], 'w+');

        fprintf(fid,['#!/bin/bash\n']);
        fprintf(fid,['# Job name:\n']);
        fprintf(fid,['#SBATCH --job-name=KurHC',num2str(k),'\n']);
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
        fprintf(fid,['#SBATCH -o KurHC',num2str(k),'.out\n']);
        fprintf(fid,['#\n']);
        fprintf(fid,['#SBATCH -e KurHC',num2str(k),'.err\n']);
        fprintf(fid,['\n']);
        %
        fprintf(fid,['\n']);
        fprintf(fid,['module load matlab/R2013a\n']);
        fprintf(fid,['\n']);
        fprintf(fid,['echo Kuramoto HandCooked Network',num2str(k),'\n']);
        fprintf(fid,['\n']);
        fprintf(fid,['matlab -nosplash -nodesktop -noFigureWindows << EOF\n']);   % was matlab -nodisplay but that stopped working on cluster.
        fprintf(fid,['\n']);
        fprintf(fid,['cd /global/home/users/cwarner/Projects\n']);
        fprintf(fid,['addpath(genpath(pwd))\n']);
        fprintf(fid,['\n']);
        fprintf(fid,['loop_main_Kuramoto_HandCookedNetwork(',loop_Gnd_Truth{i},',',num2str(loop_weak(j)),');\n']);
        fprintf(fid,['\n']);
        fprintf(fid,['exit\n']);
        fprintf(fid,['EOF\n']);
    %     fprintf(fid,['fi\n']);
        fclose(fid);
    
    end
end







% make a run file that will call sbatch script_HAWSamp##
% Note: can not just run file because I do not have permission on cluster.
%       But I can just copy and paste its contents into command line. Works.
fid = fopen(['scripts4cluster/Run_scripts_Kur_HC'], 'w+');
for i=1:k
    fprintf(fid,['sbatch script_Kuramoto_HandCookedNet',num2str(i),'\n']);
    fprintf(fid,['sleep 3 \n']);
end

fprintf(fid,['squeue\n']);
fclose(fid);