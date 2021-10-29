fileType = 'BSDS_patch'; % 'GradientBox'; %   
fileSize = '51x51_ds1';

Rmax = [1:4,inf];
sigD = 1; % set this below
sigP = [0.1, 0.2, 0.3, 0.4];
sigW = [0, 0.01*60];
Kscale = 300;
TiScale = {'1', '0p5', '0p1'};


if ~exist('scripts4cluster','dir')
   mkdir('scripts4cluster') 
end


% exp_SepVPar_avgOvrImgPatches(fileType,fileSize,Tscale,sigP,Rmax,sigD,sigW,Kscale)


k=0;
for i = 1:numel(TiScale) 
    for j = 1:numel(sigP)
        for L = 1:numel(Rmax)
            
            %set sigD = Inf or (Rmax/4)
            sigD = [Rmax(L)/4, Inf];
            sigD = unique(sigD); % for the case when Rmax = Inf.
            
            for m = 1:numel(sigD)
                
                for n = 1:numel(sigW)
                    for o = 1:numel(Kscale)

                        k=k+1;

                        fid = fopen(['scripts4cluster/script_SepVPar_avgOvrImgPatches',num2str(k)], 'w+');

                        fprintf(fid,['#!/bin/bash\n']);
                        fprintf(fid,['# Job name:\n']);
                        fprintf(fid,['#SBATCH --job-name=SvPmn',num2str(k),'\n']);
                        fprintf(fid,['#\n']);
                        fprintf(fid,['# Partition:\n']);
                        fprintf(fid,['#SBATCH --partition=cortex\n']);
                        fprintf(fid,['#\n']);
                        fprintf(fid,['# Processors:\n']);
                        fprintf(fid,['#SBATCH --nodes=1\n']);
                        fprintf(fid,['#\n']);
                        fprintf(fid,['# Memory:\n']);
                        fprintf(fid,['#SBATCH --mem-per-cpu=2G\n']);
                        fprintf(fid,['#\n']);
                        fprintf(fid,['# Wall clock limit:\n']);
                        fprintf(fid,['#SBATCH --time=2-00:0:00\n']);
                        fprintf(fid,['#\n']);
                        fprintf(fid,['#\n']);
                        fprintf(fid,['#SBATCH -o SvPmn',num2str(k),'.out\n']);
                        fprintf(fid,['#\n']);
                        fprintf(fid,['#SBATCH -e SvPmn',num2str(k),'.err\n']);
                        fprintf(fid,['\n']);
                        %
                        fprintf(fid,['\n']);
                        fprintf(fid,['module load matlab/R2013a\n']);
                        fprintf(fid,['\n']);
                        fprintf(fid,['echo SepVPar Average Over Image Patches',num2str(k),'\n']);
                        fprintf(fid,['\n']);
                        fprintf(fid,['matlab -nosplash -nodesktop -noFigureWindows << EOF\n']);   % was matlab -nodisplay but that stopped working on cluster.
                        fprintf(fid,['\n']);
                        fprintf(fid,['cd /global/home/users/cwarner/Projects\n']);
                        fprintf(fid,['addpath(genpath(pwd))\n']);
                        fprintf(fid,['\n']);
                        fprintf(fid,['exp_SepVPar_avgOvrImgPatches(''',fileType,''',''',fileSize,''',''',TiScale{i},''',',...
                                    num2str(sigP(j)),',',num2str(Rmax(L)),',',num2str(sigD(m)),',',num2str(sigW(n)),','...
                                    num2str(Kscale(o)),');\n']);
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







% make a run file that will call sbatch script_HAWSamp##
% Note: can not just run file because I do not have permission on cluster.
%       But I can just copy and paste its contents into command line. Works.
fid = fopen(['scripts4cluster/Run_scripts_SepVPar_avgOvrImgPatches'], 'w+');
for i=1:k
    fprintf(fid,['sbatch script_SepVPar_avgOvrImgPatches',num2str(i),'\n']);
    fprintf(fid,['sleep 1 \n']);
end

fprintf(fid,['squeue\n']);
fclose(fid);






% Make a run file that will display last few lines of all *.out files to look thru quickly
fid = fopen(['scripts4cluster/View_outs_SvPmn'], 'w+');
for i=1:k
    fprintf(fid,['echo --------------------------------------------------------------------------\n']);
    fprintf(fid,['echo SvPmn',num2str(i),'.out\n']);
    fprintf(fid,['tail -n 30 SvPmn',num2str(i),'.out\n']);
    fprintf(fid,['\n']);
end

fclose(fid);


% Make a run file that will display last few lines of all *.err files to look thru quickly
fid = fopen(['scripts4cluster/View_errs_SvPmn'], 'w+');
for i=1:k
    fprintf(fid,['echo --------------------------------------------------------------------------\n']);
    fprintf(fid,['echo SvPmn',num2str(i),'.err\n']);
    fprintf(fid,['tail -n 20 SvPmn',num2str(i),'.err\n']);
    fprintf(fid,['\n']);
end

fclose(fid);









