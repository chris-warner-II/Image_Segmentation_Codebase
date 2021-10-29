
% function [] = explore_Separation_vs_Parameters(fileGeneral,fileSize,fileSubset,tscale)

% Inputs to function must be of this form
fileGeneral =  '''BSDS_patch'''; % 'GradientBox'; %
fileSize = '''51x51_ds1''';
tscale = {'''0p1''','''0p5''','''1'''};           % scale for phase initialization of oscillators based on input image
% fileSubset = '_'; % '_100075_ptch7_';  %          % certain BSDS image & patch.  Can leave blank too to get all files in dir..




% For BSDS patches version, construct a vector of image reference number (ie. 100075) and patch numbers (ie. 1-10)
numPatches = 10;
numImages = 500;

% iids = imgList('all');

files = dir(['./images/',fileGeneral(2:end-1),'/',fileSize(2:end-1),'/*.mat']);


for i = 1:numImages
    for j = 1:numPatches
        fileSubset{i,j} = ['''_',files(numPatches*(i-1)+j).name(1:end-4),'_'''];  
    end
end





if ~exist('scripts4cluster','dir')
   mkdir('scripts4cluster') 
end


for i = 1:numImages

    fid = fopen(['scripts4cluster/script_explore_SepVPar',num2str(i)], 'w+');

    fprintf(fid,['#!/bin/bash\n']);
    fprintf(fid,['# Job name:\n']);
    fprintf(fid,['#SBATCH --job-name=SvP',num2str(i),'\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['# Partition:\n']);
    fprintf(fid,['#SBATCH --partition=cortex\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['# Processors:\n']);
    fprintf(fid,['#SBATCH --nodes=1\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['# Memory:\n']);
    fprintf(fid,['#SBATCH --mem-per-cpu=1500M\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['# Wall clock limit:\n']);
    fprintf(fid,['#SBATCH --time=2-00:0:00\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['#SBATCH -o SepVPar',num2str(i),'.out\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['#SBATCH -e SepVPar',num2str(i),'.err\n']);
    fprintf(fid,['\n']);
    %
    fprintf(fid,['\n']);
    fprintf(fid,['module load matlab/R2013a\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['echo Exploring Separation Vs Parameters for BSDS Image #',num2str(i),'\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['matlab -nosplash -nodesktop -noFigureWindows << EOF\n']);   % was matlab -nodisplay but that stopped working on cluster.
    fprintf(fid,['\n']);
    fprintf(fid,['cd /global/home/users/cwarner/Projects\n']);
    fprintf(fid,['addpath(genpath(pwd))\n']);
    fprintf(fid,['\n']);
    for j = 1:numPatches
        for k = 1:numel(tscale)
            fprintf(fid,['explore_Separation_vs_Parameters(',fileGeneral,',',fileSize,',',fileSubset{i,j},',',tscale{k},');\n']);
        end

        fprintf(fid,['\n']);
    end

    fprintf(fid,['exit\n']);
    fprintf(fid,['EOF\n']);
    fclose(fid);


end







% make a run file that will call sbatch script_HAWSamp##
% Note: can not just run file because I do not have permission on cluster.
%       But I can just copy and paste its contents into command line. Works.
fid = fopen(['scripts4cluster/Run_scripts_SepVPar'], 'w+');
for i=1:numImages
    fprintf(fid,['sbatch script_explore_SepVPar',num2str(i),'\n']);
    fprintf(fid,['sleep 1 \n']);
end

fprintf(fid,['squeue\n']);
fclose(fid);




% Make a run file that will display last few lines of all *.out files to look thru quickly
fid = fopen(['scripts4cluster/View_outs_SepVPar'], 'w+');
for i=1:numImages
    fprintf(fid,['echo --------------------------------------------------------------------------\n']);
    fprintf(fid,['echo SepVPar',num2str(i),'.out\n']);
    fprintf(fid,['tail -n 30 SepVPar',num2str(i),'.out\n']);
    fprintf(fid,['\n']);
end

fclose(fid);


% Make a run file that will display last few lines of all *.err files to look thru quickly
fid = fopen(['scripts4cluster/View_errs_SepVPar'], 'w+');
for i=1:numImages
    fprintf(fid,['echo --------------------------------------------------------------------------\n']);
    fprintf(fid,['echo SepVPar',num2str(i),'.err\n']);
    fprintf(fid,['tail -n 20 SepVPar',num2str(i),'.err\n']);
    fprintf(fid,['\n']);
end

fclose(fid);





