method = 5; % Mod SKHAdj
radmax = 1; %[1,4]; %[1:7,inf]; % [1:4,inf];
sigpix = 0.2; %[0.1, 0.2, 0.3, 0.4];
imdim =  {[51,51]}; % {'full'}; %
InputImages = 'BSDS_patch'; % 'GradientBox'; %  'BSDS_tile'; % 'BSDS_full'; %

% Put this here so I can now parallelize by patches extracted from single image if I like.
imgNums2procLoop(:,1) = [1:200:5000]; %1:70; % start image patch number
imgNums2procLoop(:,2) = [200:200:5000]; %1:70; % end image patch number


Kscale = 300;  
SigW = 0; % [0, 0.01*60]; 

TiScale = 1; % [1, 0.5, 0.1];


if ~exist('scripts4nersc','dir')
   mkdir('scripts4nersc') 
   
end


k=0;
for i = 1:numel(radmax)
    for j = 1:numel(sigpix)
        for L = 1:numel(imdim)
            
            if ischar(imdim{L})
                imdim_cell = ['{''',imdim{L},'''}'];
            else
                imdim_cell = ['{[',num2str(imdim{L}(1)),',',num2str(imdim{L}(2)),']}'];
            end
            
            for m = 1:numel(Kscale)
                for n = 1:numel(SigW)
                    for o = 1:numel(TiScale)
                        for p = 1:size(imgNums2procLoop,1)
                            
                            imgNums2proc = imgNums2procLoop(p,:);

                            k=k+1;

                            fid = fopen(['scripts4nersc/script_ImgSeg',num2str(k)], 'w+');

                            fprintf(fid,['#!/bin/bash -l\n']);
                            fprintf(fid,['#PBS -q debug\n']);                     % can be debug, regular, premium
                            fprintf(fid,['#PBS -l mppwidth=24\n']);               % should be some multiple of 24.  24=1core.
                            fprintf(fid,['#PBS -l walltime=00:10:00\n']);        % wall time - max time alotted for job.
                            fprintf(fid,['#PBS -N ImSeg',num2str(k),'\n']);       % Job name
                            fprintf(fid,['#PBS -o ImSeg',num2str(k),'.out\n']);   % output file name
                            fprintf(fid,['#PBS -e ImSeg',num2str(k),'.err\n']);   % error file name
                            fprintf(fid,['#PBS -V']);                             % I dunno what this does.
                            %
                            fprintf(fid,['\n']);
                            fprintf(fid,['module load matlab\n']);
                            fprintf(fid,['\n']);
                            fprintf(fid,['echo Kuramoto Image Segmentation',num2str(k),'\n']);
                            fprintf(fid,['\n']);
                            fprintf(fid,['matlab -nosplash -nodesktop << EOF\n']);
                            fprintf(fid,['\n']);
                            fprintf(fid,['cd /global/u2/c/cwarner/Projects\n']);
                            fprintf(fid,['addpath(genpath(pwd))\n']);
                            fprintf(fid,['\n']);
                            fprintf(fid,['Loop_ImgSegMethodsD(',num2str(method),',',num2str(radmax(i)),',',num2str(sigpix(j)),...
                                        ',',imdim_cell,',''',InputImages,''',[',num2str(imgNums2proc(1)),',',num2str(imgNums2proc(2)),'],',num2str(Kscale(m)),...
                                        ',',num2str(SigW(n)),',',num2str(TiScale(o)),');\n']);
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
fid = fopen(['scripts4nersc/Run_scripts_ImgSeg'], 'w+');
for i=1:k
    fprintf(fid,['qsub script_ImgSeg',num2str(i),'\n']);
    fprintf(fid,['sleep 1 \n']);
end

fprintf(fid,['qstat\n']);
fclose(fid);



% Make a run file that will display last few lines of all *.out files to look thru quickly
fid = fopen(['scripts4nersc/View_outs_ImSeg'], 'w+');
for i=1:k
    fprintf(fid,['echo --------------------------------------------------------------------------\n']);
    fprintf(fid,['echo ImSeg',num2str(i),'.out\n']);
    fprintf(fid,['tail -n 30 ImSeg',num2str(i),'.out\n']);
    fprintf(fid,['\n']);
end

fclose(fid);


% Make a run file that will display last few lines of all *.err files to look thru quickly
fid = fopen(['scripts4nersc/View_errs_ImSeg'], 'w+');
for i=1:k
    fprintf(fid,['echo --------------------------------------------------------------------------\n']);
    fprintf(fid,['echo ImSeg',num2str(i),'.err\n']);
    fprintf(fid,['tail -n 20 ImSeg',num2str(i),'.err\n']);
    fprintf(fid,['\n']);
end

fclose(fid);



