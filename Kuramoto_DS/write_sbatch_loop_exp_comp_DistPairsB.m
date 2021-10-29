MethodType = {'Mod_N&G','Mod_SKHAdj','GLnrm','AAnrm'};
radmax = [1:10,inf];
sigpix = [0.1,0.2,0.3,0.4]; % 0.1,0.2,0.3,
imdim =  {[51,51]}; % {'full'}; %
InputImages = 'BSDS_patch'; % 'GradientBox'; %  'BSDS_tile'; % 'BSDS_full'; %


ImPtchSize = [num2str(imdim{1}(1)),'x',num2str(imdim{1}(2)),'_ds1'];



jobs_per_cpu = 22; % # combinations of rM & sP
numCPUs = 4;

Ks = 300;  
SigW = 0; %[0.15*60]; % 0.01*60,0.03*60,0.05*60,0.1*60,
SigD = inf;

Ts = 1; % [1, 0.5, 0.1];




rM = zeros( 1,numel(radmax)*numel(sigpix)*numel(MethodType) );
sP = zeros( 1,numel(radmax)*numel(sigpix)*numel(MethodType) );
meths = cell( 1,numel(radmax)*numel(sigpix)*numel(MethodType) );

st=1;
st2=1;
for i = 1:numel(radmax)
    nd = st - 1 + numel(MethodType)*numel(sigpix);
    rM(st:nd) = radmax(i);
    st = nd + 1;
    
    
    for j = 1:numel(sigpix)
        nd2 = st2 - 1 + numel(MethodType);
        sP(st2:nd2) = sigpix(j);
        for B=1:numel(MethodType)
            meths{st2+B-1} = MethodType{B};
        end    
        st2 = nd2 + 1;
    end
    
end


% just a plot for visualization to see if things are lined up right.
if(0)
    figure, hold on,
    plot(1:numel(rM),rM,'b')
    plot(1:numel(sP),sP./max(sP),'r')
end






if ~exist('scripts4cluster','dir')
   mkdir('scripts4cluster') 
end


% Delete Old Versions of scripts, tasks & task files and repopulate them.
!rm scripts4cluster/script_DistPairs*
!rm scripts4cluster/taskfile_DistPairs*
!rm scripts4cluster/task_DistPairs*.sh

for k = 1:ceil( numel(rM)./(numCPUs*jobs_per_cpu) )


    %% (#1). Build up script_ImgSeg file that calls helper_ht & 
    fid = fopen(['scripts4cluster/script_DistPairs',num2str(k)], 'w+');
    fprintf(fid,['#!/bin/bash\n']);
    fprintf(fid,['# Job name:\n']);
    fprintf(fid,['#SBATCH --job-name=DistPW',num2str(k),'\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['# Partition:\n']);
    fprintf(fid,['#SBATCH --partition=cortex\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['# Processors:\n']);
    fprintf(fid,['#SBATCH --ntasks=4\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['# Exclude the Following Nodes:\n']);
    fprintf(fid,['#SBATCH -x n0000.cortex0,n0001.cortex0\n']); % 
    fprintf(fid,['#\n']);
    fprintf(fid,['# Memory:\n']);
    fprintf(fid,['#SBATCH --mem-per-cpu=3500M\n']);
    fprintf(fid,['#\n']);
    fprintf(fid,['# Wall clock limit:\n']);
    fprintf(fid,['#SBATCH --time=5-00:0:00\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['ht_helper.sh -t taskfile_DistPairs',num2str(k),' -n1 -s60 -dvk']);
    fclose(fid);
    

    
    %% (#2). Make the taskfile_ImgSeg# file that pertains to each script_ImgSeg# file
    fid = fopen(['scripts4cluster/taskfile_DistPairs',num2str(k)], 'w+');
    for i = 1:(numCPUs*jobs_per_cpu)
        if( (k-1)*(numCPUs*jobs_per_cpu) + i <= numel(rM) )
            fprintf(fid,['task_DistPairs',num2str( (k-1)*(numCPUs*jobs_per_cpu) + i  ),'.sh\n']);
        end
    end
    fclose(fid);
                            
end                          
                            
                            
                            
for i = 1: numel(rM)                            
                            
    %% (#2.) Build up each individual task_ImgSeg#.sh file
    fid = fopen(['scripts4cluster/task_DistPairs',num2str(i),'.sh'], 'w+');
    fprintf(fid,['#!/bin/bash\n']);
    fprintf(fid,['module load matlab/R2013a\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['echo Explore Compare DistPairs ',num2str(i),'\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['matlab -nosplash -nodesktop -noFigureWindows << EOF\n']);   % was matlab -nodisplay but that stopped working on cluster.
    fprintf(fid,['\n']);
    fprintf(fid,['cd /global/home/users/cwarner/Projects\n']);
    fprintf(fid,['addpath(genpath(pwd))\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['explore_compare_DistPairs(''',InputImages,''',''',ImPtchSize,''',''',meths{i},''',',num2str(rM(i)),',',num2str(SigD),...
                                        ',',num2str(sP(i)),',',num2str(SigW),',',num2str(Ts),',',num2str(Ks),');\n']);
    fprintf(fid,['\n']);
    fprintf(fid,['exit\n']);
    fprintf(fid,['EOF\n']);
    fclose(fid);
                        


end





% % Make a run file that will call sbatch script_ImgSeg##
% % Note: can not just run file because I do not have permission on cluster.
% %       But I can just copy and paste its contents into command line. Works.
% fid = fopen(['scripts4cluster/Run_scripts_ImgSeg'], 'w+');
% for i=1:( numel(rM)./(numCPUs*jobs_per_cpu) )
%     fprintf(fid,['sbatch script_ImgSeg',num2str(i),'\n']);
%     fprintf(fid,['sleep 1 \n']);
% end
% 
% fprintf(fid,['squeue\n']);
% fclose(fid);
% 
% 
% 
% % Make a run file that will display last few lines of all *.out files to look thru quickly
% fid = fopen(['scripts4cluster/View_outs_ImSeg'], 'w+');
% for i=1:k
%     fprintf(fid,['echo --------------------------------------------------------------------------\n']);
%     fprintf(fid,['echo ImSeg',num2str(i),'.out\n']);
%     fprintf(fid,['tail -n 30 ImSeg',num2str(i),'.out\n']);
%     fprintf(fid,['\n']);
% end
% 
% fclose(fid);
% 
% 
% % Make a run file that will display last few lines of all *.err files to look thru quickly
% fid = fopen(['scripts4cluster/View_errs_ImSeg'], 'w+');
% for i=1:k
%     fprintf(fid,['echo --------------------------------------------------------------------------\n']);
%     fprintf(fid,['echo ImSeg',num2str(i),'.err\n']);
%     fprintf(fid,['tail -n 20 ImSeg',num2str(i),'.err\n']);
%     fprintf(fid,['\n']);
% end
% 
% fclose(fid);



