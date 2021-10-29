% a script/function to look into a directory and plot figures of Methods
% with subplots of Parameters (to compare Parameter values for one method)
% & also, plot figures of Parameter values with subplots of Methods (to 
% compare Methods)

%%
imageInput = 'Waddle';
direcBase = ['./output/ImgSeg/simpleExamples/' imageInput];
comp_flg = 1;

foci = dir([direcBase,'/data from main']);
foci = foci(~[foci.isdir]); % get rid of '.' & '..'
if strmatch(foci(1).name,'.DS_Store')
    foci = foci(2:end); % get rid of the .DS_Store file
end
focus = foci(1).name; % starts at 4 because dir finds '.' '..' & '.Dsstore' or starts at 3 because dir finds only '.' & '..'
x = strfind(focus,'-'); 
km = 1; kp = 1;
beg = focus(1:x-1);
space = find(beg==' ');
method{km} = beg(space(1)+1:end);
params{km} = focus(x+1:end-4);

for j = 1:numel(foci)
    % 
    focus = foci(j).name;
    x = strfind(focus,'-');
    tempM = focus(space(1)+1:x-1);
    tempP = focus(x+1:end-4);
    
    % if this file is for a new setting for method, record it.
    recM=1;
    for i = 1:numel(method)
        if ~isempty(strmatch(method{i},tempM,'exact'))
            recM = 0;
        end
    end
    %
    if(recM)
        km=km+1;
        method{km} = tempM;
    end
    
    % if this file is for a new setting for params, record it.
    recP=1;
    for i = 1:numel(params)
        if ~isempty(strmatch(params{i},tempP,'exact'))
            recP = 0;
        end
    end
    %
    if(recP)
        kp=kp+1;
        params{kp} = tempP;
    end
    
end




%% Loop through methods and compare different parameter values.
comp_flg = 2;
for i = 1:numel(method)
    tag = method{i};
    plot_ImgSeg_Evectors(imageInput,tag,comp_flg)
    close
end


%% Loop through params and compare different methods.
comp_flg = 1;
for i = 1:numel(params)
    tag = params{i};
    plot_ImgSeg_Evectors(imageInput,tag,comp_flg)
    close
end