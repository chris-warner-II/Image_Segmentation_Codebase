function [] = plot_ImgSeg_Evectors(imageInput,tag,comp_flg)

% syntax: plot_ImgSeg_Evectors(tag,comp_flg);
%
% This function will plot Eigenvectors in one figure for comparison.
%
% tag = string like 'Modu_N&G' to compare params for one method
%                or 'sP0p3 sD1p5 rM2' to compare methods for one param setting
%
% comp_flg = 1 to compare Methods (for 1 parameter value combination)
%          = 2 to compare Parameters (for 1 method)


% These will need to be changed depending on what your input came from.
% That is, what images you are wanting to investigate the evectors of.
direcBase = ['./output/ImgSeg/simpleExamples/',imageInput];

vizSlope = 1e-15; % steepness parameter for the visualization of Eigenvectors


%% Organize the directories by method and parameter values
foci = dir([direcBase,'/data from main/*',tag,'*']); % tag allows you to look at subset of files in the directory.
foci = foci(~[foci.isdir]); % get rid of '.' & '..'
if strmatch(foci(1).name,'.DS_Store')
    foci = foci(2:end); % get rid of the .DS_Store file
end


%% Look at first file to have baseline to compare against
k=1;
focus = foci(k).name; % starts at 4 because dir finds '.' '..' & '.Dsstore' or starts at 3 because dir finds only '.' & '..'
[imageF,method,params] = parse_Fname(focus);
%
switch(comp_flg)
    case(1) % compare Methods
        comparison{k} = method;
        fixed = params;
        fnameOut = 'Params';
        dnameOut = 'Method';
    case(2) % compare Parameters
        comparison{k} = params;
        fixed = method;
        fnameOut = 'Method';
        dnameOut = 'Params';
end

%% Loop through files and record the value of what you are comparing.
for j = 1:numel(foci)
    % 
    focus = foci(j).name;
    [imageF,method,params] = parse_Fname(focus);
    
    switch(comp_flg)
        case(1) % compare Methods
            temp = method;
        case(2) % compare Parameters
            temp = params;
    end
    
    % if this file is for a new setting (method or param), record it.
    rec=1;
    for i = 1:numel(comparison)
        if ~isempty(strmatch(comparison{i},temp))
            rec = 0;
        end
    end
    %
    if(rec)
        k=k+1;
        comparison{k} = temp;
    end
    %
end

%% Plot Eigenvectors for varying 'comparison' values in a single figure for 'fixed' value.
h = figure; hold on
xsub = round(sqrt(numel(comparison)+1));
ysub = ceil(sqrt(numel(comparison)+1));

for j = 1:numel(comparison)

    % load file, plot and label eigenvector in subplot.
    load([direcBase,'/data from main/',foci(j).name]) % ,'/',fnameIn
    tit = [comparison{j}];
    
    % Calculate how modular or piecewise the eigenvector is by ratio
    H = hist(EvecML(:));
    PWrat = max(H)./min(H);
    
    subplot(ysub,xsub,j)
    imagesc(EvecVizF(EvecML,vizSlope))
    axis off
    title([tit,' - ',num2str(PWrat,'%0.2f')],'FontSize',12,'FontWeight','Bold')
    
end

% Plot original image for comparison
subplot(ysub,xsub,numel(comparison)+1)
imagesc(im) % note: either im or img (I want downsampled or patched one)
xlabel(fixed,'FontSize',12,'FontWeight','Bold')
set(gca,'xtick',[],'ytick',[]);
title(['Image ',imageF],'FontSize',12,'FontWeight','Bold')
set(h,'Position',get(0,'ScreenSize'))

% Save Plot
if ~exist([direcBase,'/data Evecs comp/',dnameOut],'dir')
	mkdir([direcBase,'/data Evecs comp/',dnameOut]);
end      
ftit = fixed;
saveas(h,[direcBase,'/data Evecs comp/',dnameOut,'/',fnameOut,'_',ftit],'jpg');





